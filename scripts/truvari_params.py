"""Build `truvari bench` CLI parameters from named config profiles (issue #194).

Truvari bench parameters used to be hardcoded in the run_truvari rule. They are
now selected per-analysis via the analyses-table `eval_params` column, which
names a profile under `truvari_bench_params` in resources.yml. The `default`
profile reproduces the historical hardcoded parameters, so existing analyses are
unchanged.

This module was developed with assistance from Claude (Anthropic). All code has
been reviewed and tested by the primary author.
"""

from typing import Dict

# Parameters that take a value, emitted in this order. Maps the profile key to
# the truvari bench flag.
_VALUE_FLAGS = [
    ("pick", "--pick"),
    ("refdist", "-r"),
    ("chunksize", "-C"),
    ("pctseq", "-p"),
    ("pctsize", "-P"),
    ("pctovl", "-O"),
    ("sizemin", "--sizemin"),
    ("sizemax", "--sizemax"),
    ("sizefilt", "--sizefilt"),
]

# Boolean flags (emitted when truthy in the profile).
_BOOL_FLAGS = [
    ("passonly", "--passonly"),
    ("typeignore", "--typeignore"),
    ("dup_to_ins", "--dup-to-ins"),
]

# Reproduces the parameters that were previously hardcoded in run_truvari.
DEFAULT_PROFILE = {"pick": "ac", "passonly": True, "refdist": 2000, "chunksize": 5000}


def build_truvari_bench_params(profile_name: str, config: Dict) -> str:
    """Return the truvari bench CLI fragment for a named profile.

    Order: pick, passonly, then the remaining value/bool flags. The `default`
    profile yields `--pick ac --passonly -r 2000 -C 5000`.
    """
    profiles = config.get("truvari_bench_params", {}) or {}
    if profile_name in profiles:
        params = profiles[profile_name]
    elif profile_name == "default":
        params = DEFAULT_PROFILE
    else:
        available = ", ".join(sorted(profiles)) or "(none defined)"
        raise KeyError(
            f"truvari_bench_params profile '{profile_name}' not found in "
            f"resources.yml; available: {available}"
        )

    parts = []
    # Keep pick + passonly first so `default` matches the historical command.
    if params.get("pick") is not None:
        parts.append(f"--pick {params['pick']}")
    if params.get("passonly"):
        parts.append("--passonly")
    for key, flag in _VALUE_FLAGS:
        if key == "pick":
            continue
        if params.get(key) is not None:
            parts.append(f"{flag} {params[key]}")
    for key, flag in _BOOL_FLAGS:
        if key == "passonly":
            continue
        if params.get(key):
            parts.append(flag)
    return " ".join(parts)
