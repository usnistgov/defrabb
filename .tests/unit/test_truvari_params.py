"""Tests for scripts/truvari_params.py — truvari bench param profiles (#194).

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from truvari_params import build_truvari_bench_params  # noqa: E402


def test_default_matches_historical_hardcoded_command():
    # default must reproduce the parameters previously hardcoded in run_truvari
    assert (
        build_truvari_bench_params("default", {})
        == "--pick ac --passonly -r 2000 -C 5000"
    )


def test_default_from_config_block():
    config = {
        "truvari_bench_params": {
            "default": {
                "pick": "ac",
                "passonly": True,
                "refdist": 2000,
                "chunksize": 5000,
            }
        }
    }
    assert (
        build_truvari_bench_params("default", config)
        == "--pick ac --passonly -r 2000 -C 5000"
    )


def test_custom_profile_emits_tuned_flags():
    config = {
        "truvari_bench_params": {
            "tuned": {
                "pick": "single",
                "passonly": True,
                "refdist": 3000,
                "chunksize": 8000,
                "pctseq": 0.8,
                "sizemax": 100000,
                "typeignore": True,
            }
        }
    }
    out = build_truvari_bench_params("tuned", config)
    assert out.startswith("--pick single --passonly")
    for frag in ("-r 3000", "-C 8000", "-p 0.8", "--sizemax 100000", "--typeignore"):
        assert frag in out


def test_passonly_false_omits_flag():
    config = {"truvari_bench_params": {"p": {"pick": "ac", "passonly": False}}}
    assert build_truvari_bench_params("p", config) == "--pick ac"


def test_bnddist_negative_one_emitted():
    # bnddist=-1 disables BND distance matching (Q100 SV benchmark, #194).
    # -1 is falsy-adjacent but not None, so it must still be emitted.
    config = {
        "truvari_bench_params": {
            "stvar_v5": {
                "pick": "ac",
                "passonly": True,
                "refdist": 2000,
                "chunksize": 5000,
                "bnddist": -1,
            }
        }
    }
    assert (
        build_truvari_bench_params("stvar_v5", config)
        == "--pick ac --passonly -r 2000 -C 5000 -B -1"
    )


def test_unknown_profile_raises():
    with pytest.raises(KeyError):
        build_truvari_bench_params("nope", {"truvari_bench_params": {"default": {}}})


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
