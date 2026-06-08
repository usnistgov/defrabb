"""Structural tests guarding the truvari_remap conda env fix (GitLab #200).

`truvari anno remap` needs a working bwapy with the `-a` option. bioconda's
bwapy 0.1.4 lacks `-a` and fails to load on modern glibc with
`undefined symbol: __log_finite`, so the env installs a patched fork via pip and
must carry the build toolchain to compile it. These tests fail loudly if that
wiring is dropped.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import unittest
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
REMAP_ENV = REPO_ROOT / "envs" / "truvari_remap.yml"


def _load_env():
    return yaml.safe_load(REMAP_ENV.read_text())


def _conda_deps(env):
    """Conda (non-pip) dependency strings."""
    return [d for d in env["dependencies"] if isinstance(d, str)]


def _pip_deps(env):
    """Flattened pip dependency strings, if any."""
    for d in env["dependencies"]:
        if isinstance(d, dict) and "pip" in d:
            return d["pip"]
    return []


class TruvariRemapEnvTests(unittest.TestCase):
    def test_env_file_exists(self):
        self.assertTrue(REMAP_ENV.exists(), f"missing {REMAP_ENV}")

    def test_installs_patched_bwapy_fork_via_pip(self):
        """bwapy comes from the patched git fork, pinned to an immutable ref."""
        pip_deps = _pip_deps(_load_env())
        bwapy = [d for d in pip_deps if "bwapy" in d]
        self.assertEqual(
            len(bwapy), 1, f"expected exactly one bwapy pip dep, got {pip_deps}"
        )
        spec = bwapy[0]
        self.assertIn("git+https://", spec, "bwapy must install from the git fork")
        self.assertIn("bwapy", spec)
        # Pinned to a tag or commit (contains '@'), not a floating branch tip.
        self.assertIn("@", spec, f"bwapy ref must be pinned: {spec}")

    def test_build_toolchain_present(self):
        """pip compiles bwapy (and bundled bwa) from source, so the env needs
        a compiler, make, and zlib."""
        conda_deps = " ".join(_conda_deps(_load_env()))
        for tool in ("c-compiler", "make", "zlib"):
            self.assertIn(tool, conda_deps, f"build toolchain missing: {tool}")

    def test_truvari_pinned(self):
        conda_deps = _conda_deps(_load_env())
        self.assertTrue(
            any("truvari ==5.4.0" in d for d in conda_deps),
            f"truvari must be pinned to 5.4.0: {conda_deps}",
        )

    def test_remap_rule_uses_this_env(self):
        anno_rules = (REPO_ROOT / "rules" / "bench_vcf_anno.smk").read_text()
        remap_rule = anno_rules.split("rule run_truvari_anno_remap:")[1].split(
            "rule run_truvari_anno_lcr:"
        )[0]
        self.assertIn("../envs/truvari_remap.yml", remap_rule)
        self.assertIn("truvari anno remap", remap_rule)


if __name__ == "__main__":
    unittest.main()
