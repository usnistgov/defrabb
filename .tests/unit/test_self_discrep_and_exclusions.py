"""Structural tests for #188 (Q100 SV exclusions) and #192 (self-discrep filter).

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import unittest
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
RESOURCES = REPO_ROOT / "config" / "resources.yml"
SELF_DISCREP = REPO_ROOT / "rules" / "exclusions_self_discrep.smk"


def _exclusion_sets():
    return yaml.safe_load(RESOURCES.read_text())["exclusion_set"]


class Q100StvarExclusionTests(unittest.TestCase):
    """#188: the new stvar set drops mosaic; historical sets are untouched."""

    def test_v022_set_exists_without_mosaic(self):
        sets = _exclusion_sets()
        self.assertIn("HG002Q100stvarv0.022", sets)
        s = sets["HG002Q100stvarv0.022"]
        self.assertNotIn("HG002Q100-mosaic", s)
        # still carries the other expected stvar exclusions
        for region in ("segdups", "HG002Q100-pav-inversions", "HG002Q100-errors"):
            self.assertIn(region, s)

    def test_v022_is_v017_minus_mosaic(self):
        sets = _exclusion_sets()
        expected = [r for r in sets["HG002Q100stvarv0.017"] if r != "HG002Q100-mosaic"]
        self.assertEqual(sets["HG002Q100stvarv0.022"], expected)

    def test_historical_v017_unchanged(self):
        # We add a new version rather than rewriting pinned releases.
        self.assertIn("HG002Q100-mosaic", _exclusion_sets()["HG002Q100stvarv0.017"])


class SelfDiscrepSymbolicFilterTests(unittest.TestCase):
    """#192: symbolic/breakend alleles are filtered before hap.py."""

    def setUp(self):
        self.text = SELF_DISCREP.read_text()

    def test_filter_rule_present_and_excludes_symbolic(self):
        self.assertIn("rule self_discrep_filter_symbolic:", self.text)
        # excludes symbolic (<...>) and breakend ([ ]) ALT alleles
        self.assertIn('ALT~"<"', self.text)

    def test_happy_consumes_filtered_vcf(self):
        happy = self.text.split("rule self_discrep_happy:")[1].split("rule ")[0]
        self.assertIn(".no-symbolic.vcf.gz", happy)
        # the raw standardized vcf is no longer fed directly to hap.py
        self.assertNotIn("vcf=get_standardized_vcf", happy)

    def test_filter_does_not_emit_tbi(self):
        # index is left to the generic tabix rule to avoid a rule clash
        flt = self.text.split("rule self_discrep_filter_symbolic:")[1].split(
            "rule self_discrep_happy:"
        )[0]
        self.assertNotIn(".no-symbolic.vcf.gz.tbi", flt)


if __name__ == "__main__":
    unittest.main()
