"""Pipeline-generated PAV vs dipcall exclusions (v0.023).

pav-discrep-{smvar,stvar} and pav-inv replace the HG002Q100-* beds curated
against HG002 Q100 v1.0/v1.1 so the exclusions track the assembly benchmarked.
Structural checks only (no Snakemake invocation).

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import re
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
PAV_DISCREP = REPO_ROOT / "rules" / "exclusions_pav_discrep.smk"
EXCLUSIONS_APPLY = REPO_ROOT / "rules" / "exclusions_apply.smk"
SNAKEFILE = REPO_ROOT / "Snakefile"
RESOURCES = yaml.safe_load((REPO_ROOT / "config" / "resources.yml").read_text())

GENERATED = ["pav-discrep-smvar", "pav-discrep-stvar", "pav-inv"]


def rule_block(text: str, name: str) -> str:
    block = text.split(f"rule {name}:")[1]
    return re.split(r"^rule ", block, flags=re.MULTILINE)[0]


def test_module_included():
    assert 'include: "rules/exclusions_pav_discrep.smk"' in SNAKEFILE.read_text()


def test_generated_exclusions_are_ref_agnostic():
    for excl in GENERATED:
        assert excl in RESOURCES["exclusion_ref_agnostic"]
        assert excl not in RESOURCES["exclusion_asm_agnostic"]


def test_params_defined():
    params = RESOURCES["_exclusion_params"]
    for bench_type in ["smvar", "stvar"]:
        for key in ["slop", "merge_dist"]:
            assert f"pav_discrep_{bench_type}_{key}" in params
    assert params["pav_discrep_smvar_max_indel"] == 50


def test_comparisons_shared_across_benchmarks():
    """hap.py / truvari PAV vs dipcall run once per ref + assembly, not per bench."""
    text = PAV_DISCREP.read_text()
    for name in ["pav_discrep_prep", "pav_discrep_happy", "pav_discrep_truvari"]:
        assert "{bench_id}" not in rule_block(text, name).split("log:")[0]


def test_truth_is_dipcall_query_is_pav():
    text = PAV_DISCREP.read_text()
    happy = rule_block(text, "pav_discrep_happy")
    assert re.search(r"truth=.*dipcall", happy) and re.search(r"query=.*pav", happy)
    truvari = rule_block(text, "pav_discrep_truvari")
    assert re.search(r"base=.*dipcall", truvari) and re.search(r"comp=.*pav", truvari)


def test_pav_inv_uses_pav_run_calls():
    """pav-inv must use the PAV run's calls so it works for dipcall benchmarks."""
    block = rule_block(EXCLUSIONS_APPLY.read_text(), "exclude_pav_inversions")
    assert "get_pav_discrep_inputs" in block
    assert "{input.pav_vcf}" in block
    # segdups overlapping inversions (expansion), not a union with all segdups
    assert "multiinter" not in block
    assert re.search(r"intersect -wa -a \{input\.segdups\}", block)


def test_v023_exclusion_sets():
    sets = RESOURCES["exclusion_set"]
    smvar = set(sets["HG002Q100smvarv0.023"])
    stvar = set(sets["HG002Q100stvarv0.023"])
    stvar_pav = set(sets["HG002Q100stvarv0.023pav"])
    for s in [smvar, stvar, stvar_pav]:
        assert set(GENERATED) <= s
        assert not any(
            e.startswith("HG002Q100-") and e != "HG002Q100-mosaic" for e in s
        )
    assert "HG002Q100-mosaic" in smvar
    assert "HG002Q100-mosaic" not in stvar | stvar_pav
    assert {"dipcall-bugs-T2TACE", "consecutive-svs"} <= smvar & stvar
    assert not {"dipcall-bugs-T2TACE", "consecutive-svs"} & stvar_pav
