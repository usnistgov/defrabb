"""Structural invariants for the SV-aware self-discrepancy path (#193).

For stvar benchmarks, the final ``_self-discrep.bed`` must be produced via the
truvari bench path, not the hap.py/vcfeval path.  For smvar, the existing
hap.py path is unchanged.

Pure-text/structural checks (no Snakemake invocation), matching the style of
the other rule-parsing unit tests in this suite.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SELF_DISCREP = REPO_ROOT / "rules" / "exclusions_self_discrep.smk"


def _text():
    return SELF_DISCREP.read_text()


def extract_rule_block(text: str, rule_name: str) -> str:
    """Return the source of a single ``rule <name>:`` block."""
    lines = text.splitlines()
    start = next(
        (
            i
            for i, ln in enumerate(lines)
            if re.match(rf"^rule {re.escape(rule_name)}:\s*$", ln)
        ),
        None,
    )
    assert start is not None, f"rule {rule_name} not found"
    end = len(lines)
    for j in range(start + 1, len(lines)):
        if re.match(r"^(rule|checkpoint) \w+:\s*$", lines[j]):
            end = j
            break
    return "\n".join(lines[start:end])


# --------------------------------------------------------------------------- #
# New truvari rules exist and are well-formed                                  #
# --------------------------------------------------------------------------- #


def test_self_discrep_truvari_rule_present():
    assert "rule self_discrep_truvari:" in _text(), (
        "self_discrep_truvari rule must be present in exclusions_self_discrep.smk"
    )


def test_self_discrep_truvari_extract_fpfns_rule_present():
    assert "rule self_discrep_truvari_extract_fpfns:" in _text(), (
        "self_discrep_truvari_extract_fpfns rule must be present"
    )


def test_self_discrep_truvari_has_log_directive():
    block = extract_rule_block(_text(), "self_discrep_truvari")
    assert re.search(r"^\s*log:\s*$", block, re.MULTILINE), (
        "self_discrep_truvari must declare a log: directive for debuggability"
    )


def test_self_discrep_truvari_extract_fpfns_has_log_directive():
    block = extract_rule_block(_text(), "self_discrep_truvari_extract_fpfns")
    assert re.search(r"^\s*log:\s*$", block, re.MULTILINE), (
        "self_discrep_truvari_extract_fpfns must declare a log: directive"
    )


# --------------------------------------------------------------------------- #
# stvar path: truvari bench is used, not hap.py                               #
# --------------------------------------------------------------------------- #


def test_self_discrep_truvari_constrained_to_stvar():
    """The truvari bench rule must only fire for stvar benchmarks."""
    block = extract_rule_block(_text(), "self_discrep_truvari")
    assert "stvar" in block, (
        "self_discrep_truvari must be constrained to bench_type=stvar "
        "to prevent it from matching smvar benchmarks"
    )
    assert "wildcard_constraints" in block, (
        "self_discrep_truvari must have a wildcard_constraints: directive "
        "scoping it to stvar"
    )


def test_self_discrep_truvari_extract_constrained_to_stvar():
    """The FP/FN extraction rule must also be stvar-only."""
    block = extract_rule_block(_text(), "self_discrep_truvari_extract_fpfns")
    assert "stvar" in block and "wildcard_constraints" in block, (
        "self_discrep_truvari_extract_fpfns must have wildcard_constraints "
        "scoping it to bench_type=stvar"
    )


def test_self_discrep_truvari_uses_vcf_as_both_base_and_comp():
    """Self-comparison: the same VCF is passed as both -b and -c."""
    block = extract_rule_block(_text(), "self_discrep_truvari")
    shell_section = block.split("shell:")[1] if "shell:" in block else ""
    # Both -b and -c must reference {input.vcf}
    assert shell_section.count("input.vcf") >= 2, (
        "self_discrep_truvari must pass input.vcf as both -b (base) and "
        "-c (comp) for self-comparison"
    )


def test_self_discrep_truvari_uses_includebed():
    """The benchmark BED must be passed as --includebed to restrict the comparison."""
    block = extract_rule_block(_text(), "self_discrep_truvari")
    assert "--includebed" in block, (
        "self_discrep_truvari must pass the benchmark BED via --includebed"
    )


def test_self_discrep_truvari_output_dir_contains_truvari():
    """The truvari output directory name must include '_truvari' to distinguish
    it from the hap.py intermediate files in the same parent directory."""
    block = extract_rule_block(_text(), "self_discrep_truvari")
    output_section = block.split("output:")[1].split("log:")[0]
    assert "_truvari/" in output_section, (
        "self_discrep_truvari outputs must live under a *_truvari/ subdirectory"
    )


def test_self_discrep_truvari_extract_fpfns_consumes_fn_and_fp():
    """The FP/FN extractor must read both fn.vcf.gz and fp.vcf.gz."""
    block = extract_rule_block(_text(), "self_discrep_truvari_extract_fpfns")
    assert "fn.vcf.gz" in block and "fp.vcf.gz" in block, (
        "self_discrep_truvari_extract_fpfns must consume both fn.vcf.gz "
        "and fp.vcf.gz from the truvari output directory"
    )


# --------------------------------------------------------------------------- #
# stvar final BED is routed via the truvari intermediate                       #
# --------------------------------------------------------------------------- #


def test_router_helper_present():
    """A helper function must route the FP/FN BED by bench_type."""
    assert "def get_self_discrep_fpfns_bed" in _text(), (
        "get_self_discrep_fpfns_bed helper must be defined in "
        "exclusions_self_discrep.smk to route smvar vs stvar paths"
    )


def test_router_helper_routes_stvar_to_truvari():
    """The router must return the *_truvari.fpfns.bed path for stvar."""
    text = _text()
    func_start = text.index("def get_self_discrep_fpfns_bed")
    # Grab up to the next blank line after a return statement
    func_block = text[func_start : func_start + 600]
    assert "_truvari.fpfns.bed" in func_block, (
        "get_self_discrep_fpfns_bed must return a *_truvari.fpfns.bed "
        "path when bench_type == 'stvar'"
    )
    assert "stvar" in func_block, (
        "get_self_discrep_fpfns_bed must branch on bench_type == 'stvar'"
    )


def test_router_helper_routes_smvar_to_hap_py():
    """The router must return the .fpfns.bed (no _truvari) path for smvar."""
    text = _text()
    func_start = text.index("def get_self_discrep_fpfns_bed")
    func_block = text[func_start : func_start + 600]
    # The smvar return must NOT include _truvari
    returns = [ln.strip() for ln in func_block.splitlines() if "return" in ln]
    assert any("_truvari" not in r for r in returns), (
        "get_self_discrep_fpfns_bed must return a path without '_truvari' "
        "for the smvar (hap.py) branch"
    )


def test_intersect_slop_uses_router_for_bed_input():
    """self_discrep_intersect_slop must use get_self_discrep_fpfns_bed as its
    bed input so that smvar and stvar are routed to the correct intermediate."""
    block = extract_rule_block(_text(), "self_discrep_intersect_slop")
    assert "get_self_discrep_fpfns_bed" in block, (
        "self_discrep_intersect_slop.input.bed must use get_self_discrep_fpfns_bed "
        "to route smvar vs stvar intermediate BEDs"
    )


# --------------------------------------------------------------------------- #
# smvar path is unchanged                                                      #
# --------------------------------------------------------------------------- #


def test_smvar_hap_py_rules_still_present():
    """The existing hap.py rules must not have been removed (#192 fix stays)."""
    text = _text()
    assert "rule self_discrep_filter_symbolic:" in text
    assert "rule self_discrep_happy:" in text
    assert "rule self_discrep_extract_fpfns:" in text


def test_happy_consumes_no_symbolic_vcf():
    """hap.py must still consume the .no-symbolic.vcf.gz (the #192 fix)."""
    block = extract_rule_block(_text(), "self_discrep_happy")
    assert ".no-symbolic.vcf.gz" in block, (
        "self_discrep_happy must still consume the no-symbolic VCF (#192 fix)"
    )
