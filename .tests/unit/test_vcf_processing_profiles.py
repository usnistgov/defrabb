"""Tests for VCF processing profile validation and resolution.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from validate_configs import validate_cross_references  # noqa: E402


def test_validate_vcf_profiles_accepts_valid_profile():
    """Cross-reference validation passes when profile exists in config."""
    config = {
        "assemblies": {"asm1": {}},
        "references": {"ref1": {}},
        "exclusion_set": {"ex1": []},
        "comparisons": {"ref1": {"comp1": {}}},
        "vcf_processing_profiles": {
            "trf": ["trfanno"],
            "trf_sv_lcr": ["trfanno", "svinfo", "lcr"],
        },
    }
    analyses = pd.DataFrame(
        [
            {
                "eval_id": "eval1",
                "asm_id": "asm1",
                "ref": "ref1",
                "exclusion_set": "ex1",
                "eval_comp_id": "comp1",
                "bench_vcf_processing": "trf",
            }
        ]
    ).set_index("eval_id")
    # Should not raise
    validate_cross_references(config, analyses)


def test_validate_vcf_profiles_accepts_none_sentinel():
    """Cross-reference validation allows 'none' without requiring a profile."""
    config = {
        "assemblies": {"asm1": {}},
        "references": {"ref1": {}},
        "exclusion_set": {"ex1": []},
        "comparisons": {"ref1": {"comp1": {}}},
        "vcf_processing_profiles": {},
    }
    analyses = pd.DataFrame(
        [
            {
                "eval_id": "eval1",
                "asm_id": "asm1",
                "ref": "ref1",
                "exclusion_set": "ex1",
                "eval_comp_id": "comp1",
                "bench_vcf_processing": "none",
            }
        ]
    ).set_index("eval_id")
    # Should not raise
    validate_cross_references(config, analyses)


def test_validate_vcf_profiles_rejects_missing_profile():
    """Cross-reference validation raises WorkflowError for missing profile."""
    from snakemake.exceptions import WorkflowError

    config = {
        "assemblies": {"asm1": {}},
        "references": {"ref1": {}},
        "exclusion_set": {"ex1": []},
        "comparisons": {"ref1": {"comp1": {}}},
        "vcf_processing_profiles": {"trf": ["trfanno"]},
    }
    analyses = pd.DataFrame(
        [
            {
                "eval_id": "eval1",
                "asm_id": "asm1",
                "ref": "ref1",
                "exclusion_set": "ex1",
                "eval_comp_id": "comp1",
                "bench_vcf_processing": "nonexistent_profile",
            }
        ]
    ).set_index("eval_id")
    with pytest.raises(WorkflowError) as excinfo:
        validate_cross_references(config, analyses)
    assert "Missing VCF processing profiles" in str(excinfo.value)
    assert "nonexistent_profile" in str(excinfo.value)
    assert "eval1" in str(excinfo.value)


def test_validate_vcf_profiles_groups_multiple_missing():
    """Validation groups multiple missing profiles in one error."""
    from snakemake.exceptions import WorkflowError

    config = {
        "assemblies": {"asm1": {}},
        "references": {"ref1": {}},
        "exclusion_set": {"ex1": []},
        "comparisons": {"ref1": {"comp1": {}}},
        "vcf_processing_profiles": {},
    }
    analyses = pd.DataFrame(
        [
            {
                "eval_id": "eval1",
                "asm_id": "asm1",
                "ref": "ref1",
                "exclusion_set": "ex1",
                "eval_comp_id": "comp1",
                "bench_vcf_processing": "profile_a",
            },
            {
                "eval_id": "eval2",
                "asm_id": "asm1",
                "ref": "ref1",
                "exclusion_set": "ex1",
                "eval_comp_id": "comp1",
                "bench_vcf_processing": "profile_b",
            },
        ]
    ).set_index("eval_id")
    with pytest.raises(WorkflowError) as excinfo:
        validate_cross_references(config, analyses)
    msg = str(excinfo.value)
    assert "profile_a" in msg
    assert "profile_b" in msg
    assert "eval1" in msg
    assert "eval2" in msg
