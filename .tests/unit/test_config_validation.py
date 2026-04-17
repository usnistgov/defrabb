"""Tests for cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

# Make scripts/ importable (mirrors the shim used in rules/common.smk).
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from validate_configs import validate_cross_references  # noqa: E402


def _minimal_config():
    return {
        "assemblies": {"HG002_v1.0": {}},
        "references": {"GRCh38": {}},
        "comparisons": {"GRCh38": {"HG002_v4.2.1": {}}},
        "exclusion_set": {"smvar": []},
    }


def _analyses_df(rows):
    """Return an analyses DataFrame indexed by eval_id."""
    df = pd.DataFrame(rows)
    return df.set_index("eval_id")


def test_happy_path_returns_none():
    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval1",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
            }
        ]
    )
    assert validate_cross_references(config, analyses) is None


def test_missing_asm_id_raises_with_eval_id_in_message():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_typo",
                "asm_id": "HG002_TYPO",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "HG002_TYPO" in msg
    assert "eval_typo" in msg
    assert "assemblies" in msg


def test_missing_ref_raises_with_eval_id_in_message():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_badref",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38_TYPO",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "GRCh38_TYPO" in msg
    assert "eval_badref" in msg
    assert "references" in msg


def test_missing_exclusion_set_raises():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_badexcl",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "not_a_real_set",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "not_a_real_set" in msg
    assert "eval_badexcl" in msg
    assert "exclusion" in msg


def test_exclusion_set_none_is_valid():
    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_no_excl",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "none",
            }
        ]
    )
    assert validate_cross_references(config, analyses) is None


def test_missing_comp_id_for_existing_ref_raises():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_badcomp",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "v4.2.1_TYPO",
                "exclusion_set": "smvar",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "v4.2.1_TYPO" in msg
    assert "GRCh38" in msg
    assert "eval_badcomp" in msg
    assert "comparisons" in msg


def test_missing_ref_does_not_double_report_comp_id():
    """If ref itself is missing, we don't also flag the comp_id under it."""
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_double",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38_NOTREAL",
                "eval_comp_id": "anything",
                "exclusion_set": "smvar",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "GRCh38_NOTREAL" in msg
    # comp_id error should NOT appear since the ref is missing
    assert "Missing comparisons" not in msg
