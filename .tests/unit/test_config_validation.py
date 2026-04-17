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
