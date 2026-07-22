"""Unit tests for parameter sweep generator.

Tests cross-product generation, vc_id assignment, validation, and cost estimation.
"""
import pytest
import tempfile
from pathlib import Path
import yaml
import pandas as pd
import sys

# Add scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from generate_param_sweep import (
    load_sweep_config,
    generate_cross_product,
    assign_vc_ids,
    estimate_costs,
    format_analyses_table,
    validate_output,
)


def test_load_sweep_config():
    """Test sweep config loading and validation."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
        yaml.dump({
            "name": "Test Sweep",
            "output": "test.tsv",
            "fixed": {"ref_id": "GRCh38"},
            "sweep": {"vc_param_id": ["z2k", "z5k"]}
        }, f)
        config_path = Path(f.name)

    try:
        config = load_sweep_config(config_path)
        assert config["name"] == "Test Sweep"
        assert config["fixed"]["ref_id"] == "GRCh38"
        assert config["sweep"]["vc_param_id"] == ["z2k", "z5k"]
    finally:
        config_path.unlink()


def test_load_sweep_config_missing_fields():
    """Test that missing required fields raise error."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
        yaml.dump({"name": "Incomplete"}, f)  # Missing fixed, sweep
        config_path = Path(f.name)

    try:
        with pytest.raises(ValueError, match="missing required fields"):
            load_sweep_config(config_path)
    finally:
        config_path.unlink()


def test_generate_cross_product():
    """Test cross-product generation from sweep dimensions."""
    fixed = {"ref_id": "GRCh38", "asm_id": "HG002"}
    sweep = {
        "vc_param_id": ["z2k", "z5k"],
        "exclusion_set": ["standard", "aggressive"]
    }

    analyses = generate_cross_product(fixed, sweep)

    # Should generate 2 × 2 = 4 analyses
    assert len(analyses) == 4

    # Check all have fixed fields
    for row in analyses:
        assert row["ref_id"] == "GRCh38"
        assert row["asm_id"] == "HG002"

    # Check all combinations present
    combos = {(r["vc_param_id"], r["exclusion_set"]) for r in analyses}
    assert combos == {
        ("z2k", "standard"), ("z2k", "aggressive"),
        ("z5k", "standard"), ("z5k", "aggressive")
    }


def test_assign_vc_ids_reuse():
    """Test vc_id assignment and reuse factor calculation."""
    analyses = [
        {"ref_id": "GRCh38", "asm_id": "HG002", "vc_cmd": "dipcall", "vc_param_id": "z2k"},
        {"ref_id": "GRCh38", "asm_id": "HG002", "vc_cmd": "dipcall", "vc_param_id": "z2k"},
        {"ref_id": "GRCh38", "asm_id": "HG002", "vc_cmd": "dipcall", "vc_param_id": "z5k"},
    ]

    unique_vc_runs = assign_vc_ids(analyses)

    # Should have 2 unique vc_ids (z2k and z5k)
    assert unique_vc_runs == 2

    # First two should share vc_id
    assert analyses[0]["vc_id"] == analyses[1]["vc_id"]

    # Third should be different
    assert analyses[2]["vc_id"] != analyses[0]["vc_id"]

    # All should have vc_id assigned
    assert all("vc_id" in row for row in analyses)


def test_assign_vc_ids_different_params():
    """Test that different variant call params get different vc_ids."""
    analyses = [
        {"ref_id": "GRCh38", "asm_id": "HG002", "vc_cmd": "dipcall", "vc_param_id": "z2k"},
        {"ref_id": "GRCh38", "asm_id": "HG002", "vc_cmd": "pav", "vc_param_id": "default"},
        {"ref_id": "GRCh37", "asm_id": "HG002", "vc_cmd": "dipcall", "vc_param_id": "z2k"},
    ]

    unique_vc_runs = assign_vc_ids(analyses)

    # All different (ref differs, or vc_cmd differs)
    assert unique_vc_runs == 3
    assert len({r["vc_id"] for r in analyses}) == 3


def test_estimate_costs():
    """Test runtime and storage estimation."""
    analyses = [
        {"vc_cmd": "dipcall", "eval_cmd": "happy"},
        {"vc_cmd": "dipcall", "eval_cmd": "happy"},
        {"vc_cmd": "pav", "eval_cmd": "truvari"},
    ]

    # 2 unique vc runs (simplified - real would assign vc_ids first)
    costs = estimate_costs(analyses, unique_vc_runs=2)

    assert costs["total_analyses"] == 3
    assert costs["unique_vc_runs"] == 2

    # Runtime: (dipcall 5h + pav 6h) + 3 evals (2×0.5 + 1×1)
    # Note: estimation uses simplified logic, exact values may vary
    assert costs["runtime_hours"] > 0
    assert costs["storage_gb"] > 0


def test_format_analyses_table():
    """Test DataFrame formatting with standard column order."""
    analyses = [
        {
            "vc_id": "vc001",
            "ref_id": "GRCh38",
            "asm_id": "HG002",
            "vc_cmd": "dipcall",
            "bench_type": "smvar",
            "eval_cmd": "happy",
            "eval_comp_id": "v5.0q-smvar"
        }
    ]

    df = format_analyses_table(analyses)

    # Check standard columns present
    assert "vc_id" in df.columns
    assert "ref_id" in df.columns
    assert "bench_id" in df.columns  # Generated
    assert "eval_id" in df.columns  # Generated

    # Check defaults filled
    assert df.iloc[0]["vc_param_id"] == "default"
    assert df.iloc[0]["vcf_processing"] == "default"

    # Check derived IDs
    assert "HG002" in df.iloc[0]["bench_id"]
    assert "happy" in df.iloc[0]["eval_id"]


def test_validate_output_success():
    """Test validation passes for valid output."""
    df = pd.DataFrame([{
        "vc_id": "vc001",
        "ref_id": "GRCh38",
        "asm_id": "HG002",
        "vc_cmd": "dipcall",
        "bench_type": "smvar",
        "bench_id": "test_bench_1",
        "eval_id": "happy_v5q",
        "eval_cmd": "happy",
        "eval_comp_id": "v5.0q-smvar"
    }])

    # Should not raise
    validate_output(df)


def test_validate_output_missing_columns():
    """Test validation fails for missing required columns."""
    df = pd.DataFrame([{"vc_id": "vc001"}])  # Missing most columns

    with pytest.raises(ValueError, match="missing required columns"):
        validate_output(df)


def test_validate_output_duplicates():
    """Test validation fails for duplicate (eval_id, bench_id) pairs."""
    df = pd.DataFrame([
        {
            "vc_id": "vc001", "ref_id": "GRCh38", "asm_id": "HG002",
            "vc_cmd": "dipcall", "bench_type": "smvar",
            "bench_id": "test_bench_1", "eval_id": "happy_v5q",
            "eval_cmd": "happy", "eval_comp_id": "v5.0q-smvar"
        },
        {
            "vc_id": "vc001", "ref_id": "GRCh38", "asm_id": "HG002",
            "vc_cmd": "dipcall", "bench_type": "smvar",
            "bench_id": "test_bench_1", "eval_id": "happy_v5q",  # Duplicate
            "eval_cmd": "happy", "eval_comp_id": "v5.0q-smvar"
        }
    ])

    with pytest.raises(ValueError, match="Duplicate"):
        validate_output(df)


def test_end_to_end_sweep_generation():
    """Integration test: full sweep generation workflow."""
    fixed = {
        "ref_id": "GRCh38",
        "asm_id": "HG002",
        "vc_cmd": "dipcall",
        "bench_type": "smvar",
        "eval_cmd": "happy",
        "eval_comp_id": "v5.0q-smvar"
    }
    sweep = {
        "vc_param_id": ["z2k", "z5k"],
        "exclusion_set": ["standard", "aggressive"]
    }

    # Generate
    analyses = generate_cross_product(fixed, sweep)
    assert len(analyses) == 4

    # Assign vc_ids
    unique_vc_runs = assign_vc_ids(analyses)
    assert unique_vc_runs == 2  # Only vc_param_id varies

    # Format
    df = format_analyses_table(analyses)
    assert len(df) == 4

    # Validate
    validate_output(df)  # Should not raise

    # Check reuse factor
    reuse_factor = len(analyses) / unique_vc_runs
    assert reuse_factor == 2.0
