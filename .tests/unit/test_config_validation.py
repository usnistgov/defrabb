"""Tests for cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.

Scope (TODO #19): negative/error case tests for config validation including
schema violations, missing columns, invalid enum values, and malformed inputs.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml
from snakemake.utils import validate as snakemake_validate

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
        "vcf_processing_profiles": {},
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
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
                "bench_vcf_processing": "none",
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "GRCh38_NOTREAL" in msg
    # comp_id error should NOT appear since the ref is missing
    assert "Missing comparisons" not in msg


def test_multiple_errors_grouped_in_single_raise():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_a",
                "asm_id": "HG002_TYPO",  # missing asm
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
                "bench_vcf_processing": "none",
            },
            {
                "eval_id": "eval_b",
                "asm_id": "HG002_TYPO",  # same missing asm — should group
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
                "bench_vcf_processing": "none",
            },
            {
                "eval_id": "eval_c",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38_TYPO",  # missing ref (different section)
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
                "bench_vcf_processing": "none",
            },
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)

    # Both sections present
    assert "Missing assemblies" in msg
    assert "Missing references" in msg

    # Same missing asm grouped on one line, both eval_ids listed
    assert "HG002_TYPO" in msg
    assert "eval_a" in msg
    assert "eval_b" in msg
    # Only one occurrence of the missing asm token in the asm section
    assert msg.count("HG002_TYPO") == 1


def test_real_repo_configs_validate():
    """Smoke: the actual repo config + default analyses.tsv must pass validation."""
    import yaml

    resources_path = REPO_ROOT / "config" / "resources.yml"
    analyses_path = REPO_ROOT / "config" / "analyses.tsv"

    with resources_path.open() as fh:
        config = yaml.safe_load(fh)

    analyses_df = pd.read_table(analyses_path).set_index("eval_id")

    assert validate_cross_references(config, analyses_df) is None


def test_missing_vcf_processing_profile_raises():
    """VCF processing profile referenced but not defined in resources.yml."""
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    config["vcf_processing_profiles"] = {"trf": ["trfanno"]}

    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_badprof",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "smvar",
                "bench_vcf_processing": "xy_trf",  # Not defined
            }
        ]
    )
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "xy_trf" in msg
    assert "eval_badprof" in msg
    assert "vcf_processing_profiles" in msg


def test_vcf_processing_none_is_valid():
    """bench_vcf_processing='none' should not require profile definition."""
    config = _minimal_config()
    config["vcf_processing_profiles"] = {}

    analyses = _analyses_df(
        [
            {
                "eval_id": "eval_no_prof",
                "asm_id": "HG002_v1.0",
                "ref": "GRCh38",
                "eval_comp_id": "HG002_v4.2.1",
                "exclusion_set": "none",
                "bench_vcf_processing": "none",
            }
        ]
    )
    assert validate_cross_references(config, analyses) is None


class TestAnalysesSchemaValidation:
    """Tests for analyses.tsv schema validation (missing columns, invalid enums)."""

    @pytest.fixture
    def schema_path(self):
        return REPO_ROOT / "schema" / "analyses-schema.yml"

    @pytest.fixture
    def tmp_analyses_file(self, tmp_path):
        """Create a temporary analyses.tsv file."""
        return tmp_path / "analyses.tsv"

    def test_missing_required_column_raises(self, schema_path, tmp_analyses_file):
        """Missing required column should fail schema validation."""
        # Missing 'eval_cmd' column
        incomplete_data = {
            "eval_id": ["eval1"],
            "bench_id": ["bench1"],
            "asm_id": ["HG002"],
            "ref": ["GRCh38"],
            # eval_cmd missing
            "eval_params": ["default"],
            "eval_comp_id": ["v4.2.1"],
            "eval_comp_id_is_truth": [True],
            "eval_truth_regions": [True],
            "eval_target_regions": [True],
            "bench_type": ["smvar"],
            "bench_vcf_processing": ["none"],
            "bench_bed_processing": ["none"],
            "exclusion_set": ["none"],
            "vc_id": ["vc1"],
            "vc_cmd": ["dipcall"],
            "vc_param_id": ["default"],
            "vc_params": ["default"],
        }
        df = pd.DataFrame(incomplete_data)
        df.to_csv(tmp_analyses_file, sep="\t", index=False)

        # Load and validate
        test_df = pd.read_table(tmp_analyses_file)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(test_df, str(schema_path))

        # Schema validation should catch the missing column
        assert (
            "eval_cmd" in str(exc_info.value)
            or "required" in str(exc_info.value).lower()
        )

    def test_invalid_eval_cmd_enum_raises(self, schema_path, tmp_analyses_file):
        """Invalid eval_cmd value should fail schema validation."""
        data = {
            "eval_id": ["eval1"],
            "eval_cmd": ["invalid_tool"],  # Not happy, truvari, or unhappy
            "eval_params": ["default"],
            "eval_comp_id": ["v4.2.1"],
            "eval_comp_id_is_truth": [True],
            "eval_truth_regions": [True],
            "eval_target_regions": [True],
            "bench_id": ["bench1"],
            "bench_type": ["smvar"],
            "bench_vcf_processing": ["none"],
            "bench_bed_processing": ["none"],
            "exclusion_set": ["none"],
            "vc_id": ["vc1"],
            "asm_id": ["HG002"],
            "ref": ["GRCh38"],
            "vc_cmd": ["dipcall"],
            "vc_param_id": ["default"],
            "vc_params": ["default"],
        }
        df = pd.DataFrame(data)
        df.to_csv(tmp_analyses_file, sep="\t", index=False)

        test_df = pd.read_table(tmp_analyses_file)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(test_df, str(schema_path))

        # Should mention pattern violation
        msg = str(exc_info.value).lower()
        assert "pattern" in msg or "invalid_tool" in msg

    def test_invalid_bench_type_enum_raises(self, schema_path, tmp_analyses_file):
        """Invalid bench_type value should fail schema validation."""
        data = {
            "eval_id": ["eval1"],
            "eval_cmd": ["happy"],
            "eval_params": ["default"],
            "eval_comp_id": ["v4.2.1"],
            "eval_comp_id_is_truth": [True],
            "eval_truth_regions": [True],
            "eval_target_regions": [True],
            "bench_id": ["bench1"],
            "bench_type": ["medvar"],  # Not smvar or stvar
            "bench_vcf_processing": ["none"],
            "bench_bed_processing": ["none"],
            "exclusion_set": ["none"],
            "vc_id": ["vc1"],
            "asm_id": ["HG002"],
            "ref": ["GRCh38"],
            "vc_cmd": ["dipcall"],
            "vc_param_id": ["default"],
            "vc_params": ["default"],
        }
        df = pd.DataFrame(data)
        df.to_csv(tmp_analyses_file, sep="\t", index=False)

        test_df = pd.read_table(tmp_analyses_file)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(test_df, str(schema_path))

        msg = str(exc_info.value).lower()
        assert "pattern" in msg or "medvar" in msg

    def test_invalid_vc_cmd_enum_raises(self, schema_path, tmp_analyses_file):
        """Invalid vc_cmd value should fail schema validation."""
        data = {
            "eval_id": ["eval1"],
            "eval_cmd": ["happy"],
            "eval_params": ["default"],
            "eval_comp_id": ["v4.2.1"],
            "eval_comp_id_is_truth": [True],
            "eval_truth_regions": [True],
            "eval_target_regions": [True],
            "bench_id": ["bench1"],
            "bench_type": ["smvar"],
            "bench_vcf_processing": ["none"],
            "bench_bed_processing": ["none"],
            "exclusion_set": ["none"],
            "vc_id": ["vc1"],
            "asm_id": ["HG002"],
            "ref": ["GRCh38"],
            "vc_cmd": ["delly"],  # Not dipcall or pav
            "vc_param_id": ["default"],
            "vc_params": ["default"],
        }
        df = pd.DataFrame(data)
        df.to_csv(tmp_analyses_file, sep="\t", index=False)

        test_df = pd.read_table(tmp_analyses_file)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(test_df, str(schema_path))

        msg = str(exc_info.value).lower()
        assert "pattern" in msg or "delly" in msg

    def test_invalid_ref_pattern_raises(self, schema_path, tmp_analyses_file):
        """ref value not matching the pattern should fail schema validation."""
        data = {
            "eval_id": ["eval1"],
            "eval_cmd": ["happy"],
            "eval_params": ["default"],
            "eval_comp_id": ["v4.2.1"],
            "eval_comp_id_is_truth": [True],
            "eval_truth_regions": [True],
            "eval_target_regions": [True],
            "bench_id": ["bench1"],
            "bench_type": ["smvar"],
            "bench_vcf_processing": ["none"],
            "bench_bed_processing": ["none"],
            "exclusion_set": ["none"],
            "vc_id": ["vc1"],
            "asm_id": ["HG002"],
            "ref": ["hg19"],  # Not matching GRCh3[7,8]|GRCh38_chr21|CHM13v[0-9].[0-9]
            "vc_cmd": ["dipcall"],
            "vc_param_id": ["default"],
            "vc_params": ["default"],
        }
        df = pd.DataFrame(data)
        df.to_csv(tmp_analyses_file, sep="\t", index=False)

        test_df = pd.read_table(tmp_analyses_file)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(test_df, str(schema_path))

        msg = str(exc_info.value).lower()
        assert "pattern" in msg or "hg19" in msg


class TestResourcesSchemaValidation:
    """Tests for resources.yml schema validation."""

    @pytest.fixture
    def schema_path(self):
        return REPO_ROOT / "schema" / "resources-schema.yml"

    @pytest.fixture
    def tmp_resources_file(self, tmp_path):
        """Create a temporary resources.yml file."""
        return tmp_path / "resources.yml"

    def test_missing_required_section_raises(self, schema_path, tmp_resources_file):
        """Missing required top-level section should fail schema validation."""
        incomplete_config = {
            "analyses": "config/analyses.tsv",
            "references": {},
            # assemblies missing
            "comparisons": {},
            "exclusion_set": {},
            "vcf_processing_profiles": {},
        }

        with open(tmp_resources_file, "w") as f:
            yaml.dump(incomplete_config, f)

        with open(tmp_resources_file) as f:
            config = yaml.safe_load(f)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(config, str(schema_path))

        msg = str(exc_info.value).lower()
        assert "assemblies" in msg or "required" in msg

    def test_malformed_exclusion_set_raises(self, schema_path, tmp_resources_file):
        """Exclusion set with wrong structure should fail schema validation."""
        bad_config = {
            "analyses": "config/analyses.tsv",
            "references": {"GRCh38": {}},
            "assemblies": {"HG002": {}},
            "comparisons": {},
            "exclusion_set": {"smvartest": "this_should_be_a_list"},  # Wrong type
            "vcf_processing_profiles": {},
        }

        with open(tmp_resources_file, "w") as f:
            yaml.dump(bad_config, f)

        with open(tmp_resources_file) as f:
            config = yaml.safe_load(f)

        with pytest.raises(Exception) as exc_info:
            snakemake_validate(config, str(schema_path))

        # Should mention type error
        msg = str(exc_info.value).lower()
        assert "type" in msg or "array" in msg or "list" in msg
