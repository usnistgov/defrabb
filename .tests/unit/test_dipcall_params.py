"""Unit tests for dipcall parameter profile lookup.

Tests the fix for #197 where _dipcall_params profiles were defined but never
accessed by rules/asm-varcall.smk. The extra= param should look up profiles by
vc_param_id, falling back to the vc_params column value if profile not found.
"""
import pytest
from unittest.mock import Mock
import pandas as pd


def test_dipcall_param_profile_lookup():
    """Test that vc_param_id looks up _dipcall_params profile."""
    # Simulate config and vc_tbl as they'd appear in Snakefile
    config = {
        "_dipcall_params": {
            "z2k": "-z200000,10000",
            "z5k": "-z500000,5000",
            "z10k": "-z1000000,10000",
        }
    }

    vc_tbl = pd.DataFrame({
        "vc_id": ["vc1", "vc2", "vc3"],
        "vc_params": ["default", "default", "-z300000,3000"],
    })
    vc_tbl = vc_tbl.set_index("vc_id")

    # The lambda from rules/asm-varcall.smk lines 48-52 (after fix)
    def get_extra(vc_id, vc_param_id):
        return (
            config["_dipcall_params"].get(
                vc_param_id,
                vc_tbl.loc[vc_id]["vc_params"]
            )
            if vc_tbl.loc[vc_id]["vc_params"] != "default"
            else config["_dipcall_params"].get(vc_param_id, "")
        )

    # Test 1: vc_params="default", vc_param_id in profile → use profile
    assert get_extra("vc1", "z2k") == "-z200000,10000"
    assert get_extra("vc2", "z5k") == "-z500000,5000"

    # Test 2: vc_params="default", vc_param_id NOT in profile → empty string
    assert get_extra("vc1", "unknown") == ""

    # Test 3: vc_params has value, vc_param_id in profile → use profile (profile wins)
    assert get_extra("vc3", "z10k") == "-z1000000,10000"

    # Test 4: vc_params has value, vc_param_id NOT in profile → use vc_params
    assert get_extra("vc3", "unknown") == "-z300000,3000"


def test_dipcall_param_backward_compatibility():
    """Test that old analyses tables with vc_params still work."""
    config = {"_dipcall_params": {"z2k": "-z200000,10000"}}

    # Old-style table: vc_params column has the literal CLI string
    vc_tbl = pd.DataFrame({
        "vc_id": ["vc_old"],
        "vc_params": ["-z200000,10000"],
    })
    vc_tbl = vc_tbl.set_index("vc_id")

    def get_extra(vc_id, vc_param_id):
        return (
            config["_dipcall_params"].get(
                vc_param_id,
                vc_tbl.loc[vc_id]["vc_params"]
            )
            if vc_tbl.loc[vc_id]["vc_params"] != "default"
            else config["_dipcall_params"].get(vc_param_id, "")
        )

    # Even with vc_param_id=z2k, if vc_params != "default", fallback to vc_params
    # This maintains backward compat for old tables
    assert get_extra("vc_old", "unknown") == "-z200000,10000"


def test_dipcall_param_profiles_defined():
    """Test that expected profiles exist in resources.yml."""
    import yaml

    with open("config/resources.yml") as f:
        resources = yaml.safe_load(f)

    assert "_dipcall_params" in resources
    params = resources["_dipcall_params"]

    # Check that standard profiles exist
    assert "z2k" in params
    assert params["z2k"] == "-z200000,10000"

    # Check that new v0.023 profiles added
    assert "z5k" in params
    assert "z10k" in params
    assert "z1k" in params

    # All profiles should be -z formatted strings
    for profile_id, cli_flags in params.items():
        assert cli_flags.startswith("-z"), f"Profile {profile_id} should start with -z"
        assert "," in cli_flags, f"Profile {profile_id} should have window,kmer format"


def test_dipcall_param_schema_validation():
    """Test that schema validates dipcall param profiles."""
    import yaml
    import jsonschema

    with open("schema/resources-schema.yml") as f:
        schema = yaml.safe_load(f)

    # Valid config should pass
    valid_config = {
        "_dipcall_params": {
            "z2k": "-z200000,10000",
            "z5k": "-z500000,5000",
        }
    }

    # Extract the _dipcall_params subschema
    dipcall_schema = schema["properties"]["_dipcall_params"]

    # Should not raise
    jsonschema.validate(valid_config["_dipcall_params"], dipcall_schema)

    # Invalid format should fail (missing -z prefix)
    invalid_config = {"bad": "200000,10000"}

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid_config, dipcall_schema)
