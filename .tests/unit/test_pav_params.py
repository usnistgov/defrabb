"""Unit tests for PAV parameter profile lookup.

Tests that _pav_config profiles are correctly defined and accessible.
"""
import pytest
import yaml


def test_pav_profiles_defined():
    """Test that expected PAV profiles exist in resources.yml."""
    with open("config/resources.yml") as f:
        resources = yaml.safe_load(f)

    assert "_pav_config" in resources
    pav_config = resources["_pav_config"]

    # Check that standard profiles exist
    assert "giab" in pav_config, "Default giab profile must exist"
    assert "strict" in pav_config, "Strict profile for v0.023 optimization"
    assert "lenient" in pav_config, "Lenient profile for v0.023 optimization"


def test_pav_profile_structure():
    """Test that PAV profiles have required merge parameters."""
    with open("config/resources.yml") as f:
        resources = yaml.safe_load(f)

    for profile_name, profile in resources["_pav_config"].items():
        assert "merge_ins" in profile, f"{profile_name} missing merge_ins"
        assert "merge_del" in profile, f"{profile_name} missing merge_del"
        assert "merge_inv" in profile, f"{profile_name} missing merge_inv"

        # All should be nr:: format (PAV merge algorithm syntax)
        assert profile["merge_ins"].startswith("nr::"), \
            f"{profile_name} merge_ins should use nr:: format"
        assert profile["merge_del"].startswith("nr::"), \
            f"{profile_name} merge_del should use nr:: format"
        assert profile["merge_inv"].startswith("nr::"), \
            f"{profile_name} merge_inv should use nr:: format"


def test_pav_profile_schema_validation():
    """Test that schema validates PAV profiles."""
    import jsonschema

    with open("schema/resources-schema.yml") as f:
        schema = yaml.safe_load(f)

    with open("config/resources.yml") as f:
        resources = yaml.safe_load(f)

    pav_schema = schema["properties"]["_pav_config"]

    # Valid config should pass
    jsonschema.validate(resources["_pav_config"], pav_schema)

    # Invalid format should fail (missing required field)
    invalid_profile = {
        "bad": {
            "merge_ins": "nr::exact",
            # missing merge_del and merge_inv
        }
    }

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(invalid_profile, pav_schema)


def test_pav_giab_profile_unchanged():
    """Test that giab profile maintains backward compatibility."""
    with open("config/resources.yml") as f:
        resources = yaml.safe_load(f)

    giab = resources["_pav_config"]["giab"]

    # These are the historical values - must not change
    assert giab["merge_ins"] == "nr::exact"
    assert giab["merge_del"] == "nr::exact"
    assert giab["merge_inv"] == "nr::exact:ro(0.5):szro(0.5,200):match"
