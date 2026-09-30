"""Unit tests for PAV3 parameter profile lookup.

Tests that _pav_config profiles are correctly defined and accessible.
"""

import jsonschema
import pytest
import yaml


def load_yaml(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def test_pav_profiles_defined():
    """The default giab profile must exist in resources.yml."""
    resources = load_yaml("config/resources.yml")
    assert "giab" in resources["_pav_config"], "Default giab profile must exist"


def test_pav_profiles_are_pav3_param_dicts():
    """Profiles are dicts of PAV3 params; PAV2 merge params are not valid."""
    resources = load_yaml("config/resources.yml")
    for profile_name, profile in resources["_pav_config"].items():
        assert isinstance(profile, dict), f"{profile_name} must be a mapping"
        for key in ("merge_ins", "merge_del", "merge_inv"):
            assert key not in profile, f"{profile_name}: {key} is PAV2-only"


def test_pav_profile_schema_validation():
    """Schema accepts resources.yml profiles and rejects `reference`."""
    schema = load_yaml("schema/resources-schema.yml")
    resources = load_yaml("config/resources.yml")
    pav_schema = schema["properties"]["_pav_config"]

    jsonschema.validate(resources["_pav_config"], pav_schema)

    # reference is set by the pipeline (setup_pav.py), not by profiles
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate({"bad": {"reference": "ref.fa"}}, pav_schema)


def test_pav_container_is_pav3():
    """run_pav must use a PAV3 image."""
    resources = load_yaml("config/resources.yml")
    assert "pav3" in resources["_pav_container"]
