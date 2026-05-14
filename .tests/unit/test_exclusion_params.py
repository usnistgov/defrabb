"""Tests for exclusion slop/merge resolvers and provenance script.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from write_exclusion_provenance import (  # noqa: E402
    _resolve_slop,
    _resolve_merge_dist,
    _source_type,
    _transform,
    build_provenance,
    load_bp_impact,
)


# ---------------------------------------------------------------------------
# Minimal config fixture
# ---------------------------------------------------------------------------
def _cfg(overrides=None):
    return {
        "_exclusion_params": {
            "slop": 15000,
            "slopmerge_dist": 10000,
            "overrides": overrides or {},
        },
        "exclusion_set": {
            "test_set": ["tandem-repeats", "segdups", "gaps"],
        },
        "exclusion_slop_regions": ["tandem-repeats", "gaps", "self-chains"],
        "exclusion_slopmerge_regions": ["segdups", "satellites"],
        "exclusion_asm_intersect": ["segdups", "tandem-repeats", "satellites", "self-chains"],
        "exclusion_asm_agnostic": ["gaps"],
        "exclusion_ref_agnostic": ["flanks", "svs-and-simple-repeats"],
    }


# ---------------------------------------------------------------------------
# Slop resolver
# ---------------------------------------------------------------------------
def test_resolve_slop_no_overrides_returns_global():
    val, is_pct = _resolve_slop("tandem-repeats", _cfg())
    assert val == 15000
    assert is_pct is False


def test_resolve_slop_bp_override():
    val, is_pct = _resolve_slop("tandem-repeats", _cfg({"tandem-repeats": {"slop": 20000}}))
    assert val == 20000
    assert is_pct is False


def test_resolve_slop_pct_override():
    val, is_pct = _resolve_slop("tandem-repeats", _cfg({"tandem-repeats": {"slop_pct": 0.25}}))
    assert val == 0.25
    assert is_pct is True


def test_resolve_slop_override_for_different_region_does_not_affect_others():
    val, is_pct = _resolve_slop("gaps", _cfg({"tandem-repeats": {"slop": 20000}}))
    assert val == 15000
    assert is_pct is False


# ---------------------------------------------------------------------------
# Merge dist resolver
# ---------------------------------------------------------------------------
def test_resolve_merge_dist_no_overrides_returns_global():
    assert _resolve_merge_dist("segdups", _cfg()) == 10000


def test_resolve_merge_dist_override():
    assert _resolve_merge_dist("segdups", _cfg({"segdups": {"slopmerge_dist": 15000}})) == 15000


# ---------------------------------------------------------------------------
# Source type / transform helpers
# ---------------------------------------------------------------------------
def test_source_type_asm_agnostic():
    assert _source_type("gaps", _cfg()) == "asm_agnostic"


def test_source_type_ref_agnostic():
    assert _source_type("flanks", _cfg()) == "ref_agnostic"


def test_source_type_asm_specific():
    assert _source_type("segdups", _cfg()) == "asm_specific"


def test_transform_slop():
    assert _transform("tandem-repeats", _cfg()) == "slop"


def test_transform_slopmerge():
    assert _transform("segdups", _cfg()) == "slopmerge"


def test_transform_none():
    assert _transform("flanks", _cfg()) == "none"


# ---------------------------------------------------------------------------
# build_provenance
# ---------------------------------------------------------------------------
def test_build_provenance_structure():
    cfg = _cfg()
    prov = build_provenance("test_set", cfg, bp_impact={})
    assert prov["exclusion_set"] == "test_set"
    assert prov["global_defaults"] == {"slop": 15000, "slopmerge_dist": 10000}
    names = [e["name"] for e in prov["exclusions"]]
    assert names == ["tandem-repeats", "segdups", "gaps"]


def test_build_provenance_slop_entry():
    prov = build_provenance("test_set", _cfg(), bp_impact={})
    tr = next(e for e in prov["exclusions"] if e["name"] == "tandem-repeats")
    assert tr["transform"] == "slop"
    assert tr["slop_bp"] == 15000
    assert "merge_dist" not in tr
    assert tr["asm_intersect"] is True


def test_build_provenance_slopmerge_entry():
    prov = build_provenance("test_set", _cfg(), bp_impact={})
    sd = next(e for e in prov["exclusions"] if e["name"] == "segdups")
    assert sd["transform"] == "slopmerge"
    assert sd["slop_bp"] == 15000
    assert sd["merge_dist"] == 10000


def test_build_provenance_pct_override_recorded():
    cfg = _cfg({"tandem-repeats": {"slop_pct": 0.25}})
    prov = build_provenance("test_set", cfg, bp_impact={})
    tr = next(e for e in prov["exclusions"] if e["name"] == "tandem-repeats")
    assert "slop_pct" in tr
    assert tr["slop_pct"] == 0.25
    assert "slop_bp" not in tr


def test_build_provenance_bp_impact_merged():
    bp = {"tandem-repeats": {"bp": 5000000, "pct": 2.1}}
    prov = build_provenance("test_set", _cfg(), bp_impact=bp)
    tr = next(e for e in prov["exclusions"] if e["name"] == "tandem-repeats")
    assert tr["bp_impact"] == 5000000
    assert tr["pct_of_initial"] == pytest.approx(2.1)
