"""Tests for scripts/varcall_lookup.py — cross-caller run resolution.

Regression guard for the run_pav/run_dipcall duplicate-run failure: exclusions
that need another caller's output must reuse the existing run for the same
reference + assembly, not trigger a duplicate. See
``docs/issues/run_pav_run_dipcall_failures.md``.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from varcall_lookup import resolve_asm_varcall_run  # noqa: E402


def make_vc_tbl(rows):
    """rows: list of (vc_id, ref, asm_id, vc_cmd, vc_param_id)."""
    df = pd.DataFrame(
        rows, columns=["vc_id", "ref", "asm_id", "vc_cmd", "vc_param_id"]
    ).set_index("vc_id")
    return df


# A PAV benchmark and a dipcall run for the same ref+asm coexist (the exact
# shape that caused the 20260615 failure).
MIXED = make_vc_tbl(
    [
        (
            "GRCh38_HG002-T2TQ100v1.1-dipz2k",
            "GRCh38",
            "HG2-T2TQ100-V1.1",
            "dipcall",
            "z2k",
        ),
        ("GRCh38_HG002-T2TQ100v1.1-pav", "GRCh38", "HG2-T2TQ100-V1.1", "pav", "giab"),
        (
            "CHM13_HG002-T2TQ100v1.1-dipz2k",
            "CHM13v2.0",
            "HG2-T2TQ100-V1.1",
            "dipcall",
            "z2k",
        ),
    ]
)


def test_pav_bench_resolves_to_matching_dipcall_run():
    """A PAV benchmark's consecutive-svs exclusion must reuse the dipcall run
    for the same ref+asm, not the PAV run."""
    vc_id, param = resolve_asm_varcall_run(
        MIXED, "GRCh38", "HG2-T2TQ100-V1.1", "dipcall"
    )
    assert vc_id == "GRCh38_HG002-T2TQ100v1.1-dipz2k"
    assert param == "z2k"


def test_dipcall_bench_resolves_to_itself():
    vc_id, _param = resolve_asm_varcall_run(
        MIXED, "CHM13v2.0", "HG2-T2TQ100-V1.1", "dipcall"
    )
    assert vc_id == "CHM13_HG002-T2TQ100v1.1-dipz2k"


def test_no_matching_run_raises():
    with pytest.raises(ValueError, match="No dipcall run found"):
        resolve_asm_varcall_run(MIXED, "GRCh37", "HG2-T2TQ100-V1.1", "dipcall")


def test_ambiguous_run_raises():
    """If two dipcall runs exist for one ref+asm, the appropriate params are
    ambiguous and we must fail loudly rather than guess."""
    ambiguous = make_vc_tbl(
        [
            ("GRCh38_a-dipz2k", "GRCh38", "asmA", "dipcall", "z2k"),
            ("GRCh38_a-dipdefault", "GRCh38", "asmA", "dipcall", "default"),
        ]
    )
    with pytest.raises(ValueError, match="Multiple dipcall runs"):
        resolve_asm_varcall_run(ambiguous, "GRCh38", "asmA", "dipcall")
