"""Tests for scripts/build_stratification_tsv.py (#173 / #59 shared core).

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from build_stratification_tsv import (  # noqa: E402
    format_strat_rows,
    parse_strat_args,
    write_strat_tsv,
)


def test_parse_strat_args():
    assert parse_strat_args(["a=x.bed", "b=y/z.bed"]) == [
        ("a", "x.bed"),
        ("b", "y/z.bed"),
    ]


def test_parse_strat_args_rejects_malformed():
    for bad in ["noequals", "=nopath", "noname="]:
        with pytest.raises(ValueError):
            parse_strat_args([bad])


def test_format_rows_relative_to_out_dir():
    entries = [("segdups", "/data/results/excl/segdups.bed")]
    rows = format_strat_rows(entries, "/data/results", relative=True)
    assert rows == ["segdups\texcl/segdups.bed"]


def test_format_rows_absolute_when_not_relative():
    entries = [("segdups", "/data/results/excl/segdups.bed")]
    rows = format_strat_rows(entries, "/data/results", relative=False)
    assert rows == ["segdups\t/data/results/excl/segdups.bed"]


def test_duplicate_names_raise():
    with pytest.raises(ValueError):
        format_strat_rows([("a", "1.bed"), ("a", "2.bed")], "/tmp", relative=False)


def test_write_strat_tsv_roundtrip(tmp_path):
    out = tmp_path / "strats.tsv"
    bed = tmp_path / "sub" / "ex.bed"
    bed.parent.mkdir()
    bed.write_text("chr1\t0\t10\n")
    write_strat_tsv([("ex", str(bed))], str(out), relative=True)
    assert out.read_text() == "ex\tsub/ex.bed\n"


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
