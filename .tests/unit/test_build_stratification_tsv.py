"""Tests for scripts/build_stratification_tsv.py (#173 / #59 shared core).

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from build_stratification_tsv import (  # noqa: E402
    combine_strat_tsvs,
    format_strat_rows,
    parse_strat_args,
    read_strat_tsv,
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


def test_read_strat_tsv_roundtrip(tmp_path):
    tsv = tmp_path / "s.tsv"
    tsv.write_text("a\tx.bed\nb\tsub/y.bed\n")
    assert read_strat_tsv(str(tsv)) == [("a", "x.bed"), ("b", "sub/y.bed")]


def test_combine_strat_tsvs_rewrites_extra_paths(tmp_path):
    # GIAB tsv + its beds live under giab/
    giab_dir = tmp_path / "GRCh38@all"
    giab_dir.mkdir()
    giab_tsv = giab_dir / "strats.tsv"
    giab_tsv.write_text("segdups\tunion/segdups.bed.gz\n")

    # genome-specific tsv + beds live elsewhere
    gs_dir = tmp_path / "results" / "gs"
    gs_dir.mkdir(parents=True)
    (gs_dir / "comphetsnp.bed").write_text("chr1\t0\t10\n")
    gs_tsv = gs_dir / "gs.tsv"
    gs_tsv.write_text("comphetsnp10bp\tcomphetsnp.bed\n")

    out_tsv = giab_dir / "combined.tsv"
    n = combine_strat_tsvs(str(giab_tsv), str(gs_tsv), str(out_tsv))
    assert n == 1

    rows = read_strat_tsv(str(out_tsv))
    # GIAB row copied verbatim
    assert rows[0] == ("segdups", "union/segdups.bed.gz")
    # extra row rewritten relative to the combined tsv's dir (giab_dir),
    # and must resolve back to the real bed file
    name, rel = rows[1]
    assert name == "comphetsnp10bp"
    assert (giab_dir / rel).resolve() == (gs_dir / "comphetsnp.bed").resolve()


def test_combine_strat_tsvs_requires_same_dir(tmp_path):
    giab = tmp_path / "a" / "g.tsv"
    giab.parent.mkdir()
    giab.write_text("x\ty.bed\n")
    extra = tmp_path / "e.tsv"
    extra.write_text("z\tz.bed\n")
    with pytest.raises(ValueError):
        combine_strat_tsvs(str(giab), str(extra), str(tmp_path / "b" / "out.tsv"))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
