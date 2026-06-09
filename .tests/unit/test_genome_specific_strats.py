"""Tests for scripts/genome_specific_strats.py (#59 / #173).

Validates the ported GS-pipeline variant classification and the in-module
interval merge/subtract primitives.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from genome_specific_strats import (  # noqa: E402
    build_complex_strats,
    classify_geno2haplo,
    merge_intervals,
    merge_with_count,
    subtract_intervals,
)

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
)


def _record(chrom, pos, ref, alt, gt):
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"


def _write_vcf(path, records):
    path.write_text(VCF_HEADER + "".join(records))


# --- interval primitives -----------------------------------------------------


def test_merge_intervals_within_distance():
    assert merge_intervals([(0, 10), (40, 60), (200, 210)], dist=50) == [
        (0, 60),
        (200, 210),
    ]


def test_merge_intervals_touching_only_with_zero_dist():
    assert merge_intervals([(0, 10), (40, 60)], dist=0) == [(0, 10), (40, 60)]


def test_merge_with_count():
    assert merge_with_count([(100, 101), (105, 106), (500, 501)], dist=10) == [
        (100, 106, 2),
        (500, 501, 1),
    ]


def test_subtract_intervals_punches_hole():
    assert subtract_intervals([(0, 100)], [(40, 60)]) == [(0, 40), (60, 100)]


def test_subtract_intervals_no_overlap():
    assert subtract_intervals([(0, 10)], [(50, 60)]) == [(0, 10)]


def test_subtract_intervals_full_cover():
    assert subtract_intervals([(10, 20)], [(0, 100)]) == []


# --- variant classification (rules 2-5) --------------------------------------


def test_comphet_snp_vs_indel():
    records = [
        _record("chr1", 1000, "A", "C,G", "1/2"),  # both SNP -> comphetsnp
        _record("chr1", 2000, "A", "AT,C", "1|2"),  # one indel -> comphetindel
    ]
    buckets = classify_geno2haplo(
        (
            ("chr1", 1000, "A", ["C", "G"], frozenset({"1", "2"})),
            ("chr1", 2000, "A", ["AT", "C"], frozenset({"1", "2"})),
        )
    )
    assert buckets["comphetsnp10bp"]["chr1"] == [(950, 1051)]
    assert buckets["comphetindel10bp"]["chr1"] == [(1950, 2051)]
    assert records  # header helper sanity


def test_noncomphet_complexindel_and_snpwithin():
    buckets = classify_geno2haplo(
        (
            # complex indel: REF>1, ALT>1, differing lengths
            ("chr2", 500, "ACGT", ["AC"], frozenset({"1"})),
            # MNP-like: REF==ALT length, >1bp
            ("chr2", 800, "AC", ["GT"], frozenset({"1"})),
            # plain SNP: excluded from both
            ("chr2", 900, "A", ["G"], frozenset({"1"})),
        )
    )
    assert buckets["complexindel10bp"]["chr2"] == [(450, 554)]
    assert buckets["snpswithin10bp"]["chr2"] == [(750, 852)]
    assert "chr2" not in buckets.get("comphetsnp10bp", {})


# --- end-to-end strat file production ----------------------------------------


def test_build_complex_strats_writes_all_five(tmp_path):
    g2h = tmp_path / "g2h.vcf"
    raw = tmp_path / "raw.vcf"
    _write_vcf(
        g2h,
        [
            _record("chr1", 1000, "A", "C,G", "1/2"),
            _record("chr1", 2000, "A", "AT,C", "1/2"),
        ],
    )
    # raw has two records within 10bp (a cluster) plus the classified ones
    _write_vcf(
        raw,
        [
            _record("chr1", 3000, "A", "C", "1/1"),
            _record("chr1", 3005, "G", "T", "1/1"),
        ],
    )
    paths = build_complex_strats(str(g2h), str(raw), str(tmp_path / "out"), "PFX_")
    assert set(paths) == {
        "comphetsnp10bp",
        "comphetindel10bp",
        "complexindel10bp",
        "snpswithin10bp",
        "othercomplexwithin10bp",
    }
    for p in paths.values():
        assert Path(p).is_file()
    # the raw 2-variant cluster (3000,3005) -> slop -> othercomplex region
    other = Path(paths["othercomplexwithin10bp"]).read_text()
    assert other.startswith("chr1\t2950\t")


def test_othercomplex_subtracts_primary_strata(tmp_path):
    # A cluster that coincides with a comphetsnp region must be subtracted out.
    g2h = tmp_path / "g2h.vcf"
    raw = tmp_path / "raw.vcf"
    _write_vcf(g2h, [_record("chr1", 1000, "A", "C,G", "1/2")])
    _write_vcf(
        raw,
        [
            _record("chr1", 1000, "A", "C", "1/2"),
            _record("chr1", 1002, "T", "G", "1/2"),
        ],
    )
    paths = build_complex_strats(str(g2h), str(raw), str(tmp_path / "out"), "")
    # comphetsnp covers [950,1051); the raw cluster ~[950,1053) is mostly removed
    other = Path(paths["othercomplexwithin10bp"]).read_text()
    assert "\t950\t" not in other


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
