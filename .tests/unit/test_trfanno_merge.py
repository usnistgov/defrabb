"""Test merge_trfanno_vcfs.py header merging and record translation."""

import tempfile
import subprocess
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")


def test_merge_trfanno_vcfs_header():
    """Verify merged VCF has unified header with all INFO defs."""
    # Create test annotated VCF (has TRF INFO)
    annotated_header = pysam.VariantHeader()
    annotated_header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">'
    )
    annotated_header.add_line(
        '##INFO=<ID=TRF,Number=0,Type=Flag,Description="In TR region">'
    )
    annotated_header.add_line(
        '##INFO=<ID=TRFrepeat,Number=1,Type=String,Description="Repeat motif">'
    )
    annotated_header.add_sample("sample")
    annotated_header.contigs.add("chr1", length=1000)

    # Create test oversize VCF (no TRF INFO)
    oversize_header = pysam.VariantHeader()
    oversize_header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">'
    )
    oversize_header.add_line(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">'
    )
    oversize_header.add_sample("sample")
    oversize_header.contigs.add("chr1", length=1000)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        annotated_vcf = tmpdir / "annotated.vcf"
        oversize_vcf = tmpdir / "oversize.vcf"
        merged_vcf = tmpdir / "merged.vcf"

        # Write annotated VCF
        with pysam.VariantFile(annotated_vcf, "w", header=annotated_header) as f:
            rec = annotated_header.new_record(
                contig="chr1", start=100, alleles=("A", "T")
            )
            rec.info["SVTYPE"] = "SNV"
            rec.info["TRF"] = True
            rec.samples["sample"]["GT"] = (0, 1)
            f.write(rec)

        # Write oversize VCF
        with pysam.VariantFile(oversize_vcf, "w", header=oversize_header) as f:
            rec = oversize_header.new_record(
                contig="chr1", start=500, alleles=("A", "T" * 100000)
            )
            rec.info["SVTYPE"] = "INS"
            rec.info["SVLEN"] = 99999
            rec.samples["sample"]["GT"] = (0, 1)
            f.write(rec)

        # Merge
        subprocess.run(
            [
                "python3",
                "scripts/merge_trfanno_vcfs.py",
                "--annotated",
                str(annotated_vcf),
                "--oversize",
                str(oversize_vcf),
                "--output",
                str(merged_vcf),
            ],
            check=True,
        )

        # Verify merged header
        with pysam.VariantFile(merged_vcf) as f:
            info_ids = [rec["ID"] for rec in f.header.records if rec.key == "INFO"]
            assert "SVTYPE" in info_ids
            assert "SVLEN" in info_ids
            assert "TRF" in info_ids
            assert "TRFrepeat" in info_ids

            # Verify records readable and INFO preserved
            records = list(f)
            assert len(records) == 2
            assert records[0].info["SVTYPE"] == "SNV"
            assert "TRF" in records[0].info
            assert records[1].info["SVTYPE"] == "INS"
            assert records[1].info["SVLEN"] == 99999


def test_merge_trfanno_vcfs_bcf_roundtrip():
    """Verify merged VCF can be BCF-encoded without INFO tag ID errors."""
    annotated_header = pysam.VariantHeader()
    annotated_header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">'
    )
    annotated_header.add_line(
        '##INFO=<ID=TRF,Number=0,Type=Flag,Description="In TR region">'
    )
    annotated_header.add_sample("sample")
    annotated_header.contigs.add("chr1", length=1000)

    oversize_header = pysam.VariantHeader()
    oversize_header.add_line(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">'
    )
    oversize_header.add_line(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">'
    )
    oversize_header.add_sample("sample")
    oversize_header.contigs.add("chr1", length=1000)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        annotated_vcf = tmpdir / "annotated.vcf"
        oversize_vcf = tmpdir / "oversize.vcf"
        merged_vcf = tmpdir / "merged.vcf"
        merged_bcf = tmpdir / "merged.bcf"

        # Write test VCFs
        with pysam.VariantFile(annotated_vcf, "w", header=annotated_header) as f:
            rec = annotated_header.new_record(
                contig="chr1", start=100, alleles=("A", "T")
            )
            rec.info["SVTYPE"] = "SNV"
            rec.info["TRF"] = True
            rec.samples["sample"]["GT"] = (0, 1)
            f.write(rec)

        with pysam.VariantFile(oversize_vcf, "w", header=oversize_header) as f:
            rec = oversize_header.new_record(
                contig="chr1", start=500, alleles=("A", "ATATATAT")
            )
            rec.info["SVTYPE"] = "INS"
            rec.info["SVLEN"] = 7
            rec.samples["sample"]["GT"] = (0, 1)
            f.write(rec)

        # Merge
        subprocess.run(
            [
                "python3",
                "scripts/merge_trfanno_vcfs.py",
                "--annotated",
                str(annotated_vcf),
                "--oversize",
                str(oversize_vcf),
                "--output",
                str(merged_vcf),
            ],
            check=True,
        )

        # Convert to BCF (will fail if INFO tag IDs corrupted)
        result = subprocess.run(
            ["bcftools", "view", "-Ob", "-o", str(merged_bcf), str(merged_vcf)],
            check=True,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

        # Read back BCF and verify
        with pysam.VariantFile(merged_bcf) as f:
            records = list(f)
            assert len(records) == 2
            assert records[0].info["SVTYPE"] == "SNV"
            assert records[1].info["SVTYPE"] == "INS"


def test_merge_trfanno_vcfs_preserves_phasing(tmp_path):
    """Phased GTs (1|1, 0|1) must stay phased and unphased stay unphased.

    Regression for the smvar/stvar phasing loss (619 SNP FN/FP vs v5.0q).
    """
    header = pysam.VariantHeader()
    header.add_line('##INFO=<ID=TRF,Number=0,Type=Flag,Description="In TR">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_sample("sample")
    header.contigs.add("chr1", length=1000)

    paths = {k: tmp_path / f"{k}.vcf" for k in ("annotated", "oversize", "noncanon")}
    records = {
        "annotated": [(100, (1, 1), True), (200, (0, 1), True), (300, (0, 1), False)],
        "oversize": [(500, (1, 0), True)],
        "noncanon": [],
    }
    for key, recs in records.items():
        with pysam.VariantFile(paths[key], "w", header=header) as f:
            for pos, gt, phased in recs:
                rec = header.new_record(contig="chr1", start=pos, alleles=("A", "T"))
                rec.samples["sample"]["GT"] = gt
                rec.samples["sample"].phased = phased
                f.write(rec)

    merged = tmp_path / "merged.vcf"
    subprocess.run(
        [
            "python3",
            "scripts/merge_trfanno_vcfs.py",
            "--annotated",
            str(paths["annotated"]),
            "--oversize",
            str(paths["oversize"]),
            "--noncanon",
            str(paths["noncanon"]),
            "--output",
            str(merged),
        ],
        check=True,
    )

    with pysam.VariantFile(merged) as f:
        got = {
            r.pos - 1: (r.samples["sample"]["GT"], r.samples["sample"].phased)
            for r in f
        }
    assert got == {
        100: ((1, 1), True),
        200: ((0, 1), True),
        300: ((0, 1), False),
        500: ((1, 0), True),
    }
