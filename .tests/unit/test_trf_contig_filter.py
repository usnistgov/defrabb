"""Tests for the non-canonical-contig filtering logic introduced to fix #197.

truvari anno trf errors with "invalid contig" when the input VCF contains
variants on contigs absent from the TRF database (e.g. chrM in CHM13, MT in
GRCh37).  The fix adds a filter_vcf_for_trf_contigs rule that partitions the
insize VCF into canonical (TRF contigs) and non-canonical, and threads the
non-canonical variants through merge_trfanno_vcfs.py via --noncanon so they
appear unannotated in the final output.

Tests here verify:
1. Canonical variants end up in the canonical output, non-canonical in noncanon.
2. merge_trfanno_vcfs.py --noncanon correctly includes those records in the merged output.
3. The merged VCF has a unified header covering all INFO fields from all sources.
"""

import subprocess

import pytest

pysam = pytest.importorskip("pysam")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _basic_header(contigs, info_fields=None, with_sample=True):
    """Return a VariantHeader with the given contig names, INFO fields, and sample."""
    h = pysam.VariantHeader()
    # Declare GT FORMAT so sample GT assignments succeed
    h.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if info_fields:
        for field in info_fields:
            h.add_line(field)
    if with_sample:
        h.add_sample("SAMPLE")
    for name, length in contigs:
        h.contigs.add(name, length=length)
    return h


# ---------------------------------------------------------------------------
# Test 1: bcftools -t / -t^ contig filtering gives the right partitions
# ---------------------------------------------------------------------------


def test_bcftools_contig_filter_splits_correctly(tmp_path):
    """Verify that bcftools -t canonical / -t ^canonical splits as expected.

    This exercises the shell logic used by filter_vcf_for_trf_contigs without
    running the full Snakemake rule.
    """
    # Build a VCF with chr1 (canonical) and chrM (non-canonical) variants
    h = _basic_header(
        [("chr1", 100000), ("chrM", 16569)],
        ['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">'],
    )

    input_vcf = tmp_path / "input.vcf.gz"
    with pysam.VariantFile(str(input_vcf), "wz", header=h) as f:
        rec1 = h.new_record(contig="chr1", start=100, alleles=("A", "T"))
        rec1.info["SVTYPE"] = "SNV"
        rec1.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec1)

        rec2 = h.new_record(contig="chrM", start=500, alleles=("G", "C"))
        rec2.info["SVTYPE"] = "SNV"
        rec2.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec2)

    subprocess.run(["tabix", str(input_vcf)], check=True)

    canonical_vcf = tmp_path / "canonical.vcf.gz"
    noncanon_vcf = tmp_path / "noncanon.vcf.gz"
    trf_contigs = "chr1"  # pretend TRF db only covers chr1

    subprocess.run(
        f"bcftools view -t '{trf_contigs}' -Oz -o {canonical_vcf} {input_vcf}",
        shell=True,
        check=True,
    )
    subprocess.run(
        f"bcftools view -t '^{trf_contigs}' -Oz -o {noncanon_vcf} {input_vcf}",
        shell=True,
        check=True,
    )

    with pysam.VariantFile(str(canonical_vcf)) as f:
        canonical_records = list(f)
    assert len(canonical_records) == 1
    assert canonical_records[0].contig == "chr1"

    with pysam.VariantFile(str(noncanon_vcf)) as f:
        noncanon_records = list(f)
    assert len(noncanon_records) == 1
    assert noncanon_records[0].contig == "chrM"


def test_bcftools_contig_filter_empty_noncanon(tmp_path):
    """When all variants are on canonical contigs, noncanon output is empty."""
    h = _basic_header(
        [("chr1", 100000)],
        ['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">'],
    )

    input_vcf = tmp_path / "input.vcf.gz"
    with pysam.VariantFile(str(input_vcf), "wz", header=h) as f:
        rec = h.new_record(contig="chr1", start=200, alleles=("A", "T"))
        rec.info["SVTYPE"] = "SNV"
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    subprocess.run(["tabix", str(input_vcf)], check=True)

    noncanon_vcf = tmp_path / "noncanon.vcf.gz"
    trf_contigs = "chr1"
    subprocess.run(
        f"bcftools view -t '^{trf_contigs}' -Oz -o {noncanon_vcf} {input_vcf}",
        shell=True,
        check=True,
    )

    with pysam.VariantFile(str(noncanon_vcf)) as f:
        records = list(f)
    assert (
        records == []
    ), "Expected empty noncanon output when all contigs are canonical"


# ---------------------------------------------------------------------------
# Test 2: merge_trfanno_vcfs.py --noncanon passes through non-canonical records
# ---------------------------------------------------------------------------


def _run_merge(annotated, oversize, output, noncanon=None):
    cmd = [
        "python3",
        "scripts/merge_trfanno_vcfs.py",
        "--annotated",
        str(annotated),
        "--oversize",
        str(oversize),
        "--output",
        str(output),
    ]
    if noncanon is not None:
        cmd += ["--noncanon", str(noncanon)]
    subprocess.run(cmd, check=True)


def test_merge_includes_noncanon_records(tmp_path):
    """--noncanon records appear in merged output and header is unified."""
    # annotated VCF: chr1 variant with TRF annotation
    ann_header = _basic_header(
        [("chr1", 100000)],
        [
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=TRF,Number=0,Type=Flag,Description="In TR region">',
        ],
    )
    annotated_vcf = tmp_path / "annotated.vcf"
    with pysam.VariantFile(str(annotated_vcf), "w", header=ann_header) as f:
        rec = ann_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        rec.info["SVTYPE"] = "INS"
        rec.info["TRF"] = True
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    # oversize VCF: chr1 large insertion, no TRF (plain VCF — no index needed for merge)
    over_header = _basic_header(
        [("chr1", 100000)],
        [
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        ],
    )
    oversize_vcf = tmp_path / "oversize.vcf"
    with pysam.VariantFile(str(oversize_vcf), "w", header=over_header) as f:
        rec = over_header.new_record(
            contig="chr1", start=500, alleles=("A", "ATATATAT")
        )
        rec.info["SVTYPE"] = "INS"
        rec.info["SVLEN"] = 7
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    # noncanon VCF: chrM variant, no TRF (mitochondrial); plain VCF for simplicity
    nc_header = _basic_header(
        [("chrM", 16569)],
        ['##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">'],
    )
    noncanon_vcf = tmp_path / "noncanon.vcf"
    with pysam.VariantFile(str(noncanon_vcf), "w", header=nc_header) as f:
        rec = nc_header.new_record(contig="chrM", start=1000, alleles=("G", "A"))
        rec.info["SVTYPE"] = "SNV"
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    merged_vcf = tmp_path / "merged.vcf"
    _run_merge(annotated_vcf, oversize_vcf, merged_vcf, noncanon=noncanon_vcf)

    with pysam.VariantFile(str(merged_vcf)) as f:
        info_ids = {r["ID"] for r in f.header.records if r.key == "INFO"}
        assert "SVTYPE" in info_ids
        assert "TRF" in info_ids
        assert "SVLEN" in info_ids

        records = list(f)

    assert len(records) == 3, f"Expected 3 records, got {len(records)}"
    contigs = [r.contig for r in records]
    assert "chr1" in contigs
    assert "chrM" in contigs

    # The chrM record should be present and unannotated (no TRF flag)
    chrM_records = [r for r in records if r.contig == "chrM"]
    assert len(chrM_records) == 1
    assert "TRF" not in chrM_records[0].info


def test_merge_without_noncanon_unchanged(tmp_path):
    """merge_trfanno_vcfs.py behaves identically when --noncanon is omitted."""
    ann_header = _basic_header(
        [("chr1", 100000)],
        ['##INFO=<ID=TRF,Number=0,Type=Flag,Description="In TR region">'],
    )
    annotated_vcf = tmp_path / "annotated.vcf"
    with pysam.VariantFile(str(annotated_vcf), "w", header=ann_header) as f:
        rec = ann_header.new_record(contig="chr1", start=100, alleles=("A", "T"))
        rec.info["TRF"] = True
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    over_header = _basic_header([("chr1", 100000)])
    oversize_vcf = tmp_path / "oversize.vcf"
    with pysam.VariantFile(str(oversize_vcf), "w", header=over_header) as f:
        rec = over_header.new_record(contig="chr1", start=500, alleles=("A", "ATTT"))
        rec.samples["SAMPLE"]["GT"] = (0, 1)
        f.write(rec)

    merged_vcf = tmp_path / "merged.vcf"
    _run_merge(annotated_vcf, oversize_vcf, merged_vcf)  # no --noncanon

    with pysam.VariantFile(str(merged_vcf)) as f:
        records = list(f)
    assert len(records) == 2
