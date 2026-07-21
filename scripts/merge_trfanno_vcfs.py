#!/usr/bin/env python3
"""
Merge TRF-annotated and oversized VCFs with unified header.

Reads two or three VCFs:
  - annotated: variants with TRF INFO annotations (canonical contigs, in-size)
  - oversize: variants without TRF annotations (bypassed TRF due to insertion size)
  - noncanon (optional): variants on contigs absent from the TRF database (e.g. chrM,
    MT), preserved unannotated to avoid truvari "invalid contig" errors (#197)

Outputs a single VCF with unified header containing all INFO defs.
"""
import argparse

import pysam


def _copy_record(entry, vcf_out):
    """Copy a variant record into the vcf_out header context."""
    # Create new record to avoid BCF INFO tag ID corruption across headers.
    # Don't pass stop= to avoid auto-generating END INFO field.
    new_entry = vcf_out.new_record(
        contig=entry.contig,
        start=entry.start,
        alleles=entry.alleles,
        id=entry.id,
        qual=entry.qual,
        filter=entry.filter,
    )
    # Copy INFO fields by name (not by BCF ID)
    for key in entry.info:
        new_entry.info[key] = entry.info[key]
    # Copy samples
    for sample in entry.samples:
        for fmt_key in entry.samples[sample]:
            new_entry.samples[sample][fmt_key] = entry.samples[sample][fmt_key]
    return new_entry


def _merge_info_defs(out_header, source_header):
    """Add any INFO defs from source_header not already in out_header."""
    existing_ids = {r["ID"] for r in out_header.records if r.key == "INFO"}
    for rec in source_header.records:
        if rec.key == "INFO" and rec["ID"] not in existing_ids:
            out_header.add_record(rec)
            existing_ids.add(rec["ID"])


def _merge_contigs(out_header, source_header):
    """Add any contig definitions from source_header not already in out_header."""
    existing_contigs = set(out_header.contigs)
    for contig_name, contig in source_header.contigs.items():
        if contig_name not in existing_contigs:
            length = contig.length if contig.length is not None else 0
            out_header.contigs.add(contig_name, length=length)
            existing_contigs.add(contig_name)


def merge_vcfs(annotated_path, oversize_path, output_path, noncanon_path=None):
    """Merge annotated, oversize, and optional non-canonical-contig VCFs."""
    # Open inputs
    vcf_annotated = pysam.VariantFile(annotated_path)
    vcf_oversize = pysam.VariantFile(oversize_path)
    vcf_noncanon = pysam.VariantFile(noncanon_path) if noncanon_path else None

    # Build unified header: start with annotated (has TRF defs), merge in any
    # missing INFO defs and contig definitions from the pass-through sources.
    out_header = vcf_annotated.header.copy()
    _merge_info_defs(out_header, vcf_oversize.header)
    _merge_contigs(out_header, vcf_oversize.header)
    if vcf_noncanon is not None:
        _merge_info_defs(out_header, vcf_noncanon.header)
        _merge_contigs(out_header, vcf_noncanon.header)

    # Add END INFO field if missing (pysam auto-generates it for symbolic alleles with SVLEN)
    if "END" not in {r["ID"] for r in out_header.records if r.key == "INFO"}:
        out_header.add_line(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
        )

    # Open output - mode "w" writes uncompressed VCF by default
    vcf_out = pysam.VariantFile(output_path, "w", header=out_header)

    # Write annotated records (canonical contigs, in-size, TRF-annotated)
    for entry in vcf_annotated:
        vcf_out.write(_copy_record(entry, vcf_out))

    # Write oversize pass-through records (canonical contigs, insertion too large)
    for entry in vcf_oversize:
        vcf_out.write(_copy_record(entry, vcf_out))

    # Write non-canonical-contig pass-through records (e.g. chrM / MT), unannotated
    if vcf_noncanon is not None:
        for entry in vcf_noncanon:
            vcf_out.write(_copy_record(entry, vcf_out))

    vcf_annotated.close()
    vcf_oversize.close()
    if vcf_noncanon is not None:
        vcf_noncanon.close()
    vcf_out.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotated", required=True, help="TRF-annotated VCF")
    parser.add_argument(
        "--oversize", required=True, help="Oversized (un-annotated) VCF"
    )
    parser.add_argument(
        "--noncanon",
        default=None,
        help="Non-canonical-contig VCF (contigs absent from TRF db, e.g. chrM/MT)",
    )
    parser.add_argument("--output", required=True, help="Output merged VCF")
    args = parser.parse_args()

    merge_vcfs(args.annotated, args.oversize, args.output, noncanon_path=args.noncanon)


if __name__ == "__main__":
    main()
