#!/usr/bin/env python3
"""
Merge TRF-annotated and oversized VCFs with unified header.

Reads two VCFs:
  - annotated: variants with TRF INFO annotations
  - oversize: variants without TRF annotations (bypassed TRF)

Outputs a single VCF with unified header containing all INFO defs.
"""
import argparse

import pysam


def merge_vcfs(annotated_path, oversize_path, output_path):
    """Merge annotated and oversize VCFs with unified header."""
    # Open inputs
    vcf_annotated = pysam.VariantFile(annotated_path)
    vcf_oversize = pysam.VariantFile(oversize_path)

    # Build unified header: start with annotated (has TRF defs), merge in any missing from oversize
    out_header = vcf_annotated.header.copy()
    for rec in vcf_oversize.header.records:
        if rec.key == "INFO":
            # Add INFO def if not already present
            if rec["ID"] not in [
                r["ID"] for r in out_header.records if r.key == "INFO"
            ]:
                out_header.add_record(rec)

    # Add END INFO field if missing (pysam auto-generates it for symbolic alleles with SVLEN)
    if "END" not in [r["ID"] for r in out_header.records if r.key == "INFO"]:
        out_header.add_line(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
        )

    # Open output - mode "w" writes uncompressed VCF by default
    vcf_out = pysam.VariantFile(output_path, "w", header=out_header)

    # Write annotated records
    for entry in vcf_annotated:
        # Create new record in output header context to avoid BCF ID corruption
        # Don't pass stop= to avoid auto-generating END INFO field
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
        vcf_out.write(new_entry)

    # Write oversize records
    for entry in vcf_oversize:
        new_entry = vcf_out.new_record(
            contig=entry.contig,
            start=entry.start,
            alleles=entry.alleles,
            id=entry.id,
            qual=entry.qual,
            filter=entry.filter,
        )
        for key in entry.info:
            new_entry.info[key] = entry.info[key]
        for sample in entry.samples:
            for fmt_key in entry.samples[sample]:
                new_entry.samples[sample][fmt_key] = entry.samples[sample][fmt_key]
        vcf_out.write(new_entry)

    vcf_annotated.close()
    vcf_oversize.close()
    vcf_out.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotated", required=True, help="TRF-annotated VCF")
    parser.add_argument(
        "--oversize", required=True, help="Oversized (un-annotated) VCF"
    )
    parser.add_argument("--output", required=True, help="Output merged VCF")
    args = parser.parse_args()

    merge_vcfs(args.annotated, args.oversize, args.output)


if __name__ == "__main__":
    main()
