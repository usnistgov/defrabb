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

    # Open output
    vcf_out = pysam.VariantFile(output_path, "w", header=out_header)

    # Write annotated records
    for entry in vcf_annotated:
        entry.translate(out_header)
        vcf_out.write(entry)

    # Write oversize records
    for entry in vcf_oversize:
        entry.translate(out_header)
        vcf_out.write(entry)

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
