#!/usr/bin/env python3

"""
get_sv_widen_coords.py: Create a bed file with SV coordinates from a VCF file with coordinates widened to include overlapping trandem repeat regions.

Usage: python get_sv_widen_coords.py --input input.vcf --output output.bed [--verbose] [--table] [--sort-merge --genome genome.txt --max-distance 5]

Author: ND Olson with OpenAI
Date: 2023-09-22
"""

import argparse
import logging
import sys
import os
import gzip
from pybedtools import BedTool


def extract_info_from_vcf_line(line):
    """
    Extract necessary information from a VCF line.

    Args:
    - line (str): A single line from a VCF file.

    Returns:
    - tuple: A tuple containing chromosome, position, reference allele,
             alternative allele, TRF start, and TRF end.
    """
    fields = line.strip().split("	")
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    info = fields[7]

    trf_start = None
    trf_end = None
    for item in info.split(";"):
        if item.startswith("TRFstart="):
            trf_start = int(item.split("=")[1])
        elif item.startswith("TRFend="):
            trf_end = int(item.split("=")[1])

    return chrom, pos, ref, alt, trf_start, trf_end


def main(
    input_vcf,
    output_bed,
    verbose=False,
    log_file=None,
    table_format=False,
    sort_merge=False,
    genome=None,
    max_distance=5,
):
    """
    Main function to convert VCF to BED.

    Args:
    - input_vcf (str): Path to the input VCF file.
    - output_bed (str): Path to the output BED file.
    - verbose (bool, optional): Whether to print verbose logs. Defaults to False.
    """
    logging_level = logging.DEBUG if verbose else logging.WARNING
    log_format = "%(asctime)s - %(levelname)s - %(message)s"

    if log_file:
        logging.basicConfig(filename=log_file, level=logging_level, format=log_format)
    else:
        logging.basicConfig(level=logging_level, format=log_format)

    # Determine if input is compressed based on file extension
    is_compressed = input_vcf.endswith(".gz")

    # Open the file with the appropriate method
    open_method = gzip.open if is_compressed else open
    mode = "rt"  # Read as text, even if using gzip

    # If input is '-', it's stdin, and we assume it's uncompressed
    if input_vcf == "-":
        open_method = lambda x, y: sys.stdin

    input_stream = open_method(input_vcf, mode)
    output_stream = open(output_bed, "w") if output_bed != "-" else sys.stdout

    bed_entries = []
    table_entries = []

    with input_stream, output_stream:
        for line in input_stream:
            if not line.startswith("#"):
                chrom, pos, ref, alt, trf_start, trf_end = extract_info_from_vcf_line(
                    line
                )

                # Table output format
                if table_format:
                    table_entries.append(
                        f"{chrom}	{pos}	{len(ref)}	{len(alt)}	{trf_start or '-'}	{trf_end or '-'}\n"
                    )

                ## Including TR region for trf annotated SVs
                if trf_start:
                    start = min(pos, trf_start, trf_end) - 50
                    end = max(pos + len(ref), trf_start, trf_end) + 50
                else:
                    start = pos - 50
                    end = pos + len(alt) + 50

                # Adjust for base encoding: VCF is 1-based and BED is 0-based
                start -= 1

                ## Only writing coordinates for SVs 50bp or larger
                if (len(alt) > 49) | (len(ref) > 49):
                    bed_entries.append((chrom, start, end))

        # If sort and merge
        if sort_merge:
            bed = BedTool(bed_entries).sort(g=genome).merge(d=max_distance)
            with open(output_bed, "w") as outfile:
                for interval in bed:
                    outfile.write(
                        f"{interval.chrom}\t{interval.start}\t{interval.end}\n"
                    )

        # Write table output if --table is specified
        if table_format:
            table_output_file = output_bed.rsplit(".", 1)[0] + ".tsv"
            with open(table_output_file, "w") as outfile:
                outfile.write(
                    "chrom\tpos\tlength_ref\tlength_alt\ttrf_start\ttrf_end\n"
                )
                for entry in table_entries:
                    outfile.write(entry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a VCF file with tandem repeat annotations to a BED file and optionally a table."
    )
    parser.add_argument(
        "--input", default="-", help="Path to the input VCF file. Use '-' for stdin."
    )
    parser.add_argument(
        "--output", default="-", help="Path to the output BED file. Use '-' for stdout."
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose logging."
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Path to the log file. If not provided, logs are printed to console.",
    )
    parser.add_argument(
        "--table",
        action="store_true",
        help="Output in table format in addition to BED format.",
    )
    parser.add_argument(
        "--sort-merge", action="store_true", help="Sort and merge the BED entries."
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=5,
        help="Maximum distance between two intervals to be merged. Default is 5.",
    )
    parser.add_argument(
        "--genome",
        default=None,
        help="Path to the genome index file for sorting. Required if --sort-merge is used.",
    )

    args = parser.parse_args()

    if args.sort_merge and not args.genome:
        print("Error: --genome is required when using --sort-merge.")
        sys.exit(1)

    main(
        args.input,
        args.output,
        args.verbose,
        args.log,
        args.table,
        args.sort_merge,
        args.genome,
        args.max_distance,
    )
