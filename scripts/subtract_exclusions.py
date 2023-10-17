#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Dwarshuis
import sys
from os import stat as os_stat
import pandas as pd
from pybedtools import BedTool
from itertools import accumulate
import warnings

"""
subtract_exclusions script:

Subtracts regions specified in multiple BED files from an input BED file,
prints a summary of exclusions, and writes the resulting BED to an output file.

Usage: subtract_exclusions INPUT_BED OUTPUT_BED EXCLUDED_BEDS ...
"""

VALID_CHROMOSOMES = (
    [f"chr{i}" for i in range(1, 23)] + [range(1, 23)] + ["chrX", "chrY", "X", "Y"]
)


def count_bp(bedfile):
    """
    Count the total base pairs in a BED file.

    Parameters:
    - bedfile (BedTool): A BedTool object of the BED file to count base pairs.

    Returns:
    - int: Total number of base pairs in the BED.
    """
    df = bedfile.to_dataframe(names=["chr", "start", "end"])
    ## To avoid error with empty bed files
    if len(df.index) < 1:
        return int(0)
    return int((df["end"] - df["start"]).sum())


def print_summary(after, excluded, exclusion_stats):
    """
    Print a summary of exclusions.

    Parameters:
    - after (list): List of BedTool objects after exclusions.
    - excluded (list): List of BedTool objects of excluded regions.
    - exclusion_stats (str): Path to write the summary file.
    """
    # ASSUME: 'excluded' will have one less than 'after', since 'after' will
    # have the initial bed file before any exclusions were performed
    paths = ["initial"] + [*map(lambda b: b.fn, excluded)]
    bed_lengths = [0] + [*map(count_bp, excluded)]
    after_lengths = [*map(count_bp, after)]
    pd.DataFrame(
        [*zip(paths, bed_lengths, after_lengths)],
        columns=["exclusion", "exclusion_length", "resulting_length"],
    ).to_csv(exclusion_stats, index=False, sep="\t")


def get_excluded(paths):
    """
    Load BED files from the paths and exclude empty ones.

    Parameters:
    - paths (list): List of file paths.

    Returns:
    - list: List of BedTool objects of non-empty BED files.
    """
    non_empty_paths = []
    for path in paths:
        if os_stat(path).st_size != 0:
            non_empty_paths.append(path)
        else:
            warnings.warn(f"The BED file at path {path} is empty and will be ignored.")

    return [*map(BedTool, non_empty_paths)]


def exclude_beds(input_bed, excluded_beds):
    """
    Subtract excluded regions from the input BED file.

    Parameters:
    - input_bed (BedTool): A BedTool object of the input BED file.
    - excluded_beds (list): List of BedTool objects to be subtracted.

    Returns:
    - list: List of BedTool objects after each subtraction.
    """
    after = [
        *accumulate(
            excluded_beds,
            lambda acc, to_exclude: acc.subtract(to_exclude),
            initial=input_bed,
        )
    ]
    return after


def filter_chromosomes(bedtool, valid_chromosomes):
    """
    Filter out regions from BED file not in valid chromosomes.

    Parameters:
    - bedtool (BedTool): A BedTool object to filter.
    - valid_chromosomes (list): List of valid chromosome names.

    Returns:
    - BedTool: Filtered BedTool object.
    """
    filtered = bedtool.filter(lambda b: b.chrom in valid_chromosomes)
    ## Making sure the filtered bed bedtools object is file backed,
    ## not doing so resulted in downstream errors
    filt_bed = filtered.saveas()
    return filt_bed


def main():
    input_bed = BedTool(sys.argv[1])
    # Filter only valid chromosomes
    filtered_bed = filter_chromosomes(input_bed, VALID_CHROMOSOMES)
    # Exclude regions from filtered input bed
    excluded_beds = get_excluded(sys.argv[4:])
    after_exclusions = exclude_beds(filtered_bed, excluded_beds)
    # Print summary and excluded bed
    print_summary(after_exclusions, excluded_beds, sys.argv[3])
    after_exclusions[-1].saveas(sys.argv[2])


if __name__ == "__main__":
    main()
