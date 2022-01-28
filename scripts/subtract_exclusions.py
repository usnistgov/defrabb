#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Dwarshuis
import sys
import pandas as pd
from pybedtools import BedTool
from itertools import accumulate

# usage: subtract_exclusions INPUT_BED OUTPUT_BED EXCLUDED_BEDS ...


def count_bp(bedfile):
    df = bedfile.to_dataframe(names=["chr", "start", "end"])
    return int((df["end"] - df["start"]).sum())


def print_summary(after, excluded):
    # ASSUME: 'excluded' will have one less than 'after', since 'after' will
    # have the initial bed file before any exclusions were performed
    paths = ["initial"] + [*map(lambda b: b.fn, excluded)]
    bed_lengths = [0] + [*map(count_bp, excluded)]
    after_lengths = [*map(count_bp, after)]
    pd.DataFrame(
        [*zip(paths, bed_lengths, after_lengths)],
        columns=["exclusion", "exclusion_length", "resulting_length"],
    ).to_csv(sys.stdout, index=False, sep="\t")


def get_excluded(paths):
    return [*map(BedTool, paths)]


def exclude_beds(input_bed, excluded_beds):
    after = [
        *accumulate(
            excluded_beds,
            lambda acc, to_exclude: acc.subtract(to_exclude),
            initial=input_bed,
        )
    ]
    return after


def main():
    input_bed = BedTool(sys.argv[1])
    excluded_beds = get_excluded(sys.argv[3:])
    after_exclusions = exclude_beds(input_bed, excluded_beds)
    print_summary(after_exclusions, excluded_beds)
    after_exclusions[-1].saveas(sys.argv[2])


main()
