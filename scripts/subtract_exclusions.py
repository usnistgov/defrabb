#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Dwarshuis
import sys
from pybedtools import BedTool
from itertools import accumulate

# usage: subtract_exclusions INPUT_BED OUTPUT_BED EXCLUDED_BEDS ...


def count_bp(bedfile):
    df = bedfile.to_dataframe(names=["chr", "start", "end"])
    return int((df["end"] - df["start"]).sum())


def summarize(header, beds, paths):
    print("%s:" % header)
    for b, p in zip(beds, paths):
        count = count_bp(b)
        print("%-50s %s" % (p, count))
    print("")


def get_excluded(paths):
    beds = [*map(BedTool, paths)]
    summarize("Total base pairs in each excluded bed", beds, paths)
    return beds


def exclude_beds(input_bed, excluded_beds):
    after = [
        *accumulate(
            excluded_beds,
            lambda acc, to_exclude: acc.subtract(to_exclude),
            initial=input_bed,
        )
    ]
    print("Total base pairs in initial dip.bed: %i\n" % count_bp(after[0]))
    paths = map(lambda b: b.fn, excluded_beds)
    summarize("Total base pairs in dip.bed after each exclusion", after[1:], paths)
    return after


def main():
    input_bed = BedTool(sys.argv[1])
    excluded_beds = get_excluded(sys.argv[3:])
    after_exclusions = exclude_beds(input_bed, excluded_beds)
    after_exclusions[-1].saveas(sys.argv[2])


main()
