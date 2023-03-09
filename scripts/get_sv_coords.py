#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: ND Olson
import sys
import pandas as pd
from os import path
from snakemake.shell import shell

# log = snakemake.log_fmt_shell(stdout=True, stderr=True)


def load_sv_tbl(tbl_path):
    ## Read table as pandas data frame
    df = pd.read_table(
        tbl_path,
        header=None,
        names=["CHROM", "POS", "REF", "ALT", "REFWIDENED", "REPTYPE"],
    )
    ## Parsing ref widen coords
    refwiden_df = df.REFWIDENED.str.split(r".*:(.*)-(.*)", expand=True)
    df["refwiden_start"] = refwiden_df[1].astype("int")
    df["refwiden_end"] = refwiden_df[2].astype("int")

    ## Adding SV end position coordinates
    df["REF_len"] = [len(x) for x in df["REF"]]
    df["ALT_len"] = [len(x) for x in df["ALT"]]
    df["END"] = df["POS"] + df["REF_len"] + df["ALT_len"]

    ## Filtering small variants
    df = df.loc[(df["REF_len"] > 49) | (df["ALT_len"] > 49)]

    ## Returing data frame with ref widen and SV coords
    return df


def get_sv_exclusion_coords(sv_df):
    ## SV widen coordinates does not always include original SV
    ## - using min and max of svWiden coordinates and original SV coordinates
    sv_min = []
    sv_max = []
    for index, row in sv_df.iterrows():
        row_min = min(row["POS"], row["refwiden_start"], row["refwiden_end"])
        sv_min.append(row_min)
        row_max = max(row["END"], row["refwiden_start"], row["refwiden_end"])
        sv_max.append(row_max)

    sv_df["min"] = sv_min
    sv_df["max"] = sv_max
    exclusion_coords = sv_df.loc[:, ["CHROM", "min", "max"]]

    return exclusion_coords


def write_sv_bed(sv_coords, out_path):
    ## Convert to pybedtools object and sort
    ## Print pybedtools object
    sv_coords.to_csv(out_path, header=False, index=False, sep="\t")


def main():
    ## Load table with sv and ref widen coords
    sv_tbl_path = snakemake.input[0]
    # sv_tbl_path = "GRCh38_T2T-XY-v2.7_dipcall-z2k_dip_SVs.tsv"
    sv_df = load_sv_tbl(sv_tbl_path)
    ## Find new min/max coordinates
    sv_exclusion_coords = get_sv_exclusion_coords(sv_df)

    # Output coordinates as sorted bed file
    output_path = snakemake.output[0]
    # output_path = "test.bed"
    write_sv_bed(sv_exclusion_coords, output_path)


main()
