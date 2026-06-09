#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Olson
import sys
from os import path

from snakemake.shell import shell

sys.stdout = open(snakemake.log[0], "a")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

## Optional parameters
engine = snakemake.params.get("engine", "")
if engine:
    engine = f"--engine {engine}"

truth_regions = snakemake.input.get("truth_regions", "")
if truth_regions:
    truth_regions = f"-R {truth_regions}"

target_regions = snakemake.input.get("target_regions", "")
if target_regions:
    target_regions = f"-T {target_regions}"

gender = snakemake.params.get("gender_param", "")

verbose = "--verbose" if snakemake.config.get("debug") else ""

## Extracting stratification tarball
## Can add if statement with eval_params to see if starts are used
ref_id = snakemake.wildcards.ref_id
strat_id = snakemake.config["references"][ref_id]["stratifications"]["id"]
strat_tsv = f"{snakemake.params.strat_tsv}"

print("Extracting Stratifications")
shell("tar --skip-old-files -xf {snakemake.input.strat_tb}" + log)

if path.isfile(strat_tsv):
    print("Stratification tsv file present")
else:
    print(f"stratifications file, {strat_tsv}, not present!!! help!")

## Opt-in: merge genome-specific stratifications into the GIAB strat TSV
## (#59/#173). Only present when config["genome_specific_strats"] is enabled for
## an smvar evaluation; absent otherwise, so default behavior is unchanged.
gs_strats = snakemake.input.get("genome_specific_strats", "")
if gs_strats:
    sys.path.insert(0, "scripts")
    from build_stratification_tsv import combine_strat_tsvs

    combined_tsv = path.join(path.dirname(strat_tsv), "with_genome_specific.tsv")
    n_added = combine_strat_tsvs(strat_tsv, gs_strats, combined_tsv)
    print(f"Added {n_added} genome-specific stratifications -> {combined_tsv}")
    strat_tsv = combined_tsv


## Running Happy
shell(
    "(hap.py "
    "    --threads {snakemake.params.threads} "
    "    {engine} "
    "    {snakemake.params.engine_extra} "
    "    {gender} "
    "    -r {snakemake.input.genome}  "
    "    {truth_regions} "
    "    --stratification {strat_tsv} "
    "    -o {snakemake.params.prefix} "
    "    {verbose} "
    "    {target_regions} "
    "    {snakemake.input.truth} "
    "    {snakemake.input.query} )" + log
)
