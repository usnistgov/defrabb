#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Olson
from os import path
import sys
from snakemake.shell import shell

sys.stdout = open(snakemake.log[0], "a")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

## Optional parameters
engine = snakemake.params.get("engine", "")
if engine:
    engine = f"--engine {engine}"

truth_regions = snakemake.input.get("truth_regions", "")
if truth_regions:
    truth_regions = f"-f {truth_regions}"

target_regions = snakemake.input.get("target_regions", "")
if target_regions:
    target_regions = f"-T {target_regions}"

## Extracting stratification tarball
## Can add if statement with eval_params to see if starts are used
ref_id = snakemake.wildcards.ref_id
strat_id = snakemake.config["references"][ref_id]["stratifications"]["id"]
strat_tsv = f"{snakemake.params.strat_tsv}"

print("Extracting Stratifications")
shell("tar -xf {snakemake.input.strat_tb}" + log)

if path.isfile(strat_tsv):
    print("Stratification tsv file present")
else:
    print(f"stratifications file, {strat_tsv}, not present!!! help!")


## Running Happy
shell(
    "(hap.py "
    "    --threads {snakemake.params.threads} "
    "    {engine} "
    "    {snakemake.params.engine_extra} "
    "    -r {snakemake.input.genome}  "
    "    {truth_regions} "
    "    --stratification {strat_tsv} "
    "    -o {snakemake.params.prefix} "
    "    --verbose "
    "    {target_regions} "
    "    {snakemake.input.truth} "
    "    {snakemake.input.query} )" + log
)
