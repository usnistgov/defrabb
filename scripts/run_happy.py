#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Olson
from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

## Optional parameters
engine = snakemake.params.get("engine", "")
if engine:
	engine = "--engine {}".format(engine)

truth_regions = snakemake.input.get("truth_regions", "")
if truth_regions:
	truth_regions = "-f {}".format(truth_regions)

target_regions = snakemake.input.get("target_regions","")
if target_regions:
	target_regions = "-T {}".format(target_regions)

## Extracting stratification tarball
shell(
	"(tar -xvf {snakemake.input.strat_tb})"
	"{log}"
)

## Running Happy
shell(
	"(hap.py "
	"    --threads {snakemake.params.threads} "
	"    {engine} "
	"    -r {snakemake.input.genome}  "
	"    {truth_regions} "
	"    --stratification {snakemake.params.strat_tsv} "
	"    -o {snakemake.params.prefix} "
	"    --verbose "
	"    {target_regions} "
	"    {snakemake.input.truth} "
	"    {snakemake.input.query} )"
	"{log}"
)