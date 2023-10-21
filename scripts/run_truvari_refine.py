#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Olson
from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

## Optional parameters
query_regions = snakemake.input.get("query_regions", "")
if query_regions:
    query_regions = f"--regions {query_regions}"

## Running Truvari
shell(
    "rm -rf {snakemake.params.refine_output} && "
    "truvari refine --threads {snakemake.threads} "
    "--align mafft "
    "--use-original "
    "--reference {snakemake.input.ref} "
    "{query_regions} "
    "{snakemake.params.bench_output} "
    "&> {snakemake.log}",
)
