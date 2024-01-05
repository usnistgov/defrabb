#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Nathan Olson
from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

## Optional parameters
## - using candidate.refine.bed produced by truvari bench as 
##   truvari refine hangs on bcftool consensus step when 
##   bed file with large regions is used https://github.com/ACEnglish/truvari/issues/182
# query_regions = snakemake.input.get("query_regions", "")
# if query_regions:
#     query_regions = f"--regions {query_regions}"
query_regions = f"--regions {snakemake.input.refine_bed}"

## Running Truvari
shell(
    "rm -rf {snakemake.params.bench_output}/phab_bench && "
    "truvari refine --threads {snakemake.threads} "
    "--align mafft "
    "--use-original "
    "--reference {snakemake.input.ref} "
    "{query_regions} "
    "{snakemake.params.bench_output} " + log
)
