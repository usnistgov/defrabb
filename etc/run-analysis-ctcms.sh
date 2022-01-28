#!/usr/bin/bash

## Start pipeline
snakemake \
	--use-conda -j 20 \
	--local-cores 1 \
	--profile etc/slurm  \
	--verbose \
	--printshellcmds \
	--config analyses="config/0.003-analyses.tsv"

## TODO look into mounting resdata or NAS for exporting analysis output
## -d /working/geneteam/defrabb-runs/${time}_${repotag}  \
##		--mem {cluster.mem} \
