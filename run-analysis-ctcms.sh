#!/usr/bin/bash

## Start pipeline
snakemake \
	--use-conda -j 5 \
	--profile slurm  \
	--verbose \
	--printshellcmds

## TODO look into mounting resdata or NAS for exporting analysis output
## -d /working/geneteam/defrabb-runs/${time}_${repotag}  \
##		--mem {cluster.mem} \
