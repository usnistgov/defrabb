#!/bin/bash

#script to run defrabb on workstation

# SET VARIABLES WITH EACH RUN
DRYRUN=""
# DRYRUN="-n"
JOBS=12
RUNID="defrabb_test"
ANALYSES="config/analyses_${RUNID}.tsv"
RUNDIR="../${RUNID}"
REPORTNAME="${RUNID}.report.zip"
ARCHIVENAME="${RUNID}.archive.tar.gz"

## CHOOSE WHICH SNAKEMAKE COMMNAD TO USE #######################################
## Run Pipeline 
snakemake --use-conda -p --verbose \
	--config analyses=${ANALYSES} \
	${DRYRUN} \
	--jobs ${JOBS} \
	--directory ${RUNDIR}

## Generating Report
# snakemake \
# 	--config analyses=${ANALYSES} \
# 	--directory ${RUNDIR} \
# 	--report ${REPORTNAME}

## Making snakemake archive
# snakemake \
# 	--use-conda \
# 	--config analyses=${ANALYSES} \
# 	--directory ${RUNDIR} \
# 	--archive ${ARCHIVENAME}

## Archiving run - syncing run directory with NAS
# rsync -rv \
# 	--exclude=.snakemake \
# 	--exclude=resources \
# 	${RUNDIR} \
# 	/mnt/bbdhg-nas/analysis/defrabb-runs/