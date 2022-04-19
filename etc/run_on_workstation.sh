#!/bin/bash

#script to run defrabb on workstation

# SET VARIABLES WITH EACH RUN
DRYRUN=""
# DRYRUN="-n"
ANALYSES="config/analyses.tsv"
JOBS=12
RUNID="defrabb_test"
RUNDIR="../${RUNID}"
REPORTNAME="defrabb_test.report.zip"
ARCHIVENAME="defrabb_test.archive.tar.gz"

## CHOOSE WHICH SNAKEMAKE COMMNAD TO USE
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