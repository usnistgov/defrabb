#!/bin/bash

#script to run defrabb on workstation

# SET VARIABLES WITH EACH RUN
#DRYRUN=""
DRYRUN="-n"
ANALYSES="config/analyses_20220413_v0.006-HG002-HPRC-CHM13v2.tsv"
JOBS="-j 12"
REPORTNAME="20220413_v0.006-HG002-HPRC-CHM13v2.report.zip"
ARCHIVENAME="20220413_v0.006-HG002-HPRC-CHM13v2.archive.tar.gz"

## CHOOSE WHICH SNAKEMAKE COMMNAD TO USE
snakemake --use-conda -p --verbose --config analyses=${ANALYSES} _dipcall_threads=3 ${DRYRUN} ${JOBS} --rerun-incomplete --keep-going 

#snakemake --report ${REPORTNAME} --config analyses=${ANALYSES}

#snakemake --archive ${ARCHIVENAME} --config analyses=${ANALYSES}
