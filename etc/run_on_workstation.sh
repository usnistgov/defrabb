#!/bin/bash

#script to run defrabb on workstation

# SET VARIABLES WITH EACH RUN
DRYRUN=""
#DRYRUN="-n"
ANALYSES="config/analyses_2022408_v0.005-HG002-HPRC.tsv"
JOBS="-j 12"
REPORTNAME="2022408_v0.005-HG002-HPRC.report.zip"
ARCHIVENAME="2022408_v0.005-HG002-HPRC.archive.tar.gz"

## CHOOSE WHICH SNAKEMAKE COMMNAD TO USE
#snakemake --use-conda -p --verbose --config analyses=${ANALYSES} _dipcall_threads=3 ${DRYRUN} ${JOBS} --rerun-incomplete --keep-going 

#snakemake --report ${REPORTNAME} --config analyses=${ANALYSES}

snakemake --archive ${ARCHIVENAME} --config analyses=${ANALYSES}
