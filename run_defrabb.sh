#!/usr/bin/env bash
## defrabb wrapper script for executing and archiving framework

## Script usage text
usage() { 
    local err=${1:-""};
    cat <<EOF
Usage: $0 [options] 
Required:
    -r STRING   Analysis RUN ID, please use following naming convention YYYYMMDD_milestone_brief-id

Optional:
    -a FILE     defrabb run analysis table, if not provided assumes at config/analyses_[RUN ID].tsv
    -o DIR      output directory for framework run, pipeline will create a named directory [RUN ID] at defined location, default is "../"
    -n          Run snakemake in dry run mode
EOF
>&2;
    echo -e "\n$err" >&2;
    exit 1;
}

## Parsing command line arguments
dry_run=""

while getopts "r:a:o:n" flag; do
    case "${flag}" in
        r) runid=${OPTARG};;
        a) analysis_file=${OPTARG};;
        o) out_dir=${OPTARG};;
        n) dry_run="-n";;
        *) usage;;
    esac
done
shift $((OPTIND-1))

extra_args=""
if [ -z "${runid}" ]; then
    usage "Missing required parameter -r";
fi

## Getting system resource information
source etc/common.sh
cores=$(find_core_limit)
mem_gb=$(find_mem_limit_gb)
log "Number of cores: $cores"
log "Memory limit: $mem_gb GB"

# SET VARIABLES WITH EACH RUN
### setting run analyses table path
if [ -z "${analyses_file}" ]; then
   analyses_file="config/analyses_${runid}.tsv"
fi

## Stating analysis table path
echo "Analysis Table Path: ${analysis_file}"

### setting run directory path
if [ -z "${out_dir}" ]; then
   run_dir="../${runid}"
else
   run_dir="${out_dir}/${run_id}"
fi

## setting report name
report_name="${runid}.report.zip"

## setting archive path
smk_archive_path="${run_dir}/${runid}.archive.tar.gz"

### Directory on NAS used for archiving runs
archive_dir="/mnt/bbdhg-nas/analysis/defrabb-runs/"

## Activating mamba environment
##   TODO add check to see if 

# Run Snakemake pipeline
set -euo pipefail

snakemake \
  --printshellcmds \
  --reason \
  --rerun-incomplete \
  --jobs "${cores}" \
  --resources "mem_gb=${mem_gb}" \
  --use-conda \
  --config analyses=${analyses_file} \
  --directory ${run_dir} \
  ${dry_run} \
  ${extra_args};


log "Done Executing DeFrABB"

## TODO
### - add conditional exit for dry run

## Generating Report
snakemake \
	--config analyses=${analyses_file} \
	--directory ${run_dir} \
	--report ${report_name} \
	${dry_run};

log "Done Generating Report"

## Making snakemake archive
snakemake \
  --use-conda \
  --config analyses=${analyses_file} \
  --archive ${smk_archive_path};

log "Done Making Snakemake Archive"

## Archiving run - syncing run directory with NAS
rsync -arv \
	--exclude=.snakemake \
	--exclude=resources \
	${run_dir} \
	${archive_dir};

log "Done Archiving Run"

log "DeFrABB run execution and archiving complete!"

### Resources for bash scripting and snakemake wrappers
## Example snakemake pipeline wrapper script
# https://github.com/marbl/verkko/blob/master/src/verkko.sh
# https://github.com/GooglingTheCancerGenome/sv-callers/blob/master/run.sh
# Handling command line options

## Collection of example bash scripts
# https://github.com/swoodford/aws/blob/master/vpc-sg-rename-group.sh


## Example python wrapper script
# https://github.com/CVUA-RRW/FooDMe/blob/master/foodme.py
## Python click based CLI
# https://github.com/gymrek-lab/haptools/blob/main/haptools/__main__.py

## Template repo with run wrapper script and utility functions
# https://github.com/fulcrumgenomics/python-snakemake-skeleton
## Bash script templates
# https://github.com/ralish/bash-script-template/blob/main/template.sh
