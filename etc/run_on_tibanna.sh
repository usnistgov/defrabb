#!/usr/bin/zsh
## Wrapper script for running Asm benchmarking pipeline using GIAB team tibanna

## Input parameters will convert to script arguments
## Number of jobs to run
JOBS=50
RUNDIR="asm-bench-dev"
DISKMB=50000
DRYRUN=""
#DRYRUN="--dryrun"


### Personal
PROFILE="default"
UNICORN="tibanna_unicorn_giab_test3"

## Setting Tibanna Unicorn
export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}

## Running Snakemake
time \
	snakemake \
		--use-conda -j${JOBS} -p --verbose ${DRYRUN}\
		--tibanna \
		--tibanna-config \
			subnet=subnet-083080d579317ad61 \
			log_bucket=giab-tibanna-logs \
			root_ebs_size=32 \
			spot_instance=True \
			behavior_on_capacity_limit=retry_without_spot \
		--precommand "cat etc/nist_dns.txt >> /etc/resolv.conf; cat /etc/resolv.conf" \
		--default-remote-prefix=giab-tibanna-runs/${RUNDIR} \
		--default-resources disk_mb=50000 \
		--forceall