#!/usr/bin/zsh
## Wrapper script for running Asm benchmarking pipeline using GIAB team tibanna

## Input parameters will convert to script arguments
## Number of jobs to run
JOBS=1
RUNDIR="asm-bench-dev-jenny-test"
DISKMB=50000
#DRYRUN=""
DRYRUN="--dryrun"


### Personal
PROFILE="default"
UNICORN="tibanna_unicorn_giab_test3"

### Setting AWS credentials
KEYID=$(aws configure get aws_access_key_id --profile ${PROFILE})
export AWS_ACCESS_KEY_ID=${KEYID}

SECKEY=$(aws configure get aws_secret_access_key --profile ${PROFILE})
export AWS_SECRET_ACCESS_KEY=${SECKEY}

## Setting Tibanna Unicorn
export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}

## Running Snakemake
snakemake \
	--use-conda -j${JOBS} -p --verbose ${DRYRUN}\
	--tibanna \
	--tibanna-config \
		subnet=subnet-083080d579317ad61 \
		log_bucket=giab-tibanna-logs \
		root_ebs_size=32 \
	--precommand "cat etc/nist_dns.txt >> /etc/resolv.conf; cat /etc/resolv.conf" \
	--default-remote-prefix=giab-tibanna-runs/${RUNDIR} \
	--default-resources disk_mb=50000 \
	--rerun-incomplete \
	--keep-going
#  	--config analyses="config/resources.yml"