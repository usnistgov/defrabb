#!/bin/sh -x 
## upload defrabb run to s3
RUNID="20240523_v0.016_HG002Q100"
aws s3 sync \
    . s3://giab-data/defrabb_runs/${RUNID} \
    --exclude "*" \
    --include "results/*" \
    --include "resources/comparison*" \
    --include "resources/exclusions/*" \
    --include "config/resources.yml" \
    --include "config/analyses_${RUNID}.tsv" \
    --include "${RUNID}.tar" \
    --include "environment.yml" \
    --include "run.log" \
    --include "run_README.md" \
    --acl "public-read" "$@"
    