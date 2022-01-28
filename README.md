# GIAB Assembly-Based Benchmark Set Development Framework
<!--

[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/CI/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/Linting/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)
[![Actions Status](https://github.com/mrvollger/SmkTemplate/workflows/black/badge.svg)](https://github.com/mrvollger/SmkTemplate/actions)

This is a Snakemake project template. The `Snakefile` is under `workflow`.

[Slides](https://mrvollger.github.io/SmkTemplate/slides) describing and justifying the use of this template.
-->
<!--gitlab badges-->
[![pipeline status](https://gitlab.nist.gov/gitlab/bbd-human-genomics/giab-asm-bench-whole-genome/badges/master/pipeline.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/giab-asm-bench-whole-genome/-/commits/master)
[![coverage report](https://gitlab.nist.gov/gitlab/bbd-human-genomics/giab-asm-bench-whole-genome/badges/master/coverage.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/giab-asm-bench-whole-genome/-/commits/master)

<!-- Background -->
The GIAB benchmark set development framework is a snakemake bioinformatic pipeline for the development of transparent and reproducible genome assembly based small and structural variant benchmark sets. 

<!-- Usage -->
# Usage
##  Initial Setup and Test Run
1. clone repository `git clone https://gitlab.nist.gov/gitlab/njd2/giab-asm-bench-whole-genome.git`
2. generate conda snakemake environment `mamba create -n snakemake --file envs/env.yml`
3. Activate environment `conda activate snakemake`
4. Run built in test analyses `snakemake --use-conda -j1`

## Running a defined analysis set
1. Define analysis runs by editing `config/analyses.tsv`
2. Update `config/resources.yml` as necessary.
3. Run pipeline using `snakemake --use-conda -j1`

## Development and Testing
small example datasets are included / made available for testing framework code. 
Test pipeline for runtime errors using `snakemake --use-conda -j1`

# Running Framework on CTCMS Cluster
Documentation for running pipeline on CTCMS using tmux and the ctcms profile.
CTCMS snakemake cluster deployment profile generated using ADD GIT REPO

__Steps__

1. Log into CTCMS using `ssh username@ruth.nist.gov` or `ssh ctcms` if `.ssh/config` is setup to do so (ask Nate O. if you want to do this).
1. Create `tmux` session `tmux new-session -s [session name]`.
	This will create a detachable terminal session so that the pipeline does not fail if the ssh connection dropped.
1. Switch to appropriate git branch `git checkout [branchname]`, and make sure up to date `git status` and `git pull`. 
1. Activate conda environment for running pipeline `conda activate defrabb`. 
	This environment was create by Nate O. using `mamba env create -f envs/env.yml` from the root directory of this repository.
1. Use complete `config/analyses.tsv` with required run information, update `config/resources.yml` if necessary.
1. Run pipeline using `sh etc/run-analysis-ctcms.sh`


## Notes / Gotchas for the CTCMS cluster

* snakemake is executed on headnode as job nodes are not connected to the internet.
* Need to first create conda environments using `snakemake --use-conda --conda-create-envs-only`
* Can not define job memory requirement, for jobs with high memory requirements try increasing the number of threads.

# General Execution and Documenting Analysis Runs
1. Use snakedeploy to create run directory (future work)
1. Update relevant config files
1. Run snakemake (see documentation above)
1. Create snakemake report (snakemake --archive 20220119_v0.002.tar.gz)
	activate conda environment - `conda activate defrabb`
	generate report - `snakemake --report v0.002-report.html`
	Will want to run on CTCMs headnode as it requires a network connection
1. Create snakemake archive for rerunning analyses
	In conda environment `snakemake --archive 20220119_v0.002.tar.gz`
1. Create directory and copy results files for archiving analysis run
	- Creating directory `mkdir ../defrabb-runs/20220119_v0.002`
	- Moving report and archive tarball to run archive directory `mv 20220119* ../defrabb-runs/20220119_v0.002/`
	- Copying results to `cp -r results ../defrabb-runs/20220119_v0.002/`

1. Fill out README with relevant run information - framework repo info - v0.002 tag (with some potential - hopefully minor-differences), who ran the framework and where/ how, justification / reasoning for analyses, JZ notes (what did we learn)

1. Copy run archive to local directory `   rsync -rv --progress ctcms:/working/geneteam/defrabb-runs ~/Desktop` for upload to google drive and storage on NAS or resdata (this is a temporary hack until we work out an automated process, will want to check with Andrew Reid on best way to mount, can then use NAS utility to copy files to google drive). 

Automating - copy output and config files to directory for archiving, script to automate archiving with call for report, updating analysis run log google sheet

<!-- Resources/ Citations -->
