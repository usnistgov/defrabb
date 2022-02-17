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
A small dataset with chr 21 for GRCh38 [^1] was developed for testing and development. 
The resources are either included in the repository or are hosted in a GIAB S3 bucket with appropriate urls included in the `resources.yml`
small example datasets are included / made available for testing framework code. 
Test pipeline for runtime errors using `snakemake --use-conda -j1`

# Running Framework on CTCMS Cluster
Documentation for running pipeline on CTCMS using tmux and the ctcms profile.
CTCMS snakemake cluster deployment profile generated using [Snakemake Slurm profile template](https://github.com/Snakemake-Profiles/slurm).

__Steps__

1. Log into CTCMS using `ssh username@ruth.nist.gov` or `ssh ctcms` if `.ssh/config` is setup to do so (ask Nate O. if you want to do this).
1. Create `tmux` session `tmux new-session -s [session name]`. This will create a detachable terminal session so that the pipeline does not fail if the ssh connection dropped.
1. Switch to appropriate git branch `git checkout [branchname]`, and make sure up to date `git status` and `git pull`. 
1. Activate conda environment for running pipeline `conda activate defrabb`. This environment was create by Nate O. using `mamba env create -f envs/env.yml` from the root directory of this repository.
1. Use complete `config/analyses.tsv` with required run information, update `config/resources.yml` if necessary.
1. Run pipeline using `sh etc/run-analysis-ctcms.sh`


## Notes / Gotchas for the CTCMS cluster

* snakemake is executed on headnode as job nodes are not connected to the internet, which is required for conda environments, snakemake wrappers, and downloading resource files.
* Need to first create conda environments using `snakemake --use-conda --conda-create-envs-only`
* Can not define job memory requirement, for jobs with high memory requirements try increasing the number of threads.

# Running Framework on AWS with Tibanna
Low memory jobs, defined in Snakefile, are run on "local" instance. Higher memory, compute intensive jobs will be run on instances started by Tibanna with appropriate resources. 
1. Set up AWS instance [^2] e.g. i3.large w/ 2vCPU, 15GB and 1TB storage (less storage migth be possible if fewer assemblies will be run)
2. Start instance, e.g., (ssh command for your instance can be found under "Connect" button)\
`ssh -v -i "~/.ssh/user.pem" ec2-user@10.208.44.53`
3. install dependencies
	- miniconda\
`wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh`\
`sh Miniconda3-py37_4.11.0-Linux-x86_64.sh`\
`echo $HOME (/home/ec2-user)`\
`export PATH="$HOME/miniconda3/bin:$PATH"`\
	- tmux\
`sudo yum install tmux`\
	- git\
`sudo yum install git`\
	- mamba\
`conda install -c conda-forge mamba`\
	- nano (this might already be installed)
4. initialize conda for shell interaction
conda init bash (needed to restart shell after this)
5. clone defrabb git repo and switch desired branch\
`git clone https://gitlab.nist.gov/gitlab/njd2/giab-asm-bench-whole-genome.git`\
`git checkout desired-branch`
6. use mamba to set up defrabb environment\
`mamba env create -f envs/env.yml -n defrabb`
7. Add AWS credentials\
get "DEFAULT" credentials `cat ~/.aws/credentials`\
`aws configure` entering the following from your file above when propted
	- AWS Access Key ID
	- AWS Secret Access Key
	- Default region name 
	- Default output format
8. Set up directory in S3 bucket `giab-tibanna-runs` for run output
9. Make any necessary changes to `etc/run_on_tibanna_giab.sh`
	- DRYRUN = comment/uncomment, it is suggested you start with a dry run first
	- ANALYSES = add path to file, e.g. `config/myfile.tsv` if you want to overide use of `anslyses.tsv` defined in `resources.yml` 
	- RUNDIR = directory in S3 bucket `giab-tibanna-runs` that outputs should go to
	- JOBS and DISKMB = adjust as appropriate for your run requirements
## preparing to start run on AWS instance
1. start tmux session. See [online tmux cheatsheet](https://tmuxcheatsheet.com) for helpful tmux commands [^3]\
`tmux new-session -s my-session-name`
2. activate defrabb environment\
`conda activate defrabb`
3. start run w/ tibanna\
`sh etc/run_on_tibanna_giab.sh`
4. You can look at tibanna logs for each job
	- list the jobs, in this example 10\
`tibanna stat -s tibanna_unicorn_giab_test3 -n 10`\
	- job IDs will be on far left, copy jobID\
`tibanna log -s tibanna_unicorn_giab_test3 -j <paste jobID#>`
## Notes / Gotchas for AWS instance
- Tibanna and all jobs started with tibanna are run in docker containters.  There is a limit on how many docker containers can be started in a given time.
- If you are downloading numerous large files and run into issues with downloads failing the instance might be out of storage and you might have to increase size of storage through EC2 console after stopping instance. 

# General Execution and Documenting Analysis Runs
1. Use snakedeploy to create run directory (future work)
1. Update relevant config files
1. Run snakemake (see documentation above)
1. Create snakemake report
	activate conda environment - `conda activate defrabb`
	generate report - `snakemake --report [milestone]-report.html`
	Will want to run on CTCMs headnode as it requires a network connection
1. Create snakemake archive for rerunning analyses
	In conda environment `snakemake --archive [YYYYMMDD_milestone].tar.gz`
1. Create directory and copy results, logs, and benchmark files for archiving analysis run
	- Creating directory `mkdir ../defrabb-runs/[YYYYMMDD_milestone]`
	- Moving report and archive tarball to run archive directory `mv [YYYYMMDD_milestone]* ../defrabb-runs/[YYYYMMDD_milestone]/`
	- Copying results to `cp -r results ../defrabb-runs/[YYYYMMDD_milestone]`
	- Copying results to `cp -r logs ../defrabb-runs/[YYYYMMDD_milestone]`
	- Copying results to `cp -r benchmark ../defrabb-runs/[YYYYMMDD_milestone]`

1. Fill out README with relevant run information - framework repo info - [milestone] tag (with some potential - hopefully minor-differences), who ran the framework and where/ how, justification / reasoning for analyses, JZ notes (what did we learn)

1. If run on CTCMS, copy run archive to local directory `rsync -rv --progress ctcms:/working/geneteam/defrabb-runs ~/Desktop` for upload to google drive and storage on NAS in `giab/analyses/defrabb-runs`. Directory set to automatically sync with `BBD_Human_Genomics/defrabb_runs` team google drive directory.  

Automating - copy output and config files to directory for archiving, script to automate archiving with call for report, updating analysis run log google sheet

<!-- Resources/ Citations -->

# Footnotes
[^1]: Chromosome 13 was included in the test dataset reference as dipcall incorrectly included line breaks in dip.vcf when only chr21 was included. We might want to submit an issue to the dipcall repo about this issue.
[^2]: [Team resource on setting up AWS instance](https://docs.google.com/document/d/1IdAKastyUShjVl_8msWSR-n1ReKuhXpI_fqLcU829QQ/edit)
[^3]: tmux shortcut commands like `ctrl + b s`  means press `b` while holding down `ctrl` then release both `ctrl` and `b` then press `s` alone. 
