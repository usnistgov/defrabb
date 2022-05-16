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
# Framework Diagrams
Detailed diagram by Jenny https://lucid.app/lucidchart/a1feb68c-838b-4851-8259-8289d8cd5c53/edit?invitationId=inv_977874a3-d753-4518-869d-3f0a8ca5eb2c&page=0_0#
High-level diagram by Nate -https://lucid.app/lucidchart/aea8aae0-c550-420d-80df-95b73c0cc840/edit?page=0_0#

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

## Executing and Archiving DeFrABB Analysis Runs
Steps below assume running defrabb on workstation with team NAS mounted.
1. Fill out `config/analyses_[YYYYMMDD_milestone_brief-id]` and update `config/resources.yml` if necessary.
1. Run pipeline using `./run_defrabb.sh` providing run in using the defined format (i.e. `-r [YYYYMMDD_milestone_brief-id]`) or `-r` along with `-a`.

```
Usage: ./run_defrabb.sh [options] 
Required:
    -r STRING   Analysis RUN ID, please use following naming convention YYYYMMDD_milestone_brief-id

Optional:
    -a FILE     defrabb run analysis table, if not provided assumes at config/analyses_[RUN ID].tsv
    -o DIR      output directory for framework run, pipeline will create a named directory [RUN ID] at defined location, default is "../"
    -n          Run snakemake in dry run mode
```
1. Fill out README with relevant run information - framework repo info - [milestone] tag (with some potential - hopefully minor-differences), who ran the framework and where/ how, justification / reasoning for analyses, JZ notes (what did we learn), use [defrabb run README template](https://docs.google.com/document/d/1yTXP-3OQxXfGl7kIyXWMTac-USMgiMNPhz10GXwBro0/edit?usp=sharing).
1. Add run information to the [defrabb run log spreadsheet](https://docs.google.com/spreadsheets/d/183LuUat1VCJo2dL7fu0LFMOy8CBA5FTo4WyOVsx4U6o/edit?usp=sharing)

# Running DeFrABB using NIST Compute Resources (NIST Internal Use Only)
## Running Framework on CTCMS Cluster
Documentation for running pipeline on CTCMS using tmux and the ctcms profile.
CTCMS snakemake cluster deployment profile generated using [Snakemake Slurm profile template](https://github.com/Snakemake-Profiles/slurm).

__Steps__

1. Log into CTCMS using `ssh username@ruth.nist.gov` or `ssh ctcms` if `.ssh/config` is setup to do so (ask Nate O. if you want to do this).
1. Create `tmux` session `tmux new-session -s [session name]`. This will create a detachable terminal session so that the pipeline does not fail if the ssh connection dropped.
1. Switch to appropriate git branch `git checkout [branchname]`, and make sure up to date `git status` and `git pull`. 
1. Activate conda environment for running pipeline `conda activate defrabb`. This environment was create by Nate O. using `mamba env create -f envs/env.yml` from the root directory of this repository.
1. Use complete `config/analyses.tsv` with required run information, update `config/resources.yml` if necessary.
1. Run pipeline using `sh etc/run-analysis-ctcms.sh`


### Notes / Gotchas

* snakemake is executed on headnode as job nodes are not connected to the internet, which is required for conda environments, snakemake wrappers, and downloading resource files.
* Need to first create conda environments using `snakemake --use-conda --conda-create-envs-only`
* Can not define job memory requirement, for jobs with high memory requirements try increasing the number of threads.

## Running Framework on AWS with Tibanna
### Setting up an instance on AWS to run defrabb
Low memory jobs, defined in Snakefile, are run on "local" instance. Higher memory, compute intensive jobs will be run on instances started by Tibanna with appropriate resources. 
1. Set up AWS instance [^2] e.g. i3.large w/ 2vCPU, 15GB and 1TB storage (less storage migth be possible if fewer assemblies will be run)
2. Start instance. The ssh connection command can be found under the "Connect" button in the EC2 console.  Note: you will need to add the path to your `user.pem` file\
`ssh -v -i "~/.ssh/user.pem" ec2-user@##.###.##.##`
3. install dependencies
	- miniconda, find appropriate version of miniconda [here](https://docs.conda.io/en/latest/miniconda.html)\
`wget https://repo.anaconda.com/miniconda/<your_miniconda>.sh`\
`sh <your_miniconda>.sh`\
`echo $HOME` e.g. (/home/ec2-user)\
`export PATH="$HOME/miniconda3/bin:$PATH"`
	- tmux\
`sudo yum install tmux`
	- git\
`sudo yum install git`
	- mamba\
`conda install -c conda-forge mamba`
	- nano (this might already be installed)
4. initialize conda for shell interaction
`conda init bash` (needed to restart shell after this)
5. clone defrabb git repo and switch desired branch\
`git clone https://gitlab.nist.gov/gitlab/njd2/giab-asm-bench-whole-genome.git`\
`git checkout desired-branch`
6. use mamba to set up defrabb environment\
`mamba env create -f envs/env.yml -n defrabb`
7. Add AWS credentials\
get "DEFAULT" credentials ,e.g., `cat ~/.aws/credentials`\
`aws configure` entering the following from your file above when prompted
	- AWS Access Key ID
	- AWS Secret Access Key
	- Default region name 
	- Default output format

### Starting a run on AWS instance that has been configured (see above)
1. Set up `your-directory` using name specific to run, in S3 bucket `giab-tibanna-runs` for run output
2. Setup `config/analyses.tsv`, see `schema/analyses-schema.yml` for description of fields in `analyses.tsv`
3. Setup `config/resources.yml` with specifics for the run(s) outlined in `config/analyses.tsv`.  Make sure to review dipcall and hap.py resources at the bottom of file.
4. Make any necessary changes to `etc/run_on_tibanna_giab.sh`
	- DRYRUN = comment/uncomment, it is suggested you start with a dry run first
	- ANALYSES = add path to file, e.g. `config/myfile.tsv`, if you want to overide use of `anslyses.tsv` defined in `resources.yml` 
	- RUNDIR = `your-directory` in S3 bucket `giab-tibanna-runs` that outputs should go to
	- JOBS and DISKMB = adjust as appropriate for your run requirements
	- set tibanna unicorn  
		If running on AWS from local terminal you can set/use unicorn two ways. Unicorn we are currently using for most work is "tibanna_unicorn_giab_test3"  
		
		**OPTION 1, use if going to consistenly use the same unicorn:**  
		Add the following to your ~ /.bash_profile  
			`export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=<insert unicorn name>`  
		In run_on_tibanna.sh change as follows such that default is used:  
			`PROFILE="default"`  
			`#UNICORN="<insert unicorn name>"`  
			`#export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}`  
		
		**OPTION 2, set in shell script if you tend to use different unicorns or are running from AWS instance**  
		In run_on_tibanna.sh change as follows such that unicorn is defined:  
			`#PROFILE="default"`  
			`UNICORN="<insert unicorn name>"`  
			`export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=${UNICORN}`  
	
	- review snakemake command and options used by snakemake and tibanna to ensure they are appropriate for your run
5. start tmux session. See [online tmux cheatsheet](https://tmuxcheatsheet.com) for helpful tmux commands [^3]\
`tmux new-session -s my-session-name` multiple people can log in to the session using `tmux a -t my-session-name`
6. activate defrabb environment\
`conda activate defrabb`
7. start run w/ tibanna\
`sh etc/run_on_tibanna_giab.sh`

### Monitoring run on AWS instance
1. Switch from tmux session running script to new tmux session[^3]\
`ctrl + b n` (n=next window) this switches between session windows and can be used to get back to your session.  If timeout occured and you need to log back in to EC2 you can re-attach session after login using `tmux a -t my-session-name` where `my-session-name` is the name of the session you started the run in.
2. Viewing tibanna logs for each job.  This is helpful because job information in the job seesion window goes away as new jobs are started. This will allow you to see if jobs have completed or failed and hopefully gives you information on why they failed. 
	- list the jobs, in this example 10\
`tibanna stat -s tibanna_unicorn_giab_test3 -n 10`\
	- job IDs will be on far left, copy jobID\
`tibanna log -s tibanna_unicorn_giab_test3 -j <paste jobID#>`
3. To view snakemake output while the pipeline is running view snakemake log.  
`cd less .snakemake/log/`\
`less less YYYY-MM-DD<somenumber>.snakemake.log`
4. You can go back to session running script using `ctrl + b n`
5. Monitor memory usage with Tibanna `plot_metrics`.  Note: need newer version of Tibanna than what is available in conda pkg that is installed in defrabb env.  
	a. start local terminal and run the following\
	b. `mamba install pip`\
	c. `pip install -U tibanna`\
	d. `tibanna --version` (currently using tibanna 1.9.1 for this)\
	e. `tibanna plot_metrics -j <insert tibanna jobID#>` this will pull up and html page in browser with metrics

### Notes / Gotchas
- Tibanna and all jobs started with tibanna are run in docker containters.  There is a limit on how many docker containers can be started in a given time. For anonymous login (which we are currently using), limit is 100 pull/6hr, see [Docker pull limit documenation](https://www.docker.com/blog/checking-your-current-docker-pull-rate-limits-and-status/) for mor info on this.  
- If you are downloading numerous large files and run into issues with downloads failing the instance might be out of storage and you might have to increase size of storage through EC2 console after stopping instance. 
- Had a job that failed with no info in tibanna log. The issue was the spot instance was shut down and there wasn't another "reasonable" instance to use.  Users can control how Tibanna uses instances on AWS by adjusting the snakemake command in `etc/run_on_tibanna_giab.sh`. Simply restarting run if it fails migth resolve situation. 
- If hap.py or dipcall fails its worth checking tibanna plot_metrics to make sure they are being given enough memory.  If not adjustments can be made in the `resources.yml`


# Development and Testing
A small dataset with chr 21 for GRCh38 [^1] was developed for testing and development. 
The resources are either included in the repository or are hosted in a GIAB S3 bucket with appropriate urls included in the `resources.yml`
small example datasets are included / made available for testing framework code. 
Test pipeline for runtime errors using `snakemake --use-conda -j1`

## Unit Testing Framework
Unit tests for individual python functions and snakemake rules are implemented with python unittest and pytest respectively. 
The pytest unit tests were initially generated using snakemakes `--generate-unit-tests` functionality. 
The test scripts were modifed as needed;
removing unnecessary tests, including config directory, modifying commands for appropriate inputs, and limiting the number of test data files for smaller tests.
Additional modifications were made for bam and vcf comparisons, specifically ignoring file headers as the metadata for the test and expected files are not consistent.

### Python Function Unit Tests
The functions need to be in a .py file.
1. Copy `rules/common.smk` to `test/common.py` for running tests.
2. Run tests using `python -m unittest rules/common.py test/unit/config.py`

### Pytest Snakemake Rule Unit Tests
- Tests are run using `pytest .tests`
- Tests assume `GRCh38_chr21.fa` and `GRCh38_chr21.fa.fai` are in `.tests/integration/resources/references`. 
Not including these files in the repository for now to avoid including large data files in repo, therefore these files might need to be downloaded before running tests.

## Testing Process
```
## unittest
cp rules/common.smk test/common.py
python -m unittest test/common.py test/unit/config.py 
rm test/common.py

## pytest
pytest .tests

## Snakemake pipeline
snakemake --use-conda -j 1 --forceall

## Larger snakemake analysis run set
snakemake --use-conda -j 1 --forceall --config analyses=config/analyses_fulltest.tsv
```


<!-- Resources/ Citations -->

# Footnotes
[^1]: Chromosome 13 was included in the test dataset reference as dipcall incorrectly included line breaks in dip.vcf when only chr21 was included. We might want to submit an issue to the dipcall repo about this issue.
[^2]: [Team resource on setting up AWS instance](https://docs.google.com/document/d/1IdAKastyUShjVl_8msWSR-n1ReKuhXpI_fqLcU829QQ/edit)
[^3]: tmux shortcut commands like `ctrl + b s`  means press `b` while holding down `ctrl` then release both `ctrl` and `b` then press `s` alone. 
