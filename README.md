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
The following usage documentation is for running the pipeline locally using mamba on a linux system (requirement for hap.py). 
Please use these as a starting point for running the pipeline locally. Contact Nathan Olson at nolson@nist.gov, with questions or submit an issue, if you are unable to run the pipeline locally.

##  Initial Setup and Test Run
1. clone repository `git clone https://github.com/nate-d-olson/defrabb.git`
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
1. Run pipeline using `./run_defrabb.sh` providing run in using the defined format (i.e. `-r [YYYYMMDD_milestone_brief-id]`) or `-r` along with `-a`. The wrapper script records the mamba runtime environment information and the git repo status and last commit tag

```
Usage: ./run_defrabb.sh [options] [extra arguments passed to snakemake]
Required:
    -r STRING   Analysis RUN ID, please use following naming convention YYYYMMDD_milestone_brief-id

Optional:
    -a FILE     defrabb run analyses table, if not provided assumes at config/analyses_[RUN ID].tsv
    -o DIR      output directory for framework run, pipeline will create a named directory [RUN ID] at defined location, default is "/defrabb_runs/runs_in_progress/", note this is a system specific path.
    -s all|pipe|report|archive|release  Defining which workflow steps are run
                                    all: pipe, report, and archive (default)
                                    pipe: just the snakemake pipeline
                                    report: generating the snakemake run report
                                    archive: generating snakemake archive tarball
                                    release: copy run output to NAS for upload to Google Drive
    -j          number of jobs used by snakemake, default number of system cores
    -n          Run snakemake in dry run mode, only runs pipe step
    -F          Force rerunning all steps, includes downloading resouces
    -k          keep going with independent jobs if one job fails
    -u          unlock snakeamke run directory
```
1. (For NIST internal run documentation) Fill out README with relevant run information - framework repo info - [milestone] tag (with some potential - hopefully minor-differences), who ran the framework and where/ how, justification / reasoning for analyses, JZ notes (what did we learn), use [defrabb run README template](https://docs.google.com/document/d/1yTXP-3OQxXfGl7kIyXWMTac-USMgiMNPhz10GXwBro0/edit?usp=sharing).
1. (For NIST internal run documentation) Add run information to the [defrabb run log spreadsheet](https://docs.google.com/spreadsheets/d/183LuUat1VCJo2dL7fu0LFMOy8CBA5FTo4WyOVsx4U6o/edit?usp=sharing) 



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
2. Run tests using `python -m unittest test/common.py test/unit/config.py`

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
