# Development Framework for Assembly-Based Benchmarks (DeFrABB)
<!--gitlab badges-->
[![pipeline status](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/badges/master/pipeline.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/commits/master)
[![coverage report](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/badges/master/coverage.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/commits/master)
[![Latest Release](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/badges/release.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/releases)

<!-- Background -->
The GIAB benchmark set development framework is a snakemake bioinformatic pipeline for the development of transparent and reproducible genome assembly based small and structural variant benchmark sets.

## Framework Diagrams

Detailed diagram by Jenny <https://lucid.app/lucidchart/a1feb68c-838b-4851-8259-8289d8cd5c53/edit?invitationId=inv_977874a3-d753-4518-869d-3f0a8ca5eb2c&page=0_0#>
High-level diagram by Nate -<https://lucid.app/lucidchart/aea8aae0-c550-420d-80df-95b73c0cc840/edit?page=0_0#>

<!-- Usage -->
## Usage

The following usage documentation is for running the pipeline locally using mamba on a linux system (requirement for hap.py).
Please use these instructions as a starting point for running the pipeline locally.
Contact Nathan Olson at <nolson@nist.gov>, with questions or submit an issue, if you are unable to run the pipeline locally,
or if you have other questions about the pipeline.
This pipeline was developed and maintained primarily for use in generating benchmark sets for the GIAB RMs by the NIST-GIAB team.
The code is provided for transparency in the benchmark set development process.

### Initial Setup and Test Run

1. clone repository `git clone https://github.com/usnistgov/giab-defrabb.git`
2. generate conda snakemake environment `mamba create -n snakemake --file envs/env.yml`
3. Activate environment `mamba activate snakemake`
4. Run built in test analyses `snakemake --use-conda -j1`

### Running a defined analysis set

1. Define analysis runs by editing `config/analyses.tsv`
2. Update `config/resources.yml` as necessary.
3. Run pipeline using `snakemake --use-conda -j1`

### Executing and Archiving DeFrABB Analysis Runs

1. Fill out `config/analyses_[YYYYMMDD_milestone_brief-id]` and update `config/resources.yml` if necessary.
1. Run pipeline using `./run_defrabb` providing run id using the defined format (i.e. `-r [YYYYMMDD_milestone_brief-id]`) or `-r` along with `-a`. The wrapper script records the mamba runtime environment information and the git repo status and last commit tag

```sh
usage: run_defrabb [-h] -r RUNID [-a ANALYSES] [-o OUTDIR] [-j JOBS] [--archive_dir ARCHIVE_DIR]
                   [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [-s]

DeFrABB wrapper script for executing and archiving framework

options:
  -h, --help            show this help message and exit
  -r RUNID, --runid RUNID
                        Analysis RUN ID, please use following naming convention YYYYMMDD_milestone_brief-id
  -a ANALYSES, --analyses ANALYSES
                        defrabb run analyses table
  -o OUTDIR, --outdir OUTDIR
                        Output directory
  -j JOBS, --jobs JOBS  Number of jobs used by snakemake
  --archive_dir ARCHIVE_DIR
                        Directory to copy pipeline run output to for release. Primarily intended for internal NIST use.
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

workflow steps:
  -s , --steps          Defining which workflow steps are run:
                            all: pipe, report, and archive (default)
                            pipe: just the snakemake pipeline
                            report: generating the snakemake run report
                            archive: generating snakemake archive tarball
                            release: copy run output to NAS for upload to Google Drive (internal NIST use-case)

Any additional arguments provided will be passed directly to Snakemake.
```

__For external users__ The default `output` and `release` directories in the `run_defrabb` 
script are configured specifically for internal use.
The `output` directory can be provided as a command line argument.
The  `release` directory is hard coded to copy files to the NIST-GAIB team NAS. 
You will need to modify the output directory path to a path that is appropriate for your setup.


1. (For NIST internal run documentation) Fill out README with relevant run information - framework repo info - [milestone] tag (with some potential - hopefully minor-differences), who ran the framework and where/ how, justification / reasoning for analyses, JZ notes (what did we learn), use [defrabb run README template](https://docs.google.com/document/d/1yTXP-3OQxXfGl7kIyXWMTac-USMgiMNPhz10GXwBro0/edit?usp=sharing).
1. (For NIST internal run documentation) Add run information to the [defrabb run log spreadsheet](https://docs.google.com/spreadsheets/d/183LuUat1VCJo2dL7fu0LFMOy8CBA5FTo4WyOVsx4U6o/edit?usp=sharing)

## DeFrABB Output

Output directory structure

```txt
.
├── [RUN ID].archive.tar.gz - pipeline archive with code, dependencies, and input generated by snakemake
├── [RUN ID].report.zip - snakemake run report, include interactive html with run information and some results
├── analyses_[RUN ID].tsv - config table defining snakemake pipeline run
├── resources.yml - config file with urls for external inputs and rule parameters
├── benchmark - rule run time along with cpu and memory usage information
├── logs - log file for individual rules
├── defrabb_environment.yml - conda env file with high level environment used to run pipeline
├── resources - externally sourced input files
│   ├── assemblies - assemblies used for benchmark generation
│   ├── comparison_variant_callsets - comparison callsets and regions use in benchmark evaluations
│   ├── exclusions - bed files used to define benchmark exclusion regions
│   ├── references - reference genomes assemblies compared to
│   └── strats - GIAB stratifications used to stratify evaluation results
├── results
│   ├── asm_varcalls - reference v. assembly comparisons including vcf annotations
│   ├── draft_benchmarksets - benchmark regions and variants along with intermediate files
│   ├── evaluations - benchmarking draft benchmark sets to comparison callsets
│   └── report - summary stats and metric used in snakemake report
└── run.log - run_defrabb.sh log file
```

## Development and Testing

A small dataset with chr 21 for GRCh38 [^1] was developed for testing and development.
The resources are either included in the repository or are hosted in a GIAB S3 bucket with appropriate urls included in the `resources.yml`
small example datasets are included / made available for testing framework code.
Test pipeline for runtime errors using `snakemake --use-conda -j1`

## Unit Testing Framework

Unit tests for individual python functions and snakemake rules are implemented with python unittest and pytest respectively.
The  pytest unit tests were initially generated using snakemake `--generate-unit-tests` functionality.
The test scripts were modified as needed;
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

```sh
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

## Footnotes

[^1]: Chromosome 13 was included in the test dataset reference as dipcall incorrectly included line breaks in dip.vcf when only chr21 was included. We might want to submit an issue to the dipcall repo about this issue.
