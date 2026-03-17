# Development Framework for Assembly-Based Benchmarks (DeFrABB)
<!--gitlab badges-->
[![pipeline status](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/badges/master/pipeline.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/commits/master)
[![coverage report](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/badges/master/coverage.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/commits/master)
[![Latest Release](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/badges/release.svg)](https://gitlab.nist.gov/gitlab/bbd-human-genomics/defrabb/-/releases)

DeFrABB is a Snakemake workflow used by the NIST/GIAB team to develop transparent, reproducible assembly-based small- and structural-variant benchmark sets.

## Status and intended audience

This repository is maintained primarily for internal benchmark-development work and is made public for transparency and reproducibility. It is not currently packaged as a broadly supported end-user toolkit. External issues and pull requests are welcome, but support and review are best-effort.

## Repository layout

- `Snakefile` - workflow entry point.
- `rules/` - modular Snakemake rule files for resource downloads, assembly variant calling, exclusions, evaluation, and reporting.
- `scripts/` - Python, R, and shell helpers used by rules.
- `config/` - default analysis tables, resource configuration, and release settings.
- `schema/` - schemas for `config/resources.yml` and analyses tables.
- `envs/` - per-rule Conda environments.
- `.tests/` - pytest-based rule and helper tests.
- `report/` - Snakemake report templates.
- `docs/` - NIST-specific operational notes and release procedures.

## Requirements

DeFrABB is developed for Linux environments. A local base environment should provide:

- Snakemake 8.30 or newer
- Conda or mamba
- Apptainer for rules that rely on containers
- Python dependencies needed by `run_defrabb` (notably `boto3`) if you use the wrapper

There is no single bootstrap environment file in this repository. Instead, install Snakemake in a local base environment and let the workflow create rule-specific software environments with `--use-conda`.

## Running the workflow directly

`config/resources.yml` points to `config/analyses.tsv` by default. For a local smoke run:

```sh
snakemake --use-conda --use-apptainer --cores 1
```

To run a different checked-in analysis table, override the `analyses` config value:

```sh
snakemake --use-conda --use-apptainer --cores 1 \
  --config analyses=config/analyses_20250708_v0.021_HG008TN.tsv
```

Before running new analyses, update `config/analyses.tsv` and `config/resources.yml` as needed. Field requirements are defined in `schema/analyses-schema.yml` and `schema/resources-schema.yml`.

## Running with `run_defrabb`

`run_defrabb` creates `OUTDIR/RUNID`, records git status, exports the active Conda environment to `environment.yml`, and can run the pipeline, generate a Snakemake report, build an archive, and release selected files.

Example:

```sh
./run_defrabb -r 20250708_v0.021_HG008TN -o ../ -s pipe
```

Run IDs must follow `YYYYMMDD_v#.###_brief-id`. If `-a/--analyses` is omitted, the wrapper looks for `config/analyses_<RUNID>.tsv`.

Release defaults in `run_defrabb` and `config/release.json` are NIST-specific. External users should expect to override `--outdir`, `--archive_dir`, `--s3_bucket`, `--s3_path`, and `--release_type`.

## Output structure

The wrapper writes results under `OUTDIR/RUNID/`, typically including:

```txt
OUTDIR/RUNID/
├── archive.tar.gz
├── snakemake_report_<RUNID>.zip
├── environment.yml
├── run.log
├── benchmark/
├── config/
├── logs/
├── resources/
└── results/
```

Generated reports also include `analysis.html` and supporting report inputs under `results/analysis_params.yml` and `results/report/`.

## Testing and development

For lightweight validation:

```sh
pytest .tests
snakemake --use-conda --use-apptainer --cores 1 --forceall
```

Test references for the bundled chr21 dataset are already tracked under `.tests/integration/resources/references`.

Versioned example analyses are kept in `config/analyses_*.tsv` and can be used as templates for new runs.

## Internal-only notes

Some documents in `docs/` and some release defaults in the wrapper assume NIST infrastructure (NAS paths, S3 destinations, and internal run documentation practices). Treat those as internal operational references rather than general public setup instructions.
