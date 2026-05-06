## Project Overview

DeFrABB (Development Framework for Assembly-Based Benchmarks) is a Snakemake pipeline developed by NIST/GIAB for creating reproducible assembly-based small-variant and structural-variant benchmark sets. It orchestrates diploid assembly variant calling, exclusion region processing, benchmark VCF/BED generation, and evaluation against comparison callsets.

## Build and Run

Local smoke test: `snakemake --use-conda --use-apptainer --cores 1` (add `--forceall` to rerun, or `--config analyses=config/analyses_YYYYMMDD_v0.###_*.tsv` to pick a config). For NIST-flavored runs use `./run_defrabb --help`. Set `debug: true` in `config/resources.yml` (or pass `--config debug=true`) to enable verbose logs from gated rules (`get_sv_widen_coords`, `run_happy`).

## Architecture

`Snakefile` is the entry point. It includes rule modules in this order:
1. `rules/common.smk` - Config loading, table parsing; includes `helpers_{ref,varcall,eval,bench}.smk` (split for AI-friendlier reads)
2. `rules/utils.smk` - Indexing, sorting, compression utilities
3. `rules/download_resources.smk` - Fetch assemblies, references, strats
4. `rules/asm-varcall.smk` - Assembly variant calling (dipcall, PAV)
5. `rules/exclusions_{download,self_discrep,apply}.smk` - Exclusion region processing (split by sub-domain)
6. `rules/report.smk` - Statistics and reporting
7. `rules/bench_vcf_{normalize,anno,finalize}.smk` - VCF post-processing and annotation (split by sub-domain)
8. `rules/evaluation.smk` - hap.py and Truvari evaluations

Two main config files drive the pipeline:

- **`config/resources.yml`** - Reference genome URLs, assembly URLs, exclusion region URLs, comparison callset URLs, and compute resources (threads/memory). Validated by `schema/resources-schema.yml`.
- **`config/analyses.tsv`** - Defines what the pipeline runs. Each row specifies an evaluation with: assembly ID, reference, variant caller (dipcall/pav), benchmark type (smvar/stvar), VCF processing steps, exclusion set, and evaluation tool (happy/truvari/unhappy). Validated by `schema/analyses-schema.yml`.

## Key Wildcard Patterns

- Assembly calls: `{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.*`
- Benchmarks: `{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`
- Evaluations: `{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`

## Conventions

- Helper functions go in the appropriate `rules/helpers_*.smk` (ref / varcall / eval / bench); scripts in `scripts/`. `rules/common.smk` keeps only schema/table loading + module-level config init, then `include:`s the helpers so they can reference module-level globals (`vc_tbl`, `bench_tbl`, `ref_config`, `asm_config`, `REFIDS`, …).
- File naming: `run_*.py`, `test_*.py`, `analyses_YYYYMMDD_v0.###_*.tsv`

## NIST-Specific Defaults

`run_defrabb` and `config/release.json` contain NIST-specific paths and S3 settings. External contributors should use CLI overrides: `--outdir`, `--archive_dir`, `--s3_bucket`, `--s3_path`.
