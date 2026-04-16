# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DeFrABB (Development Framework for Assembly-Based Benchmarks) is a Snakemake pipeline developed by NIST/GIAB for creating reproducible assembly-based small-variant and structural-variant benchmark sets. It orchestrates diploid assembly variant calling, exclusion region processing, benchmark VCF/BED generation, and evaluation against comparison callsets.

## Build and Run Commands

```sh
# Run default analysis (local smoke test)
snakemake --use-conda --use-apptainer --cores 1

# Full workflow smoke test (force rerun all)
snakemake --use-conda --use-apptainer --cores 1 --forceall

# Run with a specific analyses config
snakemake --use-conda --use-apptainer --cores 1 --config analyses=config/analyses_YYYYMMDD_v0.###_description.tsv

# Run via project wrapper
./run_defrabb -r 20240513_v0.016_test -a config/analyses.tsv -s pipe

# Run tests
pytest .tests
pytest .tests/unit/test_specific_feature.py

# Code formatting
black . --check          # verify Python formatting
snakefmt .               # format Snakemake files
```

## Environment Requirements

- Linux only
- Snakemake >= 8.30.0 in base environment
- Conda/Mamba for `--use-conda` (per-rule envs created automatically from `envs/*.yml`)
- Apptainer for PAV rules (uses Docker container `becklab/pav`)
- `boto3` in active environment if using `run_defrabb` wrapper

## Architecture

### Entry Point and Rule Composition

`Snakefile` is the entry point. It includes rule modules in this order:
1. `rules/common.smk` - Config loading, helper functions, table parsing
2. `rules/utils.smk` - Indexing, sorting, compression utilities
3. `rules/download_resources.smk` - Fetch assemblies, references, strats
4. `rules/asm-varcall.smk` - Assembly variant calling (dipcall, PAV)
5. `rules/exclusions.smk` - Exclusion region processing
6. `rules/report.smk` - Statistics and reporting
7. `rules/bench_vcf_processing.smk` - VCF post-processing and annotation
8. `rules/evaluation.smk` - hap.py and Truvari evaluations

### Configuration System

Two main config files drive the pipeline:

- **`config/resources.yml`** - Reference genome URLs, assembly URLs, exclusion region URLs, comparison callset URLs, and compute resources (threads/memory). Validated by `schema/resources-schema.yml`.
- **`config/analyses.tsv`** - Defines what the pipeline runs. Each row specifies an evaluation with: assembly ID, reference, variant caller (dipcall/pav), benchmark type (smvar/stvar), VCF processing steps, exclusion set, and evaluation tool (happy/truvari/unhappy). Validated by `schema/analyses-schema.yml`.

`rules/common.smk` parses the analyses TSV into filtered DataFrames (`dipcall_tbl`, `pav_tbl`, `bench_tbl`, `happy_analyses`, `truvari_analyses`, etc.) that drive `expand()` calls in `rule all`.

### Data Flow

1. **Download** references, assemblies, stratifications, comparison callsets
2. **Index** references (fai, bwa, mmi, sdf) for different tools
3. **Variant call** using dipcall or PAV against reference
4. **Post-process** VCFs (normalize, split multiallelics, fix XY genotypes, annotate TRs)
5. **Process exclusions** (download, apply slop/merge, subtract from callable regions)
6. **Generate** draft benchmark VCF + BED
7. **Evaluate** with hap.py (small variants) or Truvari (structural variants)
8. **Report** statistics (assembly stats, bcftools stats, BED summaries)

### Key Wildcard Patterns

File paths follow structured wildcard patterns:
- Assembly calls: `{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.*`
- Benchmarks: `{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`
- Evaluations: `{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`

### Output Structure

```
results/
├── asm_varcalls/{vc_id}/        # Raw variant calling outputs
├── draft_benchmarksets/{bench_id}/  # Processed benchmark VCF/BED
├── evaluations/{happy,truvari}/     # Evaluation results
├── report/                          # Assembly and stats reports
└── analysis_params.yml              # Serialized analysis config
```

## Coding Conventions

- Python: Black formatting, 4-space indentation, `snake_case`
- Type hints and short docstrings for reusable Python functions
- One logical transformation per Snakemake rule
- Helper functions go in `rules/common.smk`; scripts in `scripts/`
- Each rule uses its own conda env from `envs/*.yml`
- Logs: `logs/{rulename}/{wildcards}.log`; benchmarks: `benchmark/{rulename}/{wildcards}.tsv`
- File naming: `run_*.py`, `test_*.py`, `analyses_YYYYMMDD_v0.###_*.tsv`

## Testing

- Unit tests in `.tests/unit/test_*.py` using pytest
- Tests expect chr21 reference files under `.tests/integration/resources/references/` (already tracked)
- New rule tests: put fixtures under `.tests/unit/<rule>/data` with expected outputs in `expected/`
- If a change touches core workflow wiring, run the full smoke test: `snakemake --use-conda --use-apptainer -j1 --forceall`

## Commit and PR Conventions

- Short, lowercase commit subjects describing the change (e.g., `fixing path`, `adding pav exclusion`)
- PRs should summarize workflow impact, note config/resource changes, and mention validation performed
- Release PRs: update CHANGELOG, use milestone title format `milestone v#.###`

## NIST-Specific Defaults

`run_defrabb` and `config/release.json` contain NIST-specific paths and S3 settings. External contributors should use CLI overrides: `--outdir`, `--archive_dir`, `--s3_bucket`, `--s3_path`.
