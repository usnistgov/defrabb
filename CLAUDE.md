## Project Overview

DeFrABB (Development Framework for Assembly-Based Benchmarks) is a Snakemake pipeline developed by NIST/GIAB for creating reproducible assembly-based small-variant and structural-variant benchmark sets. It orchestrates diploid assembly variant calling, exclusion region processing, benchmark VCF/BED generation, and evaluation against comparison callsets.

## Build and Run

**Quick start:** `snakemake --use-conda --use-apptainer --cores 1` (add `--forceall` to rerun, or `--config analyses=config/analyses_YYYYMMDD_v0.###_*.tsv` to pick a config).

**NIST runs:** `./run_defrabb run -r <RUNID>` (see `./run_defrabb --help` for subcommands: run, report, archive, release, validate).

**Rerunning:** Simply rerun `./run_defrabb run -r <RUNID>` to complete a partial run.

**Unlocking:** If workflow is locked after interruption: `snakemake --unlock --directory <RUNID>`

**Debug mode:** Set `debug: true` in `config/resources.yml` (or pass `--config debug=true`) for verbose logs from gated rules.

## Workflow

**Recommended pattern (clone-into-runid):**
```sh
git clone <repo-url> 20260519_v0.021_HG002
cd 20260519_v0.021_HG002
./run_defrabb run -r 20260519_v0.021_HG002
```

**Legacy pattern:** `./run_defrabb run -r <RUNID> --outdir <path>` (deprecated but still supported).

## Environment

The `snakemake` CLI requires activation: `conda activate snakemake` (or `~/miniforge3/envs/snakemake`). Individual pipeline rules use their own conda environments defined in `envs/`.

## Testing & QA

```sh
pytest .tests          # Run unit and integration tests
snakefmt .            # Format Snakemake files
black scripts/        # Format Python scripts
scripts/run_full_pipeline_test.sh  # Whole-genome regression test (creates clone with dated config)
```

CI runs all three checks on every push.

## Architecture

`Snakefile` is the entry point. It includes rule modules in this order:
1. `rules/common.smk` - Config loading, table parsing; includes `helpers_{ref,varcall,eval,bench}.smk` (split for AI-friendlier reads)
2. `rules/utils.smk` - Indexing, sorting, compression utilities
3. `rules/download_resources.smk` - Fetch assemblies, references, strats
4. `rules/asm-varcall.smk` - Assembly variant calling (dipcall, PAV)
5. `rules/exclusions_{download,self_discrep,apply}.smk` - Exclusion region processing (split by sub-domain)
6. `rules/report.smk` - Statistics and reporting
7. `rules/bench_vcf_{normalize,anno,finalize}.smk` - VCF post-processing and annotation (split by sub-domain)
8. `rules/stratifications_genome_specific.smk` - Genome-specific (complex/overlapping-variant) hap.py stratifications; opt-in via config `genome_specific_strats`
9. `rules/evaluation.smk` - hap.py and Truvari evaluations

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

## Key Files

- **`docs/TODO-pre-development.md`** - Pre-development roadmap and completed items
- **`docs/development-roadmap.md`** - Planned refactoring phases
- **`docs/issues/`** - Known bugs and workarounds (e.g., truvari-refine-v5.4.0-bug.md)
- **`CHANGELOG`** - Release notes

## NIST-Specific Defaults

`run_defrabb` and `config/release.json` contain NIST-specific paths and S3 settings. External contributors should use CLI overrides: `--archive_dir`, `--s3_bucket`, `--s3_path`.

## Known Issues

- **PAV pysam FIPS self-test crash:** On FIPS hosts, PAV's bundled pysam pip wheel triggers `FATAL FIPS SELFTEST FAILURE`. `run_defrabb` auto-binds `/proc/sys/crypto/fips_enabled -> 0` inside apptainer containers to neutralize this. See `docs/issues/run_pav_fips_selftest.md`.
- **Truvari anno trf on PAV:** Large insertions (>100kb) trigger O(n²) edlib alignment causing multi-day stalls. Size-cap at 100kb routes oversized variants around TRF (kept in output, un-annotated). Config: `truvari_anno_max_ins_length` in `_vcf_processing_params`. See `docs/issues/run_pav_run_dipcall_failures.md` section E.
- **VCF merging with pysam:** Use `vcf_out.new_record()` and copy INFO by name (not `entry.translate()`) to avoid BCF INFO tag ID corruption when headers have different field counts. pysam auto-generates END for symbolic alleles with SVLEN—must declare in header.
- **Memory exhaustion:** Pass `--resources mem_mb=<budget>` to enforce per-rule reservations; without it, multiple heavy jobs (dipcall, PAV) can OOM concurrently. `run_defrabb` defaults to 80% of system memory.
- **Truvari refine v5.4.0:** Fails with 2700+ regions (samtools faidx bug). Workaround: use `truvari` (without refine) or downgrade. See `docs/issues/truvari-refine-v5.4.0-bug.md`.
- **Truvari conda FIPS conflicts:** Some conda envs trigger FIPS conflicts. Fixed in Truvari 4.3+. See `docs/truvari-env-debugging-reference.md`.
