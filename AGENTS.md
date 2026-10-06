# Repository Guidelines

Instructions for coding agents (Claude Code, Codex, and others) and human
contributors. `CLAUDE.md` imports this file; edit this file, not `CLAUDE.md`.

## Project Overview

DeFrABB (Development Framework for Assembly-Based Benchmarks) is a Snakemake
pipeline developed by NIST/GIAB for creating reproducible assembly-based
small-variant and structural-variant benchmark sets. It orchestrates diploid
assembly variant calling, exclusion region processing, benchmark VCF/BED
generation, and evaluation against comparison callsets. It produced the GIAB
HG002 v5.0q benchmarks (DeFrABB v0.020; preprint
<https://www.biorxiv.org/content/10.64898/2026.09.23.752440v1>).

## Build and Run

- **Quick start:** `snakemake --use-conda --use-apptainer --cores 1` (add
  `--forceall` to rerun, or `--config analyses=config/analyses_<RUNID>.tsv` to
  pick a config).
- **NIST runs:** `./run_defrabb run -r <RUNID>` (subcommands: run, report,
  archive, release, validate; see `./run_defrabb --help`). Rerun the same
  command to complete a partial run.
- **Unlocking:** `snakemake --unlock --directory <RUNID>`.
- **Debug mode:** `debug: true` in `config/resources.yml` (or
  `--config debug=true`) for verbose logs from gated rules.
- **Memory:** pass `--resources mem_mb=<budget>`; `run_defrabb` defaults to 80%
  of system memory. Without it, dipcall/PAV/hap.py jobs can OOM concurrently.

Recommended run pattern (clone into a directory named by run ID):

```sh
git clone <repo-url> 20260519_v0.023_HG002
cd 20260519_v0.023_HG002
./run_defrabb run -r 20260519_v0.023_HG002
```

`--outdir` (legacy pattern) is deprecated but still supported.

## Environment

The `snakemake` CLI requires activation: `conda activate snakemake` (or
`~/miniforge3/envs/snakemake`). Rules use their own conda environments from
`envs/`; PAV runs in an Apptainer container.

## Testing & QA

```sh
pytest .tests                      # unit and structural tests
snakefmt .                         # format Snakemake files (CI uses 0.10.2)
black scripts/                     # format Python scripts
scripts/run_full_pipeline_test.sh  # whole-genome regression test (clone + dated config)
```

CI runs pytest, snakefmt, and black on every push. Add or update tests in
`.tests/unit/test_<feature>.py` for rule or helper changes; tests are Python
unit and structural tests, not auto-generated rule tests. Some tests expect the
chr21 references under `.tests/integration/resources/references`. If core
workflow wiring changes, also run the chr21 workflow with `--forceall`.

## Architecture

`Snakefile` is the entry point. It includes rule modules in this order:

1. `rules/common.smk`: config loading and table parsing; includes
   `helpers_{ref,varcall,eval,bench}.smk`
2. `rules/utils.smk`: indexing, sorting, compression utilities
3. `rules/download_resources.smk`: fetch assemblies, references, strats
4. `rules/asm-varcall.smk`: assembly variant calling (dipcall; PAV3 via
   `pav3 batch`)
5. `rules/exclusions_{download,self_discrep,pav_discrep,apply}.smk`: exclusion
   region processing
6. `rules/report.smk`: statistics and reporting
7. `rules/bench_vcf_{normalize,anno,finalize}.smk`: VCF post-processing and
   annotation
8. `rules/stratifications_genome_specific.smk`: genome-specific
   (complex-variant) hap.py stratifications; opt-in via `genome_specific_strats`
9. `rules/evaluation.smk`: hap.py and Truvari evaluations

Two config files drive the pipeline:

- **`config/resources.yml`**: reference, assembly, exclusion, comparison
  callset, and stratification URLs; named parameter profiles; compute resources.
  Validated by `schema/resources-schema.yml`.
- **analyses table** (`config/analyses*.tsv`): one row per evaluation (assembly,
  reference, caller and params, `smvar`/`stvar`, VCF processing, exclusion set,
  evaluation tool and comparison). Validated by `schema/analyses-schema.yml`.

Wildcard patterns:

- Assembly calls: `{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.*`
- Benchmarks: `{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`
- Evaluations:
  `{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.*`

## Parameter Optimization (v0.023+)

- **Sweep generator:** `scripts/generate_param_sweep.py <config.yml>` (YAML to
  analyses table; vc_id reuse, cost estimation, dry-run mode)
- **Scoring:**
  `scripts/score_param_sweep.py --results-dir <dir> --baseline v5.0q --top-n 3`
- **Profiles** in `config/resources.yml`: dipcall `z2k` (default), `z5k`,
  `z10k`, `z1k`; PAV `giab` (PAV3 defaults; profiles are dicts of PAV3
  `pav.json` params, no merge params); exclusions `standard` (default),
  `conservative`, `aggressive`
- **Docs:** `docs/parameter-optimization.md`; example configs in
  `config/sweeps/`

## Coding Style & Conventions

- Python: 4-space indentation, Black, type hints, short docstrings for reusable
  functions; `snake_case` for helpers, TSV/YAML fields, and config keys.
- Snakemake: one logical transformation per rule; helper functions in the
  appropriate `rules/helpers_*.smk` (`common.smk` keeps only schema/table
  loading and module-level config, then includes the helpers); tool environments
  in `envs/*.yml`. Keep `ruleorder` on one line (CI snakefmt compatibility).
- File naming: `run_*.py`, `test_*.py`, `analyses_YYYYMMDD_v0.###_*.tsv`.
- Dated analyses tables (`config/analyses_YYYYMMDD_*.tsv`) record past runs; do
  not edit them. See GitLab #207 on moving run configs out of the codebase.
- Commits: conventional commit subjects, e.g. `fix(exclusions): ...`. Update
  `CHANGELOG` for user-visible changes.

## Documentation Placement

Repository docs (`README.md`, `docs/`) are public-facing: mirrored to GitHub and
linked from the v5.0q preprint. Keep them to method, configuration, outputs, and
usage. Do not add session notes, planning docs, investigations, validation
reports, or NIST runbooks there.

- **Planning, roadmap, TODOs:** GitLab issues and milestones (project ID 6652;
  `glab api projects/6652/issues`). Root `TODO.md` is a gitignored scratch file.
- **Investigations, validation reports, NIST runbooks, release and CI notes:**
  GitLab wiki (`git@gitlab.nist.gov:bbd-human-genomics/defrabb.wiki.git`). Code
  comments cite pages as `wiki:investigations/<page>`.
- **Session and design notes:** the maintainer's notes vault.

## NIST-Specific Defaults

`run_defrabb` uses generic defaults (`./defrabb_archive/`, template
`config/release.json`). NIST paths and S3 settings live in `profiles/nist/` and
load with `--profile nist` (or `DEFRABB_PROFILE=nist`). Others should override
`--archive_dir`, `--s3_bucket`, and `--s3_path` instead of editing defaults.

## Known Issues

- **PAV pysam FIPS self-test crash:** PAV's bundled pysam wheel triggers
  `FATAL FIPS SELFTEST FAILURE` on FIPS hosts. `run_defrabb` binds
  `/proc/sys/crypto/fips_enabled -> 0` inside Apptainer.
  (`wiki:investigations/run_pav_fips_selftest`)
- **Truvari anno trf on PAV:** insertions >100 kb cause multi-day edlib stalls;
  they are routed around TRF (kept, un-annotated) via
  `truvari_anno_max_ins_length` in `_vcf_processing_params`.
  (`wiki:investigations/run_pav_run_dipcall_failures`, section E)
- **VCF merging with pysam:** use `vcf_out.new_record()` and copy INFO by name
  (not `entry.translate()`) to avoid INFO tag ID corruption; pysam
  auto-generates END for symbolic alleles with SVLEN, so declare it in the
  header.
- **Truvari refine v5.4.0:** fails with 2700+ regions (samtools faidx bug); use
  `truvari` without refine, or filter to primary chromosomes.
  (`wiki:investigations/truvari-refine-v5.4.0-bug`)
- **Truvari conda FIPS conflicts:** fixed in Truvari 4.3+; only the truvari
  4.3.0 `trf` env still needs `OPENSSL_CONF=/dev/null`.
  (`wiki:operations/truvari-env-debugging-reference`)
