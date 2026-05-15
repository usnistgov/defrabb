# DeFrABB Pre-Development TODO List

_Prepared: 2026-03-24_

20 issues to address before further pipeline development to establish a robust, intuitive, and modular foundation for generating benchmarksets across additional genomes.

---

## Correctness & Reliability

### 3. Resolve truvari_remap.yml mixed conda/pip install — 🚧 In progress (user)

`truvari_remap.yml` previously installed Truvari 4.2.2 via pip while
`truvari.yml` installed 5.3.0 via conda. The pip install bypassed
conda's dependency resolution and was the source of the FIPS OpenSSL
conflict.

**Update (2026-04-16):** Truvari **4.3** fixed the FIPS OpenSSL issue,
so upgrading `envs/truvari_remap.yml` to ≥4.3 via conda resolves both
the pip-vs-conda divergence and the FIPS workaround in one shot. The
user is editing this file directly to switch to the conda install.

**File:** `envs/truvari_remap.yml`
**Impact:** FIPS compatibility, version consistency

### 4. Add FIPS OpenSSL workaround systematically

Rather than patching individual shell blocks with `OPENSSL_CONF=/dev/null`, consider adding the workaround to all affected conda environments via an `env_vars` activation script, or document it as a known limitation with a centralized fix.

**Update (2026-04-16):** Largely subsumed by item 3 — once Truvari is on
4.3+ the main offender is fixed. Audit any other env that still triggers
the FIPS issue after the Truvari upgrade.

**Impact:** All rules using conda envs with OpenSSL on FIPS systems

---

## CLI & Operator Experience

### 6. Add subcommands to run_defrabb

Currently `run_defrabb` uses `-s/--steps` with values like `all|pipe|report|archive|release`. Refactor into proper subcommands: `run_defrabb run`, `run_defrabb report`, `run_defrabb package`, `run_defrabb release`. This makes the CLI more discoverable and allows subcommand-specific arguments.

**File:** `run_defrabb`
**Impact:** Usability, maintainability

---

## Modularity & Extensibility

### 11. Implement named VCF processing profiles (Phase 1 from roadmap) — ✅ Done

Resolved 2026-05-13. Added `vcf_processing_profiles` registry in `config/resources.yml`
mapping profile names (e.g., `trf`, `trf_sv_lcr`, `norm_xy_trf_sv_lcr`) to ordered lists
of processing steps. Updated `get_processed_vcf()` in `rules/helpers_bench.smk` to resolve
profile names and join steps with dots before inserting into file paths. Extended
`validate_cross_references()` to check profile names exist. Added schema validation in
`schema/resources-schema.yml` following the `exclusion_set` pattern (patternProperties,
no step-name enum). Migrated active `config/analyses.tsv` to use profile names; archived
TSV files unchanged. Unit tests in `.tests/unit/test_vcf_processing_profiles.py`. Profiles
using `remap`/`repmask` steps deferred until Truvari envs are fixed (TODO #3).

**Files:** `config/resources.yml`, `schema/resources-schema.yml`, `schema/analyses-schema.yml`,
`rules/helpers_bench.smk`, `scripts/validate_configs.py`, `config/analyses.tsv`,
`.tests/unit/test_vcf_processing_profiles.py`
**Impact:** Config clarity, new genome onboarding, parameter optimization
**Related issues:** GitLab #177 (programmatic id definition) — folded into this scope on 2026-04-27 when stale draft MR !138 was closed.

### 12. Externalize exclusion slop/merge constants

~~Slop (15kb) and merge (10kb) values are hardcoded in rule params. Move them to resources.yml under the exclusion_set definition so they can be tuned per-genome during parameter optimization.~~

**Resolved 2026-05-14.** Global defaults (`slop: 15000`, `slopmerge_dist: 10000`) moved to `_exclusion_params` in `config/resources.yml`. Per-exclusion overrides (including `slop_pct` for percentage-based expansion via `bedtools slop -pct`) supported via `_exclusion_params.overrides`. Resolver functions `get_slop_value`, `get_slop_flags`, `get_merge_dist` in `rules/helpers_bench.smk`; `add_slop` and `add_slop_and_merge` rules updated to use them.
**Files:** `rules/exclusions_apply.smk`, `rules/helpers_bench.smk`, `config/resources.yml`, `schema/resources-schema.yml`

### 13. Add exclusion provenance output

~~Each run should emit a machine-readable exclusion provenance file listing every exclusion applied, its source, transform parameters, and base-pair impact. This is essential for comparing exclusion strategies across genomes.~~

**Resolved 2026-05-14.** `write_exclusion_provenance` rule added to `rules/exclusions_apply.smk` using `script:` directive. `scripts/write_exclusion_provenance.py` emits a per-benchmark YAML with `exclusion_set`, `global_defaults`, and per-exclusion entries (name, source_type, transform, resolved slop/merge, asm_intersect, bp_impact, pct_of_initial). Provenance files wired into `write_report_params` inputs. Unit tests in `.tests/unit/test_exclusion_params.py`.
**Files:** `rules/exclusions_apply.smk`, `scripts/write_exclusion_provenance.py`, `rules/report.smk`

### 14. Separate NIST-specific logic from core pipeline — ✅ Done

~~`run_defrabb`, `config/release.json`, and parts of resources.yml contain NIST-specific paths and policies. Isolate these into a `profiles/nist/` directory or similar mechanism so the core pipeline remains clean for external use and internal NIST defaults are preserved.~~

**Resolved 2026-05-15.** Changed hardcoded NIST paths in `run_defrabb` to generic defaults (`./defrabb_runs/`, `./defrabb_archive/`). Created `profiles/nist/` directory structure with `config.json` (containing NIST-specific defaults for outdir, archive_dir, release_config, release_type) and `release.json` (NIST S3 bucket configuration). Replaced `config/release.json` with generic template pointing users to profile system. Added `--profile` argument to both new subcommand parser and legacy parser. Implemented `load_profile_config()` and `apply_profile_defaults()` functions that apply profile defaults only when args not explicitly set by user. Profile can be specified via `--profile` flag or `DEFRABB_PROFILE` environment variable. Tests in `.tests/unit/test_run_defrabb.py` (ProfileLoadingTests class, 3 tests).

**Files:** `run_defrabb`, `config/release.json`, `profiles/nist/config.json`, `profiles/nist/release.json`, `.tests/unit/test_run_defrabb.py`
**Impact:** Maintainability, open-source clarity

---

## Testing & Quality

### 18. Add tests for run_defrabb wrapper — ✅ Done

~~The wrapper script previously had no test coverage.~~ Comprehensive test coverage added.

- ✅ Run ID format validation — landed 2026-05-11 (commit `afd3316`).
  `.tests/unit/test_run_defrabb.py` covers the contract enforced by
  `validate_and_set_defaults`: valid IDs, missing/empty runid, and
  format violations (short date, missing version, two-digit major,
  wrong minor width, missing trailing `_`, alphabetic date).
- ✅ Config path resolution — covered by the
  `test_default_analyses_path_uses_runid` case in the same file
  (default `config/analyses_{RUNID}.tsv` derivation when `--analyses`
  is omitted).
- ✅ S3 upload include/exclude logic — landed 2026-05-15 (commit `90f780d`).
  `ReleasePatternTests` covers pattern matching (exact, wildcard, RUNID
  substitution) and include/exclude precedence. `ReleaseRulesExpansionTests`
  covers "same as local/s3" expansion and validation.
- ✅ Manifest generation (`create_data_manifest`) — landed 2026-05-15 (commit `90f780d`).
  `ManifestGenerationTests` covers TSV structure, file type detection, reference
  ID extraction, and S3 pagination handling.

**Files:** `.tests/unit/test_run_defrabb.py` (26 tests, all passing)
**Impact:** Release reliability

### 19. Add negative/error case tests — ✅ Done

~~Current tests only validate happy paths.~~ Comprehensive negative/error case testing added 2026-05-15 (commit `5564253`).

- ✅ Invalid analyses.tsv entries (missing columns, invalid enum values) — 
  `TestAnalysesSchemaValidation` covers missing required columns (eval_cmd),
  invalid enums (eval_cmd, bench_type, vc_cmd), and invalid ref patterns.
- ✅ Missing assembly/reference in resources.yml — covered by existing
  cross-reference tests (`test_missing_asm_id_raises`, `test_missing_ref_raises`).
- ✅ Invalid exclusion_set names — covered by `test_missing_exclusion_set_raises`.
- ✅ Missing VCF processing profiles — covered by `test_missing_vcf_processing_profile_raises`.
- ✅ Malformed resources.yml structure — `TestResourcesSchemaValidation` covers
  missing required sections and type violations (e.g., exclusion_set as string
  instead of list).

**Files:** `.tests/unit/test_config_validation.py` (18 tests, all passing)
**Impact:** Robustness, new genome onboarding, clear error messages for config issues

### 20. Create an evaluation comparison utility

For parameter optimization across genomes, build a script or Snakemake rule that:
- Reads multiple hap.py summary.csv or truvari summary.json files
- Produces a comparison table (precision, recall, F1 by variant type)
- Highlights regressions or improvements
- Optionally generates plots

This directly supports the planned workflow of running multiple parameter configurations and choosing the best one.

**File:** `scripts/compare_evaluations.py` (new)
**Impact:** Parameter optimization workflow, evaluation automation

---

## Priority Order

Item 3 (Truvari) in progress by user. Item 4 effectively closed by item 3.
Also added (not in original list): pre-commit hook config (`.pre-commit-config.yaml`).

Remaining items, in recommended order:

1. **Items 11, 12, 13** — Processing profiles and exclusion externalization (enables parameter optimization)
2. **Items 18, 19** — Testing and quality (prevents regressions during optimization)
4. **Items 6, 14** — CLI and modularity improvements (improves operator velocity)
5. **Item 20** — Evaluation comparison utility (enables systematic optimization)
