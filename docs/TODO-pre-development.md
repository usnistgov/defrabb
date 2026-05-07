# DeFrABB Pre-Development TODO List

_Prepared: 2026-03-24_

20 issues to address before further pipeline development to establish a robust, intuitive, and modular foundation for generating benchmarksets across additional genomes.

---

## Correctness & Reliability

### 1. Fix GitHub Actions CI (stale paths) — ✅ Done

Resolved by deleting `.github/workflows/` entirely. The workflows were a
leftover from the snakemake template; internal GitLab CI is the canonical
pipeline.

**File:** `.github/workflows/main.yml` (removed)
**Impact:** CI/CD, contributor experience

### 2. Reconcile conda environment version inconsistencies — ✅ Done

Resolved across commits `b64af75`–`e4a1335` (2026-04-16):
- `bcftools`: aligned on 1.17 (truvari envs deferred to item 3)
- `rtg-tools`: aligned on 3.12.1
- `pybedtools`: aligned on 0.10.0
- `python`: sv_widen_coords bumped 3.8 → 3.12

Truvari version split is being addressed by user separately (see item 3).

**Files:** `envs/bcftools.yml`, `envs/rtgtools.yml`, `envs/bedtools.yml`, `envs/bcftools_and_bedtools.yml`, `envs/sv_widen_coords.yml`
**Impact:** Reproducibility, correctness

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

### 5. Upgrade sv_widen_coords.yml from Python 3.8 — ✅ Done

Resolved in commit `b64af75` (2026-04-16). Python 3.8 → 3.12,
pybedtools pinned to 0.10.0 to match other envs. Script
(`scripts/get_sv_widen_coords.py`) only uses argparse/gzip/logging/sys/
pybedtools — all available in 3.12.

**File:** `envs/sv_widen_coords.yml`
**Impact:** Security, compatibility

---

## CLI & Operator Experience

### 6. Add subcommands to run_defrabb

Currently `run_defrabb` uses `-s/--steps` with values like `all|pipe|report|archive|release`. Refactor into proper subcommands: `run_defrabb run`, `run_defrabb report`, `run_defrabb package`, `run_defrabb release`. This makes the CLI more discoverable and allows subcommand-specific arguments.

**File:** `run_defrabb`
**Impact:** Usability, maintainability

### 7. Fix run_defrabb bucket hardcoding — ✅ Done

Fixed in commit `55737a3` ("stabilize pav release and docs"). The
release flow now respects configured release inputs instead of a
hardcoded bucket name.

**File:** `run_defrabb`
**Impact:** External contributor usability

### 8. Add a `--validate-only` flag to run_defrabb — ✅ Done

Resolved 2026-04-17. `run_defrabb --validate-only` invokes
`snakemake --list-target-rules --quiet` so the workflow loads and
`validate_cross_references` runs (see item 9), without executing rules.

**File:** `run_defrabb`
**Impact:** Parameter optimization workflow, error prevention

### 9. Improve error messages for missing config entries — ✅ Done

Resolved 2026-04-17. `scripts/validate_configs.py` now provides
`validate_cross_references`, called from `rules/common.smk` after the
existing schema validation. Reports missing asm_id, ref, eval_comp_id,
and exclusion_set with grouped errors that name the offending IDs and
the eval_ids that reference them.

**File:** `rules/common.smk`
**Impact:** New genome onboarding, debugging time

### 10. Add dry-run summary output

After `snakemake --dryrun`, produce a human-readable summary showing: number of rules to execute, which assemblies and references are involved, estimated disk/memory requirements. This helps operators plan runs for new genomes.

**Impact:** Operator experience, resource planning

---

## Modularity & Extensibility

### 11. Implement named VCF processing profiles (Phase 1 from roadmap)

Replace the dot-separated `bench_vcf_processing` strings with a named profile registry in resources.yml. This is the single highest-impact refactor for enabling new genome work, because parameter optimization requires trying different processing combinations safely.

**Files:** `config/resources.yml`, `rules/common.smk`, `rules/bench_vcf_processing.smk`
**Impact:** Config clarity, new genome onboarding, parameter optimization
**Related issues:** GitLab #177 (programmatic id definition) — folded into this scope on 2026-04-27 when stale draft MR !138 was closed.

### 12. Externalize exclusion slop/merge constants

Slop (15kb) and merge (10kb) values are hardcoded in rule params. Move them to resources.yml under the exclusion_set definition so they can be tuned per-genome during parameter optimization.

**Files:** `rules/exclusions.smk`, `config/resources.yml`
**Impact:** Parameter optimization, exclusion strategy experimentation

### 13. Add exclusion provenance output

Each run should emit a machine-readable exclusion provenance file listing every exclusion applied, its source, transform parameters, and base-pair impact. This is essential for comparing exclusion strategies across genomes.

**Files:** `rules/exclusions.smk`, `scripts/subtract_exclusions.py`
**Impact:** Reproducibility, parameter optimization

### 14. Separate NIST-specific logic from core pipeline

`run_defrabb`, `config/release.json`, and parts of resources.yml contain NIST-specific paths and policies. Isolate these into a `profiles/nist/` directory or similar mechanism so the core pipeline remains clean for external use and internal NIST defaults are preserved.

**Files:** `run_defrabb`, `config/release.json`
**Impact:** Maintainability, open-source clarity

### 15. Unify report source paths — ✅ Done

Resolved in commit `55737a3` ("stabilize pav release and docs"). The
root `analysis.qmd` is now the canonical report source; the duplicate
under `scripts/reports/analysis.qmd` was removed.

**Files:** `analysis.qmd` (root), `rules/report.smk`
**Impact:** Maintainability, contributor confusion

---

## Testing & Quality

### 16. Add pytest to CI — ✅ Done

Resolved in commit `c995786` (2026-04-16). Added a `pytest` job to
`.gitlab-ci.yml` that installs pytest into the existing mamba env and
runs `pytest .tests`. As part of the same commit, the 29 broken
Snakemake-7-era auto-generated unit tests were removed (they used the
deprecated `--keep-target-files` flag and had been silently broken
since the Snakemake 8 migration). Only `test_stabilization_pr.py`
(5 working tests) remains.

**Files:** `.gitlab-ci.yml`, `.tests/unit/`
**Impact:** Regression prevention

### 17. Enable skipped config validation test — ✅ Resolved (deleted)

`.tests/unit/test_eval_config.skip` was deleted as part of commit
`c995786` along with the other broken Snakemake-7-era tests. The test
referenced hardcoded `eval1`/`eval6` IDs that no longer match the
current `analyses.tsv`, and reusing the pattern was not worthwhile.
Future config-validation coverage should be written as plain Python
unit tests against `rules/common.smk` helpers (see item 19).

**File:** `.tests/unit/test_eval_config.skip` (removed)
**Impact:** Config safety

### 18. Add tests for run_defrabb wrapper

The wrapper script has no test coverage. Add tests for:
- Run ID format validation
- Config path resolution
- S3 upload include/exclude logic
- Manifest generation

**Files:** `.tests/unit/test_run_defrabb.py` (new)
**Impact:** Release reliability

### 19. Add negative/error case tests

Current tests only validate happy paths. Add tests for:
- Invalid analyses.tsv entries (missing columns, invalid enum values)
- Missing assembly/reference in resources.yml
- Invalid exclusion_set names
- Malformed VCF inputs

**Files:** `.tests/unit/test_config_validation.py` (new)
**Impact:** Robustness, new genome onboarding

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

## Follow-ups from config validation (2026-04-17)

Flagged during the final code review of the config-validation feature
(commits `88473b5..a74ffdf`). Deferred deliberately — the feature
shipped clean enough to be useful as-is.

### 21. Stop `--validate-only` output from flooding streams — ✅ Done

Observed during validation-only smoke testing: Snakemake stringified
bound `functools.partial(...)` input functions in `rules/evaluation.smk`,
including the full analyses DataFrame and `resources.yml` dict. The
resulting warning emitted ~160 KB of workflow-load output and could break
API streaming with `response.failed` before the command completed.

Resolved by replacing those bound partials with lightweight wrapper
functions, capturing `validate_only()` subprocess output, and only
emitting captured output on failure.

**Files:** `rules/evaluation.smk`, `rules/common.smk`, `run_defrabb`
**Impact:** Operator UX for the validation flow

### 22. Defensive top-level config access in validator — ✅ Won't fix (2026-04-17)

The reviewer was concerned about bare `KeyError` if a future refactor
removed a top-level section. After re-checking, `schema/resources-schema.yml`
already marks `references`, `assemblies`, `comparisons`, and `exclusion_set`
as `required`, and `rules/common.smk:438` runs schema validation **before**
`validate_cross_references` (line 468). Adding defensive `.get()` would
guard a scenario the schema already prevents.

Resolved by (a) documenting the schema-first contract in
`validate_cross_references`'s docstring, and (b) dropping the now-redundant
outer `.get()` from the comparisons branch for consistency. If a future
refactor relaxes the schema, that change should re-introduce defensiveness
explicitly rather than relying on a `.get()` to mask the drift silently.

**File:** `scripts/validate_configs.py`
**Impact:** Code clarity, schema/validator contract

### 23. `--validate-only` should not require `--runid` — ✅ Done

Resolved 2026-05-05. `-r/--runid` is now optional at argument parsing
time and enforced only for non-validation runs. `./run_defrabb
--validate-only` defaults to `config/analyses.tsv`; passing both
`--validate-only` and `-r` still supports the versioned
`config/analyses_{RUNID}.tsv` default.

**File:** `run_defrabb`
**Impact:** Validation workflow ergonomics

---

## Priority Order

Foundation items (1, 2, 5, 7, 15, 16, 17) completed 2026-03–04.
Item 3 (Truvari) in progress by user. Item 4 effectively closed by item 3.
Also added (not in original list): pre-commit hook config (`.pre-commit-config.yaml`).

Remaining items, in recommended order:

1. ~~Items 8, 9~~ — Config validation completed 2026-04-17.
2. **Items 11, 12, 13** — Processing profiles and exclusion externalization (enables parameter optimization)
3. **Items 18, 19** — Testing and quality (prevents regressions during optimization)
4. **Items 6, 10, 14** — CLI and modularity improvements (improves operator velocity)
5. **Item 20** — Evaluation comparison utility (enables systematic optimization)
6. Items 21, 23 completed 2026-05-05. Item 22 closed as won't-fix on 2026-04-17.
