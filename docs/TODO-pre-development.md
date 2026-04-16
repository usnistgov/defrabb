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

### 2. Reconcile conda environment version inconsistencies

Multiple environments pin different versions of the same tools:
- `bcftools`: 1.14 (bcftools.yml) vs 1.17 (bcftools_and_bedtools.yml)
- `rtg-tools`: 3.10.1 (rtgtools.yml) vs 3.12.1 (happy.yml)
- `truvari`: 5.3.0 (truvari.yml, conda) vs 4.2.2 (truvari_remap.yml, pip)
- `pybedtools`: 0.8.2 vs 0.10.0
- `python`: 3.8 (sv_widen_coords.yml) vs 3.12+ elsewhere

Decide on canonical versions for each tool and align all envs. The Truvari 4.2.2 vs 5.3.0 split is especially risky — different major versions may produce different annotation outputs.

**Files:** `envs/*.yml`
**Impact:** Reproducibility, correctness

### 3. Resolve truvari_remap.yml mixed conda/pip install — recommended next env fix

`truvari_remap.yml` installs Truvari 4.2.2 via pip while `truvari.yml`
installs 5.3.0 via conda. The pip install bypasses conda's dependency
resolution and was the source of the FIPS OpenSSL conflict.

**Update (2026-04-16):** Truvari **4.3** fixed the FIPS OpenSSL issue,
so upgrading `envs/truvari_remap.yml` to ≥4.3 (preferably via conda)
resolves both the pip-vs-conda divergence and the FIPS workaround in one
shot. This is the highest-leverage env fix to do next.

**File:** `envs/truvari_remap.yml`
**Impact:** FIPS compatibility, version consistency

### 4. Add FIPS OpenSSL workaround systematically

Rather than patching individual shell blocks with `OPENSSL_CONF=/dev/null`, consider adding the workaround to all affected conda environments via an `env_vars` activation script, or document it as a known limitation with a centralized fix.

**Update (2026-04-16):** Largely subsumed by item 3 — once Truvari is on
4.3+ the main offender is fixed. Audit any other env that still triggers
the FIPS issue after the Truvari upgrade.

**Impact:** All rules using conda envs with OpenSSL on FIPS systems

### 5. Upgrade sv_widen_coords.yml from Python 3.8

`envs/sv_widen_coords.yml` pins Python 3.8, which is EOL (Oct 2024). This env is used by `get_sv_widen_coords` (SV exclusion logic). Upgrade to 3.12+ and verify the script still works.

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

### 8. Add a `--validate-only` flag to run_defrabb

Before running the full pipeline on a new genome, users need to validate their analyses TSV and resources.yml configuration. Add a `--validate-only` mode that loads configs, validates schemas, checks that referenced assemblies/references/comparisons exist in resources.yml, and reports issues without executing any rules.

**File:** `run_defrabb`
**Impact:** Parameter optimization workflow, error prevention

### 9. Improve error messages for missing config entries

When an asm_id, ref, or comp_id in analyses.tsv doesn't exist in resources.yml, the pipeline fails deep inside Snakemake with a cryptic KeyError. Add early validation in `rules/common.smk` that produces clear error messages listing exactly which IDs are missing and from which config section.

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

### 16. Add pytest to CI

Neither GitHub Actions nor GitLab CI runs `pytest .tests`. The unit test suite exists but isn't exercised in CI. Add a pytest job to the canonical CI pipeline.

**Files:** `.github/workflows/main.yml` or `.gitlab-ci.yml`
**Impact:** Regression prevention

### 17. Enable skipped config validation test

`.tests/unit/test_eval_config.skip` is explicitly skipped. Either fix the underlying issue and enable it, or document why it's skipped. Config validation is critical for new genome work.

**File:** `.tests/unit/test_eval_config.skip`
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

## Priority Order

For the planned new-genome development work, address in this order
(items 1, 7, 15 already done; 16 superseded by retiring `.github/workflows/`
in favor of GitLab CI):

1. **Item 3** — Upgrade Truvari to ≥4.3 in `envs/truvari_remap.yml`
   (closes 4 in passing; smallest high-value change)
2. **Items 2, 5** — Reconcile remaining conda env version
   inconsistencies (bcftools, rtg-tools, pybedtools, Python 3.8 in
   `sv_widen_coords.yml`)
3. **Items 8, 9** — Config validation (prevents wasted compute on bad configs)
4. **Items 11, 12, 13** — Processing profiles and exclusion externalization (enables parameter optimization)
5. **Items 17, 18, 19** — Testing and quality (prevents regressions during optimization)
6. **Items 6, 10, 14** — CLI and modularity improvements (improves operator velocity)
7. **Item 20** — Evaluation comparison utility (enables systematic optimization)
