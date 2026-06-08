# DeFrABB Pre-Development TODO List

_Prepared: 2026-03-24_  
_Last Updated: 2026-05-19_

20 issues to address before further pipeline development to establish a robust, intuitive, and modular foundation for generating benchmarksets across additional genomes.

---

## Correctness & Reliability

### 3. Resolve truvari_remap.yml mixed conda/pip install — ✅ Done (2026-06-08)

`truvari_remap.yml` previously installed Truvari 4.2.2 via pip while
`truvari.yml` installed 5.3.0 via conda. The pip install bypassed
conda's dependency resolution and was the source of the FIPS OpenSSL
conflict.

**Update (2026-04-16):** Truvari **4.3** fixed the FIPS OpenSSL issue,
so upgrading `envs/truvari_remap.yml` to ≥4.3 via conda resolves both
the pip-vs-conda divergence and the FIPS workaround in one shot. The
user is editing this file directly to switch to the conda install.

**Resolution (2026-06-08):** `envs/truvari_remap.yml` now installs
`bioconda::truvari ==5.4.0` via conda. The only remaining pip dependency is
`bwapy`, which `truvari anno remap` requires for the BwaAligner `-a` backend —
bioconda's bwapy 0.1.4 lacks `-a` and fails to load on modern glibc
(`undefined symbol: __log_finite`). It is installed from a patched fork
(`nate-d-olson/bwapy`, pinned to commit `5b91fe0` / tag `0.1.5-nist`) that drops
the `pkg_resources` import (so it builds under setuptools≥80 with normal build
isolation) and compiles with `-fno-finite-math-only` while rebuilding
`bwa/libbwa.a` (fixing the `__log_finite` load failure). Verified end-to-end:
`truvari anno remap` produces a `REMAP`-annotated VCF in the conda-built env.
Guarded by `.tests/unit/test_truvari_remap_env.py`. See GitLab #200.

The durable follow-up is GitLab #198: get the fix into the bioconda bwapy recipe
so the pip block can be replaced with a plain `bioconda::bwapy` dep.

**File:** `envs/truvari_remap.yml`
**Impact:** FIPS compatibility, version consistency

### 4. Add FIPS OpenSSL workaround systematically

Rather than patching individual shell blocks with `OPENSSL_CONF=/dev/null`, consider adding the workaround to all affected conda environments via an `env_vars` activation script, or document it as a known limitation with a centralized fix.

**Update (2026-04-16):** Largely subsumed by item 3 — once Truvari is on
4.3+ the main offender is fixed. Audit any other env that still triggers
the FIPS issue after the Truvari upgrade.

**Audit (2026-06-08):** Confirmed while resolving item 3 that the
`bioconda::truvari ==5.4.0` envs (`truvari_core`, `truvari_remap`,
`truvari_repmask`) run cleanly on this FIPS host **without** an
`OPENSSL_CONF=/dev/null` guard. The only rule that still needs the guard is
`run_truvari_anno_trf`, which is pinned to truvari `4.3.0` (`truvari_trf.yml`).
So the remaining FIPS exposure is isolated to that one pinned env; no systematic
`env_vars` workaround is required as long as the other envs stay on 5.4.0+.

**Impact:** All rules using conda envs with OpenSSL on FIPS systems

---

## Completed Items

### 6. Add subcommands to run_defrabb — ✅ Done (2026-05-14)

Resolved in commit b78b348. `run_defrabb` now supports subcommands (`run`, `validate`, `report`, `archive`, `release`) while maintaining backward compatibility with the legacy `-s/--steps` flag.

### 10. Add dry-run summary output — ✅ Done (2026-05-14)

Resolved in commit 3e2a1a4. Added `scripts/dryrun_summary.py` that parses `snakemake --dryrun` output and generates human-readable summaries of planned jobs and analyses.

### 11. Implement named VCF processing profiles — ✅ Done (2026-05-13)

Resolved in commit 8ae254c. Added `vcf_processing_profiles` registry in `config/resources.yml` and migrated analyses table to use profile names instead of inline step lists.

### 12. Externalize exclusion slop/merge constants — ✅ Done (2026-05-14)

Resolved in commit 95eb860. Moved slop and merge parameters from hardcoded rule params to `config/resources.yml` with per-exclusion override support.

### 13. Add exclusion provenance output — ✅ Done (2026-05-14)

Resolved in commit 95eb860. Each benchmark run now emits `.exclusion_provenance.yml` with machine-readable record of all exclusions, transforms, parameters, and base-pair impact.

### 14. Separate NIST-specific logic from core pipeline — ✅ Done (2026-05-14)

Resolved in commit 54c806a. Implemented profile system that separates NIST-specific defaults (in `profiles/nist/config.yml`) from generic pipeline configuration.

### 18. Add tests for run_defrabb wrapper — ✅ Done (2026-05-14)

Resolved in commit 0dc5f78. Added comprehensive test coverage for `run_defrabb` including runID validation, subcommand parsing, profile loading, and manifest generation.

### 19. Add negative/error case tests — ✅ Done (2026-05-14)

Resolved in commits 5564253 and 8d92c7e. Added negative test cases for config validation including missing fields, invalid enums, malformed schemas, and edge cases.

### Pipeline Stabilization — ✅ Done (2026-05-18)

Completed full pipeline stabilization including:
- Fixed `write_report_params.py` import error (removed incorrect snakemake import)
- Fixed truvari log path typos in `rules/evaluation.smk`
- Removed `workflow.source_path()` and `workflow.basedir` references throughout codebase for proper Snakemake archive support
- Updated `run_defrabb` to support clone-into-runid workflow pattern
- Documented truvari refine v5.4.0 bug (samtools faidx failure with 2700+ regions) in `docs/issues/truvari-refine-v5.4.0-bug.md`
- Changed truvari_refine evaluations back to truvari as workaround
- Simplified CI analysis table to 6 v5.0q evaluations (GRCh37/GRCh38/CHM13 × smvar/stvar) for streamlined testing (commit cb29418)

All 85 unit tests passing, lint/formatting clean. Ready for full pipeline test run.

---

## Testing & Quality

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

**Completed (2026-05-14):**
- Items 6, 10, 11, 12, 13, 14, 18, 19 — All core infrastructure and testing improvements complete

**Completed (2026-06-08):**
- Item 3 (Truvari env) — `truvari_remap.yml` on conda truvari 5.4.0 + patched bwapy
  fork; `truvari anno remap` verified working. See GitLab #200.
- Item 4 (FIPS) — audited; exposure isolated to the truvari 4.3.0 `trf` env, which
  already carries the `OPENSSL_CONF` guard. No systematic workaround needed.

**Remaining:**
- **Item 20** — Evaluation comparison utility (next priority for parameter optimization workflow)
- Follow-up (GitLab #198) — bioconda bwapy recipe fix so the pip block in
  `truvari_remap.yml` can become a plain `bioconda::bwapy` dep.
