# Config Validation Design (TODO #8 + #9)

**Date:** 2026-04-17
**Status:** Approved for implementation planning
**Scope:** Pre-development TODO items #8 (`--validate-only` flag) and #9 (early cross-reference validation in `rules/common.smk`).

## Problem

When a `config/analyses_*.tsv` row references an `asm_id`, `ref`, `eval_comp_id`, or `exclusion_set` that does not exist in `config/resources.yml`, the pipeline currently fails deep inside Snakemake with a cryptic `KeyError`. Operators onboarding a new genome waste time tracing the failure back to the offending config field.

The existing schema validation in `rules/common.smk` (via `snakemake.utils.validate`) catches *structural* errors but does not check cross-references between `analyses.tsv` and `resources.yml`.

## Goals

1. Detect missing cross-references at config-load time and report them with one clear, grouped message that names every offending ID and which `eval_id`(s) reference it.
2. Provide a CLI flag — `run_defrabb --validate-only` — that runs only the validation, exits 0 on success, and exits non-zero on failure, without executing any pipeline rules.

## Non-Goals

- Validating the dot-separated tokens inside `bench_vcf_processing` strings. That vocabulary is being reframed by TODO #11 (named processing profiles); validating the current ad-hoc tokens would be wasted work.
- Verifying that downloadable resources are reachable (no URL liveness checks).
- Replacing `snakemake.utils.validate` / jsonschema with pydantic. Worth considering separately as part of a larger configuration refactor; out of scope here.
- Revising the structure of `analyses.tsv` or `resources.yml`. Same reasoning — separate decision.

## Design

### Validation logic — `scripts/validate_configs.py`

The function lives in a plain Python module so unit tests can `import` it directly. `rules/common.smk` imports it with a small sys.path shim (3 lines) so the function runs at workflow load time.

```python
# scripts/validate_configs.py
from snakemake.exceptions import WorkflowError


def validate_cross_references(config, analyses):
    """Verify every ID referenced in analyses.tsv exists in resources.yml.

    Raises WorkflowError with a single grouped message if any are missing.
    """
    errors = {
        "assemblies": {},          # asm_id -> [eval_ids]
        "references": {},          # ref    -> [eval_ids]
        "exclusion_sets": {},      # excl   -> [eval_ids]
        "comparisons": {},         # (ref, comp_id) -> [eval_ids]
    }

    asm_ids = set(config["assemblies"])
    ref_ids = set(config["references"])
    excl_set_ids = set(config["exclusion_set"])

    for eval_id, row in analyses.iterrows():
        if row["asm_id"] not in asm_ids:
            errors["assemblies"].setdefault(row["asm_id"], []).append(eval_id)

        if row["ref"] not in ref_ids:
            errors["references"].setdefault(row["ref"], []).append(eval_id)

        if row["exclusion_set"] != "none" and row["exclusion_set"] not in excl_set_ids:
            errors["exclusion_sets"].setdefault(row["exclusion_set"], []).append(eval_id)

        # Only check comparisons[ref][comp_id] if ref is itself valid;
        # otherwise the missing-ref error already covers the situation.
        if row["ref"] in ref_ids:
            comp_for_ref = set(config.get("comparisons", {}).get(row["ref"], {}))
            if row["eval_comp_id"] not in comp_for_ref:
                errors["comparisons"].setdefault(
                    (row["ref"], row["eval_comp_id"]), []
                ).append(eval_id)

    if any(errors.values()):
        raise WorkflowError(_format_grouped_errors(errors))
```

Plus a small `_format_grouped_errors(errors)` helper that produces output of the form:

```
Config validation failed: missing cross-references in resources.yml.

Missing assemblies (resources.yml:assemblies):
  - HG002_v1.2 (used by eval_ids: eval_HG002_dipcall, eval_HG002_pav)

Missing references (resources.yml:references):
  - GRCh38_v1.0 (used by eval_ids: eval_HG002_dipcall)

Missing comparisons (resources.yml:comparisons):
  - GRCh38 / HG002_v4.2.1 (used by eval_ids: eval_HG002_pav)

Missing exclusion sets (resources.yml:exclusion_set):
  - my_excl_set (used by eval_ids: eval_HG002_dipcall)
```

Sections with no errors are omitted from the output.

**Call site:** in `rules/common.smk`, immediately after `analyses = analyses.set_index("eval_id")` (currently around line 460). This guarantees both inputs are loaded and indexed before validation runs, and runs before any rule logic touches the config.

**Import shim in common.smk:**

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from validate_configs import validate_cross_references
```

**Trade-off:** keeping the function in a `.py` adds the 3-line `sys.path` shim, but makes the function trivially unit-testable. Embedding the function in `common.smk` would avoid the shim but force tests to use `exec()` or `runpy` gymnastics to load a `.smk` file — worse on net. The shim is small, lives in one place, and follows a pattern already common in Snakemake workflows.

### CLI surface — `run_defrabb --validate-only`

Add a `--validate-only` flag to `run_defrabb`'s argparse setup. When present:

1. Skip run-directory creation, file logging, git status, and conda env export — none are needed for pure validation.
2. Build the command:
   ```
   snakemake --list-target-rules --quiet \
       --config analyses=<args.analyses>
   ```
   Note: no `--directory` because we are not creating a run directory.
3. Run via `subprocess.run`. On exit code 0, log "Config validation passed." and exit 0. On non-zero, propagate the exit code (Snakemake will have already printed the `WorkflowError`).

`--validate-only` is mutually exclusive with `--steps`. Implement via argparse `add_mutually_exclusive_group()` so passing both raises a clear CLI error rather than silently ignoring `--steps`.

**Fallback:** if testing reveals that `--list-target-rules` short-circuits before `common.smk` is loaded (and therefore before `validate_cross_references` runs), swap to `snakemake --dry-run --quiet`. Dry-run is guaranteed to load the workflow but may emit additional DAG-resolution noise.

### Tests — `.tests/unit/test_config_validation.py`

New pytest file. Imports `validate_cross_references` directly. Cases:

1. **Happy path** — minimal valid `config` dict + 1-row analyses DataFrame; function returns `None`.
2. **Missing `asm_id`** — error message names the asm_id and includes the offending eval_id.
3. **Missing `ref`** — same.
4. **Missing `exclusion_set`** — verify `"none"` is treated as valid (no error).
5. **Missing `eval_comp_id` for an existing ref** — verify the error groups under `comparisons` and references both the ref and the comp_id.
6. **Multiple errors at once** — verify all sections appear in one raised error and that one missing ID used by N rows is grouped (not duplicated N times).
7. **Smoke** — load the actual `config/resources.yml` + `config/analyses.tsv` from the repo and verify validation passes.

Plus an integration smoke step (manual or shell-test): `run_defrabb --validate-only -r 20260417_v0.022_validate-test -a config/analyses.tsv` exits 0 on the working configs and exits non-zero with a readable error when an asm_id is mistyped.

## Risk and Cleanup Notes

- The validation function is small enough (~50 LOC plus ~20 LOC helper) that if it gets replaced when TODO #11 (named processing profiles) lands or when a larger config refactor ships, throwing it away is cheap.
- The implementation plan should verify the `--list-target-rules` behavior before committing to it (Section "CLI surface" lists the dry-run fallback).
- The `sys.path` shim in `common.smk` could be replaced by Snakemake's `workflow.source_path` if a cleaner pattern emerges; a future refactor could centralize this if more `scripts/` modules need importing.

## Out-of-Scope Decisions Flagged for Future Work

These came up in brainstorming and deserve their own design pass; they should not be folded into this work:

- **Pydantic for schema validation.** Would replace `snakemake.utils.validate` + jsonschema. Doesn't subsume the cross-reference layer (JSON Schema can't express cross-key constraints either). Worth considering as part of a larger configuration refactor.
- **Config structure revision.** The `comparisons` schema in `schema/resources-schema.yml` looks under-validating, and `bench_vcf_processing` dot-strings are exactly what TODO #11 redoes. Bigger conversation.
- **Refactoring `common.smk`.** It is 508 lines and growing. Splitting it deserves its own brainstorming pass.
