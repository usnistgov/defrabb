# Config Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add cross-reference validation between `analyses.tsv` and `resources.yml` that fails fast with a clear, grouped error message, plus a `run_defrabb --validate-only` flag that runs only validation without executing pipeline rules.

**Architecture:** A pure-Python `validate_cross_references` function in `scripts/validate_configs.py` (importable for unit tests). `rules/common.smk` imports it via a 3-line `sys.path` shim and calls it after the existing schema validation, raising a `WorkflowError` with all errors grouped by section. `run_defrabb --validate-only` invokes `snakemake --list-target-rules --quiet` so the workflow loads (triggering validation) without any rule execution.

**Tech Stack:** Python 3.12, pandas, pytest, Snakemake 8, snakemake.exceptions.WorkflowError.

**Spec:** `docs/superpowers/specs/2026-04-17-config-validation-design.md`

---

## File Structure

**Create:**
- `scripts/validate_configs.py` — `validate_cross_references` and `_format_grouped_errors` (one focused module, ~80 LOC)
- `.tests/unit/test_config_validation.py` — pytest tests for the validation function

**Modify:**
- `rules/common.smk` — add `sys.path` shim + import + call to validation function (around line 460, after `analyses = analyses.set_index("eval_id")`)
- `run_defrabb` — add `--validate-only` flag in argparse; add `validate_only(args)` function; route in `main()`

---

### Task 1: Create validate_configs.py with happy-path test

**Files:**
- Create: `scripts/validate_configs.py`
- Test: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the failing test**

Create `.tests/unit/test_config_validation.py`:

```python
"""Tests for cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

# Make scripts/ importable (mirrors the shim used in rules/common.smk).
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from validate_configs import validate_cross_references  # noqa: E402


def _minimal_config():
    return {
        "assemblies": {"HG002_v1.0": {}},
        "references": {"GRCh38": {}},
        "comparisons": {"GRCh38": {"HG002_v4.2.1": {}}},
        "exclusion_set": {"smvar": []},
    }


def _analyses_df(rows):
    """Return an analyses DataFrame indexed by eval_id."""
    df = pd.DataFrame(rows)
    return df.set_index("eval_id")


def test_happy_path_returns_none():
    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval1",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        }
    ])
    assert validate_cross_references(config, analyses) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest .tests/unit/test_config_validation.py::test_happy_path_returns_none -v`

Expected: `ModuleNotFoundError: No module named 'validate_configs'`

- [ ] **Step 3: Create the minimal validate_configs.py**

Create `scripts/validate_configs.py`:

```python
"""Cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.
"""
from snakemake.exceptions import WorkflowError


def validate_cross_references(config, analyses):
    """Verify every ID referenced in analyses.tsv exists in resources.yml.

    Raises WorkflowError with a single grouped message if any IDs are missing.
    Returns None on success.
    """
    return None
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `pytest .tests/unit/test_config_validation.py::test_happy_path_returns_none -v`

Expected: `1 passed`

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_configs.py .tests/unit/test_config_validation.py
git commit -m "config validation: stub module + happy path test"
```

---

### Task 2: Validate asm_id

**Files:**
- Modify: `scripts/validate_configs.py`
- Modify: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the failing test**

Append to `.tests/unit/test_config_validation.py`:

```python
def test_missing_asm_id_raises_with_eval_id_in_message():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_typo",
            "asm_id": "HG002_TYPO",
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        }
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "HG002_TYPO" in msg
    assert "eval_typo" in msg
    assert "assemblies" in msg
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest .tests/unit/test_config_validation.py::test_missing_asm_id_raises_with_eval_id_in_message -v`

Expected: `Failed: DID NOT RAISE`

- [ ] **Step 3: Implement asm_id check + minimal formatter**

Replace `scripts/validate_configs.py` with:

```python
"""Cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.
"""
from snakemake.exceptions import WorkflowError


def validate_cross_references(config, analyses):
    """Verify every ID referenced in analyses.tsv exists in resources.yml.

    Raises WorkflowError with a single grouped message if any IDs are missing.
    Returns None on success.
    """
    errors = {
        "assemblies": {},
        "references": {},
        "exclusion_sets": {},
        "comparisons": {},
    }

    asm_ids = set(config["assemblies"])

    for eval_id, row in analyses.iterrows():
        if row["asm_id"] not in asm_ids:
            errors["assemblies"].setdefault(row["asm_id"], []).append(eval_id)

    if any(errors.values()):
        raise WorkflowError(_format_grouped_errors(errors))
    return None


def _format_grouped_errors(errors):
    sections = []
    if errors["assemblies"]:
        sections.append(_format_section(
            "Missing assemblies (resources.yml:assemblies)",
            errors["assemblies"],
        ))
    return (
        "Config validation failed: missing cross-references in resources.yml.\n\n"
        + "\n\n".join(sections)
    )


def _format_section(title, missing):
    lines = [f"{title}:"]
    for missing_id, eval_ids in missing.items():
        lines.append(f"  - {missing_id} (used by eval_ids: {', '.join(eval_ids)})")
    return "\n".join(lines)
```

- [ ] **Step 4: Run both tests to verify they pass**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `2 passed`

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_configs.py .tests/unit/test_config_validation.py
git commit -m "config validation: detect missing asm_id"
```

---

### Task 3: Validate ref

**Files:**
- Modify: `scripts/validate_configs.py`
- Modify: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the failing test**

Append:

```python
def test_missing_ref_raises_with_eval_id_in_message():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_badref",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38_TYPO",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        }
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "GRCh38_TYPO" in msg
    assert "eval_badref" in msg
    assert "references" in msg
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest .tests/unit/test_config_validation.py::test_missing_ref_raises_with_eval_id_in_message -v`

Expected: `Failed: DID NOT RAISE`

- [ ] **Step 3: Implement ref check + extend formatter**

In `scripts/validate_configs.py`, inside `validate_cross_references`, add after the `asm_ids` line and inside the loop:

```python
    ref_ids = set(config["references"])
```

Add inside the `for eval_id, row in analyses.iterrows():` loop, after the asm_id block:

```python
        if row["ref"] not in ref_ids:
            errors["references"].setdefault(row["ref"], []).append(eval_id)
```

In `_format_grouped_errors`, add a second block after the assemblies one:

```python
    if errors["references"]:
        sections.append(_format_section(
            "Missing references (resources.yml:references)",
            errors["references"],
        ))
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `3 passed`

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_configs.py .tests/unit/test_config_validation.py
git commit -m "config validation: detect missing ref"
```

---

### Task 4: Validate exclusion_set (and skip "none")

**Files:**
- Modify: `scripts/validate_configs.py`
- Modify: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the failing tests**

Append:

```python
def test_missing_exclusion_set_raises():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_badexcl",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "not_a_real_set",
        }
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "not_a_real_set" in msg
    assert "eval_badexcl" in msg
    assert "exclusion" in msg


def test_exclusion_set_none_is_valid():
    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_no_excl",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "none",
        }
    ])
    assert validate_cross_references(config, analyses) is None
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest .tests/unit/test_config_validation.py::test_missing_exclusion_set_raises .tests/unit/test_config_validation.py::test_exclusion_set_none_is_valid -v`

Expected: `test_missing_exclusion_set_raises` fails with `DID NOT RAISE`; `test_exclusion_set_none_is_valid` passes (no validation logic yet for this section).

- [ ] **Step 3: Implement exclusion_set check**

In `scripts/validate_configs.py`, add inside `validate_cross_references` after the `ref_ids` line:

```python
    excl_set_ids = set(config["exclusion_set"])
```

Add inside the for-loop after the ref block:

```python
        if row["exclusion_set"] != "none" and row["exclusion_set"] not in excl_set_ids:
            errors["exclusion_sets"].setdefault(row["exclusion_set"], []).append(eval_id)
```

Add to `_format_grouped_errors` after the references block:

```python
    if errors["exclusion_sets"]:
        sections.append(_format_section(
            "Missing exclusion sets (resources.yml:exclusion_set)",
            errors["exclusion_sets"],
        ))
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `5 passed`

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_configs.py .tests/unit/test_config_validation.py
git commit -m "config validation: detect missing exclusion_set, treat 'none' as valid"
```

---

### Task 5: Validate comparisons[ref][eval_comp_id]

**Files:**
- Modify: `scripts/validate_configs.py`
- Modify: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the failing tests**

Append:

```python
def test_missing_comp_id_for_existing_ref_raises():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_badcomp",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38",
            "eval_comp_id": "v4.2.1_TYPO",
            "exclusion_set": "smvar",
        }
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "v4.2.1_TYPO" in msg
    assert "GRCh38" in msg
    assert "eval_badcomp" in msg
    assert "comparisons" in msg


def test_missing_ref_does_not_double_report_comp_id():
    """If ref itself is missing, we don't also flag the comp_id under it."""
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_double",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38_NOTREAL",
            "eval_comp_id": "anything",
            "exclusion_set": "smvar",
        }
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)
    assert "GRCh38_NOTREAL" in msg
    # comp_id error should NOT appear since the ref is missing
    assert "Missing comparisons" not in msg
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest .tests/unit/test_config_validation.py::test_missing_comp_id_for_existing_ref_raises .tests/unit/test_config_validation.py::test_missing_ref_does_not_double_report_comp_id -v`

Expected: First test fails with `DID NOT RAISE`. Second test passes already (no comparisons logic yet → no double report possible).

- [ ] **Step 3: Implement comparisons check**

In `scripts/validate_configs.py`, add inside the for-loop after the exclusion_set block:

```python
        # Only check comparisons[ref][comp_id] if ref is itself valid;
        # otherwise the missing-ref error already covers the situation.
        if row["ref"] in ref_ids:
            comp_for_ref = set(config.get("comparisons", {}).get(row["ref"], {}))
            if row["eval_comp_id"] not in comp_for_ref:
                errors["comparisons"].setdefault(
                    (row["ref"], row["eval_comp_id"]), []
                ).append(eval_id)
```

Update `_format_grouped_errors` — add a comparisons block after exclusion_sets. The comparisons key is a `(ref, comp_id)` tuple, so use a dedicated formatter:

```python
    if errors["comparisons"]:
        lines = ["Missing comparisons (resources.yml:comparisons):"]
        for (ref, comp_id), eval_ids in errors["comparisons"].items():
            lines.append(f"  - {ref} / {comp_id} (used by eval_ids: {', '.join(eval_ids)})")
        sections.append("\n".join(lines))
```

- [ ] **Step 4: Run all tests to verify they pass**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `7 passed`

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_configs.py .tests/unit/test_config_validation.py
git commit -m "config validation: detect missing eval_comp_id under valid ref"
```

---

### Task 6: Multi-error grouping (no duplicates)

**Files:**
- Modify: `.tests/unit/test_config_validation.py`

This task verifies the formatter behavior when multiple sections fail and when one missing ID is referenced from multiple rows. No production-code changes expected.

- [ ] **Step 1: Write the failing test**

Append:

```python
def test_multiple_errors_grouped_in_single_raise():
    from snakemake.exceptions import WorkflowError

    config = _minimal_config()
    analyses = _analyses_df([
        {
            "eval_id": "eval_a",
            "asm_id": "HG002_TYPO",       # missing asm
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        },
        {
            "eval_id": "eval_b",
            "asm_id": "HG002_TYPO",       # same missing asm — should group
            "ref": "GRCh38",
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        },
        {
            "eval_id": "eval_c",
            "asm_id": "HG002_v1.0",
            "ref": "GRCh38_TYPO",          # missing ref (different section)
            "eval_comp_id": "HG002_v4.2.1",
            "exclusion_set": "smvar",
        },
    ])
    with pytest.raises(WorkflowError) as exc_info:
        validate_cross_references(config, analyses)
    msg = str(exc_info.value)

    # Both sections present
    assert "Missing assemblies" in msg
    assert "Missing references" in msg

    # Same missing asm grouped on one line, both eval_ids listed
    assert "HG002_TYPO" in msg
    assert "eval_a" in msg
    assert "eval_b" in msg
    # Only one occurrence of the missing asm token in the asm section
    assert msg.count("HG002_TYPO") == 1
```

- [ ] **Step 2: Run the test to verify it passes**

Run: `pytest .tests/unit/test_config_validation.py::test_multiple_errors_grouped_in_single_raise -v`

Expected: PASS — the existing implementation already groups by section and by missing ID.

If it fails, debug the formatter and fix before proceeding.

- [ ] **Step 3: Run all tests to confirm**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `8 passed`

- [ ] **Step 4: Commit**

```bash
git add .tests/unit/test_config_validation.py
git commit -m "config validation: test multi-error grouping"
```

---

### Task 7: Smoke test against real repo configs

**Files:**
- Modify: `.tests/unit/test_config_validation.py`

- [ ] **Step 1: Write the smoke test**

Append:

```python
def test_real_repo_configs_validate():
    """Smoke: the actual repo config + default analyses.tsv must pass validation."""
    import yaml

    resources_path = REPO_ROOT / "config" / "resources.yml"
    analyses_path = REPO_ROOT / "config" / "analyses.tsv"

    with resources_path.open() as fh:
        config = yaml.safe_load(fh)

    analyses_df = pd.read_table(analyses_path).set_index("eval_id")

    assert validate_cross_references(config, analyses_df) is None
```

- [ ] **Step 2: Run the smoke test**

Run: `pytest .tests/unit/test_config_validation.py::test_real_repo_configs_validate -v`

Expected: PASS.

If it fails, the real configs have a cross-reference issue. Triage:
- If the issue is a real config bug, fix the config (separate commit) and re-run.
- If the issue is a validator false positive (e.g., a key the validator misses), fix the validator.

- [ ] **Step 3: Run the full test file**

Run: `pytest .tests/unit/test_config_validation.py -v`

Expected: `9 passed`

- [ ] **Step 4: Commit**

```bash
git add .tests/unit/test_config_validation.py
git commit -m "config validation: smoke test against real repo configs"
```

---

### Task 8: Wire validation into common.smk

**Files:**
- Modify: `rules/common.smk`

- [ ] **Step 1: Add the import shim and call**

Open `rules/common.smk`. At the top of the file, immediately after the existing imports (line 1-3), add:

```python
import sys
from pathlib import Path
```

(If `sys` or `Path` are already imported there, do not duplicate.)

- [ ] **Step 2: Add the validation call**

Find the line in `rules/common.smk`:

```python
analyses = analyses.set_index("eval_id")
```

(currently around line 460). Immediately after that line, add:

```python

# Cross-reference validation: catches asm_id/ref/comp_id/exclusion_set
# typos in analyses.tsv before any rule expansion.
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from validate_configs import validate_cross_references

validate_cross_references(config, analyses)
```

- [ ] **Step 3: Verify validation triggers on the working configs**

Run from the repo root:

```bash
snakemake --list-target-rules --quiet --config analyses=config/analyses.tsv
```

Expected: command exits 0 and prints rule names. No `WorkflowError`.

If it fails with a `WorkflowError`, the real configs have an issue — same triage as Task 7 Step 2.

- [ ] **Step 4: Verify validation triggers on a deliberately-broken config**

Create a temporary bad analyses TSV:

```bash
sed 's/HG002/HG002_TYPO/g' config/analyses.tsv > /tmp/analyses_bad.tsv
snakemake --list-target-rules --quiet --config analyses=/tmp/analyses_bad.tsv 2>&1 | head -40
rm /tmp/analyses_bad.tsv
```

Expected: command exits non-zero. Output contains `Config validation failed` and `Missing assemblies` (or `Missing references`, depending on which IDs the sed replaced).

If `--list-target-rules` does NOT trigger validation (i.e., the bad config passes silently), fall back to `--dry-run --quiet` for both verification and the run_defrabb wiring in Task 9. Update Task 9 accordingly.

- [ ] **Step 5: Run the unit tests one more time**

Run: `pytest .tests/unit/ -v`

Expected: all tests pass (existing + new).

- [ ] **Step 6: Commit**

```bash
git add rules/common.smk
git commit -m "common.smk: call validate_cross_references after schema validation"
```

---

### Task 9: Add --validate-only to run_defrabb

**Files:**
- Modify: `run_defrabb`

- [ ] **Step 1: Add the --validate-only flag (mutually exclusive with --steps)**

In `run_defrabb`, locate `parse_arguments()`. The current `--steps` argument is added on a `group` (lines 128-144). Wrap `--steps` and a new `--validate-only` in a mutex group.

Replace lines 128-144 of `run_defrabb` with:

```python
    # Create an argument group for better formatting
    group = parser.add_argument_group("workflow steps")
    mode_mutex = group.add_mutually_exclusive_group()

    # Add the 'steps' argument with a detailed description using line breaks
    mode_mutex.add_argument(
        "-s",
        "--steps",
        type=str,
        choices=["all", "pipe", "report", "archive", "release"],
        default="all",
        metavar="",
        help="""Defining which workflow steps are run:
    all: pipe, report, and archive (default)
    pipe: just the snakemake pipeline
    report: generating the snakemake run report
    archive: generating snakemake archive tarball
    release: copy run output to NAS for upload to Google Drive (internal NIST use-case)""",
    )
    mode_mutex.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate analyses.tsv and resources.yml; do not run the pipeline.",
    )
```

- [ ] **Step 2: Add the validate_only function**

In `run_defrabb`, add this function near the other top-level workflow-step functions (e.g., right above `execute_snakemake_pipeline`):

```python
def validate_only(args: argparse.Namespace) -> int:
    """Run config validation only (no pipeline execution).

    Triggers Snakemake's workflow load (which runs validate_cross_references
    in rules/common.smk) but does not resolve the DAG or execute any rule.
    Returns the exit code.
    """
    cmd = [
        "snakemake",
        "--list-target-rules",
        "--quiet",
        "--config",
        f"analyses={args.analyses}",
    ]
    logging.info(f"Validation command: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode == 0:
        logging.info("Config validation passed.")
    return result.returncode
```

If Task 8 Step 4 forced the dry-run fallback, replace `["--list-target-rules", "--quiet"]` with `["--dry-run", "--quiet"]`.

- [ ] **Step 3: Route --validate-only in main()**

In `run_defrabb`, in `main()`, after the `args, extra_args = parse_arguments()` line and BEFORE `release_config = load_release_config(...)`, add:

```python
    # --validate-only is the simplest path: parse args, validate the analyses file
    # exists, then shell out to snakemake. No release config, no run dir, no logging.
    if args.validate_only:
        if not args.analyses:
            args.analyses = f"config/analyses_{args.runid}.tsv"
        if not os.path.exists(args.analyses):
            print(f"Error: Analyses file '{args.analyses}' does not exist.", file=sys.stderr)
            sys.exit(1)
        logging.basicConfig(level=logging.INFO, format="%(message)s")
        sys.exit(validate_only(args))
```

- [ ] **Step 4: Smoke test --validate-only on the real configs**

Run from the repo root:

```bash
./run_defrabb -r 20260417_v0.022_validate-test -a config/analyses.tsv --validate-only
```

Expected: prints `Config validation passed.` and exits 0.

Verify exit code:

```bash
./run_defrabb -r 20260417_v0.022_validate-test -a config/analyses.tsv --validate-only; echo "exit=$?"
```

Expected: `exit=0`.

- [ ] **Step 5: Smoke test --validate-only on a broken config**

```bash
sed 's/HG002/HG002_TYPO/g' config/analyses.tsv > /tmp/analyses_bad.tsv
./run_defrabb -r 20260417_v0.022_validate-test -a /tmp/analyses_bad.tsv --validate-only; echo "exit=$?"
rm /tmp/analyses_bad.tsv
```

Expected: prints `Config validation failed:` followed by grouped error output, then `exit=1` (or other non-zero).

- [ ] **Step 6: Smoke test that --validate-only and --steps are mutually exclusive**

```bash
./run_defrabb -r 20260417_v0.022_validate-test -a config/analyses.tsv --validate-only -s pipe
```

Expected: argparse error like `argument --validate-only: not allowed with argument -s/--steps` and exit code 2.

- [ ] **Step 7: Commit**

```bash
git add run_defrabb
git commit -m "run_defrabb: add --validate-only flag"
```

---

### Task 10: Update TODO doc to reflect items #8 + #9 done

**Files:**
- Modify: `docs/TODO-pre-development.md`

- [ ] **Step 1: Mark items #8 and #9 as done**

In `docs/TODO-pre-development.md`, update the heading for item #8 from:

```markdown
### 8. Add a `--validate-only` flag to run_defrabb
```

to:

```markdown
### 8. Add a `--validate-only` flag to run_defrabb — ✅ Done
```

Add a brief resolution note immediately under the heading:

```markdown
Resolved 2026-04-17. `run_defrabb --validate-only` invokes
`snakemake --list-target-rules --quiet` so the workflow loads and
`validate_cross_references` runs (see item 9), without executing rules.
```

Repeat for item #9: change heading to `### 9. Improve error messages for missing config entries — ✅ Done` and add:

```markdown
Resolved 2026-04-17. `scripts/validate_configs.py` now provides
`validate_cross_references`, called from `rules/common.smk` after the
existing schema validation. Reports missing asm_id, ref, eval_comp_id,
and exclusion_set with grouped errors that name the offending IDs and
the eval_ids that reference them.
```

In the Priority Order section at the bottom of the file, change:

```markdown
1. **Items 8, 9** — Config validation (prevents wasted compute on bad configs)
```

to:

```markdown
1. ~~Items 8, 9~~ — Config validation completed 2026-04-17.
```

- [ ] **Step 2: Commit**

```bash
git add docs/TODO-pre-development.md
git commit -m "docs: mark TODO #8 + #9 (config validation) done"
```

---

## Verification Summary

At the end of all 10 tasks, verify:

```bash
pytest .tests/unit/ -v
./run_defrabb -r 20260417_v0.022_validate-test -a config/analyses.tsv --validate-only
git log --oneline -12
```

Expected:
- All pytest tests pass (9 new + previously passing)
- `--validate-only` exits 0 with `Config validation passed.`
- Git log shows roughly one focused commit per task

## Notes for the Implementer

- The user's commit-message style is short, lowercase subjects (`docs:`, `config validation:`, etc.). Keep messages terse — they describe one change.
- Do NOT push to remote unless the user asks. The user works directly on master.
- If at any point a test fails for a reason that suggests the design is wrong (not just the code), STOP and surface the question rather than working around it. The spec is `docs/superpowers/specs/2026-04-17-config-validation-design.md`.
- The `workflow.basedir` reference inside the `common.smk` shim is provided by Snakemake at workflow-load time. It does not exist at module-import time outside Snakemake — that is why tests use `Path(__file__).resolve().parents[2]` instead.
