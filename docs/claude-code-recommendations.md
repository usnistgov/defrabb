# Claude Code Recommendations for DeFrABB Development

_Prepared: 2026-03-24_

This document recommends hooks, skills, custom skills, and resources to facilitate DeFrABB pipeline development using Claude Code.

## Recommended Hooks

Hooks are automated checks that run before/after tool calls. These prevent common mistakes during DeFrABB development.

### 1. PreToolUse: Protect analyses config format

Block writes to `config/analyses*.tsv` that don't follow the required column schema (18+ columns, tab-separated). Prevents accidentally corrupting the analysis matrix.

```
Rule: When editing config/analyses*.tsv files, verify tab-separated format and required column headers are preserved.
```

### 2. PreToolUse: Validate conda env version consistency

When editing `envs/*.yml` files, warn if a package version conflicts with the same package in another env file (e.g., bcftools 1.14 vs 1.17).

```
Rule: When modifying envs/*.yml, check for version conflicts of the same package across all env files.
```

### 3. PreToolUse: Guard resources.yml structure

Block destructive edits to `config/resources.yml` that might remove required top-level keys (references, assemblies, comparisons, exclusion_set).

```
Rule: When editing config/resources.yml, ensure required top-level keys are preserved.
```

### 4. PostToolUse: Remind to run tests after rule changes

After editing any `rules/*.smk` file, remind to run the corresponding unit test and/or smoke test.

```
Rule: After editing rules/*.smk, remind about relevant tests in .tests/unit/.
```

### 5. PreToolUse: Prevent NIST-specific hardcoding

When editing `run_defrabb` or scripts, flag hardcoded NIST paths (e.g., `/mnt/bbdhg-nas/`, `giab-data` bucket) and suggest using config/CLI overrides instead.

```
Rule: When writing code in run_defrabb or scripts/, flag hardcoded NIST-specific paths.
```

## Recommended Custom Skills

### 1. `/validate-analyses` — Validate an analyses TSV

Skill that reads an analyses TSV, cross-references against resources.yml, and reports:
- Missing assembly/reference/comparison IDs
- Invalid eval_cmd values
- Undefined exclusion_set names
- VCF processing steps that don't correspond to existing rules
- Duplicate eval_ids

This is critical for the parameter optimization workflow planned for new genomes.

### 2. `/dry-run` — Snakemake dry run with DAG summary

Skill that runs `snakemake --dryrun --use-conda --use-apptainer` and summarizes:
- Total rules to execute
- Which assemblies/references/evaluations are involved
- Estimated resource requirements
- Any missing inputs

### 3. `/compare-evaluations` — Compare evaluation results across runs

Skill that reads hap.py summary.csv or truvari summary.json files from two runs and produces a side-by-side comparison table. Essential for parameter optimization work.

### 4. `/new-genome` — Scaffold config for a new genome/assembly

Interactive skill that walks through adding a new assembly to the pipeline:
1. Adds entry to resources.yml (assemblies section)
2. Creates an analyses TSV with sensible defaults
3. Validates the config
4. Suggests a dry-run command

### 5. `/exclusion-report` — Summarize exclusion impact

Skill that reads exclusion_stats.txt and exclusion_intersection_summary.csv files and produces a human-readable report of what was excluded and why, with base-pair counts and percentages.

### 6. `/env-audit` — Audit conda environment consistency

Skill that scans all `envs/*.yml` files and reports:
- Version inconsistencies for shared packages
- Packages with loose version constraints
- Outdated Python versions
- Mixed conda/pip installs

## Recommended Plugins and MCP Servers

### 1. Context7 (already installed)

Use for looking up Snakemake 8.x API docs, pybedtools, truvari, and hap.py documentation during development.

### 2. bioRxiv / PubMed (already installed)

Useful for finding methodology papers when implementing new exclusion strategies or evaluation approaches for new genomes. Search for GIAB benchmark papers to ensure alignment with published methods.

## Documentation Improvements

### 1. Known Limitations Document

Create `docs/known-limitations.md` covering:
- FIPS OpenSSL incompatibility with conda environments
- PAV container dependency (requires Apptainer)
- Memory requirements for hap.py (64GB default)
- Chr21 test data limitations vs full-genome behavior
- GitHub Actions CI is stale (references old paths)

### 2. Parameter Optimization Guide

Create `docs/parameter-optimization.md` for the planned new-genome work:
- How to create a new analyses TSV for parameter sweeps
- Which parameters to vary (vc_params, exclusion_set, bench_vcf_processing)
- How to compare evaluation results across parameter sets
- Decision criteria for choosing final parameters

### 3. Exclusion Strategy Guide

Create `docs/exclusion-strategies.md` documenting:
- What each exclusion type does and why
- The slop/merge constants and their rationale
- How to design exclusion sets for new genomes
- How to interpret exclusion_stats.txt and intersection summaries

### 4. Conda Environment Management Guide

Create `docs/conda-environments.md` documenting:
- Why per-rule environments are used
- How to update tool versions safely
- Known FIPS workarounds
- Testing changes to environments

## CLAUDE.md Improvements

The existing CLAUDE.md is solid. Suggested additions:

1. **Common gotchas section** — FIPS OpenSSL workaround, PAV callable region bug, report source duplication
2. **Config editing guidance** — link to schema files, note which fields are validated
3. **Exclusion system overview** — brief explanation of the exclusion categories since this is a major development area
4. **Link to architecture diagram** — reference docs/architecture-diagram.md
