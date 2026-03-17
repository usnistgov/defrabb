# DeFrABB Codebase Architecture Review

_Prepared: 2026-03-17_

## Purpose

This note captures the current structure of the DeFrABB repository, the main workflow boundaries, and the highest-value maintainability risks observed during an initial codebase review.

## What DeFrABB currently does well

- Clear top-level separation between workflow rules, helper scripts, environments, schemas, configs, and report assets.
- Strong use of Snakemake modularization through `rules/*.smk`.
- Good scientific provenance intent: versioned analysis tables, benchmark outputs, logs, reports, archive generation, and release packaging.
- Practical support for both assembly-based calling and downstream evaluation.

## Main entrypoints

- `Snakefile`: workflow entry point; includes all rule modules.
- `run_defrabb`: operational wrapper for running the pipeline, generating reports, building archives, and releasing outputs.
- `config/analyses*.tsv`: run-specific analysis matrix.
- `config/resources.yml`: shared resource catalog, tool settings, exclusion presets, and reference/assembly metadata.

## Workflow layers

```text
analyses.tsv + resources.yml
        ↓
rules/common.smk
  - schema validation
  - derived tables: vc_tbl, bench_tbl, analyses
        ↓
resource acquisition
  - references, assemblies, comparison callsets, stratifications
        ↓
assembly variant calling
  - dipcall / pav
  - standardized VCF + baseline BED
        ↓
benchmark derivation
  - VCF post-processing
  - exclusion BED generation + subtraction
        ↓
evaluation
  - hap.py
  - truvari / truvari refine
        ↓
reporting and release
  - Quarto/HTML
  - Snakemake archive
  - local/S3 release packaging
```

## Current configuration model

### 1. `config/analyses*.tsv`

Each row describes an evaluation scenario and carries enough information to derive:

- the assembly-based callset (`vc_*`)
- the draft benchmark (`bench_*`)
- the comparison callset and evaluation mode (`eval_*`)

This table is then decomposed in `rules/common.smk` into:

- `vc_tbl`: assembly-calling jobs
- `bench_tbl`: benchmark derivation jobs
- `analyses`: evaluation jobs indexed by `eval_id`

### 2. `config/resources.yml`

This file currently mixes several concerns:

- reference metadata and download URLs
- assembly metadata
- comparison callset metadata
- exclusion presets
- exclusion transform rules
- tool defaults / thread / memory settings

That makes it powerful, but also large and hard to reason about.

## Key maintainability hotspots

## 1. Benchmark VCF processing is encoded as filename suffix choreography

The current benchmark VCF processing model is driven by values like:

- `none`
- `fix_XY_genotype`
- `norm.fix_XY_genotype.trfanno.svinfo.lcr.remap`

That works, but it couples:

- workflow behavior
- output naming
- provenance encoding
- user-facing API

into one string convention.

### Why this is risky

- It is hard to tell which processing combinations are valid.
- Ordering is implicit rather than validated.
- It is difficult to add new processing profiles cleanly.
- `bench_bed_processing` exists in the schema but is not a similarly strong first-class workflow abstraction.

## 2. Exclusions are scientifically rich but structurally hard to evolve

The exclusion system combines:

- reference-agnostic exclusions
- assembly-specific exclusions
- hard-coded slop / slop+merge behavior
- assembly-break intersection logic
- sample/version-specific exclusion presets

with heavy reliance on dynamic path construction in `get_exclusion_inputs()`.

### Current pain points

- Several transforms are controlled by hard-coded constants in rules/comments.
- Provenance is distributed across config naming, rule names, and derived filenames.
- Adding new samples/genomes likely means editing large config sections rather than extending a reusable exclusion registry.
- Exclusion generation and exclusion application are tightly coupled.

## 3. Reporting/UI is functional but not yet a coherent product surface

Observed issues:

- There are duplicate/overlapping report source locations: top-level `analysis.qmd` and `scripts/reports/analysis.qmd`.
- `rules/report.smk` renders from `scripts/reports/analysis.qmd`, while generated artifacts also exist at repo root.
- Several `report/*.rst` files are empty or placeholder-sized.
- The Quarto report includes TODO sections and partial wiring (for example assembly-region coverage notes).

### Effect

The report works more like an internal exploratory notebook than a curated “run dashboard”.

## 4. `run_defrabb` mixes orchestration, provenance, packaging, and release logic

The wrapper is useful, but it currently owns too many responsibilities:

- argument validation
- environment capture
- Snakemake execution
- archive generation
- S3 upload
- local copy rules
- manifest generation

### Concrete issues observed

- `main()` hardcodes `bucket_name = "giab-data"` even though CLI/config values are also accepted.
- `upload_to_s3()` checks include rules but ignores exclude rules.
- `data_manifest.tsv` is generated after upload, so it is not part of the uploaded set created by that same release pass.
- Release logic is strongly NIST-oriented, which is fine operationally, but it is mixed into the default user-facing wrapper path.

## 5. Test coverage exists, but much of it is old generated rule-smoke coverage

The `.tests/unit/` suite is valuable, but current limitations include:

- many tests are generated-style Snakemake target tests
- several tests use `--touch`, so they validate wiring more than scientific behavior
- one config-related test is explicitly skipped
- the harness still references Snakemake 7.3-era generation patterns
- the repo does not currently expose an easy bootstrap path for running tests locally

## Notable code-level issues worth addressing early

- `rules/asm-varcall.smk`:
  - `intersect_pav_callable_regions` uses the hap2 callable BED for both `h1_bed` and `h2_bed`, which looks like a real correctness bug.
- `rules/report.smk`:
  - report source paths are duplicated / ambiguous.
- `rules/common.smk`:
  - exclusion construction is powerful but brittle because many behaviors are inferred from string membership in config lists.
- `scripts/exclusion_intersection_summary.py`:
  - duplicated imports / logging setup suggest low cleanup and review pressure in this area.

## Suggested target architecture direction

## A. Treat “processing profiles” as first-class configuration

Instead of string suffix chains, define named profiles such as:

- `smvar_default`
- `stvar_trf_annotated`
- `male_smvar_xy_fix`

Each profile should declare:

- ordered steps
- required inputs
- expected output semantics
- compatible benchmark types

## B. Split exclusion definition from exclusion composition

Introduce two layers:

1. exclusion registry
   - one entry per exclusion source
   - metadata for scope, transform, parameters, provenance
2. exclusion sets
   - reusable named bundles composed from registry entries

## C. Make reporting a single-source system

Use one canonical Quarto source and build a simpler summary-first report:

- run metadata
- assemblies / references
- benchmark coverage
- exclusion impact
- evaluation highlights
- links to raw outputs

## D. Separate run execution from release packaging

Longer term, the wrapper should become:

- `run_defrabb run`
- `run_defrabb report`
- `run_defrabb package`
- `run_defrabb release`

or equivalent internal functions/modules with cleaner boundaries.

## Practical conclusion

DeFrABB already has the right major workflow components. The main opportunity is not inventing new functionality; it is reducing implicit behavior.

The highest-value theme for refactoring is:

> move from string-encoded pipeline behavior toward typed, named, and testable workflow concepts.

