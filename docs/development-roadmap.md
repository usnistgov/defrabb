# DeFrABB Development Roadmap

_Prepared: 2026-03-17_

## Goals

1. Make benchmark generation easier to reason about and safer to modify.
2. Reduce debugging time when adding new samples, references, and exclusion strategies.
3. Improve operator experience for running, reviewing, archiving, and releasing analyses.
4. Preserve scientific provenance while simplifying the code path.

## Guiding principles

- Prefer targeted refactors over large rewrites.
- Preserve current scientific outputs unless a change is explicitly intended.
- Make config semantics explicit before expanding feature surface area.
- Keep NIST-specific operations supported, but isolate them from core workflow logic.
- Add tests and docs at the same time as each refactor.

## Staged plan

## Phase 0 — Stabilize and document the current system

### Outcomes

- Shared understanding of current architecture
- Reduced risk of accidental regressions
- Clear short list of correctness fixes

### Proposed tasks

1. Add architecture and roadmap docs.  
   _This is now started in `docs/`.*
2. Create a developer bootstrap doc for:
   - base environment expectations
   - how to run unit tests
   - how to run a minimal chr21 smoke test
3. Fix obvious correctness issues:
   - PAV callable-region intersection input typo
   - wrapper bucket/release handling issues
   - duplicate/ambiguous report source usage
4. Add a small “known limitations” section to docs instead of leaving them implicit.

## Phase 1 — Refactor benchmark VCF processing into named profiles

### Current problem

`bench_vcf_processing` encodes workflow logic as dot-separated suffix strings.

### Target state

Introduce a named processing-profile registry, for example:

```yaml
vcf_processing_profiles:
  smvar_basic:
    bench_type: smvar
    steps: [fix_XY_genotype]
  stvar_annotated:
    bench_type: stvar
    steps: [normalize_vars, fix_XY_genotype, run_truvari_anno_trf, run_truvari_anno_svinfo]
```

### Benefits

- valid step order becomes explicit
- config becomes easier to review
- provenance becomes easier to explain
- new profiles can be added without inventing filename conventions

### First implementation milestone

- Keep existing suffix strings working.
- Add named-profile support beside them.
- Migrate checked-in analysis tables after validation.

## Phase 2 — Redesign exclusions around a registry + composition model

### Current problem

Exclusions combine source definition, transform behavior, and set membership in a way that is hard to extend.

### Target state

Each exclusion should declare metadata such as:

- scope: `reference`, `assembly`, `benchmark`
- source type: downloaded BED, derived from VCF, derived from baseline BED
- transforms: `sort`, `slop`, `merge`, `intersect_breaks`
- parameters
- applicability by benchmark type / caller / reference

Then named exclusion sets can compose these registry entries.

### Benefits

- easier to tune parameters for new genomes/samples
- easier to compare exclusion strategies
- easier to generate provenance and summary tables
- easier to test exclusions independently

### Recommended early deliverables

1. Move slop/merge constants into config.
2. Add a machine-readable exclusion metadata table.
3. Emit a per-run exclusion provenance file.

## Phase 3 — Simplify the run report / UI

### Current problem

The report is data-rich but not yet optimized for decision-making.

### Target state

A single canonical report should answer:

1. What run was performed?
2. What assemblies/references/comparisons were used?
3. What benchmarksets were generated?
4. How much sequence was excluded, and why?
5. What do the evaluation summaries say?
6. Which outputs should a reviewer or releaser look at next?

### Recommended redesign

- one canonical `analysis.qmd`
- summary cards near the top
- one table per benchmarkset family
- explicit exclusion-impact section
- explicit release/provenance section
- links to raw outputs and logs

## Phase 4 — Harden packaging, archive, and release workflows

### Current problem

Execution and release concerns are mixed, and the packaging logic is rule-pattern driven but brittle.

### Target state

A run should produce a clearly defined package manifest before any upload/copy step occurs.

### Recommended changes

1. Separate:
   - run execution
   - report generation
   - package assembly
   - release transport
2. Generate:
   - package manifest
   - checksums
   - release README template
3. Make S3/local/FTP transport adapters consume the same package spec.

## Phase 5 — Modernize tests around behavior, not only target existence

### Current problem

Much of the existing coverage is workflow-smoke oriented.

### Target state

A balanced test suite should include:

- pure Python helper tests
- fixture-driven script tests
- rule-level smoke tests
- at least one minimal end-to-end integration test

### Recommended priorities

1. Add direct tests for `scripts/subtract_exclusions.py`.
2. Add direct tests for `scripts/get_sv_widen_coords.py`.
3. Add tests for config/profile resolution helpers in `rules/common.smk`-adjacent logic.
4. Keep one minimal chr21 smoke DAG in CI.

## Suggested initial ticket list

These are the best first issues to open/work in order:

1. fix pav callable-region intersection inputs
2. document developer bootstrap and local test workflow
3. make report source single-sourced
4. fix `run_defrabb` bucket/exclude/manifest release behavior
5. introduce named VCF processing profiles
6. externalize exclusion slop/merge constants
7. add exclusion provenance output
8. add direct tests for exclusion helper scripts
9. define package manifest for release contents
10. redesign top section of analysis report

## Recommended working style for future development

For each new feature or refactor:

1. update docs for the user-facing behavior
2. update or add a fixture-backed test
3. record release/archive impact if outputs change
4. keep config migrations explicit
5. validate on one minimal dataset before broad rollout

## Practical recommendation for the next session

The best next step is a small stabilization PR that does **not** change scientific behavior:

- fix the obvious bugs
- unify report source paths
- add bootstrap/developer docs
- add a short list of known limitations

After that, start the VCF-processing profile refactor.

