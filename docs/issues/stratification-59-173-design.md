# Stratification features #59 / #173 — design & status

_Last updated: 2026-06-08_

Two related hap.py stratification features. The shared building block —
`scripts/build_stratification_tsv.py` — is implemented and unit-tested
(`.tests/unit/test_build_stratification_tsv.py`). It produces a hap.py
`--stratification` TSV (`<name>\t<bed_path>`, paths relative to the TSV) from
named BED inputs. **Neither feature is wired into the live `run_happy` rule yet**
— that step changes scientific evaluation output and needs hap.py end-to-end
validation (the default CI is dry-run only).

## Current hap.py wiring (for reference)

- Rule `run_happy` (`rules/evaluation.smk`) → `scripts/run_happy.py`.
- Stratifications: `--stratification {strat_tsv}` where `strat_tsv` =
  `config["references"][ref_id]["stratifications"]["tsv"]` (the GIAB tarball TSV,
  extracted at runtime).
- Confidence/target regions: `-R {truth_regions}` / `-T {target_regions}`, wired
  from `eval_truth_regions` / `eval_target_regions` in the analyses table
  (`rules/helpers_eval.smk: get_happy_inputs_inner`). Currently the draft
  benchmark `.benchmark.bed` (post-exclusion) or comparison BED — **not** the raw
  `dip.bed`.

## #173 — hap.py stratifying with exclusions

Goal: use `dip.bed` (dipcall callable regions,
`results/asm_varcalls/{vc_id}/{ref}_{asm}_{vc_cmd}-{vc_param_id}.dip.bed`) as the
confidence region, then stratify by the per-benchmark exclusion BEDs
(`results/draft_benchmarksets/{bench_id}/exclusions/..._{region}.bed`) to quantify
each exclusion's impact on precision/recall.

Implementation plan:
1. New rule: build a strat TSV from the analysis's exclusion BEDs via
   `build_stratification_tsv.py` (one `--strat <region>=<bed>` per exclusion).
2. New/parameterized hap.py run that (a) sets `-R dip.bed` as confidence and
   (b) passes the exclusion strat TSV (possibly combined with the GIAB TSV).
3. Surface the per-stratum results through `scripts/compare_evaluations.py`
   (hap.py `Subset` column already carries stratum names).

**Blockers:** needs a hap.py validation run on real data; decision on whether
this is a new eval mode (e.g. an `eval_cmd` / `eval_params` variant) vs. changing
the existing happy run; confirm dip.bed vs `.benchmark.bed` as the confidence set.

## #59 — genome-specific stratifications from dipcall (JZ priority)

Goal: generate genome-specific stratification BEDs from dipcall output and pass
them to hap.py in addition to the GIAB strat TSV.

**Blocker (needs JZ input):** the issue does not specify *which* stratifications
to derive from dipcall. Candidate sources available today: `dip.bed` (diploid
callable), and for PAV the per-haplotype callable regions
(`callable_regions_h{1,2}_500.bed.gz`) / merged `diploid_regions.bed`. Once the
desired strata are defined, the mechanics mirror #173: build a combined strat TSV
with `build_stratification_tsv.py` and pass it to hap.py.

## Shared primitive

`scripts/build_stratification_tsv.py`
- `parse_strat_args(["name=bed", ...])`, `format_strat_rows(entries, out_dir,
  relative)`, `write_strat_tsv(entries, out, relative)`.
- CLI: `--strat NAME=BED ... -o OUT [--no-relative]`.
- Writes bed paths relative to the output TSV (hap.py resolves them that way);
  rejects duplicate stratum names.
