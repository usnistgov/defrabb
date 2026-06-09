# Stratification features #59 / #173 — design & status

_Last updated: 2026-06-09_

Two related hap.py stratification features that share machinery. As of
2026-06-09 the **genome-specific complex/overlapping-variant stratifications are
implemented, validated, and wired (opt-in) into hap.py.** The SV/CNV stratum
(GS rules 7-10) is specified below but not yet wired.

## What is implemented (2026-06-09)

### Shared primitives
- `scripts/build_stratification_tsv.py` — builds a hap.py `--stratification`
  TSV (`<name>\t<bed>`, paths relative to the TSV) and, via
  `combine_strat_tsvs()`, merges genome-specific rows into the GIAB strat TSV
  with correct relative-path rewriting. Unit-tested
  (`.tests/unit/test_build_stratification_tsv.py`).
- `scripts/genome_specific_strats.py` — dependency-light, unit-tested port of
  the GIAB genome-stratifications **GenomeSpecific** pipeline's complex/
  overlapping-variant strata (GS rules 1-6). Classifies a `vcfgeno2haplo -w 10`
  benchmark VCF into `comphetsnp10bp`, `comphetindel10bp`, `complexindel10bp`,
  `snpswithin10bp`, and (from the raw VCF) `othercomplexwithin10bp`, with
  in-module interval merge/subtract. `.tests/unit/test_genome_specific_strats.py`.

### Pipeline rules — `rules/stratifications_genome_specific.smk`
- `genome_specific_geno2haplo` — `vcfgeno2haplo -w 10` on the draft benchmark VCF
  (`envs/genome_strats.yml`: vcflib + bedtools).
- `genome_specific_complex_strats` — runs the classifier → 5 strat BEDs.
- `genome_specific_strat_tsv` — assembles the per-benchmark genome-specific
  stratification TSV.

### hap.py wiring (OPT-IN)
- Config flag `genome_specific_strats` (schema: boolean, default false).
- When **true** and `bench_type == smvar`, `rules/helpers_eval.smk`
  (`get_happy_inputs_inner`) adds the genome-specific strat TSV as a `run_happy`
  input, and `scripts/run_happy.py` merges it into the GIAB strat TSV via
  `combine_strat_tsvs()` before calling hap.py.
- When **false** (default) the input is absent and behavior is byte-for-byte
  unchanged — pinned/released evaluations are unaffected.

### Validation performed
- **Classifier is byte-identical to the upstream GIAB awk/bedtools pipeline** for
  all five strata, run on the real HG002 chr21 dipcall benchmark VCF
  (`results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.dip.vcf.gz`,
  2653 records / 32 compound hets). Output beds well-formed (all chr21,
  0 ≤ start < end).
- `envs/genome_strats.yml` solves (vcflib 1.0.15 + bedtools 2.31.1; bedtools left
  unpinned to avoid a libzlib conflict with vcflib).
- Snakemake DAG resolves for the strat chain and for `run_happy` both with the
  flag off (unchanged) and on (genome-specific chain pulled in). hap.py was not
  executed end-to-end here (CI is dry-run; the full GRCh38 strat tarball + a real
  hap.py run are the maintainer's full-HG002 validation step).

## Remaining: #59 SV/CNV stratum (GS rules 7-10) — NOT yet wired

JZ's spec (issue thread): add the genome-specific strats **except** the external
SV/CNV-callset-dependent ones, and for the SV+CNV bed use "the complement of the
intermediate small-variant benchmark bed (after excluding partially-covered
repeats, gaps, and SVs, but **not** subsequent newer exclusions like dipcall bugs
and homopolymers>30bp)".

**defrabb wrinkle:** there is no incremental "intermediate bed" — `subtract_exclusions`
applies all exclusions in one shot. The intermediate bed must be constructed by
subtracting only the repeat/gap/SV exclusions from `baseline.bed`. Decided
category split (from the exclusion sets in `config/resources.yml`):

- **Subtract for the intermediate bed (repeat/gap/SV):** `segdups`,
  `tandem-repeats`, `satellites`, `gaps`, `flanks`, `svs-and-simple-repeats`,
  `consecutive-svs`, `self-discrep`, `pav-inv`, `TSPY2-segdups`, `VDJ`.
- **Do NOT subtract (dipcall-bug / error / homopolymer-class):**
  `dipcall-bugs-T2TACE`, `HG002Q100-pav-discrep-smvar`,
  `HG002Q100-pav-discrep-stvar`, `HG002Q100-pav-inversions`, `HG002Q100-errors`,
  `HG002Q100-mosaic`, `HG002Q100-delins-errors`.

Planned rules (mirror GS rules 7-10):
1. `intermediate_smvar_bed` = `baseline.bed` minus the repeat/gap/SV subset above.
2. `sv_cnv_bed` = genome complement of (1) → the "CNVsandSVs" stratum.
3. `complexandSVs` = union of the 5 complex strata + `sv_cnv_bed` (GS rule 8).
4. `complexandSVs_alldifficultregions` = union with the GIAB all-difficult-regions
   strat (GS rule 9) — **needs the version-specific all-difficult bed path inside
   the extracted GIAB tarball**.
5. `notin_complexandSVs_alldifficult` = genome minus (4) (GS rule 10).

Blockers for this part: (a) confirm the category split above with JZ; (b) the
all-difficult strat path is GIAB-strat-version specific and was not validatable
here without extracting the full GRCh38 tarball. (1)-(3) are pure defrabb bed
algebra and could be added next; (4)-(5) need the all-difficult path.

## How to enable / run
```sh
# generate just the genome-specific strat TSV for a benchmark:
snakemake --use-conda --cores 1 \
  results/draft_benchmarksets/<bench_id>/genome_specific_strats/<...>.genome_specific_strats.tsv

# run hap.py with genome-specific stratification:
snakemake --use-conda --cores 1 --config genome_specific_strats=True <happy target>
```
