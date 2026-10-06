# Outputs

All outputs are written under the run directory. File names follow the wildcard
patterns below. The IDs (`vc_id`, `bench_id`, `eval_id`) come from the analyses
table.

## Directory layout

- `results/asm_varcalls/{vc_id}/`: assembly-based variant calls
  - `{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz`: standardized calls (for PAV,
    PASS-only)
  - `{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.baseline.bed`: standardized diploid
    regions (from dipcall `*.dip.bed`, or from the intersection of the PAV
    haplotypes' callable regions)
  - `*.hap1.bam`, `*.hap2.bam` (dipcall): haplotype-to-reference alignments
  - `annotations/`: intermediate VCFs from each VCF processing step
- `results/draft_benchmarksets/{bench_id}/`: draft benchmark sets (see below)
- `results/evaluations/happy/{eval_id}_{bench_id}/`: hap.py output
  (`*.summary.csv`, `*.extended.csv`, annotated VCF)
- `results/evaluations/truvari/{eval_id}_{bench_id}/`: Truvari output
  (`summary.json`, TP/FP/FN VCFs; `refine` output when enabled)
- `results/report/`, `analysis.html`: summary statistics and the run's analysis
  report (rendered from `analysis.qmd`)
- `logs/`: per-rule logs
- `benchmark/`: per-rule runtime and memory (Snakemake `benchmark:` output)

`run_defrabb` also writes `run.log`, `environment.yml` (the active conda
environment), and, after `report` and `archive`, the Snakemake report and an
archive tarball.

## Draft benchmark set files

For each `{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}` prefix in
`results/draft_benchmarksets/{bench_id}/`:

- `*.vcf.gz`: benchmark variants, i.e. the processed and annotated
  assembly-based calls
- `*.benchmark.bed`: benchmark regions (diploid regions minus exclusions)
- `*_bench-vars.vcf.gz`: benchmark variants inside the benchmark regions (for
  `stvar`, only variants of 50 bp and larger)
- `*.exclusion_stats.txt`: bases removed by each exclusion
- `*.exclusion_provenance.yml`: the exclusion BEDs, slop and merge parameters,
  and processing applied
- `*_bench-vars_rtg_stats.txt`, bcftools stats: variant summary statistics

GIAB releases rename these files to `HG002_{ref}_{version}_{smvar|stvar}.vcf.gz`
and `.benchmark.bed`.

## Using a benchmark

Always use a benchmark VCF together with its `benchmark.bed`. Variants outside
the BED are not assessed and must not be counted as errors. For small variants,
use hap.py or vcfeval with the BED as truth regions, and stratify results with
GIAB genome stratifications. For SVs, use Truvari with the BED as
`--includebed`.

## Glossary

- **Assembly-based variant calling:** calling variants by aligning each
  haplotype of a diploid assembly to the reference (dipcall, PAV), rather than
  from sequencing reads.
- **Diploid regions:** regions where both haplotypes align 1:1 to the reference;
  the starting point for benchmark regions (`baseline.bed`).
- **Exclusion:** a set of regions removed from the diploid regions because
  benchmarking there would be unreliable. An **exclusion set** is a named list
  of exclusions.
- **Self-discrepancy:** an exclusion made by comparing a draft benchmark against
  itself with the benchmarking tool. Mismatches show representations the tool
  cannot handle.
- **smvar / stvar:** small-variant (generally under 50 bp) and
  structural-variant (50 bp and larger) benchmark types.
- **Stratification:** a set of genomic regions (e.g. tandem repeats, segmental
  duplications) used to report performance by genomic context.
- **PAR:** pseudoautosomal regions of chrX/chrY, which are diploid in males.
  Outside the PARs, male chrX/chrY genotypes are haploid.
- **vc_id, bench_id, eval_id:** analyses-table identifiers for a variant-call
  set, a draft benchmark, and an evaluation. Shared IDs mean shared outputs.
- **Truth / query:** in an evaluation, the callset treated as correct (truth)
  and the callset being assessed (query). DeFrABB can use either the draft or
  the comparison callset as truth.
- **Draft benchmark:** pipeline output before external evaluation and curation;
  GIAB releases are curated drafts.
