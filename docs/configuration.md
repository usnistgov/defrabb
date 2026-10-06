# Configuration

A DeFrABB run is defined by two files:

- **`config/resources.yml`**: what is available. It holds input URLs, named
  parameter profiles, exclusion definitions, and compute resources. Validated by
  [`schema/resources-schema.yml`](../schema/resources-schema.yml).
- **an analyses table** (TSV): what to run. Each row defines one evaluation of
  one draft benchmark. Validated by
  [`schema/analyses-schema.yml`](../schema/analyses-schema.yml).

The analyses table defaults to `config/analyses.tsv`, a chr21 test. Choose
another table with `--config analyses=<path>`, or let `run_defrabb` pick up
`config/analyses_<RUNID>.tsv`.

## Analyses table

Each row combines four groups of columns. Rows that share a `vc_id` reuse the
same variant calls. Rows that share a `bench_id` reuse the same draft benchmark.

**Variant calling** (identified by `vc_id`):

- `asm_id`: assembly key in `resources.yml` `assemblies`
- `ref`: reference key in `resources.yml` `references` (`GRCh37`, `GRCh38`,
  `CHM13v2.0`, `GRCh38_chr21`)
- `vc_cmd`: `dipcall` or `pav`
- `vc_param_id`: named caller profile (`_dipcall_params` or `_pav_config`)
- `vc_params`: dipcall arguments used when `vc_param_id` does not match a
  profile (a matching profile takes precedence); ignored for PAV

**Draft benchmark** (identified by `bench_id`):

- `bench_type`: `smvar` (small variants) or `stvar` (structural variants)
- `bench_vcf_processing`: named profile from `vcf_processing_profiles`, or
  `none`
- `bench_bed_processing`: BED post-processing (`exclude` to apply the exclusion
  set)
- `exclusion_set`: named list of exclusions from `exclusion_set`, or `none`
- `exclusion_profile` (optional): named slop/merge profile from
  `_exclusion_profiles` (default `standard`)

**Evaluation** (identified by `eval_id`):

- `eval_cmd`: `happy`, `truvari`, `truvari_refine`, or `unhappy` (no evaluation)
- `eval_comp_id`: comparison callset key in `resources.yml` `comparisons`
- `eval_comp_id_is_truth`: `true` uses the comparison as truth and the draft as
  query. `false` reverses them.
- `eval_truth_regions`, `eval_target_regions`: whether to pass truth regions and
  target regions to the benchmarking tool
- `eval_params`: `default`, or a named Truvari profile from
  `truvari_bench_params`

Example row (one column per line) from
`config/analyses_20260928_v0.023_HG002Q100v1.2.tsv`. It defines a GRCh38
small-variant draft benchmark from HG002 Q100 v1.2 with dipcall, evaluated
against v5.0q with hap.py:

```text
eval_id                 GRCh38_HG002_T~T2TQ100v1.2_Q~v5.0q-smvar_TR~v5.0q
bench_id                GRCh38_HG002-T2TQ100v1.2-dipz2k_smvar-excluded
eval_cmd                happy
eval_params             default
eval_comp_id            v5.0q-smvar
eval_comp_id_is_truth   FALSE
eval_truth_regions      TRUE
eval_target_regions     TRUE
vc_id                   GRCh38_HG002-T2TQ100v1.2-dipz2k
bench_type              smvar
bench_vcf_processing    xy_trf
bench_bed_processing    exclude
exclusion_set           HG002Q100smvarv0.023
asm_id                  HG2-T2TQ100-V1.2
ref                     GRCh38
vc_cmd                  dipcall
vc_param_id             z2k
vc_params               -z200000,10000
```

Full tables for past production runs are in `config/analyses_*.tsv`. The table
that produced HG002 v5.0q is
[`config/analyses_20250117_v0.020_HG002Q100v1.1.tsv`](../config/analyses_20250117_v0.020_HG002Q100v1.1.tsv),
to be read with DeFrABB v0.020.

## resources.yml

Main sections:

- `references`: per reference, the FASTA URL, PAR BED, TR annotation database,
  exclusion BED URLs, and the GIAB stratifications tarball
- `assemblies`: per assembly, the maternal and paternal FASTA URLs, `is_male`,
  and sample ID
- `comparisons`: per reference, comparison callset VCF, BED, and index URLs
- `exclusion_set`: named lists of exclusion IDs
- `exclusion_slop_regions`, `exclusion_slopmerge_regions`,
  `exclusion_asm_intersect`, `exclusion_asm_agnostic`, `exclusion_ref_agnostic`:
  how each exclusion ID is sourced and processed (see the
  [exclusion system guide](exclusion_system_guide.md))
- `_exclusion_params`, `_exclusion_profiles`: slop and merge distances
- `vcf_processing_profiles`, `_vcf_processing_params`: VCF processing step
  sequences and annotation thresholds
- `_dipcall_params`, `_pav_config`, `truvari_bench_params`: named caller and
  evaluation profiles
- `_*_threads`, `_*_jobs`, `_*_mem`: compute resources (below)

## Compute resources

Memory use depends on the thread and job settings, so the two have to be
configured together. Memory settings are reservations (`mem_mb`) that Snakemake
counts against the global `--resources mem_mb=<budget>`. They do not limit the
process itself. Reported peaks were measured on whole-genome HG002 Q100 v1.1
with the default settings below.

- **dipcall**
  - Settings: `_dipcall_jobs` parallel alignment jobs (default 4), each with
    `_dipcall_threads` threads (default 5); 20 threads in total.
  - Memory: reserved as `_dipcall_jobs` × `_dipcall_mem` (4 × 40 GB). Each
    parallel minimap2 job holds its own reference index and alignments, so peak
    memory grows with `_dipcall_jobs`.
  - Measured peak: 102–116 GB at 4 jobs (highest for CHM13). Reducing
    `_dipcall_jobs` lowers peak memory and increases runtime.
- **PAV**
  - Settings: `_pav_threads` (default 24).
  - Memory: reserved as `_pav_mem` (80 GB). PAV runs alignment and calling steps
    in parallel across the given cores, so peak memory grows with the thread
    count.
  - Measured peak: about 65 GB at 24 threads (GRCh38, about 11.5 h wall time).
- **hap.py**
  - Settings: `_happy_threads` (default 12).
  - Memory: reserved as `_happy_mem` (160 GB). hap.py processes genome chunks in
    parallel, one per thread, so peak memory scales with threads.
  - Measured peak: 95–151 GB at 12 threads (earlier whole-genome HG002 and HG008
    runs). Reduce `_happy_threads` on smaller machines.
- **Truvari**
  - Settings: `_truvari_refine_threads` (24) sets threads for both
    `truvari bench` and `refine`. `_truvari_anno_threads` (8) applies to
    `truvari anno trf`; `anno repmask` uses a fixed 5 threads.
  - Memory: `_truvari_mem` (32 GB) is reserved for `truvari bench`, and
    `_truvari_anno_mem` (16 GB) for the `trf`, `svinfo`, and `lcr` annotations.
    `refine` (MAFFT) has no reservation, so leave headroom in the global budget
    when it runs.

If you change a thread or job setting, scale the matching memory reservation by
roughly the same factor. Always pass `--resources mem_mb=<budget>` for
whole-genome runs, so Snakemake does not start several memory-heavy jobs at
once. `run_defrabb` sets the budget to 80% of system memory by default.

## Run configuration and provenance

Each production run currently commits its analyses table to `config/` with the
name `analyses_YYYYMMDD_v#.###_<id>.tsv`. `run_defrabb` records the git state,
the conda environment, and the run log. GIAB benchmark releases include the
analyses table, the resource configuration, and the Snakemake report in a
`defrabb_files/` folder. Run configs are planned to move from the codebase to
per-run records.
