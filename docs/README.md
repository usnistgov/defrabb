# DeFrABB documentation

Start with the [repository README](../README.md) for an overview of what
DeFrABB does, how to run it, and which benchmarks it has produced.

## Method

- [Exclusion system guide](exclusion_system_guide.md): how benchmark regions
  are derived by subtracting exclusion regions from diploid-assembled regions,
  and how exclusion sets are configured
- [Architecture diagram](architecture-diagram.md): rule-level view of the
  workflow, from configuration parsing through evaluation and reporting
- [Benchmark validation framework](benchmark-validation-framework.md): how
  differences between benchmark versions are classified and accepted or rejected

## Using the pipeline

- [Configuration overview](../config/README.md): `resources.yml` and analyses
  tables (field definitions are in [`schema/`](../schema/))
- [Parameter optimization](parameter-optimization.md) and the
  [walkthrough](examples/parameter-optimization-walkthrough.md): generating and
  scoring parameter sweeps across variant-caller, exclusion, and VCF-processing
  profiles
- [Developer quickstart](developer-quickstart.md): local setup, tests, and
  formatting checks

## Validation reports

- [v0.022 vs v5.0q](validation/v0.022-vs-v5.0q.md)

## Investigations and known issues

Write-ups of pipeline-specific failures and their workarounds. Several are
referenced from code comments to explain why a rule behaves as it does.

| Topic                                                                    | Write-up                                                                                                                                                                           |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PAV crash on FIPS-enabled hosts                                          | [run_pav_fips_selftest.md](issues/run_pav_fips_selftest.md)                                                                                                                        |
| dipcall OOM, PAV failures, `truvari anno trf` stalls on large insertions | [run_pav_run_dipcall_failures.md](issues/run_pav_run_dipcall_failures.md)                                                                                                          |
| Truvari refine v5.4.0 `samtools faidx` failure                           | [truvari-refine-v5.4.0-bug.md](issues/truvari-refine-v5.4.0-bug.md), [truvari-refine-primary-chr-filter.md](issues/truvari-refine-primary-chr-filter.md)                           |
| hap.py phasing mismatch against v5.0q                                    | [happy-phasing-mismatch-v5q.md](issues/happy-phasing-mismatch-v5q.md)                                                                                                              |
| Genome-specific stratifications (complex variants)                       | [stratification-59-173-design.md](issues/stratification-59-173-design.md), [genome-specific-geno2haplo-haploid-segfault.md](issues/genome-specific-geno2haplo-haploid-segfault.md) |
| Reuse of variant-caller outputs across analyses                          | [varcall-caching-investigation.md](issues/varcall-caching-investigation.md)                                                                                                        |
| `snakemake --archive` network timeout                                    | [snakemake-archive-network-timeout.md](issues/snakemake-archive-network-timeout.md)                                                                                                |

## Release notes

See [CHANGELOG](../CHANGELOG).
