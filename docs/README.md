# DeFrABB documentation

Start with the [repository README](../README.md) for an overview of what DeFrABB
does, how to run it, and which benchmarks it has produced.

## Method

- [Method overview](method-overview.md): the three pipeline components, VCF
  processing steps, exclusion types, and evaluation, with the reason for each
- [Exclusion system guide](exclusion_system_guide.md): how exclusion regions are
  sourced, buffered, merged, and subtracted, and how exclusion sets are
  configured
- [Architecture diagram](architecture-diagram.md): rule-level view of the
  workflow, from configuration parsing through evaluation and reporting

## Using the pipeline

- [Configuration](configuration.md): analyses-table columns, `resources.yml`
  sections, and compute resources (thread and memory settings)
- [Outputs](outputs.md): output directory layout, draft benchmark files, how to
  use a benchmark, and a glossary
- [Parameter optimization](parameter-optimization.md) and the
  [walkthrough](examples/parameter-optimization-walkthrough.md): generating and
  scoring parameter sweeps
- [Developer quickstart](developer-quickstart.md): local setup, tests, and
  formatting checks

## Release notes

See [CHANGELOG](../CHANGELOG).
