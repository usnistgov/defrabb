# Configuration options

defrabb uses two configuration files

See `schema/analyses-schema.yml` and `schema/resources-schema.yml` for detailed descriptions and field formats requirements.

## resource.yaml

used to define:

- parameters, threads, and memory for compute intensive steps
- urls for remote files: diploid assemblies, genome reference files, stratifications, and callsets used to evaluate draft benchmark
- exclusion sets and how they are applied

## Analyses Tables

Provides run specific configurations

- input diploid assembly
- version of reference genome
- assembly-based variant caller and parameters
- vcf and bed processing including what exclusions to use
- benchmarking method and comparison callset used for initial evaluation