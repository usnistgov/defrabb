# Configuration files

- `resources.yml`: input URLs, named parameter profiles, exclusion definitions,
  and compute resources
- `analyses.tsv`: default chr21 test analyses table
- `analyses_YYYYMMDD_v#.###_<id>.tsv`: analyses tables for past production runs
  (do not edit)
- `sweeps/`: parameter sweep configs for `scripts/generate_param_sweep.py`
- `release.json`: NIST-specific release defaults for `run_defrabb`

See [docs/configuration.md](../docs/configuration.md) for field definitions and
compute resource settings. The schemas are in [`schema/`](../schema/).
