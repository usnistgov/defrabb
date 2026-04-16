# developer quickstart

This note captures the smallest supported developer setup for local validation and small workflow changes.

## environment assumptions

- Linux environment
- Snakemake available in your active base environment
- Conda or mamba available for `--use-conda`
- Apptainer available for rules that use containers
- Python dependencies needed by `run_defrabb` available in the active environment if you use the wrapper (notably `boto3`)

DeFrABB does not ship a single bootstrap environment for all development tasks. The expected workflow is to activate a local Snakemake environment, then let Snakemake create per-rule environments from `envs/*.yml`.

## basic validation commands

Run the unit and helper tests:

```sh
pytest .tests
```

Run a local smoke test of the workflow:

```sh
snakemake --use-conda --use-apptainer --cores 1 --forceall
```

For a lighter local run, the default checked-in analysis can be executed directly with:

```sh
snakemake --use-conda --use-apptainer --cores 1
```

## running the wrapper locally

Example wrapper invocation:

```sh
./run_defrabb -r 20240513_v0.016_test -a config/analyses.tsv -s pipe
```

Useful overrides outside the NIST environment include:

- `--outdir`
- `--archive_dir`
- `--s3_bucket`
- `--s3_path`
- `--release_type`

## test data caveats

- The bundled pytest suite currently exercises Python helpers and release-flow logic only (see `.tests/unit/test_stabilization_pr.py`). The Snakemake-7-era auto-generated rule tests have been removed.
- The chr21 reference inputs under `.tests/integration/resources/references` are still tracked and used by the end-to-end smoke test (`snakemake --use-conda --use-apptainer --cores 1`).
- New tests should be plain Python unit tests targeting helpers in `scripts/` or `rules/common.smk`, or structural checks via `snakemake -n`, `--lint`, `--dag`, `--rulegraph`, `--filegraph`.

## pre-commit hooks

A `.pre-commit-config.yaml` is provided that runs `black`, `snakefmt`, and a few generic hygiene hooks before each commit. To enable locally:

```sh
pip install pre-commit
pre-commit install
```

Hook versions are pinned to match the CI toolchain.

## useful references

- `README.md` for repository layout and top-level workflow usage
- `docs/running_defrabb_at_NIST.md` for current internal operational defaults
