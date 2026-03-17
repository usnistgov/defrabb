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

- The bundled pytest suite expects local chr21 reference inputs under `.tests/integration/resources/references`.
- Those references are already tracked in this repository and are copied into temporary work directories by the generated rule tests.
- If you add new rule tests, prefer small fixtures under `.tests/unit/<rule>/data` and keep expected outputs under the matching `expected/` directory.

## useful references

- `README.md` for repository layout and top-level workflow usage
- `docs/running_defrabb_at_NIST.md` for current internal operational defaults
