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

## full-pipeline regression test (after large changes)

CI and the chr21 smoke test do not exercise the real whole-genome variant-calling,
exclusion, stratification, and evaluation paths, so a change can pass CI and still
break on a full genome. After a large change, run the harness:

```sh
conda activate snakemake   # or ~/miniforge3/envs/snakemake
scripts/run_full_pipeline_test.sh [RUNID] [--dest-dir DIR] [--go] [--cores N]
```

What it does (regular clone-into-runid workflow):

1. Clones the current repo (committed HEAD) into `<DEST>/<RUNID>` (default `DEST=..`).
   Commit the code changes you want to test first — the clone copies committed HEAD
   only; the script warns if the source tree has uncommitted tracked changes.
2. In the clone, generates a dated analyses table (`config/analyses_<RUNID>.tsv`),
   adds the `HG002Q100stvarv0.022-selfdiscrep` exclusion set, enables
   `genome_specific_strats: true`, and commits these on a `run/<RUNID>` branch.
3. Validates and dry-runs in the clone, then prints the `run_defrabb run` command
   (default is setup-only; pass `--go` to launch the compute-intensive run).

The generated HG002 v1.1 table is built to exercise the current dev round's features:
the `stvar_v5` truvari profile (#194), genome-specific stratifications (#59/#173), the
`HG002Q100stvarv0.022` exclusion set (#188), and the symbolic-allele self-discrepancy
filter via a PAV stvar row (#192). When updating the harness for a future dev round,
edit the table transforms and the enabled flags in the script to match the new changes.

After the run completes, aggregate results from the clone with the generated
`compare_results_<RUNID>.sh` (wraps `scripts/compare_evaluations.py`). All run output
and commits live in the clone; `rm -rf <clone>` cleans up. This run is compute- and
time-intensive — it is meant for post-large-change verification, not routine use.

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
