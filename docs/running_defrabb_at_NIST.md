# Running DeFrABB on NIST Compute Resources (NIST Internal Use Only)

This document describes the current internal workflow for launching and releasing DeFrABB runs on NIST-managed systems. It assumes access to NIST storage, AWS credentials where needed, and the internal Git hosting used by the team.

## Current scope

The repository no longer includes the older helper scripts previously used for CTCMS or Tibanna submission. Use the checked-in `run_defrabb` wrapper or run `snakemake` directly. If site-specific cluster wrappers still exist locally, treat them as separate operational tooling rather than part of this repository.

## Prerequisites

Before starting a run, make sure you have:

- a Linux environment with network access for Conda environment creation and remote resource downloads
- a working Conda or mamba installation
- an active environment containing at least:
  - `snakemake` 8.30+
  - `apptainer`
  - `boto3`
  - `rsync`
- access to the internal Git repo and, if releasing, access to the configured NAS and S3 destinations

There is no single bootstrap environment file in this repository. Use a team-managed environment if one already exists, or build one manually with the tools above.

## Recommended run layout

Use a dedicated directory per run:

```bash
RUNID=20250708_v0.021_HG008TN
mkdir -p /defrabb_runs/runs_in_progress/${RUNID}
cd /defrabb_runs/runs_in_progress/${RUNID}
git clone git@gitlab.nist.gov:bbd-human-genomics/defrabb.git .
git checkout <branch-or-tag>
tmux new-session -s ${RUNID}
```

Run IDs should follow `YYYYMMDD_v#.###_brief-id`.

## Preparing configuration

1. Start from an existing checked-in example such as `config/analyses_20250708_v0.021_HG008TN.tsv`.
2. Update `config/resources.yml` for the run.
3. Validate field expectations against:
   - `schema/analyses-schema.yml`
   - `schema/resources-schema.yml`

If you omit `-a/--analyses`, `run_defrabb` expects `config/analyses_<RUNID>.tsv`.

## Running the workflow

For most internal runs, use the wrapper:

```bash
./run_defrabb -r ${RUNID} -o ../ -s pipe
```

The wrapper:

- creates `../${RUNID}/`
- logs git status and the active Conda environment
- runs Snakemake with `--use-conda --use-apptainer`
- forwards any extra arguments to Snakemake

Common follow-up steps:

```bash
./run_defrabb -r ${RUNID} -o ../ -s report
./run_defrabb -r ${RUNID} -o ../ -s archive
./run_defrabb -r ${RUNID} -o ../ -s release
```

Or run the full sequence in one command:

```bash
./run_defrabb -r ${RUNID} -o ../ -s all
```

## Internal release defaults

Current defaults are NIST-specific:

- run output root: `/defrabb_runs/runs_in_progress/`
- local archive destination: `/mnt/bbdhg-nas/analysis/defrabb-runs/`
- S3 bucket/path: `giab-data/defrabb_runs`

Release behavior is controlled by `config/release.json`. Override defaults with `--archive_dir`, `--s3_bucket`, `--s3_path`, or `--release_type` when needed.

## Monitoring and troubleshooting

- Use `tmux` for long-running sessions.
- Review `run.log` in the run directory for wrapper-level status.
- Review `.snakemake/log/` inside the run directory for workflow-level logs.
- If Conda environment creation fails, rerun on a host with internet access or pre-create environments first.
- `run_defrabb` records uncommitted Git changes at startup; clean working trees are preferred for release runs.

## Outputs to expect

Typical wrapper outputs under `OUTDIR/RUNID/`:

- `archive.tar.gz`
- `snakemake_report_<RUNID>.zip`
- `environment.yml`
- `run.log`
- `benchmark/`
- `logs/`
- `resources/`
- `results/`

For downstream release or sharing steps, see `docs/defrabb_ftp_release.md` if applicable to the run.
