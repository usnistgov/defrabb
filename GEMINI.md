# DeFrABB (Development Framework for Assembly-Based Benchmarks)

## Project Overview
DeFrABB is a Snakemake workflow developed by the NIST/GIAB (Genome in a Bottle) team. It is used to generate transparent, reproducible assembly-based small- and structural-variant benchmark sets. The pipeline orchestrates the downloading of references and assemblies, variant calling (e.g., using `dipcall`), benchmarking against comparisons (using tools like `hap.py` and `truvari`), and generating comprehensive reports.

**Key Technologies & Tools:**
- **Workflow Manager:** Snakemake (>= 8.30.0)
- **Languages:** Python (scripts, wrappers, tests), R, Bash
- **Environment Management:** Conda / Mamba, Apptainer (for containerized steps)
- **Testing:** Pytest
- **Domain Tools:** bcftools, bedtools, truvari, hap.py, rtgtools, dipcall

## Directory Structure
- `Snakefile`: The main entry point for the Snakemake workflow.
- `rules/`: Modular Snakemake rules (`.smk` files) for variant calling, evaluation, exclusions, and reporting.
- `scripts/`: Python, R, and shell scripts utilized by the Snakemake rules.
- `config/`: Configuration files, schemas, and analyses tables (`analyses_*.tsv`, `resources.yml`).
- `envs/`: Conda environment definitions (`.yml`) for individual rules.
- `.tests/`: Pytest suite for unit testing Python scripts and integration-testing rule components.
- `docs/`: Developer and operational documentation.
- `run_defrabb`: A Python wrapper script for orchestrating run execution, environment logging, and S3/NAS releases.

## Building and Running

### Basic Setup
DeFrABB does not use a monolithic environment. You must provide a base environment containing Snakemake, Conda/mamba, Apptainer, and core Python dependencies (like `boto3` if using `run_defrabb`).

### Local Validation (Smoke Test)
To run a lightweight local test of the workflow using the bundled configuration:
```bash
snakemake --use-conda --use-apptainer --cores 1
```
*(Append `--forceall` to force a complete run from scratch)*

### Using the Wrapper Script
For full workflow execution, environment recording, and potential releases, use `run_defrabb`. The required `RUNID` must follow the pattern `YYYYMMDD_v#.###_brief-id`.
```bash
./run_defrabb -r 20240513_v0.016_test -a config/analyses.tsv -s pipe
```
*Note: External users may need to override NIST-specific arguments like `--outdir`, `--archive_dir`, `--s3_bucket`, and `--release_type`.*

## Testing
Run the Python and Snakemake rule unit/helper tests via `pytest`:
```bash
pytest .tests
```
*Note: Integration tests use tiny, localized data (like `chr21`) located in `.tests/integration/resources/`.*

## Development Conventions
- **Code Formatting:** Python code is formatted with `black`. Ensure all scripts conform to `black` before committing (`black . --check` is enforced in CI).
- **Snakemake Linting:** Changes to the workflow must pass `snakemake --lint`.
- **Contributions & Milestones:** Follow the procedure in `CONTRIBUTING.md`. Create a milestone for new functionality, and ensure PRs are merged onto the `dev` branch initially. Release tags follow `v#.###` format.
- **Rule Environments:** Do not install dependencies globally. Add required tools to specific `.yml` files in the `envs/` directory so Snakemake can instantiate them locally per-rule.
- **New Test Data:** Place small fixture datasets for rule testing in `.tests/unit/<rule>/data` and their expected outputs under `.tests/unit/<rule>/expected/`.