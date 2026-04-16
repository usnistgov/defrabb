# Repository Guidelines

## Project Structure & Module Organization
DeFrABB is a Snakemake-based benchmarking pipeline. `Snakefile` is the entry point and composes workflow modules from `rules/*.smk` (`asm-varcall`, `evaluation`, `exclusions`, `report`, etc.). Put Python, R, and shell helpers in `scripts/`; keep per-tool Conda environments in `envs/*.yml`. Run configuration lives in `config/` (`analyses*.tsv`, `resources.yml`, `release.json`) and is validated by `schema/*.yml`. Report templates are in `report/`. Small reference inputs are kept under `resources/`; generated outputs go to `results/`, `logs/`, and `benchmark/`. Rule and helper tests live in `.tests/unit/`.

## Build, Test, and Development Commands
- `snakemake --use-conda -j1` — run the default local test analysis.
- `snakemake --use-conda -j1 --forceall` — rerun the full workflow smoke test.
- `pytest .tests` — run generated Snakemake rule tests and Python helper tests.
- `./run_defrabb -r 20240513_v0.016_test -a config/analyses.tsv -s pipe` — run via the project wrapper.
- `snakefmt .` — format Snakemake files (used in CI).
- `black . --check` — verify Python formatting.

## Coding Style & Naming Conventions
Use 4-space indentation in Python and follow Black formatting. Prefer type hints and short docstrings for reusable Python functions. Use `snake_case` for Python helpers, TSV/YAML field names, and config keys. Keep Snakemake rules modular: one logical transformation per rule, helper functions in `rules/common.smk`, and tool environments referenced from `envs/*.yml`. Match existing file naming patterns such as `run_*.py`, `test_*.py`, and `analyses_YYYYMMDD_v0.###_*.tsv`.

## Testing Guidelines
Add or update `pytest` coverage in `.tests/unit/test_<feature>.py` for any rule or helper change. If a change touches core workflow wiring, also run `snakemake --use-conda -j1 --forceall`. Some tests expect local chr21 reference files under `.tests/integration/resources/references`; stage those before running the suite.

## Commit & Pull Request Guidelines
Follow the existing history: short, lowercase subjects describing the change, e.g. `fixing path` or `adding pav exclusion`. Keep commits focused. PRs should summarize workflow impact, note any config or resource-file changes, link the relevant issue or milestone, and mention validation performed (`pytest`, Snakemake smoke run, or `run_defrabb`). For release PRs, update `CHANGELOG` and use the milestone title format `milestone v#.###`.

## Configuration & Release Notes
`run_defrabb` defaults and `config/release.json` include NIST-specific paths. External contributors should prefer CLI overrides such as `--outdir`, `--archive_dir`, `--s3_bucket`, and `--s3_path` instead of hard-coding local environment changes.
