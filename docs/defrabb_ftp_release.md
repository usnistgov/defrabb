# Releasing a Draft Benchmark to FTP

Internal NIST notes for staging a DeFrABB run for FTP release. This is a manual packaging checklist, not an automated public release workflow.

## Scope

Use this after a run has completed and the final outputs have been reviewed. The current wrapper outputs expected here are:

- `archive.tar.gz`
- `snakemake_report_<RUNID>.zip`
- `environment.yml`
- `run.log`
- `results/`
- `resources/`

If a run-specific `run_README.md` was prepared, include it with the release package.

## Suggested staging layout

```sh
mkdir -p ftp_release/{benchmark_files,defrabb_files,dipcall_output}
```

## Benchmark files

Stage benchmark VCFs, indexes, benchmark beds, and supporting exclusion summaries from the run directory:

```sh
rsync -rv \
  --exclude="*bench-vars*" \
  --include="*.vcf.gz" \
  --include="*.vcf.gz.tbi" \
  --include="*.benchmark.bed" \
  --include="*.exclusion_intersection_summary.csv" \
  --exclude="*" \
  --no-relative --no-R --no-d \
  results/draft_benchmarksets/**/* \
  ftp_release/benchmark_files/
```

If release naming needs cleanup for a specific milestone, rename staged files after copying rather than modifying run outputs in place.

## DeFrABB provenance files

Stage workflow provenance and run documentation:

```sh
cp \
  archive.tar.gz \
  environment.yml \
  snakemake_report_*.zip \
  run.log \
  ftp_release/defrabb_files/
```

Also stage the exact config files used for the run. Prefer copies from the run directory if present; otherwise copy them from the checked-out repo used to launch the run:

```sh
cp config/analyses_<RUNID>.tsv config/resources.yml ftp_release/defrabb_files/
```

If present, also include:

```sh
cp run_README.md data_manifest.tsv ftp_release/defrabb_files/
```

## Dipcall outputs

Stage the core dipcall deliverables:

```sh
cp \
  results/asm_varcalls/*/*{dip.bed,dip.vcf.gz,dip.vcf.gz.tbi,hap1.bam,hap1.bam.bai,hap2.bam,hap2.bam.bai} \
  ftp_release/dipcall_output/
```

## Release README

Start from the previous public release README when available, then update:

- release name and date
- DeFrABB version or tag
- run ID
- assemblies, references, and comparison callsets used
- notable exclusions or evaluation choices
- directory contents and file descriptions

## Checksums

Generate checksums from inside the staged release directory:

```sh
cd ftp_release
find . -type f -exec md5sum '{}' \; > checksum.md5
```

## Upload

Upload using the current NIST-approved FTP transfer process. If the local Aspera helper is still the approved path, a typical command is:

```sh
bash ~/projects_ndo/giab-utils/handy_scripts/aspera_scripts/ftp_sra_transfer.sh \
  -m upload -r ftp \
  -d <FTP_RELEASE_DIR>/ \
  -l ftp_release
```

Replace `<FTP_RELEASE_DIR>` with the destination directory name for the release.
