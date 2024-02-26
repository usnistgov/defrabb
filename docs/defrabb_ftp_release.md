# Releasing Draft Benchmark on FTP

Internal NIST notes for structuring, documenting, and submitting draft benchmark for sharing on ftp site.

1. Create directory and subdirectories staging files with name of ftp location
2. Files to include
   1. small and sv benchmark regions
   2. benchmark vcfs and index files (use small variant vcfs, vcf will be for both small and structural variants)
   3. exclusion stats files
   4. defrabb files - snakemake archive (tar.gz), snakemake report, env yaml, resources yaml, and analysis tsv
   5. dipcall files - dip.bed, dip.vcf.gz, dip.vcf.gz.tbi, {hap1,hap2}.bam, and {hap1,hap2}.bam.bai
3. Using previous release README as starting point for README.
4. Creating md5s
5. Uploading to ftp

## Directory setup

```sh
mkdir -p ftp_release/{defrabb_files,dipcall_output}
```

## Benchmark files

Copying benchmark vcfs, bed, and exclusion stats files. Only vcf for small variants as small and structural variant vcfs are the same files.

```sh
rsync -v \
   --exclude="*bench-vars*" \
   --exclude="*stvar*vcf.gz*" \
   --include="*vcf.gz*" \
   --include="*benchmark.bed" \
   --include="*exclusion_stats.txt" \
   --exclude="*" \
   --no-relative --no-R --no-d \
   results/draft_benchmarksets/**/* \
   ftp_release/ 
```

Removing "smvar" from vcf`_dipcall-z2k` file names

```sh
rename 's/_smvar//' ftp_release/*vcf.gz*
rename 's/_dipcall-z2k//' ftp_release/*
```

## Defrabb files

```sh
cp *.tar.gz \
   environment.yml \
   snakemake_report_*.html \
   ftp_release/defrabb_files/
```

Current snakemake archives do not include pipeline code like expected.

Need to manually copy config files, e.g.

```sh
cp \
   /defrabb_runs/code_base/defrabb/config/analyses_20240215_v0.015_HG002.tsv \
   /defrabb_runs/code_base/defrabb/config/resources.yml \
   ftp_release/defrabb_files/
```

## dipcall files

```sh
cp \
   results/asm_varcalls/*_HG002-T2TQ100v1.0-dipz2k/*{dip.bed,dip.vcf.gz,dip.vcf.gz.tbi,hap1.bam,hap1.bam.bai,hap2.bam,hap2.bam.bai} \
   ftp_release/dipcall_output/
```

## README

Using previous version of the README as a template revise and update as appropriate. 

```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/README.md`
```

## Generating MD5s

```sh
find -type f -exec md5sum '{}' \; > checksum.md5
```

## Uploading to FTP

Using NIST wrapper script for upload via aspera.
Example command

```sh
bash ~/projects_ndo/giab-utils/handy_scripts/aspera_scripts/ftp_sra_transfer.sh \
   -m upload -r ftp \
   -d NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/ \
   -l .
```
