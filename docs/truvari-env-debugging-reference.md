# Truvari Env Debugging Reference

> **Resolution (2026-06-08):** the `truvari anno remap` failure is fixed. The env
> (`envs/truvari_remap.yml`) is now `bioconda::truvari ==5.4.0` + a patched bwapy
> fork (`nate-d-olson/bwapy` @ `5b91fe0`, tag `0.1.5-nist`) installed via pip,
> plus a build toolchain (`c-compiler`, `make`, `zlib`). The fork drops the
> `pkg_resources` import (builds under setuptools≥80, no `--no-build-isolation`)
> and compiles with `-fno-finite-math-only` while rebuilding `bwa/libbwa.a`
> (fixes `undefined symbol: __log_finite` on glibc≥2.31). FIPS: the 5.4.0 envs
> run clean without `OPENSSL_CONF=/dev/null`; only the 4.3.0-pinned `trf` env
> still needs that guard. See GitLab #200 (resolved) and #198 (bioconda
> follow-up). The notes below are retained as historical debugging context.

Last-known-successful Truvari env state, captured from branch `postproc-HG008TN`
(HG008 v6.3N / v3.2T run, January–March 2026). Use this as the working baseline
when debugging the current Truvari anno failures (issue #197 and related env
problems noted in TODO #3).

The current master Truvari envs went through subsequent changes (commits
`45c1152`, `24e1b93`, `165fe1a`) that have not produced a successful end-to-end
Truvari rule pass; this file preserves the configuration that did.

## envs/truvari.yml (last-successful)

```yaml
channels:
  - bioconda
  - conda-forge
dependencies:
  - bcftools
  - mafft
  - pip
  - pip:
    - Truvari>=5.3
```

Notes:
- Truvari installed via **pip**, not conda. Conda Truvari builds have been the
  source of repeatmasker/bwapy version mismatches.
- `bcftools` and `mafft` left unpinned, resolved by conda.
- No conda `repeatmasker` or `trf` in this env — those are only in
  `truvari_remap.yml`.

## envs/truvari_remap.yml (last-successful)

```yaml
channels:
  - conda-forge
  - bioconda
dependencies:
  - bcftools==1.20
  - trf==4.09.1
  - mafft>=7.515
  - python
  - pywfa >=0.5.1
  - repeatmasker==4.1.5
  - edlib>=1.3.9
  - pip
  - pip:
    - Truvari==4.2.2
```

Notes:
- Truvari **4.2.2** via pip here (older than `truvari.yml`'s 5.3 — the rules
  call this env specifically for `truvari anno trf` where v4.2.2 worked).
- `bcftools==1.20`, `trf==4.09.1`, `repeatmasker==4.1.5` all explicitly pinned.
- This is the env where the `repeatmasker`/`bwapy`/`trf` library completeness
  matters most — pinning these versions kept it working.

## rules/bench_vcf_processing.smk — `run_truvari_anno_trf` (last-successful)

```python
rule run_truvari_anno_trf:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
        ref=get_ref_file,
        trdb=get_ref_trdb,
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.trfanno.vcf",
    log:
        "logs/truvari_anno_trf/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari_remap.yml"
    params:
        min_length=20,
    threads: 8
    shell:
        """
        truvari anno trf \
            -i {input.vcf} \
            -o {output.vcf} \
            -r {input.trdb} \
            -f {input.ref} \
            -t {threads} \
            --min-length {params.min_length} \
            -C 1 \
            --trf-params "3 7 7 80 5 40 500 -h -l0.25 -ngs" \
            -e trf &> {log}
        """
```

Notes vs current master:
- `threads: 8` (master: 5) — bumped during HG008 run for throughput.
- `-C 1` flag added — limits TRF chunking; may have addressed an OOM/runtime
  symptom during the HG008 run.
- `--trf-params` changed `-l 20` → `-l0.25 -ngs` (the `-ngs` flag was already
  present elsewhere in the original; the `-l` value sets max TR length, in
  millions of bp — **`-l 20` = 20 Mbp** vs **`-l0.25` = 250 kbp**, smaller is
  faster and bounds memory).

## How to use this reference

When debugging the current Truvari env failures:

1. **Diff current `envs/truvari.yml` against this baseline** to identify what
   changed since the last successful run. The pip-vs-conda swap (commit
   `45c1152`) is a candidate root cause.
2. **For `truvari_remap.yml`**, confirm that `bcftools`, `trf`,
   `repeatmasker`, `pywfa`, `edlib` are all pinned to compatible versions.
   The conda `repeatmasker==4.1.5` here is a known-good pin.
3. **For the `run_truvari_anno_trf` rule**, evaluate whether the `-C 1` flag
   and the `-l0.25` trf-params bound on TR length are needed for the failing
   anno step on master.
4. **Issue #197** ("Truvari anno trf error - invalid contig chrM for CHM13
   and MT for GRCh37") may not be addressed by this env state alone — the
   contig-name mismatch is a separate input-prep concern.

## Source

Branch: `postproc-HG008TN`
Files captured at branch tip (commit reachable from `origin/postproc-HG008TN`).
Branch closes GitLab issue #195 (HG008 v3.2T / v6.3N run).
