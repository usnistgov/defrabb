# `run_pav` FIPS self-test crash on FIPS-enabled hosts

_Diagnosed 2026-06-17 from the 20260616_v0.022 full-pipeline run._

## Symptom

`run_pav` fails; the inner PAV snakemake aborts immediately with:

```
crypto/fips/fips.c:154: OpenSSL internal error: FATAL FIPS SELFTEST FAILURE
```

The other containerized/conda rules (dipcall, truvari, hap.py) are unaffected.

## Root cause

The host runs in **global FIPS mode** (`/proc/sys/crypto/fips_enabled == 1`).
PAV's `pysam` is installed as a **pip wheel**, which bundles its own
RedHat-built `libcrypto` (`OpenSSL 1.1.1k FIPS`). On `import pysam` that bundled
library reads the host FIPS flag, auto-enables FIPS mode, runs an integrity
self-test, and aborts because the wheel's copied module is not the validated
FIPS canister.

Key points confirmed during diagnosis:

- Python's own `ssl` module (the container's *system* OpenSSL 1.1.1n, non-FIPS)
  imports fine — only `import pysam` crashes.
- `OPENSSL_CONF=/dev/null` does **not** help: the bundled lib reacts to `/proc`,
  not to the OpenSSL config file.
- The conda-based rules use conda's non-FIPS OpenSSL, so they never hit this.
- The same crash occurs with both `becklab/pav:latest` and a locally built PAV
  image — it is a property of the pip `pysam` wheel, not the image.

## Fix

`run_defrabb` (`setup_fips_apptainer_bind`) detects a FIPS host and, for the
duration of the run, binds a file containing `0` over
`/proc/sys/crypto/fips_enabled` inside every apptainer container via the
`APPTAINER_BIND` / `SINGULARITY_BIND` env vars (inherited by the snakemake
subprocess). The bundled OpenSSL then sees FIPS disabled and skips the self-test.

- **No image rebuild**; works with the stock `docker://becklab/pav:latest`
  (the PAV version the pipeline targets).
- **No-op on non-FIPS hosts** (the flag is absent or `0`), so CI / other users
  are unaffected.
- Scope is narrow and safe: PAV does local genomics (no cryptography); the bind
  only prevents a non-validated bundled lib from force-enabling a broken FIPS
  self-test inside the container. It does not alter the host or any real
  security posture.

The PAV container is also configurable via `_pav_container` in
`config/resources.yml` (default `docker://becklab/pav:latest`) should a
FIPS-free image ever be preferred.

## Verification

With the bind in place, `import pysam` succeeds and the full PAV inner DAG
(554 jobs for HG002 GRCh38) resolves on a dry-run inside `becklab/pav:latest`.

## Note: locally built PAV images

Building a local image does **not** by itself fix this (the wheel still bundles
the FIPS OpenSSL). A local build also needs `git submodule update --init
--recursive` in the PAV clone (otherwise `dep/svpop`/`svpoplib` and the nested
`ply`/`kanapy` deps are missing) and is typically a newer PAV than the pipeline
targets. The bind workaround on the stock image is preferred.
