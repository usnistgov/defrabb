# Snakemake Archive Network Timeout

## Issue
`./run_defrabb archive` fails when Snakemake attempts to download conda package metadata from conda.anaconda.org to include in the archive tarball.

**Error:**
```
requests.exceptions.ConnectTimeout: HTTPSConnectionPool(host='conda.anaconda.org', port=443): 
Max retries exceeded with url: /conda-forge/linux-64/zstd-1.5.7-hb78ec9c_6.conda 
(Caused by ConnectTimeoutError(...'Connection to conda.anaconda.org timed out. (connect timeout=None)'))
```

## Root Cause
Snakemake's `--archive` command packages conda environments by downloading package metadata for each conda package used in the workflow. This operation can fail due to:
- Network instability
- Firewall restrictions blocking conda.anaconda.org
- Rate limiting on conda repositories
- Temporary conda.anaconda.org outages

This is an internal Snakemake operation and not controlled by DeFrABB pipeline code.

## Workarounds

### Option 1: Retry the Archive Command
Network timeouts are often intermittent. Simply retry:
```bash
./run_defrabb archive -r <RUNID>
```

### Option 2: Use Alternative Conda Mirrors
If conda.anaconda.org is consistently blocked/slow, configure conda to use alternative mirrors (e.g., prefix.dev):
```bash
# Create/edit ~/.condarc
cat >> ~/.condarc << 'EOF'
channels:
  - https://prefix.dev/conda-forge
  - conda-forge
  - defaults
EOF
```

### Option 3: Skip Archive Step
If archiving is not immediately needed, skip to release step:
```bash
./run_defrabb release -r <RUNID> [release options]
```

The archive is primarily useful for reproducibility verification. For internal NIST releases, the release step handles deployment independently.

## Related Issues
- Similar conda network issues documented in `condarc.ci` workaround for CI runners

## Status
Open - Snakemake internal behavior; no direct fix in DeFrABB code.

**Date:** 2026-07-21
