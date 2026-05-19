# CI Firewall Issue Diagnosis and Fix

## Problem Discovered

**Root cause:** CI runner cannot access conda package channels due to NIST firewall blocking `anaconda.org` and related domains.

**Evidence:**
1. `.condarc` file exists in repo with firewall workarounds BUT is in `.gitignore` (line 9)
2. User's `~/.condarc` contains working proxy and mirror configuration
3. CI has no access to this configuration → conda tries to reach blocked domains → 22-minute timeout

## Solution Implemented

### Created `.condarc.ci` (committed to repo)

A CI-specific conda configuration that uses **prefix.dev mirrors** instead of blocked anaconda.org:

```yaml
# Uses prefix.dev mirrors for conda-forge and bioconda
custom_multichannels:
  conda-forge:
    - https://prefix.dev/conda-forge
  bioconda:
    - https://prefix.dev/bioconda

# Disable SSL verification (for institutional proxies)
ssl_verify: false

# Block anaconda channel
disallow:
  - anaconda
```

### Updated `.gitlab-ci.yml`

1. **Set CONDARC environment variable** to point to `.condarc.ci`
2. **Added config verification** in before_script to confirm mirrors are active

```yaml
variables:
  CONDARC: "${CI_PROJECT_DIR}/.condarc.ci"

before_script:
  - conda config --show-sources  # Verify .condarc.ci is loaded
```

## Why This Approach

### ✅ Chosen: Separate CI config file

- **`.condarc`** → stays in `.gitignore` (user-specific, may contain credentials)
- **`.condarc.ci`** → committed to repo (CI-specific, no secrets)
- **Clean separation** between user environment and CI environment

### ❌ Rejected: Commit .condarc to repo

- Would expose user-specific proxy settings
- Mixes user config with CI config

### ❌ Rejected: Add channels to CI YAML directly

- Complex channel configuration doesn't fit well in YAML variables
- Less maintainable than dedicated config file

## Next Steps

### If CI still fails with network errors:

1. **Check if prefix.dev is accessible from runner:**
   ```bash
   # On the CI runner machine
   curl -I https://prefix.dev/conda-forge
   ```

2. **If blocked, enable proxy in .condarc.ci:**
   Edit `.condarc.ci` and uncomment proxy_servers section with runner's proxy

3. **Alternative: Use runner's system .condarc:**
   Remove `CONDARC` variable from `.gitlab-ci.yml` and instead configure
   `/home/gitlab-runner/.condarc` on the runner machine directly

### Expected outcome after this fix:

- **Without firewall issue:** CI should now complete successfully
- **First run:** ~20 min (builds cache + downloads from prefix.dev)
- **Subsequent runs:** ~2-5 min (uses cache)

## Monitoring

Check the CI logs for these indicators:

✅ **Success indicators:**
```
Conda config (using /builds/project/.condarc.ci):
==> /builds/project/.condarc.ci <==
custom_multichannels:
  conda-forge: ('https://prefix.dev/conda-forge',)
```

❌ **Failure indicators:**
```
CondaHTTPError: HTTP 000 CONNECTION FAILED for url <https://conda.anaconda.org/...>
Timeout or connection refused
```

If you see failure indicators, the runner may need additional proxy configuration.
