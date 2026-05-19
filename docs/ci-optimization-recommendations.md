# CI Optimization Recommendations

## Current Issues

1. **22-minute runtime before failure** - conda environment setup is the bottleneck
2. **No conda caching** - environments rebuilt from scratch every run
3. **Slow conda solver** - default conda solver much slower than mamba/libmamba
4. **Runner reconfiguration issues** - CI broken since 2026-04-21 after OS reinstall

## Implemented Fixes (in this commit)

### 1. Added conda environment caching
```yaml
cache:
  key: conda-envs-$CI_COMMIT_REF_SLUG
  paths:
    - .snakemake/conda/
```

**Impact:** After first successful run, conda envs cached and reused (~5-10x speedup)

### 2. Enabled libmamba solver
```yaml
variables:
  CONDA_SOLVER: libmamba
```
**Impact:** Faster dependency resolution (2-3x speedup on env creation)

### 3. Added mamba frontend to snakemake
```bash
snakemake --conda-frontend mamba ...
```
**Impact:** Uses mamba instead of conda for environment operations

**Expected total speedup:** First run may still be slow, but subsequent runs should drop from 22+ minutes to ~2-5 minutes for conda setup.

## Additional Optimization Options

### Option A: Pre-build conda environments on runner (Recommended)

Pre-create the two required environments on the GitLab runner:

```bash
# On the runner machine
cd /home/gitlab-runner/defrabb-conda-envs
snakemake --use-conda --conda-create-envs-only --conda-frontend mamba --cores 1
```

Then update `.gitlab-ci.yml`:
```yaml
run-test:
  script:
    - snakemake -p --use-conda --conda-frontend mamba --conda-prefix /home/gitlab-runner/defrabb-conda-envs/.snakemake/conda --use-apptainer --cores 1 --verbose --config "_happy_threads=1"
```

**Impact:** Eliminates conda env creation from CI entirely (~22 min → ~30 sec startup)

### Option B: Use containerized conda environments

Build apptainer/singularity containers with conda envs pre-installed, push to registry, pull in CI.

**Impact:** Most robust, but requires container registry and more setup

### Option C: Reduce test scope (if appropriate)

Current `config/analyses.tsv` triggers 49 jobs. Could create an even smaller CI-specific config:
- Single reference (chr21 only)
- Single evaluation
- Minimal VCF processing

**Impact:** Fewer rules = fewer conda envs needed

## Runner Configuration Checks

After OS reinstall, verify on the runner:

1. **Conda/mamba installed and configured:**
   ```bash
   conda config --show solver  # Should show 'libmamba' or be settable
   mamba --version  # Should be available
   ```

2. **Apptainer working:**
   ```bash
   apptainer --version
   apptainer run docker://alpine echo "test"
   ```

3. **defrabb-ci environment has mamba:**
   ```bash
   conda activate defrabb-ci
   mamba --version
   ```

4. **Disk space for conda caching:**
   ```bash
   df -h /home/gitlab-runner
   # Need at least 5GB free for conda cache
   ```

## Monitoring

After applying these changes:

1. **First run:** Will still be slow (~20 min) as it builds cache
2. **Second run:** Should be much faster (~2-5 min for conda setup)
3. **Check cache:** Look for "Extracting cache" in CI logs

If still slow after caching is established, consider Options A or B above.

## Expected Timeline

- **Immediate (this commit):** Caching + mamba solver enabled
- **First CI run:** ~20 min (building cache)
- **Subsequent runs:** ~2-5 min for conda setup + test runtime
- **With Option A:** ~30 sec startup + test runtime
