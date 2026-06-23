# Assembly Variant Caller Intermediate Reuse Investigation

**Date:** 2026-06-23  
**Status:** Investigation complete - Solution 3 (External Workflows) recommended

## Problem Statement

PAV and dipcall are resource-intensive operations (PAV: 45+ minutes for 549 jobs, peak RSS ~36GB for minimap2 step alone). When restarting the pipeline or running multiple benchmarks with the same assembly+reference combination, all variant calling work is redone from scratch. This wastes significant compute time and resources.

## Root Cause

Assembly variant callers run as nested workflows within the main DeFrABB pipeline. Snakemake treats each run directory independently, so even though the inputs (assembly + reference) are identical, outputs in `results/asm_varcalls/` are not reused across different run directories (e.g., `20260617_v0.022_fulltest/` vs `20260623_v0.023_test/`).

## Solutions Evaluated

### Solution 1: Snakemake Between-Workflow Caching ❌ **NOT VIABLE**

**Approach:** Use Snakemake's built-in `cache: True` directive (available since v7.8, current version 9.17.3).

**Why it doesn't work:**  
Snakemake's between-workflow caching requires rules with multiple outputs to use `multiext()` to declare them with a common prefix. Both `run_pav` and `run_dipcall` have heterogeneous outputs that cannot be refactored to share a common prefix:

```python
# run_dipcall outputs - cannot use multiext()
output:
    make=".mak",           # Different extension families
    vcf=".dip.vcf.gz",  
    bed=".dip.bed",
    bam1=".hap1.bam",
    bam2=".hap2.bam"

# run_pav outputs - in DIFFERENT DIRECTORIES
output:
    vcf="{vc_id}/{sample_id}.vcf.gz",
    h1_bed="{vc_id}/results/{sample_id}/callable/callable_regions_h1_500.bed.gz"  # nested
```

Refactoring to `multiext()` would require:
1. Restructuring dipcall outputs (medium effort, ~2-4 hours)
2. Restructuring PAV outputs (HIGH effort, PAV is an external tool - would need wrapper layer)
3. Updating all downstream rules that reference these outputs
4. Risk of breaking existing functionality

**Verdict:** Not worth the refactoring cost for this approach.

---

### Solution 2: Persistent Varcall Results Directory ⚠️ **VIABLE, MANUAL**

**Approach:** Create `/defrabb_runs/varcall_cache/{ref}_{asm}_{caller}_{params}/` as a persistent cache outside run directories. Add checkpoint rules that symlink from cache if results exist.

**Pros:**
- Full control over cache behavior
- Can add versioning (e.g., `pav_v2.4.6/`)
- Works with current Snakemake version
- No dependency on Snakemake's caching semantics

**Cons:**
- Requires manual cache management (pruning old entries)
- More code to maintain (symlinking logic, cache validation)
- Need to handle partial failures (incomplete cache entries)
- Race conditions if multiple runs access cache simultaneously

**Estimated implementation:** 4-6 hours

**When to use:** Short-term tactical solution if you need caching soon.

---

### Solution 3: External Workflow Pattern ✅ **RECOMMENDED**

**Approach:** Separate PAV and dipcall into independent Snakemake workflows that run once and produce versioned, reusable outputs. DeFrABB consumes these as pre-computed inputs.

**Architecture:**
```bash
# One-time varcall runs (or re-run when new assembly/ref versions arrive)
./run_varcallers.sh --ref GRCh38 --asm HG002v1.1 --caller pav --version v2.4.6
./run_varcallers.sh --ref GRCh38 --asm HG002v1.1 --caller dipcall --params z2k

# Outputs go to: /defrabb_runs/varcalls/{ref}_{asm}_{caller}_{version}/

# DeFrABB runs reference the pre-computed varcalls
./run_defrabb run -r 20260623 --varcalls /defrabb_runs/varcalls/
```

**Pros:**
- **Aligns with GIAB workflow:** varcalls are expensive assets that get reused across many benchmark iterations
- Maximum flexibility - run varcalls once, use in unlimited benchmarks
- Clean separation of concerns (variant calling vs benchmark generation)
- Can version varcall outputs independently
- Easier to parallelize across machines (run PAV on one host, dipcall on another)
- Natural fit for milestone work (#11, #12-13) where varcaller refactoring is planned

**Cons:**
- Requires new workflow scripts (`run_varcallers.sh` or separate Snakefile)
- Breaking change to current workflow (but can be phased in gradually)
- Need varcall "catalog" or registry to track available combinations

**Estimated implementation:** 1-2 days

**When to use:** Long-term strategic solution. Should be part of the milestone #11-13 refactoring.

---

## Recommendation

**Immediate (for current run):** Let the PAV job complete. Monitor with:
```bash
tail -f /defrabb_runs/runs/20260617_v0.022_fulltest-HG002v1.1/logs/asm_varcalls/GRCh38_HG002-T2TQ100v1.1-pav_HG002_pav.log
```

**Short term (next 1-2 weeks):** If you need to re-run frequently and can't wait for milestone refactoring, implement **Solution 2** (persistent cache directory with checkpoint rules).

**Long term (milestone #11-13):** Implement **Solution 3** (external varcall workflows). This aligns with:
- Issue #11: Refactor varcall module structure
- The GIAB operational model (varcalls as reusable assets)
- Better scalability for multi-assembly, multi-reference benchmarking

---

## Appendix: Snakemake Cache Error

When attempting Solution 1, Snakemake 9.17.3 returned:

```
WorkflowError in rule run_dipcall in file "rules/asm-varcall.smk", line 12:
Rule is marked for between workflow caching but has multiple output files. 
This is only allowed if multiext() is used to declare them (see docs on 
between workflow caching).
```

This confirms the technical limitation that prevents using Snakemake's native caching for these rules without significant refactoring.
