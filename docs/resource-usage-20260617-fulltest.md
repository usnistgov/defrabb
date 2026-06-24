# Resource Usage Analysis: 20260617_v0.022_fulltest-HG002v1.1

**Date:** 2026-06-24  
**Run:** 20260617_v0.022_fulltest-HG002v1.1 (whole-genome HG002 T2TQ100v1.1)  
**Purpose:** Measure actual resource usage for PAV and dipcall to right-size config

## Measured Resource Usage

All values from `benchmark/*.tsv` (Snakemake benchmark directive output).

### PAV (GRCh38, HG002-T2TQ100v1.1)

| Metric | Value | Notes |
|--------|-------|-------|
| Wall time | 11:23:00 (40,978 sec) | ~11.5 hours |
| Peak RSS | **64.9 GB** | Max resident set size |
| Peak VMS | 111.9 GB | Virtual memory |
| CPU time | 440,984 sec | 122 CPU-hours (24 threads) |
| Mean load | 878% | High parallelism |

**Key finding:** PAV peaked at 65 GB RSS, significantly higher than the previously estimated 64 GB in commit 73c825a (which isn't in current master).

### Dipcall (3 references, HG002-T2TQ100v1.1)

| Reference | Wall time | Peak RSS | Peak VMS | CPU time | Notes |
|-----------|-----------|----------|----------|----------|-------|
| GRCh37 | 1:08:00 (4,078 sec) | **102.0 GB** | 107.7 GB | 30,768 sec | Highest RSS |
| GRCh38 | 1:10:00 (4,194 sec) | **102.6 GB** | 112.2 GB | 30,292 sec | Highest RSS |
| CHM13v2.0 | 0:44:30 (2,671 sec) | **116.2 GB** | 120.7 GB | 21,327 sec | **Highest RSS - T2T genome** |

**Key finding:** Dipcall peaks at **102-116 GB** depending on reference. CHM13 (T2T) is highest at 116 GB, likely due to lack of gaps and longer chromosomes.

## Current Config vs Measured

| Resource | Current (`config/resources.yml`) | Measured Peak | Recommended |
|----------|----------------------------------|---------------|-------------|
| `_pav_mem` | **MISSING** | 64.9 GB | **80,000 MB** |
| `_dipcall_mem` | 32,000 MB (per job) | 116.2 GB (total) | **40,000 MB** (per job) |
| `_happy_mem` | 64,000 MB | Not measured this run | Check commit 73c825a (160,000 MB) |

### Memory Calculation Notes

**dipcall:**  
- Current: `_dipcall_jobs: 4` × `_dipcall_mem: 32000` = 128 GB total reservation
- Measured peak: 116 GB (CHM13)
- Current reservation is adequate but tight
- Recommendation: Increase per-job from 32,000 → **40,000 MB** for safety margin
  - New total: 4 × 40,000 = 160 GB (margin for future larger assemblies)

**PAV:**
- Currently NO `_pav_mem` in config (regression - was added in branch commit 73c825a but not merged)
- Measured peak: 65 GB
- Recommendation: Add **`_pav_mem: 80000`** (25% margin over measured)

**Note on branch 73c825a:**  
Commit 73c825a (on an unmerged branch) added:
- `_pav_mem: 64000` ← too low (measured 65 GB)
- `_happy_mem: 160000` ← based on prior measurements (95-151 GB)
- `_truvari_mem: 32000`
- `_truvari_anno_mem: 16000`

These should be merged, but `_pav_mem` should be increased to 80,000 based on actual measurements.

## Recommendations

1. **Immediate:** Add missing memory settings to `config/resources.yml`:
   ```yaml
   _pav_mem: 80000         # 25% margin over measured 65 GB peak
   _dipcall_mem: 40000     # up from 32000; total 160 GB for 4 jobs
   _happy_mem: 160000      # from commit 73c825a, based on prior measurements
   _truvari_mem: 32000     # from commit 73c825a
   _truvari_anno_mem: 16000 # from commit 73c825a
   ```

2. **Also update schema:** Add these keys to `schema/resources-schema.yml` if not present.

3. **Apply memory to rules:** Ensure `run_pav`, `run_happy`, and truvari rules declare `resources: mem_mb=config["_*_mem"]`.

4. **Monitor future runs:** Check `benchmark/*.tsv` after whole-genome runs to validate reservations match reality.

## Context

This analysis addresses the varcall caching investigation (docs/issues/varcall-caching-investigation.md). While we can't implement between-workflow caching yet (Solution 1 not viable, Solution 3 deferred), we CAN right-size memory to prevent OOM and make best use of available resources.
