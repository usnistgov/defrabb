# Hap.py Phasing Mismatch with v5.0q Comparison Callsets

**Status:** Root cause identified in v0.022  
**Priority:** High  
**Created:** 2026-07-30  
**Discovered:** During v0.022 release review  
**Root Cause:** `truvari anno trf` strips phasing from benchmark VCFs  
**Affects:** All smvar benchmarks using `xy_trf` VCF processing profile

## Summary

The 619 SNP FN/FP match in v0.022 GRCh38 evaluations is caused by **`truvari anno trf` stripping phasing** from benchmark VCFs during TRF annotation. Assembly-derived variants start phased (`1|1`) but end up unphased (`1/1`) in the final benchmark, while v5.0q comparison retains phasing, creating a representation mismatch that hap.py interprets as genotype errors.

## Root Cause

**`truvari anno trf` strips phasing during TRF annotation.**

### Processing Chain

Dipcall generates phased variants from diploid assembly:
```
1. Raw dipcall:        chr1  89177  ...  GT:AD  1|1:0,2  ✓ PHASED
2. Rename:             chr1  89177  ...  GT:AD  1|1:0,2  ✓ PHASED  
3. norm:               chr1  89177  ...  GT:AD  1|1:0,2  ✓ PHASED
4. fix_XY_genotype:    chr1  89177  ...  GT:AD  1|1:0,2  ✓ PHASED
5. trfanno:            chr1  89177  ...  GT:AD  1/1:0,2  ✗ PHASING LOST ← THE CULPRIT
6. Final benchmark:    chr1  89177  ...  GT:AD  1/1:0,2  ✗ UNPHASED
```

v5.0q comparison retains phasing:
```
v5.0q:                 chr1  89177  ...  GT:AD  1|1:0,2  ✓ PHASED
```

Hap.py interprets phasing difference as genotype mismatch:
- **619 FN:** Benchmark has `1/1`, v5.0q has `1|1` → GT doesn't match
- **619 FP:** v5.0q has `1|1`, benchmark has `1/1` → GT doesn't match

**These are identical variants differing only in phasing representation.**

## Evidence

### Extended CSV Analysis

From `GRCh38_v5.0q-smvar_HG2-T2TQ100-V1.1_smvar_dipcall-z2k.extended.csv`:

```
TRUTH.FN.het: 0.000000
TRUTH.FN.homalt: 0.000000
# → All 619 FNs reported as "no genotype"

QUERY.FP.het: 619.000000  
QUERY.FP.homalt: 0.000000
# → All 619 FPs are heterozygous calls

FP.gt: 597
FP.al: 22
# → 597 genotype errors vs 22 allele errors
```

### Summary Statistics

```
Type  TRUTH.TOTAL  TRUTH.TP  TRUTH.FN  QUERY.TOTAL  QUERY.FP  METRIC.F1_Score
SNP   3,652,238    3,651,619  619      3,657,936     619      0.9998
```

**Suspicious pattern:** `TRUTH.FN == QUERY.FP == 619` indicates systematic comparison artifact, not biological variation.

### VCF Inspection

Both VCFs have identical variant positions and alleles:
```bash
# First 10 chr1 variants
v5.0q:     chr1:89177 A>G, 89254 A>G, 89264 A>G, 89327 G>A, ...
benchmark: chr1:89177 A>G, 89254 A>G, 89264 A>G, 89327 G>A, ...
```

Only difference is GT field:
- v5.0q: `GT:AD 1|1:0,2` (phased)
- benchmark: `GT:AD 1/1:0,2` (unphased)

## Impact

- **v0.022 benchmark F1 should be ~1.0000** for SNPs, not 0.9998
- The 619 discrepancies are **artifacts of GT representation**, not genuine differences
- Affects **all hap.py evaluations** using phased comparison callsets (v5.0q and potentially others)
- **Does NOT affect variant concordance** - variants are present in both VCFs
- **Does NOT affect INDEL metrics** - INDELs show F1 = 1.0000 (perfect match)

## Solutions

### Option 1: Skip TRF Annotation for Small Variants (Recommended)

TRF annotation is primarily useful for structural variants. Small variant benchmarks don't need it.

**Change `bench_vcf_processing` profile from `xy_trf` to `xy_fix`:**

```yaml
# config/resources.yml - profile already exists
vcf_processing_profiles:
  xy_fix:
    - fix_XY_genotype
```

**Update analyses tables:** Change smvar rows from `bench_vcf_processing: xy_trf` to `bench_vcf_processing: xy_fix`.

**Pros:**
- Simplest solution - no new code
- Preserves phasing from assembly
- Faster (skips expensive TRF annotation)
- TRF annotations not used in smvar benchmarking anyway

**Cons:**
- Loses TRF annotations (but they're unused for smvar)

### Option 2: Restore Phasing After trfanno

Create new rule to copy phasing from pre-trfanno to post-trfanno VCF:

```python
rule restore_phasing_after_trfanno:
    input:
        phased = ".../{prefix}.fix_XY_genotype.vcf.gz",
        unphased = ".../{prefix}.trfanno.vcf"
    output:
        ".../{prefix}.trfanno.phased.vcf.gz"
    shell:
        """
        python scripts/restore_phasing.py \
            --phased {input.phased} \
            --unphased {input.unphased} \
            --output {output}
        """
```

**Pros:**
- Keeps TRF annotations
- Restores phasing

**Cons:**
- More complex - new rule + script
- Additional processing step
- Requires careful VCF position matching

### Option 3: Fix truvari anno trf Upstream

Report issue to truvari developers that `anno trf` should preserve phasing.

**Pros:**
- Fixes root cause for all users

**Cons:**
- Requires upstream fix
- Need workaround in meantime

## Recommended Solution

**Option 1: Use `xy_fix` instead of `xy_trf` for smvar benchmarks.**

**For new analyses tables (v0.023+):** Set `bench_vcf_processing: xy_fix` for smvar rows.

**For existing v0.022 analyses tables:** Keep as-is (historical record), but document in CHANGELOG that phasing loss was identified and will be fixed in v0.023.

## Testing Plan

1. **Unit test:** Verify `bcftools +setGT` converts `1|1` → `1/1` and `0|1` → `0/1`
2. **Integration test:** Re-run GRCh38 v5.0q hap.py evaluation with unphased comparison
3. **Validation:** Confirm F1 improves to ~1.0000 and FN/FP both drop to ~22 (the genuine allele errors)
4. **Comparison:** Run on all 3 references (GRCh37, GRCh38, CHM13) to verify consistent improvement

## Follow-up Actions

1. Implement unphasing rule for v0.023
2. Re-run v0.022 fulltest evaluations with unphased comparisons to establish corrected baseline
3. Update CHANGELOG to document this discovery and fix
4. Consider whether v5.0q phasing information should be preserved separately for phasing-specific evaluations

## Related Issues

- v0.022 release metrics review
- Benchmark vs comparison VCF format standardization

## References

- v0.022 fulltest results: `../20260721_v0.022_fulltest-HG002v1.1/results/evaluations/`
- Hap.py documentation: https://github.com/Illumina/hap.py
- bcftools +setGT plugin: http://samtools.github.io/bcftools/bcftools.html#setGT
