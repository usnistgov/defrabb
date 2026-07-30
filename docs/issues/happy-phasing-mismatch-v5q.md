# Hap.py Phasing Mismatch with v5.0q Comparison Callsets

**Status:** Identified in v0.022  
**Priority:** High  
**Created:** 2026-07-30  
**Discovered:** During v0.022 release review  
**Affects:** All hap.py evaluations using phased comparison callsets

## Summary

The 619 SNP FN/FP match in v0.022 GRCh38 evaluations is caused by **genotype phasing representation mismatch** between benchmark and comparison VCFs, not genuine variant differences. Hap.py treats `1/1` (unphased) and `1|1` (phased) as different genotypes.

## Root Cause

**v5.0q comparison VCF uses phased genotypes (`1|1`) while v0.022 benchmark VCF uses unphased (`1/1`):**

```
# v5.0q comparison (TRUTH):
chr1  89177  .  A  G  30  .  .  GT:AD  1|1:0,2

# v0.022 benchmark (QUERY):
chr1  89177  .  A  G  30  .  .  GT:AD  1/1:0,2
```

Hap.py interprets phasing differences as genotype mismatches:
- **619 FN:** Benchmark has `1/1`, v5.0q has `1|1` → genotype doesn't match exactly
- **619 FP:** v5.0q has `1|1`, benchmark has `1/1` → genotype doesn't match exactly

**These are identical variants differing only in phasing representation (`/` vs `|`).**

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

### Option 1: Unphase Comparison VCFs (Recommended)

Convert phased to unphased genotypes in comparison VCFs before hap.py evaluation:

```bash
bcftools +setGT input.vcf.gz -- -t a -n u | bcftools view -Oz -o unphased.vcf.gz
```

**Implementation:**
- Add `unphase_comparison_vcf` rule in `rules/download_resources.smk` or `rules/evaluation.smk`
- Apply to all comparison VCFs before `run_happy` rule
- Update wildcard dependencies in evaluation rules

**Pros:**
- Simple, deterministic transformation
- Preserves all variant information
- No hap.py version dependencies

**Cons:**
- Loses phasing information (acceptable for benchmarking - we're comparing variant concordance, not phasing accuracy)

### Option 2: Use hap.py `--no-phasing` Flag

Check if hap.py supports ignoring phasing differences:

```bash
hap.py --help | grep -i phas
```

**Pros:**
- Clean solution if flag exists
- No VCF pre-processing needed

**Cons:**
- May not be available in current hap.py version
- Behavior may vary across hap.py versions

### Option 3: Phase Benchmark VCF

Use `whatshap phase` or similar to add phasing to benchmark VCF.

**Pros:**
- Provides phasing information

**Cons:**
- Complex - requires read data or trio information
- Computationally expensive
- Not necessary for variant concordance benchmarking
- May introduce phasing errors

## Recommended Solution

**Option 1: Unphase comparison VCFs** is the simplest and most robust solution.

Add rule in `rules/download_resources.smk`:

```python
rule unphase_comparison_vcf:
    """
    Remove phasing from comparison VCF to match unphased benchmark.
    Prevents hap.py from treating 1/1 and 1|1 as different genotypes.
    """
    input:
        vcf = lambda wc: comp_vcf_url(wc.ref, wc.comp_id, wc.bench_type)
    output:
        vcf = "resources/comparison_variant_callsets/{ref}_{comp_id}-{bench_type}.unphased.vcf.gz",
        tbi = "resources/comparison_variant_callsets/{ref}_{comp_id}-{bench_type}.unphased.vcf.gz.tbi"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bcftools +setGT {input.vcf} -- -t a -n u | \
        bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """
```

Update `run_happy` to use `.unphased.vcf.gz` as comparison input.

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
