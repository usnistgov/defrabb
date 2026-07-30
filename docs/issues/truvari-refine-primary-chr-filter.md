# Fix truvari refine v5.4.0 faidx failure - Filter to Primary Chromosomes

**Status:** Planned for v0.023  
**Priority:** High  
**Created:** 2026-07-30  
**Related:** docs/issues/truvari-refine-v5.4.0-bug.md  
**GitLab Issue:** TBD

## Problem

Truvari refine v5.4.0 fails with samtools faidx error when processing 2700+ regions due to reference FASTA containing 218 sequences (GRCh38: 24 primary + 194 alt/random/Un/patches).

```
[E::fai_build_core] File truncated at line 1
pysam.utils.SamtoolsError: 'samtools returned with error 1: [faidx] Could not build fai index /tmp/XXXX.fa.fai'
```

Occurs with:
- GRCh38 Sawfish: 2665 regions
- GRCh37 sentLRSV-ontR9: 2740+ regions  
- CHM13 Sniffles2: 2740 regions

## Root Cause Analysis

The GRCh38 reference (`GRCh38_full_analysis_set_plus_decoy_hla.fa`) includes:
- **24 primary chromosomes:** chr1-22, chrX, chrY, chrM
- **194 additional sequences:**
  - Random contigs: chr*_KI*_random, chr*_GL*_random (68 sequences)
  - Unplaced contigs: chrUn_KI*, chrUn_GL* (127 sequences)
  - Patches: KMT2C_chr*, MAP2K3_chr*, KCNJ18_chr* (23 sequences)
  - chrEBV (1 sequence)

When truvari refine extracts 2700+ variant regions for local realignment across these 218 sequences, samtools faidx fails to build temporary FASTA indices, likely hitting internal limits on the number of sequences or command-line length when extracting regions.

**Key observation:** v5.0q comparison callsets contain **only 24 primary contigs** (chr1-22, chrX, chrY). Variants on alt/random/Un contigs in the benchmark VCF are **not evaluatable** because they don't exist in comparison callsets.

## Proposed Solution

**Filter benchmark VCFs to primary chromosomes after variant calling but before VCF processing.**

### Implementation

Add new rule in `rules/bench_vcf_normalize.smk`:

```python
rule filter_primary_chromosomes:
    """
    Remove variants on alt/random/unplaced contigs after variant calling.
    Keeps only primary chromosomes that exist in comparison callsets.
    Fixes truvari refine faidx failure with 2700+ regions.
    """
    input:
        vcf = rules.run_dipcall.output.vcf  # or run_pav
    output:
        vcf = "results/asm_varcalls/{vc_id}/dip.primary.vcf.gz",
        tbi = "results/asm_varcalls/{vc_id}/dip.primary.vcf.gz.tbi"
    params:
        primary_chr = lambda wc: get_primary_chromosomes(wc.ref)
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bcftools view -t {params.primary_chr} {input.vcf} | \
        bcftools view -Oz -o {output.vcf}
        bcftools index -t {output.vcf}
        """
```

Add helper function in `rules/helpers_ref.smk`:

```python
def get_primary_chromosomes(ref_id):
    """Return comma-separated list of primary chromosomes for reference."""
    # All references use chr1-22,chrX,chrY
    # Skip chrM (excluded in all benchmark sets)
    return "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
```

### Wildcard Updates

Update downstream rules to consume `{vc_id}/dip.primary.vcf.gz` instead of `{vc_id}/dip.vcf.gz`:
- `norm_vcf_samples`
- `split_smvar_stvar`
- Any other VCF processing rules

## Benefits

1. **Fixes truvari refine faidx failure** - Reduces sequences from 218 → 24, eliminating samtools limit
2. **Removes non-evaluatable variants** - Alt/unplaced contig variants can't be evaluated (comparison VCFs only have primary chr)
3. **Reduces VCF processing time** - Fewer contigs = faster normalization/annotation/filtering
4. **Aligns with v5.0q** - Matches comparison callset chromosome sets exactly
5. **Simplifies debugging** - Smaller VCF headers, clearer contig mismatches

## Testing Plan

1. **Unit test:** Verify `get_primary_chromosomes()` returns correct 24-chr list
2. **Integration test:** Run pipeline with primary-chr filter enabled
3. **Regression test:** Verify truvari refine succeeds on Sawfish/sentLRSV comparisons (2700+ regions)
4. **Comparison:** Confirm v0.023 benchmarks with filter match v0.022 benchmarks (no variants should be on alt contigs in current pipeline)

## Alternative Approaches Considered

1. **Increase samtools limits** - Not feasible; faidx limits are hardcoded
2. **Split truvari refine by chromosome** - Complex workflow change; doesn't address non-evaluatable variants
3. **Use stripped reference FASTA** - Breaks compatibility with existing reference URLs in config
4. **Filter at comparison callset** - Too late; truvari refine runs before comparison

## References

- Original bug report: `docs/issues/truvari-refine-v5.4.0-bug.md`
- Truvari refine docs: https://github.com/ACEnglish/truvari/wiki/refine
- GRCh38 reference spec: https://www.ncbi.nlm.nih.gov/grc/human/data
