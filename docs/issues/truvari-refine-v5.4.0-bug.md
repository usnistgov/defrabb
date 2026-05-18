# Truvari Refine v5.4.0 Bug - Samtools faidx Failure

**Status:** Open  
**Priority:** Medium  
**Created:** 2026-05-18  
**Affects:** truvari v5.4.0 (current version in envs/truvari_core.yml)

## Summary

Truvari refine v5.4.0 fails when processing large numbers of regions (2700+) due to a samtools faidx error when creating temporary reference FASTA files.

## Error Details

```
[E::fai_build_core] File truncated at line 1
pysam.utils.SamtoolsError: 'samtools returned with error 1: [faidx] Could not build fai index /tmp/XXXX.fa.fai'
```

**Traceback location:**
```python
File "truvari/phab.py", line 524, in run_phab
    vcf_info.set_regions(regions, buffer)
File "truvari/phab.py", line 310, in set_regions
    samtools.faidx(out_fn)
```

## Reproduction

Occurs when running `truvari refine` on benchmarks with 2700+ regions to refine:
- GRCh38 Sawfish: 2665 regions
- GRCh37 sentLRSV-ontR9: 2740+ regions  
- CHM13 Sniffles2: 2740 regions

All failed at the "Preparing reference/regions" step after successfully identifying regions.

## Known Upstream Issues

Similar issues reported in truvari GitHub:
- [Issue #302](https://github.com/ACEnglish/truvari/issues/302) - GRCh37 numeric chromosome issue (fixed in dtype patch)
- Related to faidx failures with large region sets

## Current Workaround

Changed `eval_cmd` from `truvari_refine` to `truvari` in analyses tables to skip the refine step while keeping successful truvari bench results.

**File:** `config/analyses_20260323_v0.020_HG002v1.1-CI.tsv`

Truvari refine is optional post-processing that improves precision of TP matches via local realignment. Skipping it retains all core benchmark metrics (precision/recall/F1).

## Action Items

1. **Investigate version history:**
   - Test truvari v5.3.0 (released 2025-04-21) to see if bug exists
   - Test truvari v5.2.0 (released 2025-02-16) to identify when bug was introduced
   - Check truvari v4.3.0 (previously used, had FIPS OpenSSL issues but may have working refine)

2. **Report upstream:**
   - File issue on ACEnglish/truvari GitHub with:
     - Error details and traceback
     - Number of regions that trigger the failure (~2700+)
     - Affected references (GRCh37, GRCh38, CHM13v2.0)
     - Request fix or workaround

3. **Update DeFrABB when fixed:**
   - Pin to working version once identified
   - Update envs/truvari_core.yml
   - Re-enable truvari_refine in default analyses templates
   - Document in CHANGELOG

## Environment

- **Truvari version:** 5.4.0 (bioconda::truvari)
- **Dependencies:** bcftools >=1.20, mafft >=7.515
- **Conda env:** envs/truvari_core.yml
- **System:** NIST FIPS-enabled Linux

## References

- Truvari refine docs: https://github.com/ACEnglish/truvari/wiki/refine
- Issue #182 (refine hangs): https://github.com/ACEnglish/truvari/issues/182
- Issue #302 (GRCh37 dtype): https://github.com/ACEnglish/truvari/issues/302
