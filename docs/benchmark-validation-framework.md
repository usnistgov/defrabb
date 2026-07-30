# Benchmark Validation Framework

**Purpose:** Systematic evaluation of differences between benchmark versions to ensure changes are intentional improvements, not regressions or artifacts.

**Principle:** v5.0q represents a validated, conservative baseline. All differences from v5.0q must be critically evaluated and categorized by root cause before accepting new benchmarks.

## Validation Workflow

### 1. Quantitative Comparison

**Region Coverage:**
- Total basepairs (by chromosome)
- Interval count
- Overlap with baseline (bp and %)
- Unique regions (new vs baseline, baseline vs new)

**Variant Counts:**
- SNP/INDEL counts by filter status
- Ti/Tv ratios
- Het/Hom ratios
- Variant length distributions

**Evaluation Metrics:**
- Precision, Recall, F1 vs established comparison callsets
- Breakdown by stratification (tandem repeats, homopolymers, segdups, etc.)
- Known artifact patterns (e.g., FN = FP indicating systematic error)

### 2. Root Cause Classification

For each difference, determine the category:

#### A. **Pipeline Artifacts** (REJECT)
- Tool bugs that strip/corrupt data
- Format conversion errors
- Unintended loss of information
- **Example:** `truvari anno trf` stripping phasing (v0.022 finding)

#### B. **Fixed Bugs** (ACCEPT with documentation)
- Corrections to previous errors
- Improved handling of edge cases
- **Example:** TRF contig mismatch filter (#197)

#### C. **Exclusion Changes** (EVALUATE)
- Different exclusion sets applied
- Updated slop/merge parameters
- New exclusion categories
- **Example:** chr2:114M segdup included in v5.0q, excluded in v0.022 (updated criteria)

#### D. **Assembly Version Differences** (ACCEPT with documentation)
- Different assembly version (Q100v1.0 vs v1.1)
- Assembly corrections/updates
- Different assembly methods
- **Example:** HG002 T2TQ100v1.1 vs v1.0 differences

#### E. **Parameter Changes** (EVALUATE via sweeps)
- Variant caller parameters (dipcall -z, PAV merge thresholds)
- VCF processing steps
- Normalization behavior
- **Example:** dipcall z2k vs z5k vs z1k parameter profiles

#### F. **Intended Improvements** (ACCEPT with validation)
- More sensitive variant detection
- Better variant representation
- Improved callable region estimation
- **Validation required:** Check that expansions are in high-confidence regions, not noise

### 3. Difference Investigation Checklist

For each identified difference region or metric change:

- [ ] **Quantify the difference** (bp, variant count, metric delta)
- [ ] **Locate the difference** (chromosomal regions, stratification categories)
- [ ] **Trace through pipeline stages**
  - [ ] Raw variant call output
  - [ ] After normalization
  - [ ] After each VCF processing step
  - [ ] After exclusion application
  - [ ] Final benchmark
- [ ] **Check for corresponding code/config changes**
  - [ ] Git log between versions
  - [ ] Config file diffs (resources.yml, analyses tables)
  - [ ] Exclusion set changes
- [ ] **Validate against independent evidence**
  - [ ] Assembly support (IGV inspection if needed)
  - [ ] Read support from HiFi/ONT/Illumina
  - [ ] Concordance with orthogonal callers
  - [ ] Known problematic regions (segdups, STRs, medically relevant genes)
- [ ] **Classify root cause** (A-F above)
- [ ] **Decision:** Accept, Reject, or Further Investigation
- [ ] **Document in issues/** or session notes

### 4. Acceptance Criteria

A difference is acceptable if:

1. **Root cause is identified and documented**
2. **Category is B, D, or F** (fixed bugs, assembly updates, intended improvements)
3. **Category C or E** (exclusions, parameters) AND:
   - Change is intentional (documented in config)
   - Improvement is validated (F1 increase, better stratification performance)
   - No regressions in other metrics
4. **Category A** (artifacts) → **MUST BE FIXED** before release

### 5. Red Flags Requiring Investigation

**Immediate investigation triggers:**
- **FN = FP** (systematic representation mismatch, like phasing issue)
- **Perfect metrics** (F1 = 1.0000 for all) that seem too good
- **Large region losses** (>1% of baseline regions missing)
- **Metric regressions** (F1 decrease, precision drop without recall gain)
- **Stratification anomalies** (tandem repeats F1 = 1.0, segdups F1 = 0.5)
- **Suspicious patterns** (all differences on one chromosome, chrY only, etc.)

## Parameter Optimization Integration

For v0.023+ parameter sweeps:

### Pre-Sweep Baseline
1. Establish v5.0q as the **validation baseline**
2. Run baseline parameter set (z2k, giab, standard exclusions) on new assembly
3. Fully vet differences from v5.0q using this framework
4. **Only proceed with sweeps after baseline is validated**

### During Sweep
1. Generate sweep results (multiple param combinations)
2. Score by metrics (F1, precision, recall)
3. **Top N candidates** → Full difference evaluation
4. Compare candidates not just to baseline, but to **each other**
5. Look for parameter-dependent patterns

### Post-Sweep Validation
1. Selected "winner" parameters → **Full validation cycle**
2. Multi-genome validation (HG002, HG005, HG007, etc.)
3. Check for:
   - Overfitting to HG002
   - Genome-specific regressions
   - Stratification trade-offs (e.g., gain in homopolymers, loss in segdups)
4. Independent reviewer sign-off before production use

## Documentation Requirements

For each benchmark release, maintain:

### Difference Report
- `docs/validation/v0.0XX-vs-baseline.md`:
  - Quantitative summary table
  - Root cause breakdown (categories A-F with counts/bp)
  - Known issues flagged for future work
  - Acceptance decision with rationale

### Change Log
- `CHANGELOG` entry documenting:
  - Intentional parameter/exclusion changes
  - Bug fixes affecting regions/variants
  - Assembly version updates
  - Pipeline improvements

### Session Notes
- `docs/sessions/YYYY-MM-DD-v0.0XX-validation.md`:
  - Detailed investigation notes
  - Commands run for difference analysis
  - Visualizations/screenshots if generated
  - Open questions for future investigation

## Example: v0.022 vs v5.0q Evaluation

### Quantitative Summary
- **GRCh38:** 100.00% overlap, +947 KB v0.022 (+0.034%)
- **GRCh37:** 99.98% overlap, +914 KB v0.022, -408 KB v5.0q
- **CHM13:** 100.00% overlap, +249 KB v0.022 (+0.009%)
- **SNP F1:** 0.9998 (619 FN = 619 FP) ← RED FLAG
- **INDEL F1:** 1.0000

### Root Cause Analysis

| Finding | Category | Root Cause | Decision |
|---------|----------|------------|----------|
| 619 SNP FN = FP | **A - Artifact** | `truvari anno trf` strips phasing | **REJECT - Fix in v0.023** |
| +947 KB GRCh38 regions | **F - Improvement** | Better callable region detection in v1.1 assembly | **ACCEPT** (pending spot checks) |
| -408 KB GRCh37 (chr2:114M) | **C - Exclusion** | Segdup excluded in v0.022, included in v5.0q (updated criteria) | **ACCEPT** (more conservative) |
| +249 KB CHM13 regions | **F - Improvement** | Assembly-supported regions | **ACCEPT** (pending validation) |

### Actions Taken
1. **Phasing artifact:** Root cause identified (`truvari anno trf`), solution designed (use `xy_fix` profile), documented in `docs/issues/happy-phasing-mismatch-v5q.md`
2. **Region expansions:** Categorized as likely improvements, flagged for spot-check validation in next session
3. **GRCh37 segdup:** Identified as intentional exclusion difference, acceptable
4. **v0.022 tag created** with known phasing issue documented in CHANGELOG

### Open Items for v0.023
- [ ] Fix phasing loss (use `xy_fix` instead of `xy_trf` for smvar)
- [ ] Spot-check sample v0.022 expansion regions for assembly support
- [ ] Re-run evaluations with phasing preserved, verify F1 → 1.0000
- [ ] Document phasing fix in v0.023 CHANGELOG

## Tools and Commands

### Region Comparison
```bash
# Overlap
bedtools intersect -a baseline.bed -b new.bed | awk '{sum+=$3-$2} END{print sum}'

# Unique to new
bedtools subtract -a new.bed -b baseline.bed | awk '{sum+=$3-$2} END{print sum}'

# Unique to baseline  
bedtools subtract -a baseline.bed -b new.bed | awk '{sum+=$3-$2} END{print sum}'

# By-chromosome breakdown
bedtools subtract -a new.bed -b baseline.bed | awk '{chr[$1]++; bp[$1]+=$3-$2} END{for(c in chr) print c, chr[c], bp[c]}'
```

### VCF Comparison
```bash
# Variant counts
bcftools stats new.vcf.gz | grep "^SN"

# Compare with baseline
scripts/compare_evaluations.py \
  --results-dir results/evaluations \
  --baseline v5.0q \
  --metric F1_Score \
  --output comparison.md
```

### Pipeline Stage Tracing
```bash
# Check genotypes at each stage
for vcf in annotations/*.vcf.gz; do
  echo "$vcf:"
  gunzip -c $vcf | grep "^chr1" | awk '$2==89177' | cut -f10
done
```

## References

- v5.0q release notes: [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/)
- Parameter optimization design: `docs/design/v0.023-parameter-optimization-design.md`
- Known issues: `docs/issues/`
- Session notes: `docs/sessions/`
