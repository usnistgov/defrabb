# Parameter Optimization Infrastructure Design

**Version:** v0.024  
**Status:** Draft  
**Created:** 2026-07-31

> **Note (2026-09-28):** the discrepancy-extraction and curation components
> (§1, `scripts/extract_discrepancies.py`) were removed from DeFrABB and are
> being developed in the separate variant-curation-review project.

## Overview

Comprehensive infrastructure for designing, running, analyzing, and reporting parameter sweep results with two primary use cases:

1. **Coverage Expansion**: Improve v5.0q by increasing coverage through multi-callset comparison and RIDE-based discrepancy verification
2. **Multi-Genome Optimization**: Define parameters for new assemblies via HG002 optimization → multi-genome validation

## Core Principles

### RIDE Framework
- **Reliability**: Precision/recall across multiple callsets and technologies
- **Interpretability**: Clear stratification breakdowns showing performance by genomic context
- **Discrimination**: Ability to detect errors in comparison callsets (not just concordance)
- **Evidence**: Technology-specific support from multiple sequencing platforms

### Avoid Over-Training
- Manual IGV curation determines ground truth (not automated metrics alone)
- Errors in comparison callsets = good (benchmark is discriminating)
- Errors in benchmark = learning opportunity for parameter tuning
- Multi-technology validation prevents optimization to single platform artifacts

## Architecture

### 1. Discrepancy Extraction & Tracking

**Purpose**: Extract FN/FP variants from evaluations and track manual curation verdicts.

#### Components

##### 1.1 Discrepancy Extractor (`scripts/extract_discrepancies.py`)
```bash
./scripts/extract_discrepancies.py \
  --evaluation-dir results/evaluations/happy/<eval_id> \
  --benchmark-vcf <benchmark.vcf.gz> \
  --comparison-vcf <comparison.vcf.gz> \
  --stratifications resources/stratifications/GRCh38/ \
  --output discrepancies/<eval_id>.tsv
```

**Outputs**:
- TSV with columns: `chrom`, `pos`, `ref`, `alt`, `type` (FN/FP), `variant_type` (SNP/INDEL/SV), `stratifications` (comma-separated), `benchmark_gt`, `comparison_gt`, `hap.py_classification`
- Linked to hap.py extended CSV for full context

##### 1.2 Curation Tracker (`scripts/curate_discrepancies.py`)
```bash
# Initialize curation session
./scripts/curate_discrepancies.py init \
  --discrepancies discrepancies/<eval_id>.tsv \
  --bams bams/HG002/*.bam \
  --ref resources/references/GRCh38.fa \
  --session curation/<eval_id>

# Record verdict
./scripts/curate_discrepancies.py record \
  --session curation/<eval_id> \
  --variant chr1:12345:A>G \
  --verdict benchmark_error \
  --confidence high \
  --technologies "illumina,hifi,ont" \
  --notes "Clear insertion artifact in assembly, all techs support REF"

# Generate IGV batch script
./scripts/curate_discrepancies.py igv \
  --session curation/<eval_id> \
  --filter "uncurated,low_confidence" \
  --output igv_session.xml
```

**Curation Schema** (`curation/<eval_id>/verdicts.tsv`):
- `variant_id`: chr:pos:ref>alt
- `verdict`: `benchmark_error | comparison_error | representation_diff | ambiguous`
- `confidence`: `high | medium | low`
- `technologies`: Comma-separated list of supporting technologies
- `notes`: Free text
- `curator`: Username
- `timestamp`: ISO 8601

##### 1.3 IGV Integration
- Auto-generate IGV session XMLs with:
  - Reference genome
  - All available BAM files (Illumina, HiFi, ONT, etc.)
  - Benchmark VCF
  - Comparison VCF
  - Stratification BED tracks
  - Bookmarks at discrepancy positions
- Filter by curation status, stratification, variant type

### 2. Multi-Callset Comparison Framework

**Purpose**: Evaluate draft benchmarks against multiple external callsets to assess RIDE metrics.

#### Components

##### 2.1 Callset Registry (`config/comparison_callsets.yml`)
```yaml
comparison_callsets:
  HG002:
    v5.0q-smvar:
      vcf: s3://giab/v5.0q/HG002_GRCh38_smvar.vcf.gz
      bed: s3://giab/v5.0q/HG002_GRCh38_smvar.bed
      type: benchmark
      version: v5.0q
      
    v4.2.1-smvar:
      vcf: s3://giab/v4.2.1/HG002_GRCh38_smvar.vcf.gz
      bed: s3://giab/v4.2.1/HG002_GRCh38_smvar.bed
      type: benchmark
      version: v4.2.1
      
    deepvariant-illumina:
      vcf: /defrabb_runs/callsets/HG002/deepvariant_illumina.vcf.gz
      type: caller_output
      technology: illumina
      caller: deepvariant
      version: 1.6.0
      
    deepvariant-hifi:
      vcf: /defrabb_runs/callsets/HG002/deepvariant_hifi.vcf.gz
      type: caller_output
      technology: hifi
      caller: deepvariant
      version: 1.6.0
      
    gatk-illumina:
      vcf: /defrabb_runs/callsets/HG002/gatk_haplotypecaller.vcf.gz
      type: caller_output
      technology: illumina
      caller: gatk
      version: 4.5.0
```

##### 2.2 Multi-Callset Evaluator (`scripts/evaluate_multi_callset.py`)
```bash
./scripts/evaluate_multi_callset.py \
  --benchmark results/draft_benchmarksets/<bench_id>/ \
  --callsets config/comparison_callsets.yml \
  --sample HG002 \
  --output results/multi_callset_eval/<bench_id>/
```

**Outputs**:
- One hap.py/truvari run per comparison callset
- Aggregated metrics TSV
- RIDE score summary
- Discrepancy overlap analysis (variants rejected by multiple callsets = likely benchmark errors)

##### 2.3 RIDE Metrics Calculator (`scripts/calculate_ride_metrics.py`)
```bash
./scripts/calculate_ride_metrics.py \
  --multi-callset-dir results/multi_callset_eval/<bench_id>/ \
  --curation-verdicts curation/<bench_id>/verdicts.tsv \
  --output ride_metrics/<bench_id>.json
```

**Metrics**:
- **Reliability**: 
  - Precision/recall vs each callset
  - Concordance rate across callsets
  - Technology-specific F1 scores
  
- **Interpretability**:
  - Stratification breakdown (GIAB difficult regions, repeats, etc.)
  - Performance heatmap by genomic context
  
- **Discrimination**:
  - Error detection rate = `comparison_error verdicts / total_discrepancies`
  - Should be >>0 for fit-for-purpose benchmarks
  
- **Evidence**:
  - Technology support matrix (variants supported by N technologies)
  - Cross-technology concordance

### 3. Reporting System

**Purpose**: Generate comprehensive markdown reports with plots and actionable summaries.

#### Components

##### 3.1 Parameter Sweep Report (`scripts/report_param_sweep.py`)
```bash
./scripts/report_param_sweep.py \
  --sweep-dir results/sweeps/20260730_hg002_smvar_opt/ \
  --multi-callset-eval results/multi_callset_eval/ \
  --curation-verdicts curation/ \
  --output reports/20260730_hg002_smvar_opt.md
```

**Report Sections**:
1. **Executive Summary**: Top-N params, key findings, recommendations
2. **Parameter Comparison**: Table/plot comparing all params across metrics
3. **RIDE Analysis**: Breakdown by R/I/D/E dimensions
4. **Discrepancy Breakdown**: FN/FP by stratification and curation verdict
5. **Technology Support**: Cross-platform validation results
6. **Recommendations**: Which params to adopt, which to validate further

##### 3.2 Multi-Genome Validation Report (`scripts/report_multi_genome.py`)
```bash
./scripts/report_multi_genome.py \
  --genomes HG002,HG003,HG004,HG005 \
  --param-set dipcall-z5k \
  --results-dirs results/genomes/ \
  --output reports/z5k_validation_4genome.md
```

**Report Sections**:
1. **Cross-Genome Consistency**: Do params work well across genomes?
2. **Genome-Specific Issues**: Where do params fail per genome?
3. **RIDE Metrics**: Per-genome breakdown
4. **Production Readiness**: Are params safe to deploy?

### 4. Sequencing Data Integration

**Purpose**: Link BAM files for multi-technology manual curation.

#### BAM Registry (`config/sequencing_data.yml`)
```yaml
sequencing_data:
  HG002:
    illumina_pcr_free:
      - /defrabb_runs/bams/HG002/NovaSeq_PCRfree_35x.bam
      - /defrabb_runs/bams/HG002/HiSeq_PCRfree_50x.bam
    
    hifi:
      - /defrabb_runs/bams/HG002/PacBio_HiFi_Revio_40x.bam
      - /defrabb_runs/bams/HG002/PacBio_HiFi_Sequel2_30x.bam
    
    ont_ul:
      - /defrabb_runs/bams/HG002/ONT_UltraLong_Guppy6_50x.bam
    
    tenx:
      - /defrabb_runs/bams/HG002/10X_LinkedReads.bam
  
  HG003:
    # ... similar structure
```

#### Integration Points
- IGV session generator uses this registry
- Curation tracker validates technology tags against registry
- RIDE evidence metrics aggregate support across technologies

## Workflows

### Workflow 1: Coverage Expansion (v5.0q → v6.0)

```bash
# 1. Run parameter sweep on HG002
./scripts/generate_param_sweep.py config/sweeps/coverage_expansion.yml
./run_defrabb run -r 20260801_coverage_expansion

# 2. Multi-callset evaluation
./scripts/evaluate_multi_callset.py \
  --benchmark results/draft_benchmarksets/<best_bench_id>/ \
  --callsets config/comparison_callsets.yml \
  --sample HG002

# 3. Extract discrepancies for curation
./scripts/extract_discrepancies.py \
  --evaluation-dir results/multi_callset_eval/<bench_id>/ \
  --output discrepancies/coverage_expansion/

# 4. Manual curation (IGV + sequencing data)
./scripts/curate_discrepancies.py init --discrepancies ... --bams ...
# [User curates in IGV]
./scripts/curate_discrepancies.py record --verdict ...

# 5. Calculate RIDE metrics
./scripts/calculate_ride_metrics.py \
  --multi-callset-dir results/multi_callset_eval/<bench_id>/ \
  --curation-verdicts curation/coverage_expansion/verdicts.tsv

# 6. Generate report
./scripts/report_param_sweep.py \
  --sweep-dir results/sweeps/20260801_coverage_expansion/ \
  --output reports/coverage_expansion.md

# 7. Decision: If discrimination rate >X% and RIDE scores improved → adopt
```

### Workflow 2: Multi-Genome Parameter Definition

```bash
# 1. Optimize on HG002
./scripts/generate_param_sweep.py config/sweeps/hg002_optimization.yml
./run_defrabb run -r 20260801_hg002_opt

# 2. Score and select top-3 params
./scripts/score_param_sweep.py --top-n 3

# 3. Validate top params on HG003/HG004/HG005/HG007
./scripts/generate_validation_sweep.py \
  --genomes HG003,HG004,HG005,HG007 \
  --params z5k,z10k,pav-strict \
  --output config/sweeps/4genome_validation.yml

./run_defrabb run -r 20260801_4genome_validation

# 4. Multi-callset eval per genome
for genome in HG003 HG004 HG005 HG007; do
  ./scripts/evaluate_multi_callset.py \
    --benchmark results/draft_benchmarksets/<bench_id>/ \
    --callsets config/comparison_callsets.yml \
    --sample $genome
done

# 5. Sample-based curation (stratified random sample, ~100 variants per genome)
./scripts/sample_for_curation.py \
  --discrepancies discrepancies/4genome_validation/ \
  --n-per-genome 100 \
  --stratify-by genomic_region,variant_type

# 6. Multi-genome report
./scripts/report_multi_genome.py \
  --genomes HG002,HG003,HG004,HG005,HG007 \
  --param-set z5k \
  --output reports/z5k_validation.md

# 7. Decision: If consistent across genomes + RIDE OK → production params
```

## File Organization

```
defrabb/
├── config/
│   ├── comparison_callsets.yml       # Multi-callset registry
│   ├── sequencing_data.yml           # BAM file registry
│   └── sweeps/
│       ├── coverage_expansion.yml
│       └── 4genome_validation.yml
│
├── scripts/
│   ├── extract_discrepancies.py      # FN/FP extraction
│   ├── curate_discrepancies.py       # Manual curation tracking
│   ├── evaluate_multi_callset.py     # Multi-callset comparison
│   ├── calculate_ride_metrics.py     # RIDE scoring
│   ├── report_param_sweep.py         # Sweep summary reports
│   ├── report_multi_genome.py        # Cross-genome validation
│   └── sample_for_curation.py        # Stratified sampling
│
├── curation/
│   └── <eval_id>/
│       ├── verdicts.tsv              # Curation results
│       ├── igv_session.xml           # IGV batch script
│       └── notes.md                  # Free-form notes
│
├── discrepancies/
│   └── <eval_id>/
│       ├── all.tsv                   # All FN/FP variants
│       ├── fn_by_strat.tsv          # FN breakdown
│       └── fp_by_strat.tsv          # FP breakdown
│
├── results/
│   ├── multi_callset_eval/
│   │   └── <bench_id>/
│   │       ├── v5.0q/                # vs v5.0q
│   │       ├── v4.2.1/               # vs v4.2.1
│   │       ├── deepvariant_illumina/ # vs caller outputs
│   │       └── summary_metrics.tsv
│   │
│   └── sweeps/
│       └── <sweep_id>/
│           └── ride_metrics.json
│
└── reports/
    ├── <sweep_id>.md                 # Parameter sweep reports
    └── <validation_id>.md            # Multi-genome reports
```

## Implementation Phases

### Phase 1: Discrepancy Extraction (Week 1)
- [ ] `extract_discrepancies.py` - FN/FP extraction from hap.py outputs
- [ ] `curate_discrepancies.py init` - Initialize curation session
- [ ] Basic IGV session generation

### Phase 2: Curation Tracking (Week 1-2)
- [ ] `curate_discrepancies.py record` - Verdict recording
- [ ] Curation schema and TSV storage
- [ ] Integration with sequencing data registry

### Phase 3: Multi-Callset Framework (Week 2-3)
- [ ] `comparison_callsets.yml` schema + validation
- [ ] `evaluate_multi_callset.py` - Parallel evaluation runner
- [ ] Aggregated metrics collection

### Phase 4: RIDE Metrics (Week 3)
- [ ] `calculate_ride_metrics.py` - R/I/D/E scoring
- [ ] Technology support matrix
- [ ] Discrimination rate calculation

### Phase 5: Reporting (Week 4)
- [ ] `report_param_sweep.py` - Sweep summaries
- [ ] `report_multi_genome.py` - Cross-genome validation
- [ ] Plot generation (matplotlib/seaborn)

### Phase 6: Integration & Testing (Week 4-5)
- [ ] End-to-end workflow tests
- [ ] Documentation and examples
- [ ] Obsidian vault integration

## Open Questions

1. **Callset Availability**: Which comparison callsets are already available locally vs need download?
2. **BAM Storage**: Preferred local path for sequencing data? Size constraints?
3. **Curation Interface**: CLI-only or should we build a simple web UI for verdict recording?
4. **Reporting Format**: Markdown + static plots OK, or need interactive dashboards?
5. **Automation Level**: Should multi-callset eval be automatic in pipeline, or manual script invocation?

## Next Steps

1. Await user clarifying answers (see questions in main conversation)
2. Implement Phase 1 (discrepancy extraction) once requirements clarified
3. Validate against existing v0.023 param sweep results
4. Iterate based on user feedback

---

**Related Documents**:
- `docs/parameter-optimization.md` - User guide for v0.023 framework
- `docs/design/v0.023-parameter-optimization-design.md` - Original v0.023 architecture
- `docs/benchmark-validation-framework.md` - RIDE principles and validation workflow
