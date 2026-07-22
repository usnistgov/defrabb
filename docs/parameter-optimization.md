# Parameter Optimization Guide

**Version:** v0.023+  
**Status:** Production-ready  
**Framework:** DeFrABB parameter sweep and scoring infrastructure

---

## Overview

DeFrABB's parameter optimization framework enables systematic exploration of variant calling, exclusion, and VCF processing parameters across single or multiple genomes. Key capabilities:

- **Declarative sweep generation** - YAML configs → analyses tables with cross-products
- **Output reuse** - Share expensive variant calls (3-6 hours each) across sweeps
- **Baseline comparison** - Score against reference benchmarks (e.g., HG002 v5q)
- **Multi-genome validation** - Test optimized parameters across cohorts
- **Cost estimation** - Preview runtime and storage before execution

---

## Quick Start

### 1. Define a Parameter Sweep

Create a YAML config describing fixed dimensions and sweep dimensions:

```yaml
# config/sweeps/my_sweep.yml
name: HG002 Small Variant Optimization
output: config/analyses_20260722_my_sweep.tsv

fixed:
  ref_id: GRCh38
  asm_id: HG2-T2TQ100-V1.1
  vc_cmd: dipcall
  bench_type: smvar
  eval_cmd: happy
  eval_comp_id: v5.0q-smvar

sweep:
  vc_param_id: [z2k, z5k, z10k]
  exclusion_set: [conservative, standard, aggressive]
  vcf_processing: [trf, xy_trf]
```

This generates 3 × 3 × 2 = **18 analyses** with only **3 unique variant calls** (6x reuse).

### 2. Generate Analyses Table

```bash
# Dry-run preview
scripts/generate_param_sweep.py config/sweeps/my_sweep.yml --dry-run

# Generate table
scripts/generate_param_sweep.py config/sweeps/my_sweep.yml

# Output: config/analyses_20260722_my_sweep.tsv
```

### 3. Run the Sweep

```bash
./run_defrabb run -r 20260722_my_sweep \
  --analyses config/analyses_20260722_my_sweep.tsv
```

### 4. Score Results

```bash
scripts/score_param_sweep.py \
  --results-dir 20260722_my_sweep \
  --baseline v5.0q \
  --rank-by f1 \
  --top-n 3 \
  --out sweep_summary.md
```

Output identifies top-3 parameter sets, regressions, and recommendations.

---

## Parameter Profiles

### Dipcall Parameters

Minimap2 alignment window settings (format: `-z<max_chain_skip>,<max_chain_iter>`):

| Profile | Flags | Use Case |
|---------|-------|----------|
| `z2k` (default) | `-z200000,10000` | Production default, balanced |
| `z5k` | `-z500000,5000` | Wider window, better long-range phasing |
| `z10k` | `-z1000000,10000` | Extreme sensitivity, slow |
| `z1k` | `-z100000,20000` | Narrow window, better small variant recall |

**Usage:** Set `vc_param_id` column in analyses table or sweep config.

### PAV Parameters

PAV variant merge strategies (format: `nr::<reciprocal_overlap>`):

| Profile | Insertions | Deletions | Inversions | Use Case |
|---------|-----------|-----------|------------|----------|
| `giab` (default) | `nr::exact` | `nr::exact` | `nr::exact:ro(0.5):szro(0.5,200):match` | GIAB production |
| `strict` | `nr::szro(0.9,100)` | `nr::szro(0.9,100)` | `nr::exact:ro(0.7):szro(0.7,200):match` | High precision |
| `lenient` | `nr::szro(0.5,200)` | `nr::szro(0.5,200)` | `nr::ro(0.3):szro(0.3,300):match` | High recall |

**Usage:** Add `_pav_config` profile name to `config/resources.yml`, reference in sweep.

### Exclusion Profiles

Buffer distances for exclusion region processing (basepairs):

| Profile | Slop | Merge Dist | Flank Bases | Use Case |
|---------|------|------------|-------------|----------|
| `conservative` | 20000 | 15000 | 20000 | Maximize precision |
| `standard` (default) | 15000 | 10000 | 15000 | Production default |
| `aggressive` | 10000 | 5000 | 10000 | Maximize recall |

**Usage:** Set `exclusion_set` column in analyses table or sweep config.

### VCF Processing Profiles

Post-calling VCF transformations (defined in `config/resources.yml:152-192`):

| Profile | Steps | Use Case |
|---------|-------|----------|
| `trf` | TRF annotation only | Small variants, no XY |
| `xy_trf` | Diploidize XY + TRF | Autosomes + sex chromosomes |
| `norm_xy_trf_sv_lcr` | Normalize + XY + TRF + SV + LCR | Full processing (default for stvar) |

---

## Sweep Generator Reference

### Command-Line Options

```bash
scripts/generate_param_sweep.py <config.yml> [options]

Options:
  -o, --output PATH       Override output path (default: from config.output)
  -n, --dry-run           Preview without writing
  --max-analyses N        Safety limit (default: 200)
  --preview-rows N        Rows to show in dry-run (default: 5)
```

### YAML Config Format

```yaml
name: <Descriptive sweep name>
output: <Path to output analyses TSV>

fixed:
  # Fields common to all analyses
  ref_id: GRCh38
  asm_id: HG002
  # ... any analyses.tsv column

sweep:
  # Fields to cross-product
  vc_param_id: [z2k, z5k, z10k]
  exclusion_set: [standard, aggressive]
  # ... any analyses.tsv column

# Optional constraints (future enhancement)
constraints:
  - if: {bench_type: stvar}
    then: {vcf_processing: [norm_xy_trf_sv_lcr]}
```

### Output Reuse Intelligence

The generator assigns shared `vc_id` to analyses differing only in downstream parameters:

```
Input:
  sweep:
    vc_param_id: [z2k, z5k]
    exclusion_set: [standard, aggressive]
    
Cross-product: 2 × 2 = 4 analyses
Unique vc_ids: 2 (z2k, z5k share across exclusion sets)
Reuse factor: 4 / 2 = 2x
Time saved: ~2 × 5 hours = 10 hours
```

Snakemake automatically reuses outputs when `vc_id` matches.

### Cost Estimation

Runtime and storage projections based on historical averages:

| Operation | Runtime (hours) | Storage (GB) |
|-----------|----------------|--------------|
| dipcall | 5.0 | 50 |
| PAV | 6.0 | 80 |
| hap.py | 0.5 | 5 |
| Truvari | 1.0 | 10 |

Example output:
```
Sweep: HG002 Small Variant Optimization
  Analyses: 18
  Unique variant call runs: 3
  Reuse factor: 6.0x
  Estimated runtime: 24 hours
  Estimated storage: 340 GB
```

---

## Scoring and Comparison

### Basic Comparison

```bash
scripts/compare_evaluations.py \
  --results-dir 20260722_sweep \
  --baseline v5.0q \
  --metric f1 \
  --threshold 0.0001 \
  --out comparison.tsv
```

Output columns:
- `analysis_id`, `variant_type`, `precision`, `recall`, `f1`
- `delta_precision`, `delta_recall`, `delta_f1`
- `status`: `baseline` | `improved` | `regressed` | `same`

### Parameter Optimization Workflow

```bash
scripts/score_param_sweep.py \
  --results-dir 20260722_sweep \
  --baseline v5.0q \
  --rank-by f1 \
  --top-n 3 \
  --out summary.md
```

Generates:
- Top-N parameter sets by metric
- Delta vs baseline for each metric
- Regression warnings
- Validation recommendations

---

## Multi-Genome Workflow

### Scenario: 5-Genome Parameter Optimization

Optimize parameters on HG002 (has v5q baseline), then validate on 4 other genomes.

#### Stage 1: HG002 Parameter Sweep

```yaml
# config/sweeps/hg002_optimization.yml
name: HG002 Optimization
output: config/analyses_hg002_opt.tsv

fixed:
  ref_id: GRCh38
  asm_id: HG2-T2TQ100-V1.1
  eval_comp_id: v5.0q-smvar

sweep:
  vc_param_id: [z2k, z5k, z10k]
  exclusion_set: [conservative, standard, aggressive]
  vcf_processing: [trf, xy_trf]
```

Generate and run:
```bash
scripts/generate_param_sweep.py config/sweeps/hg002_optimization.yml
./run_defrabb run -r hg002_opt --analyses config/analyses_hg002_opt.tsv
```

#### Stage 2: Score and Select Winners

```bash
scripts/score_param_sweep.py \
  --results-dir hg002_opt \
  --baseline v5.0q \
  --rank-by f1 \
  --top-n 3 \
  --out hg002_winners.md
```

Suppose top-3 are:
1. `z5k` + `aggressive` + `xy_trf`
2. `z2k` + `standard` + `trf`
3. `z10k` + `conservative` + `xy_trf`

#### Stage 3: 4-Genome Validation

```yaml
# config/sweeps/validation_4genome.yml
name: 4-Genome Validation
output: config/analyses_4genome_val.tsv

fixed:
  ref_id: GRCh38
  bench_type: smvar
  eval_cmd: happy
  vcf_processing: xy_trf  # From winner

sweep:
  asm_id:
    - HG008N-T2TQ100
    - HG008T-T2TQ100
    - NA12878-T2T
    - HG005-T2T

  vc_param_id: [z5k, z2k, z10k]  # Top-3 from HG002
  exclusion_set: [aggressive, standard, conservative]  # Top-3 from HG002

  eval_comp_id:
    - HG008N-comparison
    - HG008T-comparison
    - NA12878-comparison
    - HG005-comparison
```

**Note:** This cross-product generates 4×3×3×4 = 144 analyses. Post-process to match `asm_id` with corresponding `eval_comp_id` (4 genomes × 3 vc_params × 3 exclusions = 36 analyses).

#### Stage 4: Cross-Genome Analysis

```bash
scripts/compare_evaluations.py \
  --results-dir 4genome_val \
  --regions \
  --out cross_genome.tsv

# Group by parameter set, identify consistent winner
awk -F'\t' 'NR>1 {print $1, $4, $7}' cross_genome.tsv | sort
```

Look for:
- Parameter set with best cross-genome F1
- Genome-specific outliers
- Stratifications with high discrepancy rates

---

## Tips and Best Practices

### Dry-Run Always

```bash
scripts/generate_param_sweep.py config/sweeps/my_sweep.yml --dry-run
```

Catches:
- Exploding cross-products (e.g., 1000+ analyses)
- Missing required fields
- Cost estimation surprises

### Safety Limits

Default max is 200 analyses. Override with:
```bash
scripts/generate_param_sweep.py ... --max-analyses 500
```

But consider splitting into phases:
- Phase 1: Broad sweep (3-5 levels per dimension)
- Phase 2: Narrow sweep around winners

### Profile Naming Conventions

- Dipcall: `z<window_size>` (e.g., `z5k`, `z10k`)
- PAV: Descriptive (`strict`, `lenient`, `giab`)
- Exclusions: Descriptive (`conservative`, `standard`, `aggressive`)

### Reuse Factor Targets

- Good: 3-5x (typical for 2-3 sweep dimensions)
- Great: 6-10x (downstream-heavy sweeps)
- Overkill: 1x (each analysis unique variant call)

### Baseline Selection

For HG002: `--baseline v5.0q` (matches `v5.0q-smvar` and `v5.0q-stvar`)

For other genomes: Use published benchmark or high-confidence comparison callset.

### Metric Selection

- **F1** - Balanced precision/recall (default for most)
- **Precision** - When false positives are costly (clinical)
- **Recall** - When false negatives are costly (research)

---

## Example Workflows

### Single-Genome Dipcall Optimization

```bash
# 1. Create sweep config
cat > config/sweeps/dipcall_opt.yml <<EOF
name: Dipcall Window Optimization
output: config/analyses_dipcall_opt.tsv
fixed: {ref_id: GRCh38, asm_id: HG002, vc_cmd: dipcall}
sweep: {vc_param_id: [z1k, z2k, z5k, z10k]}
EOF

# 2. Generate
scripts/generate_param_sweep.py config/sweeps/dipcall_opt.yml --dry-run
scripts/generate_param_sweep.py config/sweeps/dipcall_opt.yml

# 3. Run
./run_defrabb run -r dipcall_opt --analyses config/analyses_dipcall_opt.tsv

# 4. Score
scripts/score_param_sweep.py --results-dir dipcall_opt --baseline v5.0q
```

### Exclusion Buffer Sweep

```bash
# Sweep exclusion profiles only
cat > config/sweeps/exclusion_opt.yml <<EOF
name: Exclusion Buffer Optimization
output: config/analyses_exclusion_opt.tsv
fixed:
  ref_id: GRCh38
  asm_id: HG002
  vc_cmd: dipcall
  vc_param_id: z2k
sweep:
  exclusion_set: [conservative, standard, aggressive]
  vcf_processing: [trf, xy_trf]
EOF
```

### Full Cross-Product Sweep

```bash
# WARNING: 3×3×2 = 18 analyses
cat > config/sweeps/full_sweep.yml <<EOF
name: Full Parameter Sweep
output: config/analyses_full_sweep.tsv
fixed: {ref_id: GRCh38, asm_id: HG002}
sweep:
  vc_param_id: [z2k, z5k, z10k]
  exclusion_set: [conservative, standard, aggressive]
  vcf_processing: [trf, xy_trf]
EOF
```

---

## Troubleshooting

### "Cross-product exceeds limit"

Reduce sweep dimensions or increase `--max-analyses`:
```bash
scripts/generate_param_sweep.py ... --max-analyses 500
```

### "Baseline not found"

Check available analysis_ids:
```bash
ls results/evaluations/happy/*/
# or
scripts/compare_evaluations.py --results-dir results --out tmp.tsv
cut -f1 tmp.tsv | sort -u
```

### Unexpected Reuse Factor

Check vc_id assignments in generated TSV:
```bash
cut -f1 config/analyses_my_sweep.tsv | sort | uniq -c
```

Each vc_id should appear multiple times if reuse is working.

### Slow Sweep Generation

Normal - cost estimation reads configs and does math. Typical: <2 seconds for 100 analyses.

---

## Implementation Notes

### Output Reuse Mechanism

Snakemake reuses outputs when `vc_id` matches:
```
vc_id = f"{ref}_{asm}_{vc_cmd}-{vc_param_id}"

Same vc_id → same outputs → no re-run
```

The sweep generator assigns shared vc_ids to analyses differing only in:
- `exclusion_set`
- `vcf_processing`
- `eval_cmd`, `eval_comp_id`, `eval_params`

### Schema Validation

All profiles validated by `schema/resources-schema.yml`:
- Dipcall: Must match `-z\d+,\d+` pattern
- PAV: Required `merge_ins`, `merge_del`, `merge_inv` fields
- Exclusions: Integer values ≥ 0 for all buffer/merge params

### Design Documentation

Full design (448 lines) in `docs/design/v0.023-parameter-optimization-design.md`:
- Adversarial review (over/under-engineering critique)
- 4-component architecture
- Multi-genome workflow
- Implementation estimates

---

## Future Enhancements

- **Constraint enforcement** - Filter invalid combinations in YAML
- **Multi-reference sweeps** - GRCh37/38/T2T cross-product
- **Auto-tuning** - Bayesian optimization over param space
- **WGS integration** - Platform data for error characterization
- **Region refinement** - Discrepancy hotspots → exclusion updates

---

## References

- **Design:** `docs/design/v0.023-parameter-optimization-design.md`
- **Sweep Generator:** `scripts/generate_param_sweep.py --help`
- **Scoring Script:** `scripts/score_param_sweep.py --help`
- **Comparison Tool:** `scripts/compare_evaluations.py --help`
- **Schema:** `schema/resources-schema.yml`

---

*Document version: v0.023 (2026-07-22)*  
*Framework status: Production-ready*
