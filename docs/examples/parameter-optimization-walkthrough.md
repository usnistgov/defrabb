# Parameter Optimization Walkthrough

This document walks through a complete parameter optimization example using the v0.023 optimization framework.

## Overview

The optimization workflow has three phases:

1. **Generate sweep** - Create analyses table with parameter cross-products
2. **Run pipeline** - Execute the sweep
3. **Score results** - Rank parameters and identify winners

## Example: Optimizing dipcall window size for HG002

### Step 1: Create sweep configuration

Create `config/sweeps/hg002_dipcall_opt.yml`:

```yaml
name: HG002 Dipcall Window Optimization
output: config/analyses_20260723_v0.023_hg002_dipcall_opt.tsv

# Fixed across all analyses
fixed:
  ref_id: GRCh38
  asm_id: HG2-T2TQ100-V1.1
  vc_cmd: dipcall
  bench_type: smvar
  vcf_processing: trf
  exclusion_set: standard
  eval_cmd: happy
  eval_comp_id: v5.0q-smvar
  eval_params: default

# Sweep dimensions (cross-product)
sweep:
  vc_param_id: [z2k, z5k, z10k]

# 3 analyses, 3 unique variant calls
# Estimated runtime: ~16.5 hours (3 × 5h dipcall + 3 × 0.5h happy)
```

### Step 2: Generate sweep

```bash
./scripts/generate_param_sweep.py config/sweeps/hg002_dipcall_opt.yml
```

Output:
```
Sweep: HG002 Dipcall Window Optimization
  Analyses: 3
  Unique variant call runs: 3
  Reuse factor: 1.0x
  Estimated runtime: 16.5 hours
  Estimated storage: 165.0 GB

Wrote 3 analyses to: config/analyses_20260723_v0.023_hg002_dipcall_opt.tsv
```

### Step 3: Review generated table

```bash
head -5 config/analyses_20260723_v0.023_hg002_dipcall_opt.tsv
```

```
vc_id	ref_id	asm_id	vc_cmd	vc_param_id	bench_type	exclusion_set	vcf_processing	eval_cmd	eval_comp_id	eval_params	vc_params	bench_id	eval_id
vc001	GRCh38	HG2-T2TQ100-V1.1	dipcall	z2k	smvar	standard	trf	happy	v5.0q-smvar	default		GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k_standard_trf	happy_v5.0q-smvar
vc002	GRCh38	HG2-T2TQ100-V1.1	dipcall	z5k	smvar	standard	trf	happy	v5.0q-smvar	default		GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z5k_standard_trf	happy_v5.0q-smvar
vc003	GRCh38	HG2-T2TQ100-V1.1	dipcall	z10k	smvar	standard	trf	happy	v5.0q-smvar	default		GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z10k_standard_trf	happy_v5.0q-smvar
```

### Step 4: Run pipeline

```bash
./run_defrabb run -r 20260723_v0.023_hg002_dipcall_opt \
  --config analyses=config/analyses_20260723_v0.023_hg002_dipcall_opt.tsv
```

Or using the clone-into-runid pattern (recommended):

```bash
cd ..
git clone code_base/defrabb 20260723_v0.023_hg002_dipcall_opt
cd 20260723_v0.023_hg002_dipcall_opt
./run_defrabb run -r 20260723_v0.023_hg002_dipcall_opt \
  --config analyses=config/analyses_20260723_v0.023_hg002_dipcall_opt.tsv
```

### Step 5: Score results

After the pipeline completes:

```bash
./scripts/score_param_sweep.py \
  --results-dir 20260723_v0.023_hg002_dipcall_opt/results \
  --baseline "dipcall-z2k" \
  --rank-by f1 \
  --top-n 3
```

Output:
```markdown
# Parameter Sweep Scoring Summary

**Baseline:** dipcall-z2k
**Ranked by:** f1

## Top-3 Parameter Sets (SNP f1)

### 1. dipcall-z5k
- **F1:** 0.999870 (+0.000009)
- **Precision:** 0.999872 (+0.000010)
- **Recall:** 0.999868 (+0.000008)
- **Status:** improved

### 2. dipcall-z10k
- **F1:** 0.999865 (+0.000004)
- **Precision:** 0.999866 (+0.000004)
- **Recall:** 0.999864 (+0.000003)
- **Status:** improved

### 3. dipcall-z2k (baseline)
- **F1:** 0.999861 (baseline)
- **Status:** baseline

## Recommendations

**Winner:** dipcall-z5k shows best F1 performance (+0.9 bp per 100k)

For validation, test z5k across 4 genomes (HG003, HG004, HG005, HG007).
```

### Step 6: Multi-genome validation (optional)

Create validation sweep:

```yaml
name: Dipcall z5k Validation
output: config/analyses_20260723_v0.023_z5k_validation.tsv

fixed:
  ref_id: GRCh38
  vc_cmd: dipcall
  vc_param_id: z5k  # Winner from optimization
  bench_type: smvar
  vcf_processing: trf
  exclusion_set: standard
  eval_cmd: happy
  eval_comp_id: v5.0q-smvar
  eval_params: default

sweep:
  asm_id: [HG2-T2TQ100-V1.1, HG003v1.1, HG004v1.1, HG005v1.1]
```

## Advanced: Multi-dimensional sweeps

### Example: Optimize dipcall + exclusions + VCF processing

```yaml
name: HG002 Multi-Parameter Optimization
output: config/analyses_20260723_v0.023_hg002_multi_opt.tsv

fixed:
  ref_id: GRCh38
  asm_id: HG2-T2TQ100-V1.1
  vc_cmd: dipcall
  bench_type: smvar
  eval_cmd: happy
  eval_comp_id: v5.0q-smvar
  eval_params: default

sweep:
  vc_param_id: [z2k, z5k, z10k]
  exclusion_set: [conservative, standard, aggressive]
  vcf_processing: [trf, xy_trf]

# Cross-product: 3 × 3 × 2 = 18 analyses
# Unique variant calls: 3 (reuse across exclusion/vcf combos)
# Reuse factor: 6x
# Estimated runtime: ~24 hours
```

Key benefit: **6x reuse factor** - dipcall runs only 3 times, outputs reused across 6 different benchmark/eval combinations each.

## Cost estimation

The sweep generator provides runtime and storage estimates:

```
Sweep: HG002 Multi-Parameter Optimization
  Analyses: 18
  Unique variant call runs: 3
  Reuse factor: 6.0x
  Estimated runtime: 24.0 hours
  Estimated storage: 240.0 GB
```

**Without reuse:** 18 separate dipcall runs = 90 hours  
**With reuse:** 3 dipcall runs = 15 hours + 9 hours eval = 24 hours  
**Speedup:** 3.75x

## Parameter profiles reference

### Dipcall (`vc_param_id`)

- `z2k` - minimap2 `-z 2000` (default, balanced)
- `z5k` - minimap2 `-z 5000` (wider window, more sensitive)
- `z10k` - minimap2 `-z 10000` (widest, most sensitive, slower)
- `z1k` - minimap2 `-z 1000` (narrower, faster, less sensitive)

### PAV (`vc_param_id`)

- `giab` - Default GIAB merge strategy
- `strict` - Stricter merge thresholds
- `lenient` - More permissive merging

### Exclusions (`exclusion_set`)

- `standard` - Default buffer distances
- `conservative` - Larger buffers (more exclusions)
- `aggressive` - Smaller buffers (fewer exclusions)

### VCF Processing (`vcf_processing`)

- `trf` - TRF annotation only
- `xy_trf` - XY GT fix + TRF annotation
- `default` - No processing

See `config/resources.yml` for complete definitions.

## Troubleshooting

### "Duplicate (eval_id, bench_id) combinations"

This means your sweep creates non-unique analyses. Common causes:

- Sweeping a dimension that doesn't affect `bench_id` (fixed in v0.023)
- Multiple identical parameter combinations
- Check the dry-run output to see what's being generated

### "Baseline is ambiguous"

When using `score_param_sweep.py`, provide a specific baseline substring that uniquely identifies one analysis:

```bash
# Too broad
--baseline "v5.0q"  # Matches all comparisons

# Specific
--baseline "GRCh38_v5.0q-smvar_HG2-T2TQ100-V1.1_smvar_dipcall-z2k"
```

### Reuse not working as expected

Check `vc_id` assignments in the generated table. Same `(ref_id, asm_id, vc_cmd, vc_param_id)` → same `vc_id` → reuse.

Different exclusions or VCF processing profiles do NOT require re-running variant calling.
