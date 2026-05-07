# DeFrABB Exclusion System Guide

## Overview

The exclusion system defines high-confidence benchmark regions by **subtracting**
problematic genomic areas from the diploid-assembled region BED file produced by
the variant caller (dipcall or PAV). The remaining regions form the benchmark set.

## How It Works

### Data Flow

```
Diploid assembled regions (dipcall .dip.bed or PAV callable_regions)
    ↓
Subtract exclusion BEDs (one at a time, via scripts/subtract_exclusions.py)
    ↓
Final benchmark BED (.benchmark.bed)
```

### Configuration

Exclusions are configured in `config/resources.yml` through several interrelated sections:

#### 1. Exclusion Sets (`exclusion_set`)

Named lists of exclusion region IDs. Each benchmark analysis row in `analyses.tsv`
references one set by name via the `exclusion_set` column. Use `none` for no exclusions.

```yaml
exclusion_set:
  HG008smvarv0.021:
    - segdups
    - tandem-repeats
    - satellites
    - gaps
    - flanks
    - svs-and-simple-repeats
    - VDJ
    - consecutive-svs
    - TSPY2-segdups
    - self-discrep
```

#### 2. Region Categories

Each exclusion region ID falls into one of three categories that determine
how it is sourced:

| Category | Config Key | Source | Example |
|----------|-----------|--------|---------|
| Assembly-agnostic | `exclusion_asm_agnostic` | Downloaded BED from URL in `references.{ref}.exclusions` | gaps, VDJ, HG002Q100-errors |
| Assembly-specific (ref-agnostic) | `exclusion_ref_agnostic` | Derived from variant calls/alignments | flanks, svs-and-simple-repeats, consecutive-svs, self-discrep, pav-inv |
| Assembly-intersect | `exclusion_asm_intersect` | Downloaded BED, but only intervals overlapping assembly breaks are excluded | segdups, tandem-repeats, satellites, self-chains |

#### 3. Region Processing

Before subtraction, each exclusion BED can be modified:

| Processing | Config Key | Effect | Default Values |
|-----------|-----------|--------|----------------|
| Slop | `exclusion_slop_regions` | Adds buffer around each interval | `_exclusion_params.slop` (15000 bp) |
| Slop + Merge | `exclusion_slopmerge_regions` | Adds buffer then merges nearby intervals | `_exclusion_params.slop` + `_exclusion_params.slopmerge_dist` (15000 + 10000 bp) |
| Assembly intersect | `exclusion_asm_intersect` | Only excludes intervals overlapping assembly alignment breaks (start/end positions) | N/A |

#### 4. Tunable Parameters (`_exclusion_params`)

All numeric parameters are in `config/resources.yml` under `_exclusion_params`:

| Parameter | Default | Used By | Purpose |
|-----------|---------|---------|---------|
| `slop` | 15000 | `add_slop`, `add_slop_and_merge` | Buffer around reference exclusion regions |
| `slopmerge_dist` | 10000 | `add_slop_and_merge` | Merge distance after slop |
| `flank_bases` | 15000 | `get_flanks` | Buffer around assembly alignment breaks |
| `sv_repeat_slop` | 50 | `intersect_SVs_and_simple_repeats` | Buffer for SV+repeat intersection |
| `sv_repeat_merge_dist` | 1000 | `intersect_SVs_and_simple_repeats` | Merge distance for SV+repeat regions |
| `self_discrep_max_indel` | 50 | `self_discrep_extract_fpfns` | Max indel size for self-comparison FP/FN |
| `self_discrep_slop` | 50 | `self_discrep_intersect_slop` | Buffer for self-discrepancy regions |
| `self_discrep_merge_dist` | 1000 | `self_discrep_intersect_slop` | Merge distance for self-discrepancy |
| `pav_inv_slop` | 50 | `exclude_pav_inversions` | Buffer for PAV inversion exclusions |
| `pav_inv_merge_dist` | 1000 | `exclude_pav_inversions` | Merge distance for PAV inversions |
| `consecutive_sv_min_bp` | 10 | `get_consecutive_svs` | Min size for consecutive DEL/INS events |
| `gap_min_ns` | 50 | `make_gaps_bed` | Min N-stretch length for gap detection |

## Creating a New Exclusion Set for a New Genome

### Step 1: Determine which exclusion regions apply

Start with a base set. For small variants, a typical starting point:

```yaml
  NewGenome_smvar:
    - segdups
    - tandem-repeats
    - satellites
    - gaps
    - flanks
    - svs-and-simple-repeats
    - VDJ
    - consecutive-svs
    - TSPY2-segdups
    - self-discrep
```

For structural variants, similar to small variant except svs-and-simple-repeats are not excluded:

```yaml
  NewGenome_stvar:
    - segdups
    - tandem-repeats
    - satellites
    - gaps
    - flanks
    - VDJ
    - consecutive-svs
    - TSPY2-segdups
    - self-discrep
```

### Step 2: Add genome-specific exclusions (if needed)

If the assembly has known issues (errors, mosaicism, tool-specific artifacts),
create BED files, host them, and add URLs to `references.{ref}.exclusions`
in `resources.yml`. Then add the region IDs to your exclusion set and
to `exclusion_asm_agnostic`.

### Step 3: Configure region categories

Ensure each new exclusion ID appears in exactly one category:
- `exclusion_asm_agnostic` — for pre-computed BEDs
- `exclusion_ref_agnostic` — for assembly-derived regions
- `exclusion_asm_intersect` — for BEDs filtered by assembly breaks

### Step 4: Test with a subset

Run the pipeline with the new exclusion set on a small region or chr21 subset
to verify BED paths resolve correctly:

```bash
snakemake --dryrun --cores 1 -a config/your_analyses.tsv
```

### Step 5: Tune parameters

Override `_exclusion_params` values on the command line for experimentation:

```bash
snakemake --cores 4 --config '_exclusion_params={slop: 20000, slopmerge_dist: 15000, ...}'
```

Or modify `resources.yml` directly for persistent changes.

## Exclusion Rules Reference

| Rule | Input | Output | Purpose |
|------|-------|--------|---------|
| `download_bed_gz` | URL from config | `resources/exclusions/{ref}/region.bed` | Download pre-computed exclusion BEDs |
| `make_gaps_bed` | Reference FASTA | `resources/exclusions/{ref}/gaps.bed` | N-stretches in reference |
| `add_slop` | Exclusion BED | `{region}_slop.bed` | Add buffer around intervals |
| `add_slop_and_merge` | Exclusion BED | `{region}_slopmerge.bed` | Add buffer + merge nearby |
| `intersect_start_and_end` | Baseline BED + exclusion BED | `{region}_start_sorted.bed`, `{region}_end_sorted.bed` | Assembly-break-only exclusion |
| `get_flanks` | Baseline BED | `_flanks.bed` | Buffer around non-diploid regions |
| `get_sv_widen_coords` | Processed VCF | `_SVs.bed` | SV coordinates with TR widening |
| `intersect_SVs_and_simple_repeats` | SV BED + TR BED | `_svs-and-simple-repeats.bed` | Combined SV+repeat exclusion |
| `get_consecutive_svs` | Hap1/Hap2 BAMs | `_consecutive-svs.bed` | DEL+INS artifacts from alignment |
| `self_discrep_happy` | VCF vs itself | FP/FN VCF | Small variant self-comparison discrepancies (hap.py) |
| `self_discrep_extract_fpfns` | Self-discrep VCF | `.fpfns.bed` | Extract discrepant regions (indels ≤ 50bp) |
| `self_discrep_intersect_slop` | FP/FN BED + TR BED | `_self-discrep.bed` | Expand discrepancies in repeats |
| `exclude_pav_inversions` | Standardized VCF + segdups | `_pav-inv.bed` | PAV inversion + segdup overlap |
| `subtract_exclusions` | Baseline BED + all exclusion BEDs | `.benchmark.bed` | Final subtraction |
| `generate_intersection_summary` | Baseline BED + exclusion BEDs | `.exclusion_intersection_summary.csv` | Overlap statistics |

## Known Limitations

### Self-discrepancy is small-variant only

The current `self-discrep` exclusion uses hap.py for self-comparison, which only
handles small variants. The `self_discrep_extract_fpfns` rule further filters to
indels ≤ 50bp. Ideally, self-discrepancies should be variant-type specific:

- **smvar benchmarks**: Use hap.py-based self-discrepancy (current implementation)
- **stvar benchmarks**: Would need a truvari-based self-comparison to detect SV
  self-discrepancies — **not yet implemented**

Despite this, `self-discrep` currently appears in both smvar and stvar exclusion
sets. For stvar benchmarks, the hap.py-based self-discrepancy may still catch some
relevant regions, but a truvari-based approach would be more appropriate.
