# Method overview

This page describes what DeFrABB does at each stage and why. It is aimed at
readers who want to understand how a DeFrABB benchmark (for example GIAB HG002
v5.0q) was built. For the scientific context and evaluation of the HG002
benchmark, see the preprint:
<https://www.biorxiv.org/content/10.64898/2026.09.23.752440v1>.

DeFrABB has three components:

1. assembly-based variant calling
2. draft benchmark generation
3. evaluation

Each run is driven by an analyses table. Each row of the table defines one
combination of assembly, reference, variant caller, benchmark type, VCF
processing, exclusion set, and evaluation (see
[configuration.md](configuration.md)).

## Inputs

- **Diploid assembly:** maternal and paternal haplotype FASTAs (for HG002 v5.0q,
  the T2T HG002 Q100 v1.1 assembly). The `is_male` setting determines how chrX
  and chrY are handled.
- **Reference genome:** GRCh37, GRCh38, or T2T-CHM13v2.0, with a
  pseudoautosomal-region (PAR) BED.
- **Exclusion region BEDs:** mostly GIAB genome stratifications (segmental
  duplications, tandem repeats, satellites, VDJ), plus genome-specific BEDs such
  as known assembly errors.
- **Comparison callsets:** existing benchmarks or high-quality callsets used to
  evaluate the drafts.

All inputs are referenced by URL in `config/resources.yml` and downloaded by the
workflow.

## 1. Assembly-based variant calling

Each haplotype is aligned to the reference with an assembly-to-assembly aligner
(minimap2). Variants are called from the alignments. The caller also reports
**diploid regions**: regions where both haplotypes have a single, 1:1 alignment
to the reference. Only in these regions can the assembly support confident
homozygous-reference, heterozygous, and homozygous-variant calls. For a male
sample, chrX and chrY outside the PARs are called from the single haplotype that
carries them.

Supported callers:

- **[dipcall](https://github.com/lh3/dipcall)**: used for HG002 v5.0q. The
  minimap2 `-z` (Z-drop) parameter is configurable through named profiles. The
  production profile `z2k` (`-z200000,10000`) improves alignment through SVs and
  the MHC.
- **[PAV](https://github.com/BeckLaboratory/pav3)** (PAV3): run in a container.
  Only `FILTER=PASS` calls are kept. PAV calls are also used to derive
  exclusions (see below).

Both callers' outputs are standardized to a common VCF and BED layout, so the
downstream steps do not depend on the caller.

## 2. Draft benchmark generation

### VCF processing

Each benchmark applies a named, ordered sequence of processing steps
(`vcf_processing_profiles` in `config/resources.yml`):

- `norm`: split multiallelic records, remove exact duplicates, and
  left-normalize against the reference (bcftools).
- `fix_XY_genotype`: convert genotypes in non-PAR chrX/chrY of male samples to
  haploid representation, so benchmarking tools compare them correctly. The
  change is recorded in the VCF header.
- `trfanno`: annotate tandem-repeat context (`truvari anno trf`, using the
  adotto TR catalog). Very large insertions (>100 kb by default) skip this step,
  because they cause extreme runtimes. They stay in the VCF without the
  annotation.
- `svinfo`, `repmask`, `lcr`, `remap`: further Truvari annotations. They cover
  SV type and length, RepeatMasker classification, low-complexity content, and
  whether an inserted or deleted sequence maps elsewhere in the reference.
- `end_info`, `no_inv`: add `END` tags; remove symbolic inversion records.

Annotations do not change which variants are in the benchmark. They let users
stratify results and help interpret errors.

Small-variant (`smvar`) and structural-variant (`stvar`) benchmarks are built
from the same calls. They differ in VCF processing profile and exclusion set.
For `stvar` benchmarks, the `_bench-vars.vcf.gz` output keeps only variants of
50 bp and larger.

### Benchmark regions and exclusions

The benchmark regions are the diploid regions minus a configured **exclusion
set**. Exclusions remove genomic contexts where the benchmark would not be
reliable. That happens when the assembly, the assembly-based variant calls, the
reference, the VCF representation, or the benchmarking tools are likely to
produce incorrect comparisons. Exclusion types:

- **Gaps and flanks:** reference gaps (N-stretches detected in the reference
  FASTA), and a buffer (15 kb by default) around breaks in the
  assembly-to-reference alignment.
- **Large repeats with alignment breaks:** segmental duplications, satellites,
  and large tandem repeats are excluded only where an alignment break falls
  inside them. Nearby segmental duplications and satellites are merged first, so
  a break anywhere in a cluster removes the whole cluster.
- **SVs in small-variant benchmarks:** regions containing SVs, together with
  overlapping simple repeats, are removed from `smvar` benchmark regions. Small
  variants near SVs have ambiguous representations.
- **Consecutive SVs:** adjacent deletions and insertions that typically reflect
  alignment artifacts.
- **Caller discrepancies:** regions where dipcall and PAV disagree. Small
  variants are compared with hap.py/vcfeval and SVs with `truvari bench`. (The
  HG002 v5.0q exclusion BEDs were generated with `truvari bench` followed by
  MAFFT-based `refine`.) This also covers inversions called by PAV, which
  dipcall represents inconsistently.
- **Self-discrepancies:** the draft benchmark is compared against itself with
  the benchmarking tool. Any false positives or negatives come from
  representations the tool cannot match, and those regions are excluded.
- **Genome-specific exclusions:** user-provided BEDs, e.g. known HG002 assembly
  errors, putative mosaic variants, manually identified dipcall bugs, and VDJ
  and TSPY2 regions.

Slop and merge distances are configurable through named exclusion profiles.
Every benchmark records how much sequence each exclusion removed
(`exclusion_stats.txt`) and the exact inputs and parameters used
(`exclusion_provenance.yml`). For configuration details, see the
[exclusion system guide](exclusion_system_guide.md).

By design, the remaining regions omit very large or complex SVs, CNVs, and
regions without a 1:1 assembly alignment. No standards currently exist for
representing and comparing variants in these regions.

## 3. Evaluation

Each draft benchmark is compared against one or more comparison callsets:

- **Small variants:** [hap.py](https://github.com/Illumina/hap.py) with the
  vcfeval engine and GIAB genome stratifications. Optional genome-specific
  stratifications cover complex and overlapping variants.
- **Structural variants:** [Truvari](https://github.com/ACEnglish/truvari)
  `bench`, optionally followed by `refine`. Matching parameters come from named
  profiles.

Either the draft or the comparison callset can be the truth set. Truth and
target regions are configurable per evaluation. Results are summarized along
with assembly, variant, and region statistics, and can be rendered as an
analysis report (`snakemake analysis.html`). This serves as an initial QC step
for iterating on parameters and exclusions.

## From draft to release

DeFrABB produces _draft_ benchmarks. Before a GIAB release, drafts are shared
publicly and evaluated by external groups. Putative errors are curated, and any
problems are fixed by changing parameters or exclusions and rerunning DeFrABB.
That review cycle (Fig. 1 of the preprint) happens outside the pipeline. The
released benchmark files include the analyses table, the resource configuration,
and the Snakemake report for the run that produced them.
