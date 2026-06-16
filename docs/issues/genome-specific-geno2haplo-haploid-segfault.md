# `genome_specific_geno2haplo` segfault on haploid genotypes (#59 / #173)

_Diagnosed 2026-06-16 from the 20260615_v0.022 full-pipeline regression run._

## Symptom

`rule genome_specific_geno2haplo` failed with exit 139 (SIGSEGV) on the
whole-genome GRCh38 HG002 smvar draft benchmark:

```
no alternate alleles remain at chrY:11328327 after haplotype validation
...
Segmentation fault (core dumped) vcfgeno2haplo -w 10 -r resources/references/GRCh38.fa "$tmp"
```

Run on a smaller window the crash surfaces its real cause:

```
terminate called after throwing an instance of 'std::logic_error'
  what():  basic_string: construction from null is not valid
```

## Root cause

`vcfgeno2haplo` (vcflib 1.0.15) mishandles **haploid genotypes**. dipcall emits
hemizygous chrX / chrY calls for male samples as single-token GTs (`GT=1`, and
multiallelic `2`/`3`/`4`); vcfgeno2haplo assumes a diploid GT, indexes a second
allele that does not exist, and constructs a `std::string` from the resulting
null pointer -> crash. The feature (#59/#173) was validated only on HG002
**chr21** (all-diploid, 2653 records), so no haploid GT was ever exercised until
the whole-genome run hit chrX/chrY.

Confirmed by bisection: the crash is cumulative window state on chrY haploid
records; diploidizing the GT (`1` -> `1|1`) makes it disappear (exit 0).

## Fix

`scripts/diploidize_gt.py` — a stdin->stdout VCF filter — rewrites single-token
GTs to homozygous diploid (`1` -> `1|1`, `.` -> `.|.`) before `vcfgeno2haplo`.
The `genome_specific_geno2haplo` rule pipes `zcat {vcf}` through it.

Why this is safe / classification-neutral:

- **No-op on diploid records.** Autosomes and the pseudoautosomal regions keep
  their `/`- or `|`-separated GT untouched (verified: the filter produces a
  byte-identical VCF on chr21).
- **Hemizygous == homozygous.** Representing a single-haplotype call as `1|1`
  is the standard diploid encoding of a hemizygous call.
- **Cannot fabricate a compound het.** `scripts/genome_specific_strats.py` flags
  a compound het only when the GT alleles are exactly `{1, 2}` (i.e. `1/2`),
  which a haploid region can never carry. Diploidizing `1` -> `1|1` yields
  allele set `{1}`, never a comphet.

## Validation

- Reproduced and fixed against the prior-run input
  `20260323_v0.020_HG002v1.1-CI/.../GRCh38_HG2-T2TQ100-V1.1_smvar_dipcall-z2k.vcf.gz`
  (same asm/ref/caller), vcflib 1.0.15 from `envs/genome_strats.yml`.
- Whole-genome `vcfgeno2haplo` now completes (exit 0), all 24 contigs present in
  output including chrX (140,500) and chrY (13,282) haplotype records.
- `scripts/diploidize_gt.py` is unit-tested (`.tests/unit/test_diploidize_gt.py`).

## Note (separate, harmless)

`vcfgeno2haplo` writes an uninitialized QUAL value (garbage float vs `0`) that
varies between runs. The genome-specific classifier ignores QUAL, so this
pre-existing non-determinism does not affect stratification output.
