"""Diploidize haploid (single-token) genotypes in a VCF stream (#59 / #173).

`vcfgeno2haplo` (vcflib 1.0.15) segfaults on haploid genotypes -- the
hemizygous chrX / chrY calls dipcall emits for male samples crash it with
`basic_string: construction from null`. This stdin->stdout filter rewrites any
single-token GT (no `/` or `|` separator, e.g. `1`, `2`, `.`) to the equivalent
homozygous diploid form (`1` -> `1|1`, `.` -> `.|.`) so the call no longer
crashes.

The transform is a no-op on already-diploid records (autosomes and the
pseudoautosomal regions keep their separated GT untouched) and is
classification-neutral for the downstream genome-specific strats: hemizygous and
homozygous are equivalent, and a compound het requires GT `1/2`, which a haploid
region can never carry. See
wiki:investigations/genome-specific-geno2haplo-haploid-segfault.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys


def diploidize_gt_field(gt: str) -> str:
    """Return a diploid GT: single-token GTs become homozygous, others pass through.

    The first ``:``-delimited token of a sample column is the GT. A GT with no
    ``/`` or ``|`` separator is haploid; duplicate it as ``<allele>|<allele>``.
    """
    if "/" in gt or "|" in gt:
        return gt
    return f"{gt}|{gt}"


def diploidize_line(line: str) -> str:
    """Diploidize the GT of a VCF data line; header/blank lines pass through."""
    if line.startswith("#") or not line.strip():
        return line
    fields = line.rstrip("\n").split("\t")
    if len(fields) <= 9:  # no sample column -> nothing to do
        return line
    sample = fields[9].split(":")
    sample[0] = diploidize_gt_field(sample[0])
    fields[9] = ":".join(sample)
    return "\t".join(fields) + "\n"


def main() -> None:
    for line in sys.stdin:
        sys.stdout.write(diploidize_line(line))


if __name__ == "__main__":
    main()
