"""Generate genome-specific stratification BEDs from a benchmark VCF (#59/#173).

This ports the complex-/overlapping-variant stratifications from the GIAB
genome-stratifications "GenomeSpecific" Snakemake pipeline
(genome-in-a-bottle/genome-stratifications, GRCh38/GenomeSpecific/
GS-snakemake-pipeline) into a single, dependency-light, unit-testable module.

The upstream pipeline classifies variants in a benchmark VCF that has first been
processed with ``vcfgeno2haplo -w 10`` (vcflib) to merge variants within 10 bp
into combined haplotype records. From that it derives five strata:

  comphetsnp10bp      compound hets where both alleles are length-preserving
                      (SNP/MNP only) within 10 bp
  comphetindel10bp    compound hets where at least one allele changes length
  complexindel10bp    non-comphet complex indels (REF>1 and ALT>1, len differ)
  snpswithin10bp      non-comphet MNP-like (REF==ALT length, >1 bp)
  othercomplexwithin  clusters of >1 variant within 10 bp not captured above

Each variant contributes the interval ``[POS-50, POS+len(REF)+50)`` (mirroring
the upstream ±50 bp slop, applied to the 1-based POS), and intervals are merged
within 50 bp. ``othercomplexwithin`` is computed from the raw (pre-geno2haplo)
benchmark VCF: clusters of more than one record within 10 bp, slopped ±50 and
merged, then with the other four strata subtracted.

Interval merge/subtract are implemented in-module so the whole classification is
deterministic and testable without bedtools; the surrounding Snakemake rules
handle the genome-wide bed algebra (intermediate bed, SV/CNV bed, unions) where
bedtools is the natural tool.

This module was developed with assistance from Claude (Anthropic). All code has
been reviewed and tested by the primary author.
"""

import argparse
import gzip
import os
import sys
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

SLOP = 50
Interval = Tuple[int, int]  # (start, end), 0-based half-open after construction

# Genotype allele pairs that mark a compound heterozygote (two distinct ALTs).
_COMPHET_ALLELES = frozenset({"1", "2"})


def _open_vcf(path: str):
    """Open a plain or gzipped VCF for text reading."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _gt_alleles(sample_field: str) -> frozenset:
    """Return the set of allele indices in the first sample's GT."""
    gt = sample_field.split(":", 1)[0]
    return frozenset(a for a in gt.replace("|", "/").split("/") if a not in ("", "."))


def iter_records(path: str) -> Iterator[Tuple[str, int, str, List[str], frozenset]]:
    """Yield (chrom, pos, ref, alts, gt_alleles) for each VCF data record."""
    with _open_vcf(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            alts = alt.split(",")
            gt = _gt_alleles(f[9]) if len(f) > 9 else frozenset()
            yield chrom, pos, ref, alts, gt


def _is_comphet(gt: frozenset, alts: Sequence[str]) -> bool:
    return len(alts) >= 2 and gt == _COMPHET_ALLELES


def _slop_interval(pos: int, ref: str) -> Interval:
    """Interval [POS-50, POS+len(REF)+50), clamped at 0 (mirrors upstream awk)."""
    return (max(0, pos - SLOP), pos + len(ref) + SLOP)


def classify_geno2haplo(
    records: Iterable[Tuple[str, int, str, List[str], frozenset]],
) -> Dict[str, Dict[str, List[Interval]]]:
    """Bucket geno2haplo records into the four primary strata (rules 2-5).

    Returns ``{stratum: {chrom: [interval, ...]}}`` with unmerged intervals.
    """
    out: Dict[str, Dict[str, List[Interval]]] = {
        "comphetsnp10bp": {},
        "comphetindel10bp": {},
        "complexindel10bp": {},
        "snpswithin10bp": {},
    }

    def add(stratum: str, chrom: str, iv: Interval) -> None:
        out[stratum].setdefault(chrom, []).append(iv)

    for chrom, pos, ref, alts, gt in records:
        iv = _slop_interval(pos, ref)
        if _is_comphet(gt, alts):
            a1, a2 = alts[0], alts[1]
            same_len = len(ref) == len(a1) and len(ref) == len(a2)
            if same_len:
                add("comphetsnp10bp", chrom, iv)  # rule 2
            else:
                add("comphetindel10bp", chrom, iv)  # rule 3
        else:
            alt = alts[0]
            if len(ref) != len(alt) and len(ref) > 1 and len(alt) > 1:
                add("complexindel10bp", chrom, iv)  # rule 4
            elif len(ref) == len(alt) and len(ref) > 1:
                add("snpswithin10bp", chrom, iv)  # rule 5
    return out


def merge_intervals(intervals: Sequence[Interval], dist: int = SLOP) -> List[Interval]:
    """Merge intervals that are within ``dist`` of each other (bedtools -d)."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + dist:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def merge_with_count(
    intervals: Sequence[Interval], dist: int
) -> List[Tuple[int, int, int]]:
    """Merge within ``dist`` and count how many source intervals each cluster has."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    s, e = ordered[0]
    clusters = [[s, e, 1]]
    for start, end in ordered[1:]:
        if start <= clusters[-1][1] + dist:
            clusters[-1][1] = max(clusters[-1][1], end)
            clusters[-1][2] += 1
        else:
            clusters.append([start, end, 1])
    return [(s, e, c) for s, e, c in clusters]


def subtract_intervals(a: Sequence[Interval], b: Sequence[Interval]) -> List[Interval]:
    """Return ``a`` with all regions overlapping ``b`` removed (bedtools subtract)."""
    if not a:
        return []
    b_merged = merge_intervals(b, dist=0)
    result: List[Interval] = []
    bi = 0
    for start, end in sorted(a):
        cur = start
        while bi > 0 and b_merged[bi - 1][1] > cur:
            bi -= 1
        while bi < len(b_merged) and b_merged[bi][1] <= cur:
            bi += 1
        j = bi
        while j < len(b_merged) and b_merged[j][0] < end:
            bs, be = b_merged[j]
            if bs > cur:
                result.append((cur, min(bs, end)))
            cur = max(cur, be)
            if cur >= end:
                break
            j += 1
        if cur < end:
            result.append((cur, end))
    return result


def othercomplex_intervals(
    raw_records: Iterable[Tuple[str, int, str, List[str], frozenset]],
    primary: Dict[str, Dict[str, List[Interval]]],
) -> Dict[str, List[Interval]]:
    """Rule 6: clusters of >1 raw variant within 10 bp, minus the four strata."""
    by_chrom: Dict[str, List[Interval]] = {}
    for chrom, pos, ref, _alts, _gt in raw_records:
        by_chrom.setdefault(chrom, []).append((pos, pos + len(ref)))

    out: Dict[str, List[Interval]] = {}
    for chrom, ivs in by_chrom.items():
        clusters = merge_with_count(ivs, dist=10)
        slopped = [
            (max(0, s - SLOP), e + SLOP) for s, e, count in clusters if count > 1
        ]
        merged = merge_intervals(slopped, dist=0)
        for stratum in (
            "comphetindel10bp",
            "comphetsnp10bp",
            "complexindel10bp",
            "snpswithin10bp",
        ):
            merged = subtract_intervals(
                merged, merge_intervals(primary[stratum].get(chrom, []))
            )
        if merged:
            out[chrom] = merged
    return out


def _write_bed(path: str, by_chrom: Dict[str, List[Interval]]) -> int:
    """Write a sorted BED; return the number of intervals written."""
    n = 0
    with open(path, "w") as fh:
        for chrom in sorted(by_chrom):
            for start, end in by_chrom[chrom]:
                fh.write(f"{chrom}\t{start}\t{end}\n")
                n += 1
    return n


STRATA_ORDER = [
    "comphetsnp10bp",
    "comphetindel10bp",
    "complexindel10bp",
    "snpswithin10bp",
    "othercomplexwithin10bp",
]


def build_complex_strats(
    geno2haplo_vcf: str, raw_vcf: str, outdir: str, prefix: str
) -> Dict[str, str]:
    """Write the five complex-variant strat BEDs; return {stratum: path}."""
    primary = classify_geno2haplo(iter_records(geno2haplo_vcf))
    other = othercomplex_intervals(iter_records(raw_vcf), primary)

    merged = {
        s: {c: merge_intervals(ivs) for c, ivs in d.items()} for s, d in primary.items()
    }
    merged["othercomplexwithin10bp"] = other

    os.makedirs(outdir, exist_ok=True)
    paths: Dict[str, str] = {}
    for stratum in STRATA_ORDER:
        path = os.path.join(outdir, f"{prefix}{stratum}_slop50.bed")
        _write_bed(path, merged.get(stratum, {}))
        paths[stratum] = path
    return paths


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geno2haplo", required=True, help="vcfgeno2haplo -w 10 VCF")
    p.add_argument("--raw", required=True, help="raw benchmark VCF (pre-geno2haplo)")
    p.add_argument("--outdir", required=True, help="output directory for strat BEDs")
    p.add_argument("--prefix", default="", help="filename prefix for outputs")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    paths = build_complex_strats(args.geno2haplo, args.raw, args.outdir, args.prefix)
    for stratum, path in paths.items():
        print(f"{stratum}\t{path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
