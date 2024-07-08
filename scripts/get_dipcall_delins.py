import argparse
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

import pybedtools
import pysam

cigar_dict = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
    9: "B",
}


def process_bam(aln_file: str) -> Dict[str, List[Tuple[int, str]]]:
    """
    Process a BAM file to identify consecutive SVs (DELINS and INSDEL).

    Args:
        aln_file (str): Path to the input BAM file.

    Returns:
        Dict[str, List[Tuple[int, str]]]: A dictionary containing consecutive SVs,
            where the keys are reference names and the values are lists of tuples
            containing the location and SV type.
    """
    aln_bam = pysam.AlignmentFile(aln_file, "rb")
    consecutive_svs = defaultdict(list)
    for r in aln_bam.fetch():
        cigar_tuples = r.cigartuples
        loc = r.reference_start
        for i, cigar in enumerate(cigar_tuples[:-1]):
            if (
                cigar_dict[cigar[0]] == "D"
                and cigar_dict[cigar_tuples[i + 1][0]] == "I"
                and cigar[1] >= 35
                and cigar_tuples[i + 1][1] >= 35
            ):
                consecutive_svs[r.reference_name].append((loc, "DELINS"))
            if (
                cigar_dict[cigar[0]] == "I"
                and cigar_dict[cigar_tuples[i + 1][0]] == "D"
                and cigar[1] >= 35
                and cigar_tuples[i + 1][1] >= 35
            ):
                consecutive_svs[r.reference_name].append((loc, "INSDEL"))
            if cigar_dict[cigar[0]] not in ["I", "S", "H"]:
                loc += cigar[1]

    return consecutive_svs


def write_intervals_to_bed(
    intervals: Dict[str, List[Tuple[int, str]]], bed_file: str
) -> None:
    """
    Write the consecutive SV intervals to a BED file.

    Args:
        intervals (Dict[str, List[Tuple[int, str]]]): A dictionary containing consecutive SVs,
            where the keys are reference names and the values are lists of tuples
            containing the location and SV type.
        bed_file (str): Path to the output BED file.
    """
    bed_lines = []
    for c in intervals:
        for loc, sv_type in intervals[c]:
            bed_lines.append(f"{c}\t{loc-100}\t{loc+100}\t{sv_type}")
    print("Num intervals:", len(bed_lines), file=sys.stderr)

    bed = pybedtools.BedTool("\n".join(bed_lines), from_string=True)
    if len(bed_lines) > 0:
        bed = bed.sort().merge(c=3, o="distinct")
    bed.saveas(bed_file)


def main():
    parser = argparse.ArgumentParser(
        description="Identify consecutive SVs from BAM files."
    )
    parser.add_argument("--hap1_bam", help="Path to the haplotype 1 BAM file")
    parser.add_argument("--hap2_bam", help="Path to the haplotype 2 BAM file")
    parser.add_argument("--output_bed", help="Path to the output BED file")
    args = parser.parse_args()

    print("Processing haplotype 1 BAM file")
    hap1_svs = process_bam(args.hap1_bam)
    print("Processing haplotype 2 BAM file")
    hap2_svs = process_bam(args.hap2_bam)

    print("Combining haplotype 1 and haplotype 2 intervals and writing to BED file")
    combined_intervals = defaultdict(list)
    for c in hap1_svs:
        combined_intervals[c].extend(hap1_svs[c])
    for c in hap2_svs:
        combined_intervals[c].extend(hap2_svs[c])

    write_intervals_to_bed(combined_intervals, args.output_bed)


if __name__ == "__main__":
    main()
