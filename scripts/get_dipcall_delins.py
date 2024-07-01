import sys
from collections import defaultdict

import pybedtools
import pysam

cigar_dict={
    0:"M",
    1:"I",
    2:"D",
    3:"N",
    4:"S",
    5:"H",
    6:"P",
    7:"=",
    8:"X",
    9:"B"
}
 
def process_bam(aln_file): 
    aln_bam= pysam.AlignmentFile(aln_file, "rb")
    consecutive_insdels = defaultdict(list)
    consecutive_delins = defaultdict(list)
    
    for r in aln_bam.fetch():
        cigar_tuples = r.cigartuples
        loc = r.reference_start
        for i, cigar in enumerate(cigar_tuples[:-1]):
            if cigar_dict[cigar[0]] == "D" and cigar_dict[cigar_tuples[i+1][0]] == "I" and cigar[1] >= 35 and cigar_tuples[i+1][1] >= 35:
                consecutive_delins[r.reference_name].append(loc)
            if cigar_dict[cigar[0]] == "I" and cigar_dict[cigar_tuples[i+1][0]] == "D" and cigar[1] >= 35 and cigar_tuples[i+1][1] >= 35:
                consecutive_insdels[r.reference_name].append(loc)
            if cigar_dict[cigar[0]] != "I" and cigar_dict[cigar[0]] != "S" and cigar_dict[cigar[0]] != "H":
                loc += cigar[1]
    print(f"{aln_file} insdels: {len(consecutive_insdels)}  delins: {len(consecutive_delins)}")
    return consecutive_insdels, consecutive_delins

def write_intervals_to_bed(intervals, bed_file):
    bed_lines = []
    for c in intervals:
        for loc in intervals[c]:
            bed_lines.append(f"{c}\t{loc-100}\t{loc+100}")

    bed = pybedtools.BedTool("\n".join(bed_lines), from_string=True)
    bed = bed.sort().merge()
    bed.saveas(bed_file)

hap1_bam = sys.argv[1]
hap2_bam = sys.argv[2]
output_bed = sys.argv[3]

hap1_insdels, hap1_delins = process_bam(hap1_bam)
hap2_insdels, hap2_delins = process_bam(hap2_bam)

combined_intervals = defaultdict(list)

svs_count = 0
for c in hap1_insdels:
    svs_count += len(hap1_insdels[c])
    combined_intervals[c].extend(hap1_insdels[c])
print(f"hap1 insdels: {svs_count}; total: {len(combined_intervals)}")
svs_count = 0
for c in hap2_insdels:
    svs_count += len(hap2_insdels[c])
    combined_intervals[c].extend(hap2_insdels[c])
print(f"hap2 insdels: {svs_count}; total: {len(combined_intervals)}")
svs_count = 0
for c in hap1_delins:
    svs_count += len(hap1_delins[c])
    combined_intervals[c].extend(hap1_delins[c])
print(f"hap1 delins: {svs_count}; total: {len(combined_intervals)}")
svs_count = 0
for c in hap2_delins:
    svs_count += len(hap2_delins[c])
    combined_intervals[c].extend(hap2_delins[c])
print(f"hap2 delins: {svs_count}; total: {len(combined_intervals)}")

write_intervals_to_bed(combined_intervals, output_bed)