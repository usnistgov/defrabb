#!/usr/bin/env python3
"""Extract discrepancies (FN/FP variants) from hap.py evaluation outputs.

Reads hap.py output VCF (with BD:TP/FN/FP annotations) and extracts all false
negative and false positive variants with genomic context for manual curation.

Outputs TSV compatible with GIAB curation table format.

Part of v0.024 parameter optimization infrastructure.
"""

import argparse
import csv
import gzip
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set


class DiscrepancyExtractor:
    """Extract FN/FP variants from hap.py outputs with genomic context."""

    def __init__(
        self,
        evaluation_dir: Path,
        stratifications_dir: Optional[Path] = None,
        output: Optional[Path] = None
    ):
        self.eval_dir = evaluation_dir
        self.strat_dir = stratifications_dir
        self.output: Path = output if output else (evaluation_dir / "discrepancies.tsv")

        # Auto-detect evaluation VCF
        self.vcf_path = self._find_vcf()

    def _find_vcf(self) -> Path:
        """Find the hap.py output VCF in evaluation directory."""
        candidates = list(self.eval_dir.glob("*.vcf.gz"))
        if not candidates:
            raise FileNotFoundError(
                f"No *.vcf.gz found in {self.eval_dir}"
            )
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple *.vcf.gz found in {self.eval_dir}: {candidates}"
            )
        return candidates[0]

    def load_stratifications(self) -> Dict[str, Set[tuple]]:
        """Load stratification BED files into position sets.

        Returns dict mapping stratification name to set of (chrom, start, end) tuples.
        """
        if not self.strat_dir or not self.strat_dir.exists():
            return {}

        strats = {}
        bed_files = list(self.strat_dir.glob("*.bed.gz"))
        bed_files += list(self.strat_dir.glob("*.bed"))

        for bed_file in bed_files:
            strat_name = bed_file.stem.replace(".bed", "")
            positions = set()

            opener = gzip.open if bed_file.suffix == ".gz" else open
            with opener(bed_file, "rt") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("track"):
                        continue
                    fields = line.strip().split("\t")
                    if len(fields) < 3:
                        continue
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    positions.add((chrom, start, end))

            strats[strat_name] = positions

        return strats

    def annotate_stratifications(
        self, chrom: str, pos: int, strats: Dict[str, Set[tuple]]
    ) -> List[str]:
        """Return list of stratification names overlapping this position."""
        overlapping = []
        for strat_name, positions in strats.items():
            for s_chrom, s_start, s_end in positions:
                if s_chrom == chrom and s_start <= pos < s_end:
                    overlapping.append(strat_name)
                    break
        return overlapping

    def extract_from_vcf(self) -> List[Dict]:
        """Extract FN/FP variants from hap.py output VCF using bcftools.

        Hap.py annotates variants with BD (benchmark decision) field:
        - BD:TP = true positive
        - BD:FN = false negative
        - BD:FP = false positive

        Returns list of discrepancy records compatible with GIAB curation table format.
        """
        discrepancies = []
        strats = self.load_stratifications()

        # Use bcftools to avoid pysam FIPS issues
        cmd = [
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%BD\t]\n",
            str(self.vcf_path)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"bcftools query failed: {result.stderr}")

        for line in result.stdout.strip().split("\n"):
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2]
            alt = fields[3]
            qual = fields[4] if fields[4] != "." else "."
            filter_str = fields[5] if fields[5] != "." else "PASS"

            # Hap.py has two samples: TRUTH and QUERY
            # Format: GT_truth, BD_truth, GT_query, BD_query
            gt_truth = fields[6] if len(fields) > 6 else "."
            bd_truth = fields[7] if len(fields) > 7 else "."
            gt_query = fields[8] if len(fields) > 8 else "."
            bd_query = fields[9] if len(fields) > 9 else "."

            # Determine FN vs FP vs TP
            if bd_truth == "FN" or bd_query == "FN":
                disc_type = "fn"
            elif bd_truth == "FP" or bd_query == "FP":
                disc_type = "fp"
            else:
                # TP - skip
                continue

            # Determine variant type
            alt_alleles = alt.split(",")
            if len(ref) == len(alt_alleles[0]) == 1:
                var_type = "SNP"
            else:
                var_type = "INDEL"

            # Annotate with stratifications
            strat_list = self.annotate_stratifications(chrom, pos, strats)

            discrepancies.append({
                "chrom": chrom,
                "chromStart": pos - 1,  # 0-based for GIAB format
                "chromEnd": pos,  # 1-based
                "var_type": var_type,
                "label": disc_type,
                "VCF_REF": ref,
                "VCF_ALT": alt,
                "VCF_QUAL": qual,
                "VCF_FILTER": filter_str,
                "VCF_GT": gt_query,
                "stratifications": ",".join(strat_list) if strat_list else "none",
                "gt_truth": gt_truth,
                "gt_query": gt_query,
                "bd_truth": bd_truth,
                "bd_query": bd_query,
            })

        return discrepancies

    def write_tsv(self, discrepancies: List[Dict]):
        """Write discrepancies to TSV file (GIAB curation table format)."""
        if not discrepancies:
            print(f"No discrepancies found in {self.vcf_path}", file=sys.stderr)
            return

        # GIAB curation table core fields + extras
        fieldnames = [
            "chrom", "chromStart", "chromEnd", "var_type", "label",
            "VCF_REF", "VCF_ALT", "VCF_QUAL", "VCF_FILTER", "VCF_GT",
            "stratifications", "gt_truth", "gt_query", "bd_truth", "bd_query"
        ]

        self.output.parent.mkdir(parents=True, exist_ok=True)
        with open(self.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(discrepancies)

        # Summary stats
        fn_count = sum(1 for d in discrepancies if d["label"] == "fn")
        fp_count = sum(1 for d in discrepancies if d["label"] == "fp")
        snp_count = sum(1 for d in discrepancies if d["var_type"] == "SNP")
        indel_count = sum(1 for d in discrepancies if d["var_type"] == "INDEL")

        print(f"Extracted {len(discrepancies)} discrepancies:")
        print(f"  FN (missing from benchmark): {fn_count}")
        print(f"  FP (extra in benchmark): {fp_count}")
        print(f"  SNPs: {snp_count}, INDELs: {indel_count}")
        print(f"Written to: {self.output}")

    def run(self):
        """Execute discrepancy extraction."""
        print(f"Extracting discrepancies from: {self.eval_dir}")
        print(f"  VCF: {self.vcf_path}")
        if self.strat_dir:
            print(f"  Stratifications: {self.strat_dir}")

        discrepancies = self.extract_from_vcf()
        self.write_tsv(discrepancies)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--evaluation-dir",
        type=Path,
        required=True,
        help="Hap.py evaluation directory containing *.vcf.gz output"
    )
    parser.add_argument(
        "--stratifications",
        type=Path,
        help="Directory containing stratification BED/BED.gz files (GIAB format)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV file (default: <evaluation-dir>/discrepancies.tsv)"
    )

    args = parser.parse_args()

    extractor = DiscrepancyExtractor(
        evaluation_dir=args.evaluation_dir,
        stratifications_dir=args.stratifications,
        output=args.output
    )

    extractor.run()


if __name__ == "__main__":
    main()
