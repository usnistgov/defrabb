#!/usr/bin/env python3
"""Evaluate draft benchmark against multiple comparison callsets.

Runs hap.py evaluation in parallel against multiple comparison callsets
(v5.0q, v4.2.1, external caller outputs) and aggregates results.

Part of v0.024 parameter optimization infrastructure - RIDE discrimination analysis.
"""

import argparse
import concurrent.futures
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional


class MultiCallsetEvaluator:
    """Evaluate benchmark against multiple comparison callsets."""

    def __init__(
        self,
        benchmark_vcf: Path,
        benchmark_bed: Path,
        reference: Path,
        comparison_callsets: Dict[str, Dict],
        output_dir: Path,
        threads: int = 4,
        stratifications: Optional[Path] = None
    ):
        self.benchmark_vcf = benchmark_vcf
        self.benchmark_bed = benchmark_bed
        self.reference = reference
        self.callsets = comparison_callsets
        self.output_dir = output_dir
        self.threads = threads
        self.strat_dir = stratifications

    def run_happy_evaluation(self, callset_id: str, callset_info: Dict) -> Dict:
        """Run hap.py evaluation for one comparison callset.

        Returns dict with callset_id, status, and output_dir.
        """
        output_subdir = self.output_dir / callset_id
        output_subdir.mkdir(parents=True, exist_ok=True)

        output_prefix = output_subdir / f"{callset_id}_vs_benchmark"

        # Determine which is query vs truth based on callset type
        # For benchmarks (v5.0q, v4.2.1), we compare as: truth=comparison, query=our_benchmark
        # This makes FN = "missing from our benchmark" and FP = "extra in our benchmark"
        is_benchmark = callset_info.get("type") == "benchmark"

        if is_benchmark:
            truth_vcf = callset_info["vcf"]
            truth_bed = callset_info.get("bed")
            query_vcf = self.benchmark_vcf
            query_bed = self.benchmark_bed
        else:
            # For caller outputs, flip: truth=our_benchmark, query=caller
            # This makes FN = "caller missed" and FP = "caller spurious"
            truth_vcf = self.benchmark_vcf
            truth_bed = self.benchmark_bed
            query_vcf = callset_info["vcf"]
            query_bed = callset_info.get("bed")

        cmd = [
            "hap.py",
            str(truth_vcf),
            str(query_vcf),
            "-f", str(truth_bed) if truth_bed else str(self.benchmark_bed),
            "-r", str(self.reference),
            "-o", str(output_prefix),
            "--pass-only",
            "--engine=vcfeval",
            "--threads", str(self.threads)
        ]

        # Add stratifications if provided
        if self.strat_dir:
            cmd.extend(["--stratification", str(self.strat_dir)])

        # Add query bed if available
        if not is_benchmark and query_bed:
            cmd.extend(["-T", str(query_bed)])

        print(f"[{callset_id}] Running hap.py evaluation...")
        print(f"  Truth: {truth_vcf}")
        print(f"  Query: {query_vcf}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )

            if result.returncode != 0:
                print(f"[{callset_id}] FAILED: {result.stderr}", file=sys.stderr)
                return {
                    "callset_id": callset_id,
                    "status": "failed",
                    "error": result.stderr,
                    "output_dir": str(output_subdir)
                }

            print(f"[{callset_id}] Completed successfully")
            return {
                "callset_id": callset_id,
                "status": "success",
                "output_dir": str(output_subdir),
                "summary_csv": str(output_prefix) + ".summary.csv",
                "extended_csv": str(output_prefix) + ".extended.csv"
            }

        except subprocess.TimeoutExpired:
            print(f"[{callset_id}] TIMEOUT after 2 hours", file=sys.stderr)
            return {
                "callset_id": callset_id,
                "status": "timeout",
                "output_dir": str(output_subdir)
            }
        except Exception as e:
            print(f"[{callset_id}] ERROR: {e}", file=sys.stderr)
            return {
                "callset_id": callset_id,
                "status": "error",
                "error": str(e),
                "output_dir": str(output_subdir)
            }

    def run_all_evaluations(self) -> List[Dict]:
        """Run evaluations in parallel across all comparison callsets."""
        print(f"Evaluating benchmark against {len(self.callsets)} comparison callsets...")
        print(f"  Benchmark VCF: {self.benchmark_vcf}")
        print(f"  Benchmark BED: {self.benchmark_bed}")
        print(f"  Output: {self.output_dir}")
        print()

        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, len(self.callsets))) as executor:
            futures = {
                executor.submit(self.run_happy_evaluation, cid, cinfo): cid
                for cid, cinfo in self.callsets.items()
            }

            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)

        return results

    def write_summary(self, results: List[Dict]):
        """Write summary JSON of all evaluations."""
        summary_path = self.output_dir / "multi_callset_summary.json"

        summary = {
            "benchmark_vcf": str(self.benchmark_vcf),
            "benchmark_bed": str(self.benchmark_bed),
            "reference": str(self.reference),
            "num_callsets": len(self.callsets),
            "results": results,
            "success_count": sum(1 for r in results if r["status"] == "success"),
            "failed_count": sum(1 for r in results if r["status"] != "success")
        }

        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)

        print(f"\nSummary written to: {summary_path}")
        print(f"  Successful: {summary['success_count']}/{len(results)}")
        print(f"  Failed: {summary['failed_count']}/{len(results)}")

    def run(self):
        """Execute multi-callset evaluation."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        results = self.run_all_evaluations()
        self.write_summary(results)
        return results


def load_callsets_from_resources(
    resources_yml: Path,
    sample_id: str,
    ref_id: str,
    callset_filter: Optional[List[str]] = None
) -> Dict[str, Dict]:
    """Load comparison callsets from resources.yml.

    Args:
        resources_yml: Path to config/resources.yml
        sample_id: Sample ID (e.g., HG002)
        ref_id: Reference ID (e.g., GRCh38)
        callset_filter: Optional list of callset IDs to include (default: all)

    Returns:
        Dict mapping callset_id to callset info (vcf, bed, type, etc.)
    """
    import yaml

    with open(resources_yml) as f:
        resources = yaml.safe_load(f)

    # Extract comparison callsets for this sample + reference
    comparison_config = resources.get("comparison_callsets", {})
    sample_callsets = comparison_config.get(sample_id, {})

    callsets = {}
    for callset_id, callset_data in sample_callsets.items():
        # Filter by reference if callset has ref-specific URLs
        if ref_id not in callset_data:
            continue

        ref_data = callset_data[ref_id]

        # Apply filter if provided
        if callset_filter and callset_id not in callset_filter:
            continue

        callsets[callset_id] = {
            "vcf": Path(ref_data["vcf"]),
            "bed": Path(ref_data.get("bed")) if ref_data.get("bed") else None,
            "type": callset_data.get("type", "benchmark"),
            "version": callset_data.get("version"),
            "technology": callset_data.get("technology"),
            "caller": callset_data.get("caller")
        }

    return callsets


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--benchmark-vcf",
        type=Path,
        required=True,
        help="Draft benchmark VCF.gz to evaluate"
    )
    parser.add_argument(
        "--benchmark-bed",
        type=Path,
        required=True,
        help="Draft benchmark BED (high-confidence regions)"
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Reference genome FASTA"
    )
    parser.add_argument(
        "--callsets-json",
        type=Path,
        help="JSON file with comparison callsets (alternative to --resources-yml)"
    )
    parser.add_argument(
        "--resources-yml",
        type=Path,
        help="DeFrABB resources.yml (auto-load comparison_callsets)"
    )
    parser.add_argument(
        "--sample",
        help="Sample ID for resources.yml lookup (e.g., HG002)"
    )
    parser.add_argument(
        "--ref-id",
        help="Reference ID for resources.yml lookup (e.g., GRCh38)"
    )
    parser.add_argument(
        "--callset-filter",
        nargs="+",
        help="Only evaluate these callset IDs (default: all)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for all evaluations"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads per hap.py run (default: 4)"
    )
    parser.add_argument(
        "--stratifications",
        type=Path,
        help="GIAB stratifications directory (optional)"
    )

    args = parser.parse_args()

    # Load comparison callsets
    if args.callsets_json:
        with open(args.callsets_json) as f:
            callsets = json.load(f)
    elif args.resources_yml:
        if not args.sample or not args.ref_id:
            parser.error("--sample and --ref-id required with --resources-yml")
        callsets = load_callsets_from_resources(
            args.resources_yml,
            args.sample,
            args.ref_id,
            args.callset_filter
        )
    else:
        parser.error("Either --callsets-json or --resources-yml required")

    if not callsets:
        print("ERROR: No comparison callsets found", file=sys.stderr)
        return 1

    evaluator = MultiCallsetEvaluator(
        benchmark_vcf=args.benchmark_vcf,
        benchmark_bed=args.benchmark_bed,
        reference=args.reference,
        comparison_callsets=callsets,
        output_dir=args.output_dir,
        threads=args.threads,
        stratifications=args.stratifications
    )

    results = evaluator.run()

    # Exit code: 0 if all succeeded, 1 if any failed
    failed = sum(1 for r in results if r["status"] != "success")
    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
