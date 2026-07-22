#!/usr/bin/env python3
"""Score parameter sweep results against v5q baseline.

Wrapper around compare_evaluations.py optimized for parameter optimization
workflow: identifies baseline, ranks by metric, reports top-N param sets.

Part of v0.023 parameter optimization framework.
"""

import argparse
import subprocess
import sys
from pathlib import Path
import csv


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Pipeline results directory containing evaluations/"
    )
    parser.add_argument(
        "--baseline",
        default="v5.0q",
        help="Baseline analysis_id substring (default: v5.0q)"
    )
    parser.add_argument(
        "--rank-by",
        choices=["f1", "precision", "recall"],
        default="f1",
        help="Metric to rank by (default: f1)"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0001,
        help="Min delta to flag regression/improvement (default: 0.0001)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=3,
        help="Report top-N parameter sets (default: 3)"
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output markdown summary (default: stdout)"
    )

    args = parser.parse_args()

    # Generate comparison TSV
    tsv_path = args.results_dir / "param_sweep_comparison.tsv"

    cmd = [
        "scripts/compare_evaluations.py",
        "--results-dir", str(args.results_dir),
        "--baseline", args.baseline,
        "--metric", args.rank_by,
        "--threshold", str(args.threshold),
        "--out", str(tsv_path),
        "--regions",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: compare_evaluations.py failed:\n{result.stderr}", file=sys.stderr)
        return 1

    # Load TSV and analyze
    records = []
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            records.append(row)

    # Filter to non-baseline, non-regressed
    candidates = [
        r for r in records
        if r.get("status") in ["improved", "same"]
        and r.get("variant_type") == "SNP"  # Rank on SNP performance
    ]

    # Sort by metric (descending)
    metric_key = args.rank_by
    candidates.sort(key=lambda r: float(r.get(metric_key, 0)), reverse=True)

    # Extract top-N unique analysis_ids
    seen_ids = set()
    top_analyses = []
    for r in candidates:
        aid = r["analysis_id"]
        if aid not in seen_ids:
            seen_ids.add(aid)
            top_analyses.append(r)
        if len(top_analyses) >= args.top_n:
            break

    # Generate summary report
    lines = [
        f"# Parameter Sweep Scoring Summary",
        f"",
        f"**Results directory:** {args.results_dir}",
        f"**Baseline:** {args.baseline}",
        f"**Ranked by:** {args.rank_by}",
        f"**Threshold:** {args.threshold}",
        f"",
        f"## Top-{args.top_n} Parameter Sets (SNP {args.rank_by})",
        f"",
    ]

    for i, r in enumerate(top_analyses, 1):
        delta_key = f"delta_{args.rank_by}"
        delta = float(r.get(delta_key, 0))
        delta_str = f"{delta:+.6f}" if delta != 0 else "baseline"

        lines.append(f"### {i}. {r['analysis_id']}")
        lines.append(f"")
        lines.append(f"- **{args.rank_by.capitalize()}:** {r[metric_key]} ({delta_str})")
        lines.append(f"- **Precision:** {r['precision']} ({r.get('delta_precision', 'N/A')})")
        lines.append(f"- **Recall:** {r['recall']} ({r.get('delta_recall', 'N/A')})")
        lines.append(f"- **F1:** {r['f1']} ({r.get('delta_f1', 'N/A')})")
        lines.append(f"- **Status:** {r['status']}")
        lines.append(f"")

    # Count regressions
    regressions = [r for r in records if r.get("status") == "regressed"]
    if regressions:
        lines.append(f"## ⚠️ Regressions Detected")
        lines.append(f"")
        lines.append(f"**{len(regressions)} analysis/variant combinations regressed:**")
        lines.append(f"")
        for r in regressions[:10]:  # Limit to first 10
            delta_key = f"delta_{args.rank_by}"
            lines.append(
                f"- {r['analysis_id']} ({r['variant_type']}): "
                f"{args.rank_by}={r[metric_key]} ({r.get(delta_key, 'N/A')})"
            )
        if len(regressions) > 10:
            lines.append(f"- ... ({len(regressions) - 10} more)")
        lines.append(f"")

    # Recommendations
    lines.append(f"## Recommendations")
    lines.append(f"")
    if top_analyses:
        top_id = top_analyses[0]["analysis_id"]
        lines.append(f"**Recommended parameter set:** {top_id}")
        lines.append(f"")
        lines.append(f"For validation across multiple genomes, test these top-{args.top_n}:")
        for r in top_analyses:
            lines.append(f"- {r['analysis_id']}")
    else:
        lines.append(f"**No improved or same-performing parameter sets found.**")
        lines.append(f"All sweeps regressed against baseline.")
    lines.append(f"")

    lines.append(f"---")
    lines.append(f"*Generated by scripts/score_param_sweep.py*")
    lines.append(f"*Full comparison: {tsv_path}*")

    report = "\n".join(lines)

    # Write or print
    if args.out:
        args.out.write_text(report)
        print(f"Summary written to: {args.out}")
    else:
        print(report)

    return 0


if __name__ == "__main__":
    sys.exit(main())
