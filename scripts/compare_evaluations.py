#!/usr/bin/env python
"""Compare DeFrABB evaluation results across analyses.

Aggregates hap.py (`*.summary.csv`) and Truvari (`summary.json`) outputs under a
results directory into one tidy table of precision/recall/F1 by variant type,
and (optionally) computes deltas versus a baseline analysis to highlight
regressions and improvements. Writes a machine-readable TSV and prints a
human-readable Markdown table.

Plotting is handled separately by `scripts/plot_evaluations.R`, which consumes
the TSV written by this script (pass --out, then run the R script on it; or use
--plot to invoke it automatically when Rscript is available).

This module was developed with assistance from Claude (Anthropic). All code has
been reviewed and tested by the primary author.
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

# hap.py summary.csv: keep the genome-wide PASS rows for the two core types.
HAPPY_TYPES = ("SNP", "INDEL")
HAPPY_FILTER = "PASS"

# Metrics reported in the comparison table, in display order.
METRICS = ("precision", "recall", "f1")

_BENCH_TYPE = re.compile(r"_(smvar|stvar)_")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group()
    src.add_argument(
        "--results-dir",
        default="results",
        help="Pipeline results dir to search for evaluation outputs "
        "(default: results)",
    )
    src.add_argument(
        "--inputs",
        nargs="+",
        help="Explicit evaluation files (hap.py *.summary.csv and/or Truvari "
        "summary.json) instead of discovering them under --results-dir",
    )
    p.add_argument(
        "--baseline",
        help="analysis_id (or unique substring) to use as the baseline; deltas "
        "and regression/improvement flags are computed against it",
    )
    p.add_argument(
        "--metric",
        choices=METRICS,
        default="f1",
        help="Metric used for sorting and regression flagging (default: f1)",
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Minimum absolute delta in --metric to flag a "
        "regression/improvement (default: 0.0)",
    )
    p.add_argument("--out", help="Write the tidy comparison table to this TSV path")
    p.add_argument(
        "--plot",
        metavar="DIR",
        help="After writing the TSV, render plots into DIR via "
        "scripts/plot_evaluations.R (requires --out and Rscript)",
    )
    p.add_argument(
        "--regions",
        action="store_true",
        help="Include benchmark region (BED) statistics in output",
    )
    return p.parse_args(argv)


def discover_eval_files(results_dir: str) -> List[Path]:
    """Find hap.py summary.csv and Truvari summary JSONs under a results dir.

    Truvari outputs live two levels deep: ``truvari/{group}/{callset}/summary.json``
    (and ``refine.variant_summary.json`` when refine ran).
    """
    root = Path(results_dir)
    happy = sorted(root.glob("evaluations/happy/*/*.summary.csv"))
    truvari = sorted(root.glob("evaluations/truvari/*/*/summary.json"))
    refine = sorted(root.glob("evaluations/truvari/*/*/refine.variant_summary.json"))
    return happy + truvari + refine


def discover_benchmark_beds(results_dir: str) -> List[Path]:
    """Find benchmark BED files under draft_benchmarksets."""
    root = Path(results_dir)
    return sorted(root.glob("draft_benchmarksets/*/*.benchmark.bed"))


def parse_eval_path(path: Path) -> Dict[str, str]:
    """Derive a stable analysis identity from an evaluation output path.

    Paths look like
      results/evaluations/happy/{eval_id}_{bench_id}/{callset}.summary.csv
      results/evaluations/truvari/{eval_id}_{bench_id}/{callset}/summary.json
    where both {eval_id}_{bench_id} and {callset} may contain underscores, so we
    keep them as opaque labels rather than splitting fragilely.
    """
    parts = path.parts
    eval_tool = (
        "happy" if "happy" in parts else "truvari" if "truvari" in parts else "?"
    )
    if path.suffix == ".json":
        # Truvari: <group>/<callset>/{summary,refine.variant_summary}.json
        callset = path.parent.name
        group = path.parent.parent.name
    else:
        # hap.py: <group>/<callset>.summary.csv
        callset = path.name[: -len(".summary.csv")]
        group = path.parent.name
    bench = _BENCH_TYPE.search(callset)
    return {
        "analysis_id": f"{group}/{callset}",
        "eval_tool": eval_tool,
        "bench_type": bench.group(1) if bench else "",
    }


def _to_float(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.lower() in {"na", "nan", "none"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def load_happy_summary(path: Path) -> List[Dict]:
    """Records for the genome-wide PASS SNP and INDEL rows of a hap.py summary."""
    meta = parse_eval_path(path)
    records = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("Type") not in HAPPY_TYPES:
                continue
            if row.get("Filter") != HAPPY_FILTER:
                continue
            records.append(
                {
                    **meta,
                    "variant_type": row["Type"],
                    "precision": _to_float(row.get("METRIC.Precision")),
                    "recall": _to_float(row.get("METRIC.Recall")),
                    "f1": _to_float(row.get("METRIC.F1_Score")),
                    "tp": _to_float(row.get("TRUTH.TP")),
                    "fp": _to_float(row.get("QUERY.FP")),
                    "fn": _to_float(row.get("TRUTH.FN")),
                    "gt_concordance": None,
                    "source": str(path),
                }
            )
    return records


def load_truvari_summary(path: Path) -> List[Dict]:
    """Single SV record from a Truvari bench/refine summary JSON.

    ``refine.variant_summary.json`` is labelled ``SV-refine`` so it doesn't
    collide with the bench ``summary.json`` (``SV``) for the same analysis.
    """
    meta = parse_eval_path(path)
    variant_type = "SV-refine" if "refine" in path.name else "SV"
    with open(path) as fh:
        data = json.load(fh)
    return [
        {
            **meta,
            "variant_type": variant_type,
            "precision": _to_float(data.get("precision")),
            "recall": _to_float(data.get("recall")),
            "f1": _to_float(data.get("f1")),
            "tp": _to_float(data.get("TP-comp", data.get("TP-base"))),
            "fp": _to_float(data.get("FP")),
            "fn": _to_float(data.get("FN")),
            "gt_concordance": _to_float(data.get("gt_concordance")),
            "source": str(path),
        }
    ]


def load_records(files: List[Path]) -> List[Dict]:
    records: List[Dict] = []
    for path in files:
        if path.name.endswith(".summary.csv"):
            records.extend(load_happy_summary(path))
        elif path.suffix == ".json" and path.name.endswith("summary.json"):
            records.extend(load_truvari_summary(path))
        else:
            sys.stderr.write(f"Skipping unrecognized evaluation file: {path}\n")
    return records


def resolve_baseline(records: List[Dict], baseline: str) -> str:
    """Resolve a baseline argument to exactly one analysis_id."""
    ids = sorted({r["analysis_id"] for r in records})
    if baseline in ids:
        return baseline
    matches = [i for i in ids if baseline in i]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise SystemExit(f"Error: baseline '{baseline}' matched no analysis_id")
    raise SystemExit(
        "Error: baseline '{}' is ambiguous; matched:\n  {}".format(
            baseline, "\n  ".join(matches)
        )
    )


def add_deltas(
    records: List[Dict], baseline_id: str, metric: str, threshold: float
) -> List[Dict]:
    """Add delta_* columns and a status flag versus the baseline analysis.

    Baseline metrics are matched per variant_type. Status is based on --metric:
    ``improved`` / ``regressed`` when |delta| > threshold, else ``same``;
    baseline rows are marked ``baseline``.
    """
    base_by_type = {
        r["variant_type"]: r for r in records if r["analysis_id"] == baseline_id
    }
    for r in records:
        base = base_by_type.get(r["variant_type"])
        for m in METRICS:
            ref = base.get(m) if base else None
            cur = r.get(m)
            r[f"delta_{m}"] = (
                round(cur - ref, 6) if (ref is not None and cur is not None) else None
            )
        if r["analysis_id"] == baseline_id:
            r["status"] = "baseline"
            continue
        d = r.get(f"delta_{metric}")
        if d is None:
            r["status"] = "n/a"
        elif d > threshold:
            r["status"] = "improved"
        elif d < -threshold:
            r["status"] = "regressed"
        else:
            r["status"] = "same"
    return records


def sort_records(records: List[Dict], metric: str) -> List[Dict]:
    def sort_key(r: Dict):
        value = r.get(metric)
        return (
            r["variant_type"],
            -(value if value is not None else -1.0),
            r["analysis_id"],
        )

    return sorted(records, key=sort_key)


def _fmt(value, ndigits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{ndigits}f}"
    return str(value)


def format_markdown(records: List[Dict], with_deltas: bool, metric: str) -> str:
    cols = ["analysis_id", "eval_tool", "variant_type", "precision", "recall", "f1"]
    if with_deltas:
        cols += [f"delta_{metric}", "status"]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for r in records:
        lines.append("| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |")
    return "\n".join(lines)


def write_tsv(records: List[Dict], path: str, with_deltas: bool, metric: str) -> None:
    cols = [
        "analysis_id",
        "eval_tool",
        "bench_type",
        "variant_type",
        "precision",
        "recall",
        "f1",
        "tp",
        "fp",
        "fn",
        "gt_concordance",
    ]
    if with_deltas:
        cols += [f"delta_{m}" for m in METRICS] + ["status"]
    cols += ["source"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=cols, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for r in records:
            writer.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def compute_bed_stats(bed_path: Path) -> Dict[str, float]:
    """Compute total basepairs and interval count from a BED file."""
    total_bp = 0
    interval_count = 0
    try:
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start, end = int(parts[1]), int(parts[2])
                    total_bp += end - start
                    interval_count += 1
                except (ValueError, IndexError):
                    continue
    except OSError:
        pass
    return {"total_bp": total_bp, "interval_count": interval_count}


def load_bed_records(bed_files: List[Path]) -> List[Dict]:
    """Load benchmark region statistics from BED files."""
    records = []
    for bed in bed_files:
        # Parse bench_id from path: draft_benchmarksets/{bench_id}/{name}.benchmark.bed
        bench_id = bed.parent.name
        bench = _BENCH_TYPE.search(bench_id)
        stats = compute_bed_stats(bed)
        records.append(
            {
                "analysis_id": bench_id,
                "bench_type": bench.group(1) if bench else "",
                "total_bp": stats["total_bp"],
                "interval_count": stats["interval_count"],
                "source": str(bed),
            }
        )
    return records


def format_bed_markdown(records: List[Dict], baseline: Optional[str] = None) -> str:
    """Format BED statistics as Markdown table."""
    cols = ["analysis_id", "bench_type", "total_bp", "interval_count"]
    if baseline:
        cols += ["delta_bp", "delta_pct"]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for r in records:
        lines.append("| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |")
    return "\n".join(lines)


def add_bed_deltas(records: List[Dict], baseline_id: str) -> List[Dict]:
    """Add delta columns to BED records versus baseline."""
    base_by_type = {
        r["bench_type"]: r for r in records if r["analysis_id"] == baseline_id
    }
    for r in records:
        base = base_by_type.get(r["bench_type"])
        if base and r["analysis_id"] != baseline_id:
            base_bp = base.get("total_bp", 0)
            cur_bp = r.get("total_bp", 0)
            if base_bp > 0:
                r["delta_bp"] = cur_bp - base_bp
                r["delta_pct"] = round(100.0 * (cur_bp - base_bp) / base_bp, 2)
            else:
                r["delta_bp"] = None
                r["delta_pct"] = None
        else:
            r["delta_bp"] = None
            r["delta_pct"] = None
    return records


def run_plot(tsv_path: str, out_dir: str) -> int:
    rscript = shutil.which("Rscript")
    if not rscript:
        sys.stderr.write("Warning: Rscript not found; skipping --plot\n")
        return 1
    script = Path(__file__).resolve().parent / "plot_evaluations.R"
    cmd = [rscript, str(script), tsv_path, out_dir]
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    return subprocess.run(cmd).returncode


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    files = (
        [Path(p) for p in args.inputs]
        if args.inputs
        else discover_eval_files(args.results_dir)
    )
    if not files:
        sys.exit(
            "Error: no evaluation outputs found "
            f"(looked under {args.results_dir}/evaluations/)"
        )
    records = load_records(files)
    if not records:
        sys.exit("Error: no parseable evaluation records found")

    baseline_id: Optional[str] = None
    with_deltas = bool(args.baseline)
    if with_deltas:
        baseline_id = resolve_baseline(records, args.baseline)
        add_deltas(records, baseline_id, args.metric, args.threshold)
        print(f"Baseline: {baseline_id}", file=sys.stderr)

    records = sort_records(records, args.metric)
    print(format_markdown(records, with_deltas, args.metric))

    # Add benchmark region comparison if requested
    if args.regions:
        bed_files = discover_benchmark_beds(args.results_dir)
        if bed_files:
            print("\n## Benchmark Regions\n")
            bed_records = load_bed_records(bed_files)
            bed_baseline_str = None
            if with_deltas and baseline_id:
                # Match baseline_id to BED records (may need fuzzy match)
                bed_baseline = None
                for br in bed_records:
                    if baseline_id in br["analysis_id"] or br["analysis_id"] in baseline_id:
                        bed_baseline = br["analysis_id"]
                        break
                if bed_baseline:
                    add_bed_deltas(bed_records, bed_baseline)
                    bed_baseline_str = bed_baseline
            print(format_bed_markdown(bed_records, bed_baseline_str))
        else:
            print("\nNo benchmark BED files found", file=sys.stderr)

    if args.out:
        write_tsv(records, args.out, with_deltas, args.metric)
        print(f"Wrote {len(records)} rows to {args.out}", file=sys.stderr)

    if args.plot:
        if not args.out:
            sys.exit("Error: --plot requires --out (the R script reads the TSV)")
        return run_plot(args.out, args.plot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
