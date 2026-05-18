#!/usr/bin/env python
"""Summarize a Snakemake dry-run for a DeFrABB analyses table.

Runs `snakemake --dryrun --quiet=rules` for the given analyses TSV, parses the
"Job stats" table for rule counts, and reads the analyses TSV directly to
report which assemblies, references, comparison callsets, exclusion sets, and
evaluations are touched. Self-contained: no rule or schema changes.
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-a", "--analyses", required=True, help="Path to analyses TSV")
    p.add_argument(
        "--snakefile",
        default="Snakefile",
        help="Snakefile path (default: Snakefile in cwd)",
    )
    p.add_argument(
        "--top",
        type=int,
        default=10,
        help="Show the N rules with the largest job counts (default: 10)",
    )
    p.add_argument(
        "--all-rules",
        action="store_true",
        help="List every rule, not just the top N",
    )
    return p.parse_args()


def run_dryrun(analyses: str, snakefile: str) -> str:
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--dryrun",
        "--quiet=rules",
        "--config",
        f"analyses={analyses}",
    ]
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        sys.stderr.write(result.stdout)
        raise SystemExit(f"snakemake --dryrun failed (exit {result.returncode})")
    return result.stdout


_JOB_LINE = re.compile(r"^(\S+)\s+(\d+)\s*$")


def parse_job_stats(dryrun_output: str) -> Tuple[Dict[str, int], int]:
    """Return (rule_name -> count, total). Skips header/separator/total."""
    counts: Dict[str, int] = {}
    total = 0
    in_table = False
    for line in dryrun_output.splitlines():
        if line.startswith("Job stats:"):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("---") or line.startswith("job "):
            continue
        m = _JOB_LINE.match(line)
        if not m:
            if line.strip() == "":
                if counts:
                    break
                continue
            break
        name, n = m.group(1), int(m.group(2))
        if name == "total":
            total = n
            break
        counts[name] = n
    return counts, total


def read_analyses(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row for row in reader]


def unique_ordered(values) -> List[str]:
    return list(OrderedDict.fromkeys(v for v in values if v))


def summarize_analyses(rows: List[Dict[str, str]]) -> Dict[str, List[str]]:
    cols = [
        "ref",
        "asm_id",
        "vc_cmd",
        "bench_type",
        "exclusion_set",
        "eval_comp_id",
        "eval_cmd",
    ]
    return {c: unique_ordered(r.get(c, "") for r in rows) for c in cols}


def fmt_list(values: List[str]) -> str:
    return ", ".join(values) if values else "(none)"


def print_summary(
    analyses_path: str,
    rule_counts: Dict[str, int],
    total: int,
    rows: List[Dict[str, str]],
    top: int,
    all_rules: bool,
) -> None:
    summary = summarize_analyses(rows)
    print("=" * 60)
    print(f"DeFrABB dry-run summary: {analyses_path}")
    print("=" * 60)
    print(f"Analyses rows         : {len(rows)}")
    print(f"References            : {fmt_list(summary['ref'])}")
    print(f"Assemblies            : {fmt_list(summary['asm_id'])}")
    print(f"Variant callers       : {fmt_list(summary['vc_cmd'])}")
    print(f"Benchmark types       : {fmt_list(summary['bench_type'])}")
    print(f"Exclusion sets        : {fmt_list(summary['exclusion_set'])}")
    print(f"Comparison callsets   : {fmt_list(summary['eval_comp_id'])}")
    print(f"Evaluation tools      : {fmt_list(summary['eval_cmd'])}")
    print()
    print(f"Total jobs to run     : {total}")
    print(f"Distinct rules        : {len(rule_counts)}")
    print()
    sorted_rules = sorted(rule_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if all_rules:
        shown = sorted_rules
        header = "All rules (by job count)"
    else:
        shown = sorted_rules[:top]
        header = f"Top {min(top, len(sorted_rules))} rules (by job count)"
    print(header)
    print("-" * len(header))
    width = max((len(name) for name, _ in shown), default=0)
    for name, n in shown:
        print(f"  {name.ljust(width)}  {n:>5}")


def main() -> int:
    args = parse_args()
    if not Path(args.analyses).is_file():
        sys.exit(f"Error: analyses file not found: {args.analyses}")
    output = run_dryrun(args.analyses, args.snakefile)
    rule_counts, total = parse_job_stats(output)
    if not rule_counts:
        sys.stderr.write(output)
        sys.exit("Error: could not parse a 'Job stats' table from snakemake output")
    rows = read_analyses(args.analyses)
    print_summary(args.analyses, rule_counts, total, rows, args.top, args.all_rules)
    return 0


if __name__ == "__main__":
    sys.exit(main())
