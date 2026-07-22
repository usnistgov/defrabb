#!/usr/bin/env python3
"""Generate parameter sweep analyses tables from declarative YAML configs.

Takes a sweep configuration (fixed dimensions + sweep dimensions) and generates
the cross-product as a DeFrABB analyses table. Includes validation, cost
estimation, dry-run mode, and vc_id reuse intelligence.

Part of v0.023 parameter optimization framework.
"""

import argparse
import sys
import yaml
from pathlib import Path
from itertools import product
from typing import Dict, List, Any
import pandas as pd


# Runtime estimates per operation (hours)
RUNTIME_ESTIMATES = {
    "dipcall": 5.0,
    "pav": 6.0,
    "happy": 0.5,
    "truvari": 1.0,
    "unhappy": 0.5,
}

# Storage estimates per output (GB)
STORAGE_ESTIMATES = {
    "dipcall": 50,
    "pav": 80,
    "happy": 5,
    "truvari": 10,
    "unhappy": 5,
}


def load_sweep_config(path: Path) -> Dict[str, Any]:
    """Load and validate sweep YAML config."""
    if not path.exists():
        raise FileNotFoundError(f"Sweep config not found: {path}")

    with open(path) as f:
        config = yaml.safe_load(f)

    # Required fields
    required = {"name", "fixed", "sweep"}
    missing = required - set(config.keys())
    if missing:
        raise ValueError(f"Sweep config missing required fields: {missing}")

    # Validate sweep dimensions are lists
    for dim, values in config["sweep"].items():
        if not isinstance(values, list):
            raise ValueError(f"Sweep dimension '{dim}' must be a list, got {type(values).__name__}")
        if not values:
            raise ValueError(f"Sweep dimension '{dim}' is empty")

    return config


def generate_cross_product(fixed: Dict[str, Any], sweep: Dict[str, List[Any]]) -> List[Dict[str, Any]]:
    """Generate all combinations of sweep dimensions combined with fixed fields."""
    sweep_dims = sorted(sweep.keys())  # Deterministic order
    sweep_values = [sweep[dim] for dim in sweep_dims]

    analyses = []
    for combo in product(*sweep_values):
        row = fixed.copy()
        for dim, value in zip(sweep_dims, combo):
            row[dim] = value
        analyses.append(row)

    return analyses


def assign_vc_ids(analyses: List[Dict[str, Any]]) -> int:
    """Assign shared vc_ids to analyses with matching variant call params.

    Returns number of unique variant call runs needed.
    """
    vc_groups = {}

    for row in analyses:
        # Key = (ref_id, asm_id, vc_cmd, vc_param_id)
        # Same key → same variant calls → share vc_id
        vc_param = row.get("vc_param_id", "default")
        key = (row["ref_id"], row["asm_id"], row["vc_cmd"], vc_param)

        if key not in vc_groups:
            vc_id_num = len(vc_groups) + 1
            vc_groups[key] = f"vc{vc_id_num:03d}"

        row["vc_id"] = vc_groups[key]

    return len(vc_groups)


def estimate_costs(analyses: List[Dict[str, Any]], unique_vc_runs: int) -> Dict[str, float]:
    """Estimate runtime (hours) and storage (GB) for sweep."""
    total_runtime = 0.0
    total_storage = 0.0

    # Variant calling cost (shared across analyses with same vc_id)
    vc_cmd_counts = {}
    for row in analyses:
        vc_cmd = row["vc_cmd"]
        vc_cmd_counts[vc_cmd] = vc_cmd_counts.get(vc_cmd, 0) + 1

    # Only count unique variant calls
    for vc_cmd, count in vc_cmd_counts.items():
        unique_count = unique_vc_runs if vc_cmd in ["dipcall", "pav"] else count
        total_runtime += RUNTIME_ESTIMATES.get(vc_cmd, 2.0) * unique_count
        total_storage += STORAGE_ESTIMATES.get(vc_cmd, 30) * unique_count

    # Evaluation cost (per analysis)
    for row in analyses:
        eval_cmd = row.get("eval_cmd", "happy")
        total_runtime += RUNTIME_ESTIMATES.get(eval_cmd, 0.5)
        total_storage += STORAGE_ESTIMATES.get(eval_cmd, 5)

    return {
        "runtime_hours": round(total_runtime, 1),
        "storage_gb": round(total_storage, 1),
        "unique_vc_runs": unique_vc_runs,
        "total_analyses": len(analyses),
    }


def format_analyses_table(analyses: List[Dict[str, Any]]) -> pd.DataFrame:
    """Convert analyses list to DataFrame with standard column order."""
    df = pd.DataFrame(analyses)

    # Standard column order (based on analyses.tsv schema)
    standard_cols = [
        "vc_id", "ref_id", "asm_id", "vc_cmd", "vc_param_id", "vc_params",
        "bench_type", "bench_id", "exclusion_set", "vcf_processing",
        "eval_id", "eval_cmd", "eval_comp_id", "eval_params"
    ]

    # Use standard order for columns that exist, append any extras
    present_cols = [c for c in standard_cols if c in df.columns]
    extra_cols = [c for c in df.columns if c not in standard_cols]
    df = df[present_cols + extra_cols]

    # Fill missing standard columns with defaults
    defaults = {
        "vc_params": "",
        "vc_param_id": "default",
        "exclusion_set": "default",
        "vcf_processing": "default",
        "eval_params": "default",
    }
    for col, default in defaults.items():
        if col not in df.columns:
            df[col] = default

    # Generate derived IDs if not provided
    # bench_id must be unique per evaluation, includes exclusion_set
    if "bench_id" not in df.columns:
        exclusion_suffix = df["exclusion_set"].apply(
            lambda x: f"_{x}" if x and x != "default" else ""
        )
        df["bench_id"] = (
            df["ref_id"] + "_" + df["asm_id"] + "_" +
            df["bench_type"] + "_" + df["vc_cmd"] + "-" + df["vc_param_id"] +
            exclusion_suffix
        )

    if "eval_id" not in df.columns:
        df["eval_id"] = df["eval_cmd"] + "_" + df["eval_comp_id"]

    return df


def validate_output(df: pd.DataFrame) -> None:
    """Validate output table against schema expectations."""
    required_cols = [
        "vc_id", "ref_id", "asm_id", "vc_cmd", "bench_type",
        "bench_id", "eval_id", "eval_cmd", "eval_comp_id"
    ]

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Output table missing required columns: {missing}")

    # Check for duplicates in primary key
    pk_cols = ["eval_id", "bench_id"]
    duplicates = df[df.duplicated(subset=pk_cols, keep=False)]
    if not duplicates.empty:
        raise ValueError(
            f"Duplicate (eval_id, bench_id) combinations found:\n"
            f"{duplicates[pk_cols].to_string()}"
        )

    # Warn if any cells are null
    null_counts = df.isnull().sum()
    null_cols = null_counts[null_counts > 0]
    if not null_cols.empty:
        print(f"WARNING: Null values found in columns: {null_cols.to_dict()}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Sweep configuration YAML file"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output TSV path (default: from config.output field)"
    )
    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Print summary and preview without writing output"
    )
    parser.add_argument(
        "--max-analyses",
        type=int,
        default=200,
        help="Maximum number of analyses to generate (safety limit, default: 200)"
    )
    parser.add_argument(
        "--preview-rows",
        type=int,
        default=5,
        help="Number of rows to show in dry-run preview (default: 5)"
    )

    args = parser.parse_args()

    # Load config
    try:
        config = load_sweep_config(args.config)
    except Exception as e:
        print(f"ERROR loading config: {e}", file=sys.stderr)
        return 1

    # Generate cross-product
    analyses = generate_cross_product(config["fixed"], config["sweep"])

    # Check max analyses limit
    if len(analyses) > args.max_analyses:
        print(
            f"ERROR: Cross-product would generate {len(analyses)} analyses, "
            f"exceeding limit of {args.max_analyses}.",
            file=sys.stderr
        )
        print(
            f"Increase --max-analyses or reduce sweep dimensions.",
            file=sys.stderr
        )
        return 1

    # Assign vc_ids for output reuse
    unique_vc_runs = assign_vc_ids(analyses)

    # Estimate costs
    costs = estimate_costs(analyses, unique_vc_runs)

    # Format as DataFrame
    df = format_analyses_table(analyses)

    # Validate
    try:
        validate_output(df)
    except ValueError as e:
        print(f"ERROR: Output validation failed: {e}", file=sys.stderr)
        return 1

    # Print summary
    print(f"Sweep: {config['name']}")
    print(f"  Analyses: {costs['total_analyses']}")
    print(f"  Unique variant call runs: {costs['unique_vc_runs']}")
    reuse_factor = costs['total_analyses'] / costs['unique_vc_runs']
    print(f"  Reuse factor: {reuse_factor:.1f}x")
    print(f"  Estimated runtime: {costs['runtime_hours']} hours")
    print(f"  Estimated storage: {costs['storage_gb']} GB")

    if args.dry_run:
        print(f"\nPreview (first {args.preview_rows} rows):")
        print(df.head(args.preview_rows).to_string(index=False))
        print(f"\nDry-run mode: output NOT written.")
        return 0

    # Determine output path
    output_path = args.output or config.get("output")
    if not output_path:
        print(
            "ERROR: No output path specified (use --output or config.output)",
            file=sys.stderr
        )
        return 1

    output_path = Path(output_path)

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"\nWrote {len(df)} analyses to: {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
