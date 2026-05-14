"""Write a YAML provenance file documenting every exclusion applied to a benchmark set.

Invoked via Snakemake script: directive. Reads snakemake.config directly.

Inputs:
    snakemake.input.intersection_summary  -- exclusion_intersection_summary.csv path
    snakemake.params.exclusion_set_id     -- exclusion set name from analyses table
    snakemake.output[0]                   -- output YAML path
    snakemake.log[0]                      -- log path
"""

import logging

import pandas as pd
import yaml


def _resolve_slop(region, config):
    """Return (value, is_pct) for the given exclusion region."""
    overrides = (
        config["_exclusion_params"].get("overrides", {}).get(region, {})
    )
    if "slop_pct" in overrides:
        return overrides["slop_pct"], True
    return overrides.get("slop", config["_exclusion_params"]["slop"]), False


def _resolve_merge_dist(region, config):
    overrides = (
        config["_exclusion_params"].get("overrides", {}).get(region, {})
    )
    return overrides.get("slopmerge_dist", config["_exclusion_params"]["slopmerge_dist"])


def _source_type(region, config):
    if region in config["exclusion_asm_agnostic"]:
        return "asm_agnostic"
    if region in config["exclusion_ref_agnostic"]:
        return "ref_agnostic"
    return "asm_specific"


def _transform(region, config):
    if region in config["exclusion_slop_regions"]:
        return "slop"
    if region in config["exclusion_slopmerge_regions"]:
        return "slopmerge"
    return "none"


def build_provenance(exclusion_set_id, config, bp_impact):
    """Build the full provenance dict."""
    exclusion_list = config["exclusion_set"][exclusion_set_id]
    excl_params = config["_exclusion_params"]

    entries = []
    for region in exclusion_list:
        transform = _transform(region, config)
        slop_val, is_pct = _resolve_slop(region, config)
        entry = {
            "name": region,
            "source_type": _source_type(region, config),
            "transform": transform,
            "asm_intersect": region in config["exclusion_asm_intersect"],
            "bp_impact": bp_impact.get(region, {}).get("bp", None),
            "pct_of_initial": bp_impact.get(region, {}).get("pct", None),
        }
        if transform in ("slop", "slopmerge"):
            if is_pct:
                entry["slop_pct"] = slop_val
            else:
                entry["slop_bp"] = slop_val
        if transform == "slopmerge":
            entry["merge_dist"] = _resolve_merge_dist(region, config)
        entries.append(entry)

    return {
        "exclusion_set": exclusion_set_id,
        "global_defaults": {
            "slop": excl_params["slop"],
            "slopmerge_dist": excl_params["slopmerge_dist"],
        },
        "exclusions": entries,
    }


def load_bp_impact(intersection_summary_path):
    """Load per-exclusion bp impact from intersection summary CSV.

    The summary CSV uses file basenames as region names; strip the path suffix
    pattern to recover the exclusion name (e.g. 'segdups_slopmerge_sorted_start_sorted').
    We match on prefix so that both start and end halves map to the same region.
    """
    df = pd.read_csv(intersection_summary_path)
    # Skip header/total rows
    df = df[~df["genomic_region"].isin(["initial", "benchmark_regions", "total_excluded"])]
    bp_by_region = {}
    for _, row in df.iterrows():
        name = row["genomic_region"]
        # Strip known suffixes to recover the base exclusion name
        for suffix in (
            "_slopmerge_sorted_start_sorted.bed",
            "_slopmerge_sorted_end_sorted.bed",
            "_slopmerge_sorted.bed",
            "_slop_sorted.bed",
            "_sorted.bed",
            ".bed",
        ):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        # The remaining name may still have a path prefix for asm-specific
        # exclusions; take just the final segment after the last underscore group
        # that matches a known exclusion region name.
        base = name.split("_")[-1] if "_" in name else name
        entry = bp_by_region.setdefault(base, {"bp": 0, "pct": 0.0})
        entry["bp"] += int(row["bp"])
        entry["pct"] += float(row["pct_of_initial"])
    return bp_by_region


def main():
    logging.basicConfig(
        filename=snakemake.log[0],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    config = snakemake.config
    exclusion_set_id = snakemake.params.exclusion_set_id
    intersection_summary_path = snakemake.input.intersection_summary
    output_path = snakemake.output[0]

    logging.info(f"Building provenance for exclusion set: {exclusion_set_id}")
    bp_impact = load_bp_impact(intersection_summary_path)
    provenance = build_provenance(exclusion_set_id, config, bp_impact)

    with open(output_path, "w") as fh:
        yaml.dump(provenance, fh, default_flow_style=False, sort_keys=False)
    logging.info(f"Provenance written to {output_path}")


if "snakemake" in dir():
    main()
