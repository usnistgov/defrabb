## Exclusion slop/merge resolvers — named functions for use in rule params:
def get_slop_value(wildcards):
    overrides = (
        config["_exclusion_params"]
        .get("overrides", {})
        .get(wildcards.genomic_region, {})
    )
    if "slop_pct" in overrides:
        return overrides["slop_pct"]
    return overrides.get("slop", config["_exclusion_params"]["slop"])


def get_slop_flags(wildcards):
    overrides = (
        config["_exclusion_params"]
        .get("overrides", {})
        .get(wildcards.genomic_region, {})
    )
    return "-pct" if "slop_pct" in overrides else ""


def get_merge_dist(wildcards):
    overrides = (
        config["_exclusion_params"]
        .get("overrides", {})
        .get(wildcards.genomic_region, {})
    )
    return overrides.get(
        "slopmerge_dist", config["_exclusion_params"]["slopmerge_dist"]
    )


def get_bench_exclusion_set_id(wildcards):
    return bench_tbl.loc[wildcards.bench_id, "exclusion_set"]


def get_bench_exclusion_profile(wildcards):
    """Return the exclusion profile name for a benchmark, defaulting to 'standard'."""
    if "exclusion_profile" in bench_tbl.columns:
        return bench_tbl.loc[wildcards.bench_id, "exclusion_profile"]
    return "standard"


def _bench_profile_param(wildcards, param_key: str):
    """Look up a profile parameter for a benchmark, falling back to _exclusion_params."""
    profile_name = get_bench_exclusion_profile(wildcards)
    profile = config.get("_exclusion_profiles", {}).get(profile_name, {})
    return profile.get(param_key, config["_exclusion_params"].get(param_key))


## Benchmark VCF / BED standardization + exclusion-input helpers
def get_processed_vcf(wildcards):
    # Filter rows based on bench_type using the query method
    filtered_df = bench_tbl.query(f'bench_type == "{wildcards.bench_type}"')
    # Further filter the DataFrame based on bench_id using the loc method
    subset_df = filtered_df.loc[[wildcards.bench_id]]

    # Remove duplicate entries
    subset_df = subset_df.drop_duplicates()

    # Ensure that only one unique row remains after removing duplicates
    assert (
        subset_df.shape[0] == 1
    ), f"Error: Multiple entries found for bench_id {wildcards.bench_id} and bench_type {wildcards.bench_type}"

    # Now, you can grab the value of vc_id and bench_vcf_processing from the first row
    vc_id = subset_df.iloc[0]["vc_id"]
    profile_name = subset_df.iloc[0]["bench_vcf_processing"]

    if profile_name == "none":
        return f"results/asm_varcalls/{vc_id}/annotations/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.vcf.gz"
    else:
        # Resolve profile name to ordered list of steps, then join with dots
        steps = config["vcf_processing_profiles"][profile_name]
        vcf_suffix = ".".join(steps)
        return f"results/asm_varcalls/{vc_id}/annotations/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.{vcf_suffix}.vcf.gz"


def get_std_base(wildcards):
    if wildcards.get("vc_id", ""):
        vc_id = wildcards.vc_id
    else:
        vc_id = bench_tbl.loc[wildcards.bench_id, "vc_id"]
    ref_id = wildcards.ref_id
    asm_id = wildcards.asm_id
    vc_cmd = wildcards.vc_cmd
    vc_param_id = wildcards.vc_param_id
    return f"results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}"


def get_standardized_vcf(wildcards):
    basename = get_std_base(wildcards)
    return f"{basename}.vcf.gz"


def get_standardized_vcfidx(wildcards):
    basename = get_std_base(wildcards)
    return f"{basename}.vcf.gz.tbi"


def get_standardized_bed(wildcards):
    basename = get_std_base(wildcards)
    return f"{basename}.baseline.bed"


# Update draft benchmark generation to use standardized outputs
def get_draft_benchmark_inputs(wildcards):
    return {
        "vcf": get_standardized_vcf(wildcards),
        "vcfidx": get_standardized_vcfidx(wildcards),
        "bed": get_standardized_bed(wildcards),
    }


## Cross-caller exclusion inputs
def get_consecutive_svs_bams(wildcards):
    """Hap BAMs for the consecutive-svs exclusion.

    consecutive-svs is computed from dipcall's assembly-to-reference hap BAMs.
    Resolve to the *existing* dipcall run for this reference + assembly (via
    ``get_asm_varcall_run``) rather than the current benchmark's own vc_id, so a
    PAV benchmark does not trigger a redundant dipcall run inside its directory.
    For a dipcall benchmark this resolves back to its own run.

    In parameter sweeps with multiple dipcall runs, prefer the one matching this
    benchmark's vc_param_id.
    """
    dip_vc_id, dip_param_id = get_asm_varcall_run(
        wildcards.ref_id,
        wildcards.asm_id,
        "dipcall",
        prefer_vc_param_id=wildcards.vc_param_id,
    )
    base = (
        f"results/asm_varcalls/{dip_vc_id}/"
        f"{wildcards.ref_id}_{wildcards.asm_id}_dipcall-{dip_param_id}"
    )
    return {
        "hap1_bam": f"{base}.hap1.bam",
        "hap1_bai": f"{base}.hap1.bam.bai",
        "hap2_bam": f"{base}.hap2.bam",
        "hap2_bai": f"{base}.hap2.bam.bai",
    }


def get_pav_discrep_runs(wildcards):
    """Resolve the (vc_id, vc_param_id) dipcall and PAV runs for a ref + assembly.

    The pav-discrep / pav-inv exclusions compare dipcall and PAV calls for the
    same reference + assembly, so both runs must be declared in the analyses
    table. Resolved via ``get_asm_varcall_run`` so a benchmark of either caller
    reuses the existing runs instead of triggering duplicate caller runs.
    """
    return {
        vc_cmd: get_asm_varcall_run(
            wildcards.ref_id,
            wildcards.asm_id,
            vc_cmd,
            prefer_vc_param_id=wildcards.get(
                f"{vc_cmd}_param", wildcards.get("vc_param_id")
            ),
        )
        for vc_cmd in ["dipcall", "pav"]
    }


def get_pav_discrep_inputs(wildcards):
    """Standardized dipcall and PAV calls/beds (keys dip_* / pav_*)."""
    inputs = {}
    for vc_cmd, (vc_id, vc_param_id) in get_pav_discrep_runs(wildcards).items():
        base = (
            f"results/asm_varcalls/{vc_id}/"
            f"{wildcards.ref_id}_{wildcards.asm_id}_{vc_cmd}-{vc_param_id}"
        )
        key = "dip" if vc_cmd == "dipcall" else "pav"
        inputs[f"{key}_vcf"] = f"{base}.vcf.gz"
        inputs[f"{key}_vcfidx"] = f"{base}.vcf.gz.tbi"
        inputs[f"{key}_bed"] = f"{base}.baseline.bed"
    return inputs


def get_pav_discrep_fpfns_bed(wildcards):
    """Shared (per ref + assembly + runs) PAV vs dipcall FP/FN bed for a benchmark."""
    runs = get_pav_discrep_runs(wildcards)
    method = "happy" if wildcards.discrep_type == "smvar" else "truvari"
    return (
        f"results/pav_discrep/{wildcards.ref_id}_{wildcards.asm_id}"
        f"_dipcall-{runs['dipcall'][1]}_pav-{runs['pav'][1]}/{method}.fpfns.bed"
    )


## Exclusions
def get_exclusion_inputs(wildcards):
    ## Getting list of excluded regions
    exclusion_set_id = bench_tbl.loc[wildcards.bench_id, "exclusion_set"]
    if exclusion_set_id == "none":
        return []
    try:
        exclusion_set = config["exclusion_set"][exclusion_set_id]
    except KeyError:
        print(f"{exclusion_set_id} is not defined in resources yaml")

    ## Determine exclusion profile for this benchmark.
    ## Non-standard profiles use a subdirectory in resources/exclusions so that
    ## agnostic beds (gaps, etc.) with different slop values get distinct paths.
    excl_profile = get_bench_exclusion_profile(wildcards)
    use_profile_dir = excl_profile and excl_profile != "standard"

    ## Initiating empty list for storing paths for beds to excluded from
    ## diploid assembled regions
    exc_paths = []
    for exclusion in exclusion_set:
        ## Determining path for asm specific exclusions and asm agnostic exclusions.
        ## Agnostic exclusions are profile-scoped when profile is non-standard so
        ## different slop values produce distinct resource files.
        if exclusion in config["exclusion_asm_agnostic"]:
            if use_profile_dir:
                exc_path = f"resources/exclusions/{{ref_id}}/{excl_profile}/{exclusion}"
            else:
                exc_path = f"resources/exclusions/{{ref_id}}/{exclusion}"
        else:
            exc_path = f"results/draft_benchmarksets/{{bench_id}}/exclusions/{{ref_id}}_{{asm_id}}_{{bench_type}}_{{vc_cmd}}-{{vc_param_id}}_{exclusion}"

        ## Adding slop around excluded regions
        if exclusion in config["exclusion_slop_regions"]:
            exc_path = f"{exc_path}_slop"
        ## Adding slop then merging adjacent intervals
        elif exclusion in config["exclusion_slopmerge_regions"]:
            exc_path = f"{exc_path}_slopmerge"

        ## Ensuring bed files are sorted before intersect
        if exclusion in config["exclusion_asm_intersect"]:
            exc_path = f"{exc_path}_sorted"

        ## Defining which regions are excluded based on diploid assembly breaks
        if exclusion in config["exclusion_asm_intersect"]:
            exc_paths += [f"{exc_path}_start", f"{exc_path}_end"]
        else:
            exc_paths = exc_paths + [exc_path]

    ## Adding to exc_paths list and ensuring all beds are sorted
    ## prior to exclusion from dip assembled regions
    exclusion_paths = [f"{exc}_sorted.bed" for exc in exc_paths]

    ## Returning list of bed paths for exclusion
    return exclusion_paths
