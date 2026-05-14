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

    ## Initiating empty list for storing paths for beds to excluded from
    ## diploid assembled regions
    exc_paths = []
    for exclusion in exclusion_set:
        ## Determining path for asm specific exclusions and asm agnostic exclusions
        if exclusion in config["exclusion_asm_agnostic"]:
            exc_path = f"resources/exclusions/{{ref_id}}/{exclusion}"
        else:
            exc_path = f"results/draft_benchmarksets/{{bench_id}}/exclusions/{{ref_id}}_{{asm_id}}_{{bench_type}}_{{vc_cmd}}-{{vc_param_id}}_{exclusion}"

        ## Adding slop - currently a 15kb hard coded buffer around excluded repeat regions
        if exclusion in config["exclusion_slop_regions"]:
            exc_path = f"{exc_path}_slop"
        ## Adding slop then merging - hard coded 15kb slop then merging with 10kb hard coded dist
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
