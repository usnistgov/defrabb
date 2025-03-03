from pathlib import Path
from snakemake.utils import validate
import pandas as pd


################################################################################
## Config processing functions


def load_df(path, schema):
    df = pd.read_table(path, dtype={"eval_target_regions": bool})
    validate(df, schema)
    return df


def load_analyses(path, schema):
    return load_df(path, schema).astype(dtype={"eval_target_regions": bool})


def _filter_subtable(df, filter_re, id_cols, new_index):
    params = df.filter(regex=filter_re).drop_duplicates()
    ids = df[[new_index] + id_cols].drop_duplicates()
    tbl = pd.merge(ids, params, how="inner", on=new_index).set_index(new_index)
    return (params, tbl)


def analyses_to_vc_tbl(analyses):
    return _filter_subtable(analyses, "vc_", ["asm_id", "ref"], "vc_id")


def analyses_to_bench_tbls(analyses):
    id_cols = [
        "asm_id",
        "vc_id",
        "vc_cmd",
        "vc_param_id",
        "ref",
        "exclusion_set",
    ]
    params, tbl = _filter_subtable(analyses, "bench_", id_cols, "bench_id")
    excluded_tbl = tbl[tbl.exclusion_set != "none"]
    return (params, tbl, excluded_tbl)


################################################################################
## Rule parameters
def get_ref_id(wildcards):
    ref = wildcards.get("ref", "")
    if ref:
        ref_id = ref
    else:
        ref_id = wildcards.get("ref_id", "")
        if not ref_id:
            prefix = wildcards.get("prefix", "")
            for id in REFIDS:
                if id in prefix:
                    ref_id = id
    if not ref_id:
        print(f"ref_id could not be determined from wildcards or {prefix}")

    return ref_id


def get_ref_file(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.fa"


def get_ref_index(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.fa.fai"


def get_ref_bwaindex(wildcards):
    ref_id = get_ref_id(wildcards)
    return [
        f"resources/references/{ref_id}.fa.{ext}"
        for ext in ["amb", "ann", "bwt", "pac", "sa"]
    ]


def get_genome_file(wildcards):
    ref_id = get_ref_id(wildcards)
    return workflow.source_path(f"../resources/references/{ref_id}.genome")


def get_ref_sdf(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.sdf"


def get_ref_trdb(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}_adotto_trf.bed.gz"


def get_addoto_tr_anno_db_url(wildcards):
    return ref_config[wildcards.ref_id]["adotto_db"]


def get_trf_db_url(wildcards):
    return ref_config[wildcards.ref_id]["trf_db"]


def get_sample_id(wildcards):
    return asm_config[wildcards.asm_id]["sample_id"]


def get_par_bed(wildcards):
    root = config["_par_bed_root"]
    ref_id = get_ref_id(wildcards)
    filename = ref_config[ref_id]["par_bed"]
    return Path(workflow.basedir) / root / filename


def get_dipcall_par_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    par_path = get_par_bed(wildcards)
    return f"-x {par_path}" if is_male else ""


## Happy Inputs and Parameters
def get_happy_gender_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    sex = "male" if is_male else "female"
    return f"--gender {sex}"


def get_happy_inputs(analyses, config, wildcards):
    return get_happy_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


def get_happy_inputs_inner(ref_id, eval_id, analyses, config):
    ## Creating empty dictionary for storing inputs
    inputs = {}

    ## Reference genome and stratifications
    inputs["genome"] = f"resources/references/{ref_id}.fa"
    inputs["genome_index"] = f"resources/references/{ref_id}.fa.fai"
    strat_tb = config["references"][ref_id]["stratifications"]["tarball"]
    inputs["strat_tb"] = f"resources/strats/{ref_id}/{strat_tb}"

    ## draft benchmark variant calls
    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz"
    draft_bench_vcfidx = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi"

    ## draft benchmark regions
    if analyses.loc[eval_id, "exclusion_set"] == "none":
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed"
    else:
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed"

    ## comparison variant call paths
    comp_vcf = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz"
    comp_vcfidx = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi"
    comp_bed = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    ## Determining which callsets and regions are used as truth
    if analyses.loc[eval_id, "eval_comp_id_is_truth"]:
        query = "draft_bench"
    else:
        query = "comp"

    ## Defining truth calls and regions along with query calls
    if query == "draft_bench":
        inputs["query"] = draft_bench_vcf
        inputs["query_vcfidx"] = draft_bench_vcfidx
        inputs["truth"] = comp_vcf
        inputs["truth_vcfidx"] = comp_vcfidx
        inputs["truth_regions"] = comp_bed
    else:
        inputs["query"] = comp_vcf
        inputs["query_vcfidx"] = comp_vcfidx
        inputs["truth"] = draft_bench_vcf
        inputs["truth_vcfidx"] = draft_bench_vcfidx
        inputs["truth_regions"] = draft_bench_bed

    ## Determining Target regions
    if analyses.loc[eval_id, "eval_target_regions"]:
        if query == "draft_bench":
            inputs["target_regions"] = draft_bench_bed
        else:
            inputs["target_regions"] = comp_bed

    ## Returning happy inputs
    return inputs


def get_eval_beds(analyses, wildcards):
    bench_id = analyses.loc[wildcards.eval_id, "bench_id"]
    ref_id = analyses.loc[wildcards.eval_id, "ref"]
    asm_id = analyses.loc[wildcards.eval_id, "asm_id"]
    bench_type = analyses.loc[wildcards.eval_id, "bench_type"]
    vc_cmd = analyses.loc[wildcards.eval_id, "vc_cmd"]
    vc_param_id = analyses.loc[wildcards.eval_id, "vc_param_id"]
    comp_id = analyses.loc[wildcards.eval_id, "eval_comp_id"]

    ## Getting draft benchmark bed
    if analyses.loc[wildcards.eval_id, "exclusion_set"] == "none":
        bench_bed = f"results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed"
    else:
        bench_bed = f"results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed"

    ## Getting comparison bed
    comp_bed = f"resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    return {"left": bench_bed, "right": comp_bed}


def get_truvari_inputs(analyses, config, wildcards):
    return get_truvari_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


def get_truvari_inputs_inner(ref_id, eval_id, analyses, config):
    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz"
    draft_bench_vcfidx = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi"

    if analyses.loc[eval_id, "exclusion_set"] == "none":
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed"
    else:
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed"

    comp_vcf = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz"
    comp_vcfidx = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi"
    comp_bed = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    draft_is_query = bool(analyses.loc[eval_id, "eval_comp_id_is_truth"])
    truth_regions = bool(analyses.loc[eval_id, "eval_truth_regions"])
    query_regions = bool(analyses.loc[eval_id, "eval_target_regions"])

    inputs = {
        "query": draft_bench_vcf if draft_is_query else comp_vcf,
        "query_vcfidx": draft_bench_vcfidx if draft_is_query else comp_vcfidx,
        "truth": comp_vcf if draft_is_query else draft_bench_vcf,
        "truth_vcfidx": comp_vcfidx if draft_is_query else draft_bench_vcfidx,
        "genome": f"resources/references/{ref_id}.fa",
        "genome_index": f"resources/references/{ref_id}.fa.fai",
    }

    ## Defining Evaluation Regions
    if truth_regions:
        if query_regions:
            inputs["eval_regions"] = (
                "results/evaluations/truvari/{eval_id}_{bench_id}/eval_regions.bed"
            )
        else:
            inputs["eval_regions"] = comp_bed if draft_is_query else draft_bench_bed
    elif query_regions:
        inputs["eval_regions"] = draft_bench_bed if draft_is_query else comp_bed
    else:
        inputs["eval_regions"] = "NOT DEFINED ERROR"

    return inputs


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


## Benchmark VCF generation
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
    vcf_suffix = subset_df.iloc[0]["bench_vcf_processing"]

    if vcf_suffix == "none":
        return f"results/asm_varcalls/{vc_id}/annotations/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.vcf.gz"
    else:
        return f"results/asm_varcalls/{vc_id}/annotations/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.{vcf_suffix}.vcf.gz"


################################################################################
# load config


configfile: workflow.source_path("../config/resources.yml")


validate(config, "../schema/resources-schema.yml")

asm_config = config["assemblies"]
comp_config = config["comparisons"]
ref_config = config["references"]

################################################################################
# init analyses

## Loading analysis table with run information
analyses = load_analyses(
    workflow.source_path(f"../{config['analyses']}"), "../schema/analyses-schema.yml"
)


## Generating seperate tables for individual framework components
## asm variant calls
vc_params, vc_tbl = analyses_to_vc_tbl(analyses)

## draft benchmark set generation
bench_params, bench_tbl, bench_excluded_tbl = analyses_to_bench_tbls(analyses)

## Setting index for analysis run lookup
analyses = analyses.set_index("eval_id")

################################################################################
# init wildcard constraints

## Wildcard variables and ids

## Variables for assembly based variant calling
VCIDS = set(vc_tbl.index.tolist())
REFIDS = set(vc_tbl["ref"].tolist())
ASMIDS = set(vc_tbl["asm_id"].tolist())
VCCMDS = set(vc_tbl["vc_cmd"].tolist())
VCPARAMIDS = set(vc_tbl["vc_param_id"].tolist())

## Draft benchmark set generation variables
BENCHIDS = set(bench_tbl.index.tolist())
BENCHTYPS = set(bench_tbl["bench_type"].tolist())


## Evaluations
EVALIDS = set(analyses.index.tolist())
EVALCOMPIDS = set(analyses["eval_comp_id"].tolist())


# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_id="|".join(ASMIDS),
    comp_id="|".join(EVALCOMPIDS),
    ref_id="|".join(REFIDS),
    bench_id="|".join(BENCHIDS),
    bench_type="|".join(BENCHTYPS),
    eval_id="|".join(EVALIDS),
    vc_id="|".join(VCIDS),
    vc_cmd="|".join(VCCMDS),
    vc_param_id="|".join(VCPARAMIDS),


## Using zip in rule all to get config sets by config table rows

# defining variables for cleaner rule all
happy_analyses = analyses[analyses["eval_cmd"] == "happy"]
truvari_analyses = analyses[analyses["eval_cmd"] == "truvari"]
truvari_refine_analyses = analyses[analyses["eval_cmd"] == "truvari_refine"]
dipcall_tbl = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]
