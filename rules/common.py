import pandas as pd
from pathlib import Path
from snakemake.utils import validate

################################################################################
## Config processing functions


def load_df(path, schema):
    df = pd.read_table(path)
    validate(df, schema)
    return df


def load_analyses(path, schema):
    return load_df(path, schema).astype(dtype={"eval_target_regions": str})


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
## Utility Functions

## Loading and validating analysis tables
def get_analyses(path):
    # target_regions must be a string even though it might only contain
    # 'boolean' values
    analyses = pd.read_table(path, dtype={"eval_target_regions": str})
    validate(analyses, "schema/analyses-schema.yml")

    try:
        return analyses.set_index("eval_id", verify_integrity=True)
    except ValueError:
        print("All keys in column 'eval_id' must by unique")


## Generating concatenated string for wildcard constraints
def format_constraint(xs):
    return "|".join(set(xs))


################################################################################
## Rule parameters


def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    root = config["_par_bed_root"]
    par_path = Path(root) / ref_config[wildcards.ref_id]["par_bed"]
    return f"-x {str(par_path)}" if is_male else ""


## Happy Inputs and Parameters
def get_happy_inputs(analyses, config, wildcards):
    return get_happy_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


def get_happy_inputs_inner(ref_id, eval_id, analyses, config):
    _analyses = analyses.set_index("eval_id")
    ref_id = _analyses.loc[eval_id, "ref"]
    strat_tb = config["stratifications"][ref_id]["tarball"]

    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz"
    draft_bench_bed = (
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed"
        if _analyses.loc[(eval_id, "exclusion_set")] == "none"
        else "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed"
    )

    comp_vcf = "resources/comparison_variant_callsets/{comp_id}.vcf.gz"
    comp_bed = "resources/comparison_variant_callsets/{comp_id}.bed"

    comp_is_truth = _analyses.loc[eval_id, "eval_comp_id_is_truth"] == True

    def get_target_regions():
        trs = _analyses.loc[eval_id, "eval_target_regions"]
        if trs.lower() != "false":
            if trs.lower() == "true":
                return draft_bench_bed if comp_is_truth else comp_bed
            else:
                return f"resources/manual/target_regions/{trs}"
        else:
            return None

    return {
        "genome": f"resources/references/{ref_id}.fa",
        "genome_index": f"resources/references/{ref_id}.fa.fai",
        "strat_tb": f"resources/strats/{ref_id}/{strat_tb}",
        "comp_idx": "resources/comparison_variant_callsets/{comp_id}.vcf.gz.tbi",
        "query": draft_bench_vcf if comp_is_truth else comp_vcf,
        "truth": comp_vcf if comp_is_truth else draft_bench_vcf,
        "truth_regions": comp_bed if comp_is_truth else draft_bench_bed,
        "target_regions": get_target_regions(),
    }


## Exclusions
def get_exclusion_inputs(wildcards):

    ## Getting list of excluded regions
    exclusion_set_id = bench_tbl.loc[wildcards.bench_id, "exclusion_set"]
    exclusion_set = config["exclusion_set"][exclusion_set_id]

    ## Initiating empty list for storing paths for beds to excluded from
    ## diploid assembled regions
    exc_paths = []
    for exclusion in exclusion_set:

        ## Determining path for asm specific exclusions and asm agnostic exclusions
        if exclusion in config["exclusion_asm_agnostic"]:
            exc_path = f"resources/exclusions/{{ref_id}}/{exclusion}"
        else:
            exc_path = f"results/draft_benchmarksets/{{bench_id}}/exclusions/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}_{exclusion}"

        ## Adding slop - currently a 15kb hard coded buffer around excluded repeat regions
        if exclusion in config["exclusion_slop_regions"]:
            exc_path = f"{exc_path}_slop"
        elif exclusion in config["exclusion_asm_intersect"]:
            ## Ensuring bed files are sorted before intersect
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
