import pandas as pd
from pathlib import Path
from snakemake.utils import validate

################################################################################
## Config processing functions


def load_df(path, schema):
    df = pd.read_table(path, dtype={"eval_target_regions": str})
    validate(df, schema)
    return df


def load_analyses(path, schema):
    return load_df(path, schema).astype(dtype={"eval_target_regions": str})


def _filter_subtable(df, filter_re, id_cols):
    params = df.filter(regex=filter_re)
    tbl = params.join(df[id_cols]).drop_duplicates()
    return (params, tbl)


def analyses_to_vc_tbl(analyses):
    return _filter_subtable(analyses, "vc_", ["asm_id", "ref"])


def analyses_to_bench_tbls(analyses):
    id_cols = [
        "asm_id",
        "vc_cmd",
        "vc_param_id",
        "ref",
        "bench_exclusion_set",
    ]
    params, tbl = _filter_subtable(analyses, "bench_", id_cols)
    excluded_tbl = tbl[tbl.bench_exclusion_set != "none"]
    return (params, tbl, excluded_tbl)


################################################################################
## Rule parameters
def get_genome_file(wildcards):
    ## Getting genome path from wildcard.ref_id
    ref_id = wildcards.get("ref", "")
    if ref_id:
        return workflow.source_path(f"../resources/references/{ref_id}.genome")
    ## get ref_id from prefix
    prefix = wildcards.get("prefix", "")
    for id in REFIDS:
        if id in prefix:
            return workflow.source_path(f"../resources/references/{id}.genome")
    print("ref_id not found in bed file prefix")


def get_male_bed(wildcards):
    root = config["_par_bed_root"]
    filename = ref_config[wildcards.ref]["par_bed"]
    return workflow.source_path(f"../{root}/{filename}")


def get_dipcall_par_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    par_path = get_male_bed(wildcards)
    return f"-x {str(par_path)}" if is_male else ""


## Evaluation Inputs and Parameters (Generic method for happy and truvari)
def get_eval_inputs(analyses, config, wildcards):
    return get_eval_inputs_inner(
        eval_cmd =  wildcards.eval_cmd,
        ref_id = wildcards.ref,
        eval_query = wildcards.eval_query,
        eval_truth = wildcards.eval_truth,
        eval_truth_regions = wildcards.eval_truth_regions,
        eval_target_regions = wildcards.eval_target_regions,
        exclusion_set = wildcards.bench_exclusion_set,
        bench_space = bench_space,
        config = config,
    )

def get_eval_inputs_inner(
    eval_cmd,
    ref_id,
    eval_query,
    eval_truth,
    eval_truth_regions,
    eval_target_regions,
    exclusion_set,
    bench_space,
    config,
):
    print(eval_query)
    print("HERE")
    ## Creating empty dictionary for storing inputs
    inputs = {}
    # print(analyses["comparisons"][ref_id][eval_query])
    ## Reference genome and stratifications
    inputs["genome"] = f"resources/references/{ref_id}.fa"
    inputs["genome_index"] = f"resources/references/{ref_id}.fa.fai"
    if eval_cmd == "happy":
        ## Only defining stratifications for happy
        strat_tb = config["references"][ref_id]["stratifications"]["tarball"]
        inputs["strat_tb"] = f"resources/strats/{ref_id}/{strat_tb}"

    ## Getting input vcf, vcf idx, and bed files for query and  truth
    def get_comp_files(eval_id, ref_id, eval_regions, exclusion_set):
        ## Getting VCFs
        ## ~~~~~ using draft bench vcf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if eval_id == "this_row":
            vcf = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.vcf.gz"
            vcfidx = (
                f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.vcf.gz.tbi"
            )
        ## ~~~~~ using comparison callset defined in resources.yml ~~~~~~~~~~~~~~~~~~~~~~~~~~
        else:
            try:
                eval_id in analyses["comparisons"][ref_id][eval_query]
            except:
                print(
                    "eval_query in analyses table should be this_row, a comparison id in resources,",
                    "or reference another benchmark (WIP)",
                )

            ## defining vcf paths for comparisons listed in resources
            vcf = f"resources/comparison_variant_callsets/{ref_id}_{eval_id}.vcf.gz"
            vcfidx = (
                f"resources/comparison_variant_callsets/{ref_id}_{eval_id}.vcf.gz.tbi"
            )

        ## Getting evaluation regions paths
        ## ~~~~~~~~~~ target and truth region schema values should be true, false, or a bed file

        ## Using predefined bed file
        if eval_regions.lower().endswith("bed"):
            tr_dir = "resources/manual/target_regions"
            bed = workflow.source_path(f"../{tr_dir}/{eval_regions}")
        elif eval_regions == "yes":
            ### using draft benchmark regions
            if eval_id == "this_row":
                if exclusion_set == "none":
                    bed = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.bed"
                else:
                    bed = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.excluded.bed"
            ### using comparison regions
            else:
                bed = f"resources/comparison_variant_callsets/{ref_id}_{eval_id}.bed"
        elif eval_regions == "none":
            bed = ""
        else:
            print(
                f"regions not properly defined {eval_regions}, should be true, false, or a bed file"
            )
            
        return {"vcf": vcf, "vcfidx": vcfidx, "bed": bed}

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~ QUERY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    query_inputs = get_comp_files(eval_query, ref_id, eval_target_regions, exclusion_set)
    inputs["query"] = query_inputs["vcf"]
    inputs["query_vcfidx"] = query_inputs["vcfidx"]
    inputs["target_regions"] = query_inputs["bed"]

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~ TRUTH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    truth_inputs = get_comp_files(eval_truth, ref_id, eval_truth_regions, exclusion_set)

    inputs["truth"] = truth_inputs["vcf"]
    inputs["truth_vcfidx"] = truth_inputs["vcfidx"]
    inputs["truth_regions"] = truth_inputs["bed"]

    ## Returning evaluation inputs
    return inputs


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
            exc_path = f"resources/exclusions/{ref_id}/{exclusion}"
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


## Benchmark VCF generation


def get_processed_vcf(wildcards):
    vcf_suffix = bench_tbl.loc[wildcards.bench_id, "bench_vcf_processing"]
    if vcf_suffix == "none":
        return f"results/draft_benchmarksets/{bench_space.wildcard_params}.raw.vcf.gz"
    else:
        return f"results/draft_benchmarksets/{bench_space.wildcard_params}.{vcf_suffix}.vcf.gz"
