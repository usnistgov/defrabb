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
# init resources


configfile: workflow.source_path("../config/resources.yml")


validate(config, "../schema/resources-schema.yml")

asm_config = config["assemblies"]
comp_config = config["comparisons"]
ref_config = config["references"]

################################################################################
# init analyses

## Loading analysis table with run information
analyses = load_analyses(
    workflow.source_path(f"{config['analyses']}"), "../schema/analyses-schema.yml"
)

vc_params, vc_tbl = analyses_to_vc_tbl(analyses)

ASMIDS = set(vc_tbl["asm_id"])
REFIDS = set(vc_tbl["ref"])

## Wildcard variables and ids

## Variables for assembly based variant calling
REFIDS = set(vc_tbl["ref"].tolist())
ASMIDS = set(vc_tbl["asm_id"].tolist())
VCCMDS = set(vc_tbl["vc_cmd"].tolist())
BENCHVCFPROC = set(analyses["bench_vcf_processing"])
BENCHBEDPROC = set(analyses["bench_bed_processing"])
BENCHEXCLUSIONSET = set(analyses["bench_exclusion_set"])
COMPIDS = set(analyses["eval_query"].tolist() + analyses["eval_truth"].tolist())
EVALIDS = set(
    analyses["eval_query"].tolist() + analyses["eval_truth"].tolist() + ["this_row"]
)
EVALTRUTHREGIONS = set(analyses["eval_truth_regions"].tolist())
EVALTARGETREGIONS = set(analyses["eval_target_regions"].tolist())
GENOMICREGIONS=["hifi-pacbioDV-XY-discrep", "imperfecthomopol-gt30","segdups","tandem-repeats","satellites","gaps","self-chains","flanks","svs-and-simple-repeats"]
# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).


wildcard_constraints:
    asm_id="|".join(ASMIDS),
    ref_id="|".join(REFIDS),
    bench_vcf_processing="|".join(BENCHVCFPROC),
    bench_bed_processing="|".join(BENCHBEDPROC),
    bench_exclusion_set="|".join(BENCHEXCLUSIONSET),
    comp_dir="asm_varcalls|draft_benchmarksets|evaluations|report",
    comp_id="|".join(COMPIDS),
    comp_ext="vcf.gz|vcf|bed|bed.gz",
    eval_truth="|".join(EVALIDS),
    eval_query="|".join(EVALIDS),
    eval_truth_regions="|".join(EVALTRUTHREGIONS),
    eval_target_regions="|".join(EVALTARGETREGIONS),
    genomic_region="|".join(GENOMICREGIONS)


## Using Paramspace for file paths
## - preparing tables for output naming
run_tbl = analyses
run_tbl["run_id"] = [f"r{idx}" for idx in run_tbl.index]

bench_cols = [
    "ref",
    "asm_id",
    "vc_cmd",
    "vc_param_id",
    "bench_type",
    "bench_vcf_processing",
    "bench_bed_processing",
    "bench_exclusion_set",
]

bench_tbl = analyses[bench_cols].drop_duplicates()

bench_tbl["bench_id"] = [f"b{idx}" for idx in bench_tbl.index]
## Creating analyses table with run and bench ids.
analyses_wids = run_tbl.merge(bench_tbl)

happy_tbl = analyses_wids[analyses_wids["eval_cmd"] == "happy"]

happy_space = Paramspace(
    happy_tbl,
    filename_params=[
        "run_id",
        "bench_id",
        "asm_id",
        "ref",
    ],
)

truvari_tbl = analyses_wids[analyses_wids["eval_cmd"] == "truvari"]
truvari_space = Paramspace(
    truvari_tbl,
    filename_params=[
        "run_id",
        "bench_id",
        "asm_id",
        "ref",
    ],
)

dipcall_space = Paramspace(
    analyses.loc[
        analyses["vc_cmd"] == "dipcall",
        ["asm_id", "ref", "vc_cmd", "vc_param_id"],
    ].drop_duplicates(),
    filename_params=["asm_id", "ref", "vc_cmd"],
)


bench_space = Paramspace(
    bench_tbl, filename_params=["bench_id", "asm_id", "ref", "vc_cmd", "bench_type"]
)

excluded_bench_space = Paramspace(
    bench_tbl.loc[
        bench_tbl["bench_exclusion_set"] != "none",
    ],
    filename_params=["bench_id", "asm_id", "ref", "vc_cmd", "bench_type"],
)

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
    return workflow.source_path(f"../resources/references/{ref_id}.fa")


def get_genome_file(wildcards):
    ref_id = get_ref_id(wildcards)
    return workflow.source_path(f"../resources/references/{ref_id}.genome")


def get_male_bed(wildcards):
    root = config["_par_bed_root"]
    filename = ref_config[wildcards.ref_id]["par_bed"]
    return workflow.source_path(f"../{root}/{filename}")


def get_dipcall_par_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    par_path = get_male_bed(wildcards)
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
    ## Creating empty dictionary for storing inputs
    inputs = {}

    ## Reference genome and stratifications
    inputs["genome"] = f"resources/references/{ref_id}.fa"
    inputs["genome_index"] = f"resources/references/{ref_id}.fa.fai"
    strat_tb = config["references"][ref_id]["stratifications"]["tarball"]
    inputs["strat_tb"] = f"resources/strats/{ref_id}/{strat_tb}"

    ## draft benchmark variant calls
    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz"
    draft_bench_vcfidx = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi"

    ## draft benchmark regions
    if analyses.loc[eval_id, "exclusion_set"] == "none":
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed"
    else:
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.benchmark.bed"

    ## comparison variant call paths
    comp_vcf = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz"
    comp_vcfidx = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi"
    comp_bed = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    ## Determining which callsets and regions are used as truth
    if analyses.loc[eval_id, "eval_comp_id_is_truth"] == True:
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
    trs = analyses.loc[eval_id, "eval_target_regions"]
    if trs.lower() != "false":
        if trs.lower() == "true":
            if query == "draft_bench":
                inputs["target_regions"] = draft_bench_bed
            else:
                inputs["target_regions"] = comp_bed
        else:
            tr_dir = "resources/manual/target_regions"
            inputs["target_regions"] = workflow.source_path(f"../{tr_dir}/{trs}")

    ## Returning happy inputs
    return inputs


def get_truvari_inputs(analyses, config, wildcards):
    return get_truvari_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


def get_truvari_inputs_inner(ref_id, eval_id, analyses, config):
    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz"
    draft_bench_vcfidx = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi"

    if analyses.loc[eval_id, "exclusion_set"] == "none":
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed"
    else:
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed"

    comp_vcf = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz"
    comp_vcfidx = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi"
    comp_bed = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    draft_is_query = analyses.loc[eval_id, "eval_comp_id_is_truth"] == True

    return {
        "query": draft_bench_vcf if draft_is_query else comp_vcf,
        "query_vcfidx": draft_bench_vcfidx if draft_is_query else comp_vcfidx,
        "truth": comp_vcf if draft_is_query else draft_bench_vcf,
        "truth_vcfidx": comp_vcfidx if draft_is_query else draft_bench_vcfidx,
        "truth_regions": comp_bed if draft_is_query else draft_bench_bed,
        "genome": f"resources/references/{ref_id}.fa",
        "genome_index": f"resources/references/{ref_id}.fa.fai",
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


## Benchmark VCF generation


def get_processed_vcf(wildcards):
    vcf_suffix = bench_tbl.loc[wildcards.bench_id, "bench_vcf_processing"]
    if vcf_suffix == "none":
        return f"results/draft_benchmarksets/{{bench_id}}/intermediates/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.vcf.gz"
    else:
        return f"results/draft_benchmarksets/{{bench_id}}/intermediates/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.{vcf_suffix}.vcf.gz"
