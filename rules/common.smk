import pandas as pd
from pathlib import Path
import pprint
from warnings import WarningMessage
from snakemake.common import parse_uri
from snakemake.utils import validate
from snakemake.remote import AUTO
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider(keep_local=True)
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
# init analyses

## Loading analysis table with run information
analyses = load_analyses(
    workflow.source_path(f"../{config['analyses']}"), "../schema/analyses-schema.yml"
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
    ref_id = wildcards.get("ref", "")
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


## Helper functions for downloading and validating remote files
def get_remote_key(remote, key):
    try:
        return remote[key]
    except KeyError:
        raise ValueError(f"No {key} provided for external source {remote}")


## ~~~~~ reference genomes
def get_ref_config(ref_id):
    return ref_config[ref_id]["ref"]


def get_ref_uri(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "uri")


def get_ref_checksum(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum")


def get_ref_checksum_algo(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum_algo")


## ~~~~~ assemblies
def get_asm_config(asm_id, haplotype):
    return asm_config[asm_id][haplotype]


def get_asm_uri(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "uri")


def get_asm_checksum(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "checksum")


def get_asm_checksum_algo(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "checksum_algo")


## ~~~~~ stratifications
def get_strats_config(ref_id):
    return config["references"][ref_id]["stratifications"]


def get_strats_uri(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "uri")


def get_strats_checksum(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum")


def get_strats_checksum_algo(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum_algo")


## ~~~~~ comparison callsets
def get_comp_config(ref_id, comp_id, comp_ext):
    if comp_ext == "vcf.gz":
        comp_ext = "vcf"
    if comp_ext == "bed.gz":
        comp_ext = "bed"
    return comp_config[ref_id][comp_id][comp_ext]


def get_comp_uri(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "uri")


def get_comp_checksum(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "checksum")


def get_comp_checksum_algo(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "checksum_algo")


def get_comp_outfmt(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)

    comp_ext = wildcards.comp_ext
    if comp_ext == "vcf.gz" or comp_ext == "bed.gz":
        return "bgzip"
    elif comp_ext == "vcf" or comp_ext == "bed" or comp_ext == "fa":
        return "decompressed"
    else:
        return "comp outfmt not found"


## ~~~~~ bed files used for exclusion
def get_exclusion_config(ref_id, genomic_region):
    return ref_config[ref_id]["exclusions"][genomic_region]


def get_exclusion_uri(wildcards):
    remote = get_exclusion_config(wildcards.ref_id, wildcards.genomic_region)
    return get_remote_key(remote, "uri")


def get_exclusion_checksum(wildcards):
    remote = get_exclusion_config(wildcards.ref_id, wildcards.genomic_region)
    return get_remote_key(remote, "checksum")


def get_exclusion_checksum_algo(wildcards):
    remote = get_exclusion_config(wildcards.ref_id, wildcards.genomic_region)
    return get_remote_key(remote, "checksum_algo")


## Helper functions for other resource files
def get_genome_file(wildcards):
    ## Getting genome path from wildcard.ref_id
    ref_id = get_ref_id(wildcards)
    return workflow.source_path(f"../resources/references/{ref_id}.genome")


def get_male_bed(wildcards):
    filepath = f"{config['_par_bed_root']}/{ref_config[wildcards.ref]['par_bed']}"
    return filepath


def get_dipcall_par_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    par_path = get_male_bed(wildcards)
    return f"-x {str(par_path)}" if is_male else ""


## Evaluation Inputs and Parameters (Generic method for happy and truvari)
def get_eval_inputs(analyses, config, bench_space, wildcards):
    return get_eval_inputs_inner(
        eval_cmd=wildcards.eval_cmd,
        ref_id=wildcards.ref,
        eval_query=wildcards.eval_query,
        eval_truth=wildcards.eval_truth,
        eval_truth_regions=wildcards.eval_truth_regions,
        eval_target_regions=wildcards.eval_target_regions,
        exclusion_set=wildcards.bench_exclusion_set,
        bench_bed_processing=wildcards.bench_bed_processing,
        bench_space=bench_space,
        config=config,
    )


def get_eval_inputs_inner(
    eval_cmd,
    ref_id,
    eval_query,
    eval_truth,
    eval_truth_regions,
    eval_target_regions,
    exclusion_set,
    bench_bed_processing,
    bench_space,
    config,
):
    ## Creating empty dictionary for storing inputs
    inputs = {}
    ## Reference genome and stratifications
    inputs["genome"] = f"resources/references/{ref_id}.fa"
    inputs["genome_index"] = f"resources/references/{ref_id}.fa.fai"
    if eval_cmd == "happy":
        ## Only defining stratifications for happy
        strat_tb = config["references"][ref_id]["stratifications"]["tarball"]
        inputs["strat_tb"] = f"resources/strats/{ref_id}/{strat_tb}"

    ## Getting input vcf, vcf idx, and bed files for query and  truth
    def get_comp_files(eval_id, ref_id, eval_regions, exclusion_set):
        comp_file_dict = {}
        ## Getting VCFs
        ## ~~~~~ using draft bench vcf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if eval_id == "this_row":
            comp_file_dict[
                "vcf"
            ] = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.processed.vcf.gz"
        ## ~~~~~ using comparison callset defined in resources.yml ~~~~~~~~~~~~~~~~~~~~~~~~~~
        elif eval_id in config["comparisons"][ref_id]:
            ## defining vcf paths for comparisons listed in resources
            comp_file_dict[
                "vcf"
            ] = f"resources/comparison_variant_callsets/{ref_id}_{eval_id}.vcf.gz"
        else:
            print(
                "eval_query in analyses table should be this_row, a comparison id in resources,",
                "or reference another benchmark (WIP)",
            )

        ## Getting evaluation regions paths
        ## ~~~~~~~~~~ target and truth region schema values should be true, false, or a bed file

        ## Using predefined bed file
        if eval_regions.lower().endswith("bed"):
            tr_dir = "resources/manual/target_regions"
            comp_file_dict["bed"] = f"{tr_dir}/{eval_regions}"

        elif eval_regions == "yes":
            ### using draft benchmark regions
            if eval_id == "this_row":
                if bench_bed_processing == "manual":
                    comp_file_dict[
                        "bed"
                    ] = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.manual.bed"
                elif bench_bed_processing == "none":
                    comp_file_dict[
                        "bed"
                    ] = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.bed"
                elif bench_bed_processing == "exclude":
                    comp_file_dict[
                        "bed"
                    ] = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.excluded.bed"
                else:
                    print(
                        "bench_bed_processing should be either none, manual, or exclude"
                    )
                    comp_file_dict[
                        "bed"
                    ] = f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.FIXME.bed"
            ### using comparison regions
            else:
                comp_file_dict[
                    "bed"
                ] = f"resources/comparison_variant_callsets/{ref_id}_{eval_id}.bed"
        elif eval_regions != "none":
            print(
                f"regions not properly defined {eval_regions}, should be yes, none, or a bed file"
            )

        return comp_file_dict

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~ QUERY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    query_inputs = get_comp_files(
        eval_query, ref_id, eval_target_regions, exclusion_set
    )
    inputs["query"] = query_inputs["vcf"]
    inputs["query_vcfidx"] = f"{query_inputs['vcf']}.tbi"
    if "bed" in query_inputs.keys():
        inputs["target_regions"] = query_inputs["bed"]

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~ TRUTH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    truth_inputs = get_comp_files(eval_truth, ref_id, eval_truth_regions, exclusion_set)

    inputs["truth"] = truth_inputs["vcf"]
    inputs["truth_vcfidx"] = f"{truth_inputs['vcf']}.tbi"
    if "bed" in truth_inputs.keys():
        inputs["truth_regions"] = truth_inputs["bed"]

    ## Returning evaluation inputs
    return inputs


## Exclusions
def get_exclusion_inputs(wildcards):

    ## Getting list of excluded regions
    exclusion_set = config["exclusion_set"][wildcards.bench_exclusion_set]

    ## Initiating empty list for storing paths for beds to excluded from
    ## diploid assembled regions
    exc_paths = []
    for exclusion in exclusion_set:
        ## Determining path for asm specific exclusions and asm agnostic exclusions
        if exclusion in config["exclusion_asm_agnostic"]:
            ref_id = get_ref_id(wildcards)
            exc_path = f"resources/exclusions/{ref_id}/{exclusion}"
        else:
            exc_path = f"results/draft_benchmarksets/intermediates/exclusions/{bench_space.wildcard_pattern}"
        ## Adding slop - currently a 15kb hard coded buffer around excluded repeat regions
        if exclusion in config["exclusion_slop_regions"]:
            exc_path = f"{exc_path}_slop"
        elif exclusion in config["exclusion_asm_intersect"]:
            ## Ensuring bed files are sorted before intersect
            exc_path = f"{exc_path}.sorted"

        ## Defining which regions are excluded based on diploid assembly breaks
        if exclusion in config["exclusion_asm_intersect"]:
            exc_paths += [f"{exc_path}_start", f"{exc_path}_end"]
        else:
            exc_paths = exc_paths + [exc_path]

    ## Adding to exc_paths list and ensuring all beds are sorted
    ## prior to exclusion from dip assembled regions
    exclusion_paths = [f"{exc}.sorted.bed" for exc in exc_paths]

    ## Returning list of bed paths for exclusion
    return exclusion_paths


## Benchmark VCF generation


def get_processed_vcf(wildcards):
    vcf_suffix = wildcards.bench_vcf_processing
    if vcf_suffix == "none":
        return f"results/draft_benchmarksets/intermediates/{bench_space.wildcard_pattern}.vcf.gz"
    else:
        return f"results/draft_benchmarksets/intermediates/{bench_space.wildcard_pattern}.{vcf_suffix}.vcf.gz"
