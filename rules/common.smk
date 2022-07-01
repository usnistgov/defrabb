import pandas as pd
from pathlib import Path
from warnings import WarningMessage
from snakemake.common import parse_uri
from snakemake.utils import validate
from snakemake.remote import AUTO
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(keep_local = True)

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
## Rule parameters
# def get_remotes(wildcards):
#     uri="ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
#     ## Downloading file
#     parsed_uri = parse_uri(uri)
#     print
#     scheme = parsed_uri.scheme

#     if scheme == "s3":
#         return(S3.remote(uri))
#     elif scheme in ["https","ftp"]:
#         return(AUTO.remote(uri))
#     else:
#         WarningMessage(f"Unexpected remote file protocol {scheme}, expected values are s3, https, and ftp")
def get_remote_key(remote, key):
    try:
        return remote[key]
    except KeyError:
        raise ValueError(
            f"No {key} provided for external source {remote}"
        )

## Helper functions for downloading and validating reference genomes
def get_ref_config(ref_id):
    return(ref_config[ref_id]["ref"])

def get_ref_uri(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "uri")

def get_ref_checksum(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum")

def get_ref_checksum_algo(wildcards):
    remote = get_ref_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum_algo")

## Helper functions for downloading and validating input assemblies
def get_asm_config(asm_id, haplotype):
    return(asm_config[asm_id][haplotype])

def get_asm_uri(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "uri")

def get_asm_checksum(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "checksum")

def get_asm_checksum_algo(wildcards):
    remote = get_asm_config(wildcards.asm_id, wildcards.haplotype)
    return get_remote_key(remote, "checksum_algo")

## Helper functions for downloading and validating stratifications
def get_strats_config(ref_id):
    return(config['references'][ref_id]['stratifications'])

def get_strats_uri(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "uri")

def get_strats_checksum(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum")

def get_strats_checksum_algo(wildcards):
    remote = get_strats_config(wildcards.ref_id)
    return get_remote_key(remote, "checksum_algo")

## Helper functions for downloading and validating comparison callsets
def get_comp_config(ref_id, comp_id, comp_ext):
    if comp_ext == "vcf.gz":
        comp_ext = "vcf"
    if comp_ext == "bed.gz":
        comp_ext = "bed"
    return(comp_config[ref_id][comp_id][comp_ext])

def get_comp_uri(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "uri")

def get_comp_checksum(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "checksum")

def get_comp_checksum_algo(wildcards):
    remote = get_comp_config(wildcards.ref_id, wildcards.comp_id, wildcards.comp_ext)
    return get_remote_key(remote, "checksum_algo")

## Helper functions for other resource files
def get_genome_file(wildcards):
    ## Getting genome path from wildcard.ref_id
    ref_id = wildcards.get("ref_id", "")
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
