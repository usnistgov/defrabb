## hap.py / Truvari evaluation input helpers
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


def get_happy_inputs_for_rule(wildcards):
    return get_happy_inputs_inner(wildcards.ref_id, wildcards.eval_id, analyses, config)


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


def get_eval_beds_for_rule(wildcards):
    return get_eval_beds(analyses, wildcards)


def get_truvari_inputs(analyses, config, wildcards):
    return get_truvari_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


def get_truvari_inputs_for_rule(wildcards):
    return get_truvari_inputs_inner(
        wildcards.ref_id,
        wildcards.eval_id,
        analyses,
        config,
    )


# scripts/ is on sys.path via rules/common.smk
from truvari_params import build_truvari_bench_params


def get_truvari_bench_params(wildcards):
    """truvari bench CLI flags for this analysis, from its eval_params profile."""
    profile = analyses.loc[wildcards.eval_id, "eval_params"]
    return build_truvari_bench_params(profile, config)


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
