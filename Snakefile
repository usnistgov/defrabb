import pandas as pd
from pathlib import Path
from snakemake.utils import min_version, validate


include: "rules/common.smk"
include: "rules/exclusions.smk"
include: "rules/report.smk"
include: "rules/bench_vcf_processing.smk"


# include: "rules/bench_vcf_processing.smk"


min_version("7.3.0")


## Rule ordering for ambiguous rules
ruleorder: download_bed_gz > sort_bed


################################################################################
# init resources


configfile: "config/resources.yml"


validate(config, "schema/resources-schema.yml")

asm_config = config["assemblies"]
comp_config = config["comparisons"]
ref_config = config["references"]

################################################################################
# init analyses
## TODO add checks for setting indecies - maybe move to function
analyses = analyses = pd.read_table(
    config["analyses"], dtype={"eval_target_regions": str}
)
validate(analyses, "schema/analyses-schema.yml")

## Generating seperate tables for individual framework components
## asm variant calls
vc_params = analyses.filter(regex="vc_").drop_duplicates()
vc_ids = analyses[["vc_id", "asm_id", "ref"]].drop_duplicates()
vc_tbl = pd.merge(vc_ids, vc_params, how="inner", on="vc_id").set_index("vc_id")

## draft benchmark set generation
bench_params = analyses.filter(regex="bench_").drop_duplicates()
bench_ids = analyses[
    ["bench_id", "asm_id", "vc_id", "vc_cmd", "vc_param_id", "ref", "exclusion_set"]
].drop_duplicates()
bench_tbl = pd.merge(bench_ids, bench_params, how="inner", on="bench_id").set_index(
    "bench_id"
)
bench_excluded_tbl = bench_tbl[bench_tbl.exclusion_set != "none"]

## Evaluation Runs
eval_params = analyses.filter(regex="eval_").drop_duplicates()
eval_ids = analyses[["eval_id", "bench_id", "ref"]].drop_duplicates()
eval_tbl = pd.merge(eval_ids, eval_params, how="inner", on="eval_id").set_index(
    "eval_id"
)

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


## Evaluations
EVALIDS = set(eval_tbl.index.tolist())
EVALCOMPIDS = set(eval_tbl["eval_comp_id"].tolist())


# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_id="|".join(ASMIDS),
    comp_id="|".join(EVALCOMPIDS),
    ref_id="|".join(REFIDS),
    bench_id="|".join(BENCHIDS),
    eval_id="|".join(EVALIDS),
    vc_id="|".join(VCIDS),
    vc_param_id="|".join(VCPARAMIDS),


################################################################################
# main rule
#


## Rules to run locally
localrules:
    get_ref,
    get_assemblies,
    get_comparison_vcf,
    get_comparison_bed,
    get_strats,
    download_bed_gz,
    get_SVs_from_vcf,
    subtract_exclusions,
    add_flanks,
    intersect_start_and_end,
    intersect_SVs_and_homopolymers,
    get_SVs_from_vcf,
    postprocess_bed,
    sort_bed,
    add_slop,
    tabix,
    move_asm_vcf_to_draft_bench,


## Snakemake Report
report: "report/workflow.rst"


## Using zip in rule all to get config sets by config table rows

# defining variables for cleaner rule all
happy_analyses = analyses[analyses["eval_cmd"] == "happy"]
dipcall_tbl = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]


rule all:
    input:
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        ## Bench VCF Processing
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.hap1.bam.bai",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.hap2.bam.bai",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        ## rules for report
        expand(
            "results/report/assemblies/{asm_id}_{haplotype}_stats.txt",
            asm_id=ASMIDS,
            haplotype=["maternal", "paternal"],
        ),
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip_bcftools_stats.txt",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip_rtg_stats.txt",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.summary.csv",
            zip,
            eval_id=happy_analyses.index.tolist(),
            bench_id=happy_analyses["bench_id"].tolist(),
            ref_id=happy_analyses["ref"].tolist(),
            comp_id=happy_analyses["eval_comp_id"].tolist(),
            asm_id=happy_analyses["asm_id"].tolist(),
            vc_cmd=happy_analyses["vc_cmd"].tolist(),
            vc_param_id=happy_analyses["vc_param_id"].tolist(),
        ),
        expand(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.extended.csv",
        zip,
        eval_id=analyses[analyses["eval_cmd"] == "happy"].index.tolist(),
        bench_id=analyses[analyses["eval_cmd"] == "happy"]["bench_id"].tolist(),
        ref_id=analyses[analyses["eval_cmd"] == "happy"]["ref"].tolist(),
        comp_id=analyses[analyses["eval_cmd"] == "happy"]["eval_comp_id"].tolist(),
        asm_id=analyses[analyses["eval_cmd"] == "happy"]["asm_id"].tolist(),
        vc_cmd=analyses[analyses["eval_cmd"] == "happy"]["vc_cmd"].tolist(),
        vc_param_id=analyses[analyses["eval_cmd"] == "happy"][
        "vc_param_id"
            ].tolist(),
        ),


#       expand("results/bench/truvari/{tvi_bench}.extended.csv", tvi_bench = analyses[analyses["bench_cmd"] == "truvari"].index.tolist()), ## Not yet used

################################################################################
################################################################################
##
## Downloading Resource Data Files
##
################################################################################
################################################################################


# Get and prepare assemblies
rule get_assemblies:
    output:
        "resources/assemblies/{asm_id}/{haplotype}.fa",
    params:
        url=lambda wildcards: asm_config[wildcards.asm_id][wildcards.haplotype],
    log:
        "logs/get_assemblies/{asm_id}_{haplotype}.log",
    shell:
        "curl -f -L {params.url} 2> {log} | gunzip -c 1> {output} 2>> {log}"


# Get and prepare reference
rule get_ref:
    output:
        "resources/references/{ref_id}.fa",
    params:
        url=lambda wildcards: ref_config[wildcards.ref_id]["ref_url"],
    log:
        "logs/get_ref/{ref_id}.log",
    shell:
        "curl -f --connect-timeout 120 -L {params.url} 2> {log} | gunzip -c 1> {output} 2>> {log}"


rule index_ref:
    input:
        "resources/references/{ref_id}.fa",
    output:
        "resources/references/{ref_id}.fa.fai",
    log:
        "logs/index_ref/{ref_id}.log",
    resources:
        mem_mb=16000,
    wrapper:
        "0.79.0/bio/samtools/faidx"


rule index_ref_mmi:
    input:
        "resources/references/{ref_id}.fa",
    output:
        "resources/references/{ref_id}.mmi",
    log:
        "logs/index_ref_mmi/{ref_id}.log",
    threads: 4
    resources:
        mem_mb=2400,
    conda:
        "envs/dipcall.yml"
    shell:
        "minimap2 -x asm5 -d {output} {input}"


rule index_ref_sdf:
    input:
        "resources/references/{ref_id}.fa",
    output:
        directory("resources/references/{ref_id}.sdf"),
    log:
        "logs/index_ref_sdf/{ref_id}.log",
    conda:
        "envs/rtgtools.yml"
    shell:
        "rtg format -o {output} {input}"


################################################################################
# Get stratifications


rule get_strats:
    output:
        "resources/strats/{ref_id}/{strat_id}.tar.gz",
    params:
        url=lambda wildcards: f"{config['stratifications'][wildcards.ref_id]['url']}",
    log:
        "logs/get_strats/{ref_id}_{strat_id}.log",
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations


rule get_comparison_vcf:
    output:
        "resources/comparison_variant_callsets/{comp_id}.vcf.gz",
    params:
        url=lambda wildcards: comp_config[wildcards.comp_id]["vcf_url"],
    log:
        "logs/get_comparisons/{comp_id}_vcf.log",
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


use rule get_comparison_vcf as get_comparison_bed with:
    output:
        "resources/comparison_variant_callsets/{comp_id}.bed",
    params:
        url=lambda wildcards: comp_config[wildcards.comp_id]["bed_url"],
    log:
        "logs/get_comparisons/{comp_id}_bed.log",


## General indexing rule for vcfs
rule tabix:
    input:
        "{filename}.vcf.gz",
    output:
        "{filename}.vcf.gz.tbi",
    params:
        extra="-t",
    log:
        "logs/tabix/{filename}.log",
    wrapper:
        "v1.0.0/bio/bcftools/index"


################################################################################
################################################################################
##
## Component: Assembly Based Variant Calling
##
################################################################################
################################################################################

## Run Dipcall


rule run_dipcall:
    input:
        h1=ancient("resources/assemblies/{asm_id}/paternal.fa"),
        h2=ancient("resources/assemblies/{asm_id}/maternal.fa"),
        ref=ancient("resources/references/{ref_id}.fa"),
        ref_idx=ancient("resources/references/{ref_id}.fa.fai"),
        ref_mmi=ancient("resources/references/{ref_id}.mmi"),
    output:
        make="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.mak",
        vcf="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
        bed="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.bed",
        bam1="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.hap1.bam",
        bam2="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.hap2.bam",
    conda:
        "envs/dipcall.yml"
    params:
        prefix="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
        male_bed=get_male_bed,
        ts=config["_dipcall_threads"],
        make_jobs=config["_dipcall_jobs"],
        extra=lambda wildcards: ""
        if vc_tbl.loc[wildcards.vc_id]["vc_params"] == "default"
        else vc_tbl.loc[wildcards.vc_id]["vc_params"],
    log:
        "logs/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    resources:
        mem_mb=config["_dipcall_jobs"] * config["_dipcall_mem"],  ## GB per make job run in parallel
    threads: config["_dipcall_threads"] * config["_dipcall_jobs"]
    shell:
        """
        echo "Writing Makefile defining dipcall pipeline"
        run-dipcall \
            -t {params.ts} \
            -d {input.ref_mmi} \
            {params.extra} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            1> {output.make}

        echo "Running dipcall pipeline"
        make -j{params.ts} -f {output.make}
        """


rule sort_bed:
    input:
        in_file="{prefix}.bed",
        ## TODO remove hardcoding for genome file
        genome="resources/exclusions/GRCh38.genome",
    output:
        "{prefix}_sorted.bed",
    log:
        "logs/sort_bed/{prefix}.log",
    wrapper:
        "0.74.0/bio/bedtools/sort"


rule index_dip_bam:
    input:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam",
    output:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam.bai",
    log:
        "logs/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam.bai.log",
    wrapper:
        "v1.1.0/bio/samtools/index"


################################################################################
################################################################################
##
## Component: Generating Draft Benchmarkset from Assembly-Based Variant Calls
##
################################################################################
################################################################################

## Moved to rules/bench_vcf_processing
# rule postprocess_vcf:
#     input:
#         lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
#     output:
#         "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
#     log:
#         "logs/process_benchmark_vcf/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
#     shell:
#         "cp {input} {output} &> {log}"


rule postprocess_bed:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed",
    log:
        "logs/process_benchmark_bed/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    shell:
        "cp {input} {output} &> {log}"


################################################################################
################################################################################
##
## Component: Evaluating Draft Benchmarksets
##
################################################################################
################################################################################

## Run happy


rule run_happy:
    input:
        unpack(get_happy_inputs),
    output:
        multiext(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
            ".runinfo.json",
            ".vcf.gz",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
        report(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.summary.csv",
            caption="report/happy_summary.rst",
            category="Happy",
        ),
    params:
        prefix="results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
        strat_tsv=lambda wildcards: f"{wildcards.ref_id}/{config['stratifications'][wildcards.ref_id]['tsv']}",
        threads=config["_happy_threads"],
        engine="vcfeval",
        engine_extra=lambda wildcards: f"--engine-vcfeval-template resources/references/{wildcards.ref_id}.sdf",
    resources:
        mem_mb=config["_happy_mem"],
    threads: config["_happy_threads"]
    log:
        "logs/run_happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/run_happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "envs/happy.yml"
    script:
        "scripts/run_happy.py"


################################################################################
## Run Truvari
# rule run_truvari:
#     input:
#         query="results/dipcall/{bench_id}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}_dipcall.dip.split_multi.vcf.gz",
#         truth=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bmk_prefix, 'compare_var_id']}.vcf.gz",
#         truth_regions=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bmk_prefix, 'compare_var_id']}.bed",
#         truth_tbi=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bmk_prefix, 'compare_var_id']}.vcf.gz.tbi",
#         genome="resources/references/{ref_prefix}.fa",
#         genome_index="resources/references/{ref_prefix}.fa.fai",
#     output:
#         "results/bench/truvari/{bmk_prefix}/summary.txt",
#     log: "logs/run_truvari_{comp_prefix}/truvari.log",
#     params:
#         extra=lambda wildcards: analyses.loc[(wildcards.bmk_prefix, "bench_params")],
#         prefix="results/bench/truvari/{comp_prefix}/",
#         tmpdir="tmp/truvari",
#     conda:
#         "envs/truvari.yml"
#     # TODO this tmp thing is a workaround for the fact that snakemake
#     # over-zealously makes output directories when tools like truvari expect
#     # them to not exist. Also, /tmp is only a thing on Linux (if that matters)
#     shell:
#         """
#         truvari bench \
#             -b {input.truth} \
#             -c {input.query} \
#             -o {params.tmpdir} \
#             -f {input.genome} \
#             --includebed {input.truth_regions} \
#             {params.extra}
#         mv {params.tmpdir}/* {params.prefix}
#         rm -r {params.tmpdir}
#         """
