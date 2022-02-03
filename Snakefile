import pandas as pd
from pathlib import Path
from snakemake.utils import min_version, validate


include: "rules/common.smk"
include: "rules/exclusions.smk"
include: "rules/report.smk"


min_version("6.0")

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
bench_ids = analyses[["bench_id", "vc_id", "ref", "exclusion_set"]].drop_duplicates()
bench_tbl = pd.merge(bench_ids, bench_params, how="inner", on="bench_id").set_index(
    "bench_id"
)

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


## Benchmark Evaluations
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
    vc_param_id="|".join(VCPARAMIDS)


################################################################################
# main rule
#
# Define what files we want hap.py to make, and these paths will contain the
# definitions for the assemblies, variant caller, etc to use in upstream rules.


## Rules to run locally
localrules:
    get_ref,
    get_assemblies,
    get_comparison_vcf,
    get_comparison_bed,
    get_comparison_tbi,
    get_strats,
    index_ref,
    download_bed_gz,
    link_gaps,
    get_SVs_from_vcf,


## Snakemake Report
report: "report/workflow.rst"


## Using zip in rule all to get config sets by config table rows
rule all:
    input:
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
            zip,
            vc_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"].index.tolist(),
            ref=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["ref"].tolist(),
            asm_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["asm_id"].tolist(),
            vc_cmd=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_cmd"].tolist(),
            vc_param_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_param_id"].tolist(),
        ),
        expand(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.summary.csv",
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
        ## rules for report
        expand(
            "results/report/assemblies/{asm_id}_{haplotype}_stats.txt",
            asm_id=ASMIDS,
            haplotype=["maternal", "paternal"],
        ),
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}_stats.txt",
            zip,
            vc_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"].index.tolist(),
            ref=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["ref"].tolist(),
            asm_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["asm_id"].tolist(),
            vc_cmd=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_cmd"].tolist(),
            vc_param_id=vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_param_id"].tolist(),
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
        protected("resources/references/{ref_id}.fa"),
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
    wrapper:
        "0.79.0/bio/samtools/faidx"


################################################################################
# Get stratifications


rule get_strats:
    output:
        protected("resources/strat_{ref_id}/{strat_id}.tar.gz"),
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


use rule get_comparison_vcf as get_comparison_tbi with:
    output:
        "resources/comparison_variant_callsets/{comp_id}.vcf.gz.tbi",
    params:
        url=lambda wildcards: comp_config[wildcards.comp_id]["tbi_url"],
    log:
        "logs/get_comparisons/{comp_id}_vcfidx.log",

## TODO - fix for when tbi url not provided
# rule tabix:
#     input:
#         "{filename}.vcf.gz",
#     output:
#         "{filename}.vcf.gz.tbi",
#     params:
#         extra="-t",
#     log:
#         "logs/tabix/{filename}.log",
#     wrapper:
#         "v1.0.0/bio/bcftools/index"


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
        h1="resources/assemblies/{asm_id}/paternal.fa",
        h2="resources/assemblies/{asm_id}/maternal.fa",
        ref="resources/references/{ref_id}.fa",
        ref_idx="resources/references/{ref_id}.fa.fai",
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
        mem_mb=45000,  ## GB per thread - 16 Gb per job for sorting and estimating 30 max for alignment steps
    threads: config["_dipcall_threads"] * config["_dipcall_jobs"]  ## For diploid
    shell:
        """
        echo "Writing Makefile defining dipcall pipeline"
        run-dipcall \
            -t {params.ts} \
            {params.extra} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            1> {output.make} 2> {log}

        echo "Running dipcall pipeline"
        make -j {params.make_jobs} -f {output.make} &>> {log}
        """


################################################################################
################################################################################
##
## Component: Generating Draft Benchmarkset from Assembly-Based Variant Calls
##
################################################################################
################################################################################

# rule dip_gap2homvarbutfiltered:
#     input: "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz"
#     output: "draft_benchmark_sets/{bench_id}/intermediate/{vc_id}.dip.gap2homvarbutfiltered.vcf.gz"
#     # bgzip is part of samtools, which is part of the diptcall env
#     conda: "envs/dipcall.yml"
#     shell: """
#         gunzip -c {input} |\
#         sed 's/1|\./1|1/' |\
#         grep -v 'HET\|GAP1\|DIP' |\
#         bgzip -c > {output}
#     """


## Using vc_name as a shortened wildcard to avoid writing out full variant call fileaname with multiple wildcards
## TODO workout identification and naming for final draft benchmark variant callset
# rule split_multiallelic_sites:
#     input: lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, "vc_id")]}/{{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}}.dip.vcf.gz",
#     output:
#         vcf="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.split_multi.vcf.gz",
#         vcf_tbi="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.split_multi.vcf.gz.tbi",
#     conda: "envs/bcftools.yml"
#     log: "logs/split_multiallelic_sites/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
#     shell:
#         """
#         bcftools norm -m - {input} -Oz -o {output.vcf} 2> {log}
#         tabix -p vcf {output.vcf} 2>> {log}
#         """


rule postprocess_vcf:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/process_benchmark_vcf/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    shell:
        "cp {input} {output} &> {log}"


rule postprocess_bed:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.bed",
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
        threads=config["happy_threads"],
        engine="vcfeval",
    resources:
        mem_mb=config["happy_mem"],
    threads: config["happy_threads"]
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
