import pandas as pd
from pathlib import Path
from itertools import product
# from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
# from snakemake.io import apply_wildcards
# from functools import partial


include: "rules/common.smk"
include: "rules/exclusions.smk"

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

# analyses = get_analyses("config/analyses.tsv")
#bench_tbl = pd.read_table("config/bench_tbl.tsv").set_index("bench_id", verify_integrity=True)
#eval_tbl = pd.read_table("config/eval_tbl.tsv").set_index("eval_id", verify_integrity=True)
#vc_tbl = pd.read_table("config/vc_tbl.tsv").set_index("vc_id", verify_integrity=True)

## Making a combined data frame for easier variable look-up
#run_tbl = pd.merge(eval_tbl, bench_tbl, how = "outer", on = "bench_id")
#run_tbl = pd.merge(run_tbl, vc_tbl, how = "outer", on = "vc_id")

# analyses = get_analyses("config/analyses.tsv")
vc_tbl = pd.read_table("config/vc_tbl.tsv")
bench_tbl = pd.read_table("config/bench_tbl.tsv")
eval_tbl = pd.read_table("config/eval_tbl.tsv")


## Making a combined data frame
analyses = pd.merge(eval_tbl, bench_tbl, how = "outer", on = "bench_id")
analyses = pd.merge(analyses, vc_tbl, how = "outer", on = "vc_id").set_index("eval_id")

## Setting indices
vc_tbl = vc_tbl.set_index("vc_id")
bench_tbl = bench_tbl.set_index("bench_id")
eval_tbl = eval_tbl.set_index("eval_id")

################################################################################
# init wildcard constraints

## Wildcard variables and ids

## Variables for assembly based variant calling
VCIDS=set(vc_tbl.index.tolist())
REFIDS=set(vc_tbl["ref"].tolist())
ASMIDS=set(vc_tbl["asm_id"].tolist())
VCCMDS=set(vc_tbl["vc_cmd"].tolist())
VCPARAMIDS=set(vc_tbl["vc_param_id"].tolist())

## Draft benchmark set generation variables
BENCHIDS=set(bench_tbl.index.tolist())


## Benchmark Evaluations
EVALIDS=set(eval_tbl.index.tolist())
EVALCOMPIDS=set(eval_tbl["eval_comp_id"].tolist())

# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_id="|".join(ASMIDS),
    comp_id="|".join(EVALCOMPIDS),
    ref_id="|".join(REFIDS),
    bench_id="|".join(BENCHIDS),
    eval_id ="|".join(EVALIDS),


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
    get_genome,
    download_bed_gz,
    link_gaps,
    get_satellites,
    get_SVs_from_vcf

## Using zip in rule all to get config sets by config table rows
rule all:
    input:
        expand("results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
            zip,
            vc_id       = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"].index.tolist(),
            ref         = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["ref"].tolist(),
            asm_id      = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["asm_id"].tolist(),
            vc_cmd      = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_cmd"].tolist(),
            vc_param_id = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]["vc_param_id"].tolist()),
        expand("results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.extended.csv",
            zip,
            eval_id = analyses[analyses["eval_cmd"] == "happy"].index.tolist(),
            bench_id = analyses[analyses["eval_cmd"] == "happy"]["bench_id"].tolist(),
            ref_id = analyses[analyses["eval_cmd"] == "happy"]["ref"].tolist(),
            comp_id = analyses[analyses["eval_cmd"] == "happy"]["eval_comp_id"].tolist(),
            asm_id = analyses[analyses["eval_cmd"] == "happy"]["asm_id"].tolist(),
            vc_cmd = analyses[analyses["eval_cmd"] == "happy"]["vc_cmd"].tolist(),
            vc_param_id = analyses[analyses["eval_cmd"] == "happy"]["vc_param_id"].tolist(),
        )
#       expand("results/bench/truvari/{tvi_bench}.extended.csv", tvi_bench = analyses[analyses["bench_cmd"] == "truvari"].index.tolist()), ## Not yet used

################################################################################
# Get and prepare assemblies


rule get_assemblies:
    output:
        "resources/assemblies/{asm_id}/{haplotype}.fa",
    params:
        url=lambda wildcards: asm_config[wildcards.asm_id][wildcards.haplotype],
    shell:
        "curl -f -L {params.url} | gunzip -c > {output}"


################################################################################
# Get and prepare reference


rule get_ref:
    output: "resources/references/{ref_id}.fa",
    params:
        url=lambda wildcards: ref_config[wildcards.ref_id]["ref_url"],
    shell:
        "curl -f --connect-timeout 120 -L {params.url} | gunzip -c > {output}"


rule index_ref:
    input:
        "resources/references/{ref_id}.fa",
    output:
        "resources/references/{ref_id}.fa.fai",
    wrapper:
        "0.79.0/bio/samtools/faidx"


################################################################################
# Get stratifications


# rule get_strats:
#     output:
#         tsv_full_path,
#     params:
#         root=config["_strats_root"],
#         target=strats_full_path,
#     shell:
#         """
#         curl -L \
#             {params.root}/v3.0/v3.0-stratifications-{wildcards.ref_prefix}.tar.gz | \
#             gunzip -c | \
#             tar x -C {params.target}
#         """


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations

rule get_comparison_vcf:
    output:
        "resources/comparison_variant_callsets/{comp_id}.vcf.gz",
    params:
        url=lambda wildcards: comp_config[wildcards.comp_id]["vcf_url"],
    shell:
        "curl -f -L -o {output} {params.url}"


use rule get_comparison_vcf as get_comparison_bed with:
    output:
        "resources/comparison_variant_callsets/{comp_id}.bed",
    params:
        url=lambda wildcards: comp_config[wildcards.comp_id]["bed_url"],


use rule get_comparison_vcf as get_comparison_tbi with:
    output: 
        "resources/comparison_variant_callsets/{comp_id}.vcf.gz.tbi",
    params: url=lambda wildcards: comp_config[wildcards.comp_id]["tbi_url"],


################################################################################
# Run Dipcall

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
        extra= lambda wildcards: "" if  vc_tbl.loc[wildcards.vc_id]['vc_params'] == "default" else vc_tbl.loc[wildcards.vc_id]['vc_params'],
    log: "logs/asm_varcalls-{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    resources:
        mem_mb=config["_dipcall_threads"] * 2 * 4000,  ## GB per thread
    threads: config["_dipcall_threads"] * 2  ## For diploid
    shell:
        """
        echo "Writing Makefile defining dipcall pipeline"
        run-dipcall \
            {params.extra} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            1> {output.make} 2> {log}

        echo "Running dipcall pipeline"
        make -j{params.ts} -f {output.make} &>> {log}
        """

################################################################################
# Postprocess variant caller output

# rule dip_gap2homvarbutfiltered:
#     input: rules.run_dipcall.output.vcf
#     output: "{}.dip.gap2homvarbutfiltered.vcf.gz".format(vcr_full_prefix)
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
    input: lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz"
    output: "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz"
    shell: "cp {input} {output}"

rule postprocess_bed:
    input: lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.bed"
    output: "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.bed"
    shell: "cp {input} {output}"

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
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
    params:
        prefix="results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
        strat_tsv="resources/stratifications/v3.0/{ref_id}/v3.0-{ref_id}-all-stratifications.tsv",
        threads=config["happy_threads"],
        engine="vcfeval",
    resources:
        mem_mb=config["happy_mem"],
    threads: config["happy_threads"]
    log: "logs/run_happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "envs/happy.yml"
    script: "scripts/run_happy.py"


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
