import pandas as pd
from pathlib import Path
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.io import apply_wildcards
from functools import partial


include: "rules/common.smk"
include: "rules/exclusions.smk"

min_version("6.0")

################################################################################
# init resources


configfile: "config/resources.yml"


validate(config, "schema/resources-schema.yml")

asm_config = config["assemblies"]
bmk_config = config["benchmarks"]
ref_config = config["references"]

################################################################################
# init analyses

analyses = get_analyses("config/analyses.tsv")

################################################################################
# init wildcard constraints


# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_prefix=format_constraint(asm_config),
    bench_id=format_constraint(bmk_config),
    rench_idprefix=format_constraint(ref_config),


################################################################################
# main rule
#
# Define what files we want hap.py to make, and these paths will contain the
# definitions for the assemblies, variant caller, etc to use in upstream rules.

## Rules to run locally
localrules:
    get_ref,
    get_assemblies,
    get_benchmark_vcf,
    get_benchmark_bed,
    get_genome,
    download_bed_gz,
    link_gaps,
    get_satellites,
    get_SVs_from_vcf


rule all:
    input:
        expand("results/bench/happy/{bench}_{ref}_{asm}_{vcr}-{vcrp}.extended.csv",
            zip,
            bench = analyses[analyses["bench_cmd"] == "happy"].index.tolist(),
            ref = analyses.loc[analyses["bench_cmd"] == "happy"]["ref"].tolist(),
            asm = analyses.loc[analyses["bench_cmd"] == "happy"]["asm_id"].tolist(),
            vcr = analyses.loc[analyses["bench_cmd"] == "happy"]["varcaller"].tolist(),
            vcrp = analyses.loc[analyses["bench_cmd"] == "happy"]["vc_param_id"].tolist()),
#       expand("results/bench/truvari/{tvi_bench}.extended.csv", tvi_bench = analyses[analyses["bench_cmd"] == "truvari"].index.tolist()), ## Not yet used

################################################################################
# Get and prepare assemblies


rule get_assemblies:
    output:
        "resources/assemblies/{asm_prefix}/{haplotype}.fa",
    params:
        url=lambda wildcards: asm_config[wildcards.asm_prefix][wildcards.haplotype],
    shell:
        "curl -f -L {params.url} | gunzip -c > {output}"


################################################################################
# Get and prepare reference


rule get_ref:
    output:
        "resources/references/{ref_prefix}.fa",
    params:
        url=lambda wildcards: bmk_config[wildcards.ref_prefix]["ref_url"],
    shell:
        "curl -f --connect-timeout 120 -L {params.url} | gunzip -c > {output}"


rule index_ref:
    input:
        "resources/references/{ref_prefix}.fa",
    output:
        "resources/references/{ref_prefix}.fai",
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
# Get benchmark vcf.gz and .bed

rule get_benchmark_vcf:
    output:
        "resources/benchmarks/{bmk_prefix}.vcf.gz",
    params:
        url=lambda wildcards: ref_config[wildcards.bmk_prefix]["vcf_url"],
    shell:
        "curl -f -L -o {output} {params.url}"


use rule get_benchmark_vcf as get_benchmark_bed with:
    output:
        "resources/benchmarks/{bmk_prefix}.bed",
    params:
        url=lambda wildcards: ref_config[wildcards.bmk_prefix]["bed_url"],


use rule get_benchmark_vcf as get_benchmark_tbi with:
    output:
        "resources/benchmarks/{bmk_prefix}.vcf.gz.tbi",
    params:
        url=lambda wildcards: ref_config[wildcards.bmk_prefix]["tbi_url"],


################################################################################
# Run Dipcall

rule run_dipcall:
    input:
        h1="resources/assemblies/{asm_prefix}/paternal.fa",
        h2="resources/assemblies/{asm_prefix}/maternal.fa",
        ref="resources/references/{ref_prefix}.fa",
        ref_idx="resources/references/{ref_prefix}.fa.fai",
    output:
        make="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.mak",
        vcf="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.vcf.gz",
        bed="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed",
        bam1="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.hap1.bam",
        bam2="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.hap2.bam",
    conda:
        "envs/dipcall.yml"
    params:
        prefix="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall",
        male_bed=get_male_bed,
        ts=config["_dipcall_threads"],
        extra= lambda wildcards: analyses.loc[wildcards.vcr_param_id]["vc_params"],
    log: "log/dipcall_{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.log",
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
            > {output.make} &> {log}

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


rule split_multiallelic_sites:
    input: "results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.vcf.gz",
    output:
        vcf="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.split_multi.vcf.gz",
        vcf_tbi="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.split_multi.vcf.gz.tbi",
    conda: "envs/bcftools.yml"
    log: "logs/split_multiallelic_sites/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}.log",
    shell:
        """
        bcftools norm -m - {input} -Oz -o {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


################################################################################
## Run happy

rule run_happy:
    input:
        unpack(get_happy_inputs),
        strat_tb="resources/v3.0-stratifications-{ref_prefix}.tar.gz",
        genome="resources/references/{ref_prefix}.fa",
        genome_index="resources/references/{ref_prefix}.fa.fai",
    output:
        multiext(
            "results/bench/happy/{bench_id}_{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}",
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
        prefix="results/bench/happy/{bench_id}",
        strat_tsv="resources/stratifications/v3.0/{ref_prefix}/v3.0-{ref_prefix}-all-stratifications.tsv",
        threads=config["happy_threads"],
        engine="vcfeval",
        extra=lambda wildcards: "" if analyses.loc[(wildcards.bench_id, "target_regions")] == "false" else "-T {{input.target_regions}}",
    resources:
        mem_mb=config["happy_mem"],
    threads: config["happy_threads"]
    log: "logs/run_happy_{bench_id}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}.log",
    conda:
        "envs/happy.yml"
    shell:
        """
        echo "Extracting Stratifications tarball"
        tar -xvf {input.strat_tb}

        echo "Starting Happy Run"

        hap.py \
            --threads {params.threads} \
            --engine {params.engine} \
            -r {input.genome}  \
            -f {input.truth_regions} \
            --stratification {params.strat_tsv} \
            -o {params.prefix} \
            --verbose \
            {params.extra} \
            {input.truth} \
            {input.query}
        """


################################################################################
## Run Truvari

rule run_truvari:
    input:
        query="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.split_multi.vcf.gz",
        truth=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bench_id, 'compare_var_id']}.vcf.gz",
        truth_regions=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bench_id, 'compare_var_id']}.bed",
        truth_tbi=lambda wildcards: f"resources/benchmarks/{analyses.loc[wildcards.bench_id, 'compare_var_id']}.vcf.gz.tbi",
        genome="resources/references/{ref_prefix}.fa",
        genome_index="resources/references/{ref_prefix}.fa.fai",
    output:
        "results/bench/truvari/{bench_id}/summary.txt",
    log: "logs/run_truvari_{bench_id}/truvari.log",
    params:
        extra=lambda wildcards: analyses.loc[(wildcards.bench_id, "bench_params")],
        prefix="results/bench/truvari/{bench_id}/",
        tmpdir="tmp/truvari",
    conda:
        "envs/truvari.yml"
    # TODO this tmp thing is a workaround for the fact that snakemake
    # over-zealously makes output directories when tools like truvari expect
    # them to not exist. Also, /tmp is only a thing on Linux (if that matters)
    shell:
        """
        truvari bench \
            -b {input.truth} \
            -c {input.query} \
            -o {params.tmpdir} \
            -f {input.genome} \
            --includebed {input.truth_regions} \
            {params.extra}
        mv {params.tmpdir}/* {params.prefix}
        rm -r {params.tmpdir}
        """
