import pandas as pd
from pathlib import Path
import hashlib
from snakemake.utils import min_version, validate
from snakemake.utils import Paramspace


min_version("7.3.0")


## Rule ordering for ambiguous rules
ruleorder: download_bed_gz > sort_bed >postprocess_bed > normalize_for_svwiden > run_svwiden > fix_XY_genotype > move_asm_vcf_to_draft_bench


## Loading external rules
include: "rules/common.smk"
include: "rules/exclusions.smk"
include: "rules/report.smk"
include: "rules/bench_vcf_processing.smk"


################################################################################
# init resources


configfile: workflow.source_path("config/resources.yml")


validate(config, "schema/resources-schema.yml")

asm_config = config["assemblies"]
comp_config = config["comparisons"]
ref_config = config["references"]

################################################################################
# init analyses

## Loading analysis table with run information
analyses = load_analyses(
    workflow.source_path(config["analyses"]), "schema/analyses-schema.yml"
)

vc_params, vc_tbl = analyses_to_vc_tbl(analyses)

ASMIDS = set(vc_tbl["asm_id"])
REFIDS = set(vc_tbl["ref"])

## Wildcard variables and ids

## Variables for assembly based variant calling
REFIDS = set(vc_tbl["ref"].tolist())
ASMIDS = set(vc_tbl["asm_id"].tolist())
VCCMDS = set(vc_tbl["vc_cmd"].tolist())
VCPARAMIDS = set(vc_tbl["vc_param_id"].tolist())
BENCHVCFPROC = set(analyses["bench_vcf_processing"])
BENCHBEDPROC = set(analyses["bench_bed_processing"])
COMPIDS = set(analyses["eval_query"].tolist() + analyses["eval_truth"].tolist())

# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_id="|".join(ASMIDS),
    ref_id="|".join(REFIDS),
    bench_vcf_processing="|".join(BENCHVCFPROC),
    bench_bed_processing="|".join(BENCHBEDPROC),
    comp_dir="asm_varcalls|draft_benchmarksets|evaluations|report",
    comp_id="|".join(COMPIDS)


## Using Paramspace for file paths
happy_space = Paramspace(analyses[analyses["eval_cmd"] == "happy"])
truvari_space = Paramspace(analyses[analyses["eval_cmd"] == "truvari"])

dipcall_space = Paramspace(
    analyses.loc[
        analyses["vc_cmd"] == "dipcall",
        ["asm_id", "ref", "vc_cmd", "vc_param_id", "vc_params"],
    ].drop_duplicates()
)

bench_space = Paramspace(
    analyses[
        [
            "ref",
            "asm_id",
            "vc_cmd",
            "vc_params",
            "vc_param_id",
            "bench_type",
            "bench_vcf_processing",
            "bench_bed_processing",
            "bench_exclusion_set",
        ]
    ].drop_duplicates()
)


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


################################################################################
# main rule
#


rule all:
    input:
        ## Asm based variant calls
        expand(
            "results/asm_varcalls/{params}.dip.vcf.gz",
            params=dipcall_space.instance_patterns,
        ),
        expand(
            "results/asm_varcalls/{params}.hap1.bam.bai",
            params=dipcall_space.instance_patterns,
        ),
        expand(
            "results/asm_varcalls/{params}.hap2.bam.bai",
            params=dipcall_space.instance_patterns,
        ),
        ## Bench VCF Processing
        expand(
            "results/draft_benchmarksets/{params}.processed.vcf.gz",
            params=bench_space.instance_patterns,
        ),
        # expand(
        #     "results/draft_benchmarksets/{params}.processed.vcf.gz.tbi",
        #     params=bench_space.instance_patterns,
        # ),
        expand(
            "results/draft_benchmarksets/{params}.bed",
            params=bench_space.instance_patterns,
        ),
        ## Evaluations
        expand(
            "results/evaluations/happy/{params}.summary.csv",
            params=happy_space.instance_patterns,
        ),
        expand(
            "results/evaluations/happy/{params}.extended.csv",
            params=happy_space.instance_patterns,
        ),
        expand(
            "results/evaluations/truvari/{params}/summary.txt",
            params=truvari_space.instance_patterns,
        ),
        ## Summary Stats/ Reports
        expand(
            "results/report/assemblies/{asm_id}_{haplotype}_stats.txt",
            asm_id=ASMIDS,
            haplotype=["maternal", "paternal"],
        ),
        expand(
            "results/asm_varcalls/{params}.dip_bcftools_stats.txt",
            params=dipcall_space.instance_patterns,
        ),
        expand(
            "results/asm_varcalls/{params}.dip_rtg_stats.txt",
            params=dipcall_space.instance_patterns,
        ),


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
        url=lambda wildcards: f"{config['references'][wildcards.ref_id]['stratifications']['url']}",
    log:
        "logs/get_strats/{ref_id}_{strat_id}.log",
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations


rule get_comparison_vcf:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "vcf_url"
        ],
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_vcf.log",
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


use rule get_comparison_vcf as get_comparison_bed with:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "bed_url"
        ],
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_bed.log",


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
        h1="resources/assemblies/{asm_id}/paternal.fa",
        h2="resources/assemblies/{asm_id}/maternal.fa",
        ref="resources/references/{ref}.fa",
        ref_idx="resources/references/{ref}.fa.fai",
        ref_mmi="resources/references/{ref}.mmi",
        par=get_male_bed,
    output:
        make=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.mak",
        vcf=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.dip.vcf.gz",
        bed=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.dip.bed",
        bam1=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.hap1.bam",
        bam2=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.hap2.bam",
    conda:
        "envs/dipcall.yml"
    params:
        prefix=f"results/asm_varcalls/{dipcall_space.wildcard_pattern}",
        male_bed=get_dipcall_par_param,
        ts=config["_dipcall_threads"],
        make_jobs=config["_dipcall_jobs"],
        extra=lambda wildcards: ""
        if wildcards.vc_params == "default"
        else wildcards.vc_params,
    log:
        multiext(
            f"results/asm_varcalls/{dipcall_space.wildcard_pattern}",
            ".hap1.paf.gz.log",
            ".hap2.paf.gz.log",
            ".hap1.sam.gz.log",
            ".hap2.sam.gz.log",
        ),
        rulelog=f"logs/asm_varcalls/{dipcall_space.wildcard_pattern}.log",
    benchmark:
        f"benchmark/asm_varcalls/{dipcall_space.wildcard_pattern}.tsv"
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
        in_file="results/{comp_dir}/{prefix}.bed",
        genome=get_genome_file,
    output:
        "results/{comp_dir}/{prefix}.sorted.bed",
    log:
        "logs/sort_bed/{comp_dir}/{prefix}.log",
    wrapper:
        "0.74.0/bio/bedtools/sort"


rule index_dip_bam:
    input:
        f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.{{hap}}.bam",
    output:
        f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.{{hap}}.bam.bai",
    log:
        f"logs/asm_varcalls/{dipcall_space.wildcard_pattern}.{{hap}}.bam.bai.log",
    wrapper:
        "v1.1.0/bio/samtools/index"


################################################################################
################################################################################
##
## Component: Generating Draft Benchmarkset from Assembly-Based Variant Calls
##
################################################################################
################################################################################


rule move_asm_vcf_to_draft_bench:
    input:
        f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.dip.vcf.gz",
    output:
        f"results/draft_benchmarksets/intermediates/{bench_space.wildcard_pattern}.vcf.gz",
    log:
        f"logs/process_benchmark_vcf/{bench_space.wildcard_pattern}.log",
    shell:
        "cp {input} {output} &> {log}"


rule move_processed_draft_bench_vcf:
    input:
        get_processed_vcf,
    output:
        f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.processed.vcf.gz",
    log:
        f"logs/move_processed_draft_bench_vcf/{bench_space.wildcard_pattern}.log",
    shell:
        "cp {input} {output}"


rule postprocess_bed:
    input:
        f"results/asm_varcalls/{dipcall_space.wildcard_pattern}.dip.bed",
    output:
        f"results/draft_benchmarksets/{bench_space.wildcard_pattern}.bed",
    log:
        f"logs/process_benchmark_bed/{bench_space.wildcard_pattern}.log",
    shell:
        "cp {input} {output} &> {log}"


################################################################################
################################################################################
##
## Component: Evaluating Draft Benchmarksets
##
################################################################################
################################################################################

## hap.py small variant benchmarking tool


rule run_happy:
    input:
        unpack(partial(get_eval_inputs, analyses, config)),
    output:
        multiext(
            f"results/evaluations/happy/{happy_space.wildcard_pattern}",
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
            f"results/evaluations/happy/{happy_space.wildcard_pattern}.summary.csv",
            caption="report/happy_summary.rst",
            category="Happy",
        ),
    params:
        prefix=f"results/evaluations/happy/{happy_space.wildcard_pattern}",
        strat_tsv=lambda wildcards: f"{wildcards.ref}/{ref_config[wildcards.ref]['stratifications']['tsv']}",
        threads=config["_happy_threads"],
        engine="vcfeval",
        engine_extra=lambda wildcards: "--engine-vcfeval-template resources/references/{wildcards.ref}.sdf",
    resources:
        mem_mb=config["_happy_mem"],
    threads: config["_happy_threads"]
    log:
        f"logs/run_happy/{happy_space.wildcard_pattern}.log",
    benchmark:
        f"benchmark/run_happy/{happy_space.wildcard_pattern}.tsv"
    conda:
        "envs/happy.yml"
    script:
        "scripts/run_happy.py"


################################################################################
## Run Truvari


rule run_truvari:
    input:
        unpack(partial(get_eval_inputs, analyses, config)),
    output:
        f"results/evaluations/truvari/{truvari_space.wildcard_pattern}/summary.txt",
    log:
        f"logs/run_travari/{truvari_space.wildcard_pattern}/truvari.log",
    # TODO this tmp thing is a workaround for the fact that snakemake
    # over-zealously makes output directories when tools like truvari expect
    # them to not exist. Also, /tmp is only a thing on Linux (if that matters).
    # Also^2, certain cluster admins (such as those that run Nisaba) don't like
    # it when we use /tmp
    params:
        dir=lambda wildcards, output: Path(output[0]).parent,
        tmpdir=f"tmp_truvari",
    conda:
        "envs/truvari.yml"
    shell:
        """
        ## Removing temp directory if present before run
        rm -rf {params.dir}

        truvari bench \
            -b {input.truth} \
            -c {input.query} \
            -o {params.dir} \
            -f {input.genome} \
            --passonly \
            --includebed {input.truth_regions} \
        2> {log}

        echo {params.dir}
        # mv {params.tmpdir}/* {params.dir}
        # rm -r {params.tmpdir}
        """
