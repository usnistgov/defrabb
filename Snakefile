import pandas as pd
from pathlib import Path
import hashlib
from snakemake.utils import min_version, validate
from snakemake.utils import Paramspace


min_version("7.3.0")


## Rule ordering for ambiguous rules
ruleorder: get_beds_for_exclusions > sort_bed > postprocess_bed > normalize_for_svwiden > run_svwiden > fix_XY_genotype > move_asm_vcf_to_draft_bench


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
BENCHVCFPROC = set(analyses["bench_vcf_processing"])
BENCHBEDPROC = set(analyses["bench_bed_processing"])
COMPIDS = set(analyses["eval_query"].tolist() + analyses["eval_truth"].tolist())
EVALIDS = set(
    analyses["eval_query"].tolist() + analyses["eval_truth"].tolist() + ["this_row"]
)
EVALTRUTHREGIONS = set(analyses['eval_truth_regions'].tolist())
EVALTARGETREGIONS = set(analyses['eval_target_regions'].tolist())

# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).


wildcard_constraints:
    asm_id="|".join(ASMIDS),
    ref_id="|".join(REFIDS),
    bench_vcf_processing="|".join(BENCHVCFPROC),
    bench_bed_processing="|".join(BENCHBEDPROC),
    comp_dir="asm_varcalls|draft_benchmarksets|evaluations|report",
    comp_id="|".join(COMPIDS),
    comp_ext="vcf.gz|vcf|bed|bed.gz",
    eval_truth="|".join(EVALIDS),
    eval_query="|".join(EVALIDS),
    eval_truth_regions="|".join(EVALTRUTHREGIONS),
    eval_target_regions="|".join(EVALTARGETREGIONS)


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


## Rules to run locally
localrules:
    get_ref,
    get_assemblies,
    get_comparisons,
    get_strats,
    get_beds_for_exclusions,
    get_SVs_from_vcf,
    subtract_exclusions,
    add_flanks,
    intersect_start_and_end,
    intersect_SVs_and_simple_repeats,
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
## Variables for eval renaming
defrabb_version = 0.008
sample_id = "HG002"


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
        expand(
            "results/draft_benchmarksets/{params}.bed",
            params=bench_space.instance_patterns,
        ),
        ## Evaluations
        expand(
            "results/evaluations/{params}.extended.csv",
            params=happy_space.instance_patterns,
        ),
        expand(
            "results/evaluations/{params}/summary.txt",
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
    conda:
        "envs/get_remotes.yml"
    params:
        source_uri=get_asm_uri,
        source_hash=get_asm_checksum,
        hash_algo=get_asm_checksum_algo,
        outfmt="decompressed",
    log:
        "logs/get_assemblies/{asm_id}_{haplotype}.log",
    script:
        "scripts/download_resources.py"


# Get and prepare reference
rule get_ref:
    output:
        "resources/references/{ref_id}.fa",
    conda:
        "envs/get_remotes.yml"
    params:
        source_uri=get_ref_uri,
        source_hash=get_ref_checksum,
        hash_algo=get_ref_checksum_algo,
        outfmt="decompressed",
    log:
        "logs/get_ref/{ref_id}.log",
    script:
        "scripts/download_resources.py"


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
    conda:
        "envs/get_remotes.yml"
    params:
        source_uri=get_strats_uri,
        source_hash=get_strats_checksum,
        hash_algo=get_strats_checksum_algo,
        outfmt="gzip",
    log:
        "logs/get_strats/{ref_id}_{strat_id}.log",
    script:
        "scripts/download_resources.py"


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations


rule get_comparisons:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.{comp_ext}",
    conda:
        "envs/get_remotes.yml"
    params:
        source_uri=get_comp_uri,
        source_hash=get_comp_checksum,
        hash_algo=get_comp_checksum_algo,
        outfmt=get_comp_outfmt,
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_{comp_ext}.log",
    script:
        "scripts/download_resources.py"


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
        if wildcards.vc_param_id == "default"
        else config["_dipcall_params"][wildcards.vc_param_id],
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


# downloading beds used for exclusion
rule get_beds_for_exclusions:
    output:
        "resources/exclusions/{ref_id}/{genomic_region}.bed",
    conda:
        "envs/get_remotes.yml"
    log:
        "logs/download_bed_gz/{ref_id}-{genomic_region}.log",
    params:
        source_uri=get_exclusion_uri,
        source_hash=get_exclusion_checksum,
        hash_algo=get_exclusion_checksum_algo,
        outfmt="decompressed",
    script:
        "scripts/download_resources.py"


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
            f"results/evaluations/{happy_space.wildcard_pattern}",
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
            f"results/evaluations/{happy_space.wildcard_pattern}.summary.csv",
            caption="report/happy_summary.rst",
            category="Happy",
        ),
    params:
        prefix=f"results/evaluations/{happy_space.wildcard_pattern}",
        strat_tsv=lambda wildcards: f"{wildcards.ref}/{ref_config[wildcards.ref]['stratifications']['tsv']}",
        threads=config["_happy_threads"],
        engine="vcfeval",
        engine_extra=lambda wildcards: "--engine-vcfeval-template resources/references/{wildcards.ref}.sdf",
    resources:
        mem_mb=config["_happy_mem"],
    threads: config["_happy_threads"]
#    log:
 #       f"logs/run_happy/{happy_space.wildcard_pattern}.log",
        # "logs/run_happy/defrabb-0.008_eval_HG002_{ref}_{eval_cmd}_T~{eval_truth}_TR~{eval_truth_regions}_Q~{eval_query}_QR~{eval_target_regions}.log",
#    benchmark:
 #       f"benchmark/run_happy/{happy_space.wildcard_pattern}.tsv"
        # "benchmark/run_happy/defrabb-0.008_eval_HG002_{ref}_{eval_cmd}_T~{eval_truth}_TR~{eval_truth_regions}_Q~{eval_query}_QR~{eval_target_regions}.tsv"
    conda:
        "envs/happy.yml"
    script:
        "scripts/run_happy.py"


################################################################################
## Run Truvari


rule run_truvari:
    input:
        unpack(partial(get_eval_inputs, analyses_wids, config)),
    output:
        f"results/evaluations/{truvari_space.wildcard_pattern}/summary.txt",
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
