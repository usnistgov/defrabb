import pandas as pd
from pathlib import Path
from snakemake.utils import min_version, validate


min_version("7.26.0")


## Rule ordering for ambiguous rules
ruleorder: make_gaps_bed > download_bed_gz > sort_bed


## Loading external rules
include: "rules/common.smk"
include: "rules/utils.smk"
include: "rules/download_resources.smk"
include: "rules/asm-varcall.smk"
include: "rules/exclusions.smk"
include: "rules/report.smk"
include: "rules/bench_vcf_processing.smk"
include: "rules/evaluation.smk"


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
    copy_asm_vcf_to_annotations,


## Snakemake Report
report: "report/workflow.rst"


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
        expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz.tbi",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        ## Bench VCF Processing
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            bench_type=bench_tbl["bench_type"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            bench_type=bench_tbl["bench_type"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
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
        ## Bench BED Processing
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark_bed-summary.tsv",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        ## rules for report
        "analysis.html",
        "results/analysis_params.yml",
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
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip_bed-summary.tsv",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_bench-vars_rtg_stats.txt",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            bench_type=bench_tbl["bench_type"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_intersection_summary.csv",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        expand(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_smvar_{vc_cmd}-{vc_param_id}.summary.csv",
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
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_smvar_{vc_cmd}-{vc_param_id}.extended.csv",
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
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_stvar_{vc_cmd}-{vc_param_id}/summary.json",
            zip,
            eval_id=truvari_analyses.index.tolist(),
            bench_id=truvari_analyses["bench_id"].tolist(),
            ref_id=truvari_analyses["ref"].tolist(),
            comp_id=truvari_analyses["eval_comp_id"].tolist(),
            asm_id=truvari_analyses["asm_id"].tolist(),
            vc_cmd=truvari_analyses["vc_cmd"].tolist(),
            vc_param_id=truvari_analyses["vc_param_id"].tolist(),
        ),
        expand(
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_stvar_{vc_cmd}-{vc_param_id}/refine.variant_summary.json",
            zip,
            eval_id=truvari_refine_analyses.index.tolist(),
            bench_id=truvari_refine_analyses["bench_id"].tolist(),
            ref_id=truvari_refine_analyses["ref"].tolist(),
            comp_id=truvari_refine_analyses["eval_comp_id"].tolist(),
            asm_id=truvari_refine_analyses["asm_id"].tolist(),
            vc_cmd=truvari_refine_analyses["vc_cmd"].tolist(),
            vc_param_id=truvari_refine_analyses["vc_param_id"].tolist(),
        ),
