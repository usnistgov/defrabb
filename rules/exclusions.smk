import pandas as pd

# from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate

# results paths
dip_bed_path = (
    "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.bed"
)

wildcard_constraints:
    ref_id="GRCh38|GRCh37|GRCh38_chr21",
    genomic_regions="homopolymers|segdups|tandem_repeat|gaps|self_chains|satellites"

# downloading stuff


rule download_bed_gz:
    output:
        "resources/exclusions/{ref_id}/{genomic_region}.bed",
    log:
        "logs/download_bed_gz/{ref_id}-{genomic_region}.log",
    params:
        url=lambda wildcards: config["exclusion_beds"][wildcards.genomic_region],
    group: "download"
    shell:
        "curl -L {params.url} 2> {log} | gunzip -c 1> {output} 2>> {log}"


# TODO hack since this is the only bed file that isn't processed according to
# the output from dipcall
rule link_gaps:
    input:
        "resources/exclusions/{ref_id}/gaps.bed",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/gaps.bed",
    log:
        "logs/exclusions/{bench_id}_gaps_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    group: "postprocess"
    shell:
        "cp {input} {output} &> {log}"


# structural variants


rule get_SVs_from_vcf:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/dip_SVs.bed",
    log:
        "logs/exclusions/{bench_id}_div_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    group: "postprocess"
    shell:
        """
        gunzip -c {input} | \
            awk 'length($4)>49 || length($5)>49' | \
            awk '{{FS=OFS="\\t"}} {{print $1,$2-1,$2+length($4)}}' \
            1> {output} 2> {log}
        """


rule intersect_SVs_and_homopolymers:
    input:
        sv_bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/dip_SVs.bed",
        homopoly_bed="resources/exclusions/{ref_id}/homopolymers.bed",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/structural_variants.bed",
    log:
        "logs/exclusions/{bench_id}_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    group: "postprocess"
    shell:
        """
        intersectBed -wa \
        -a {input.homopoly_bed} \
        -b {input.sv_bed} | \
        multiIntersectBed -i stdin {input.sv_bed} | \
        awk '{{FS=OFS="\\t"}} {{print $1,$2-50,$3+50}}' | \
        mergeBed -i stdin -d 1000 \
        1> {output} 2>{log}
        """


# segdups/tandemreps/self chains


rule add_slop:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_regions}.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/slop/{genomic_regions}.bed",
    log:
        "logs/exclusions/{bench_id}_slop_{genomic_regions}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed",
    conda:
        "../envs/bedtools.yml"
    group: "postprocess"
    shell:
        "bedtools slop -i {input.bed} -g {input.genome} -b 15000 1> {output} 2> {log}"


rule intersect_start_and_end:
    input:
        dip=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.bed",
        xregions="resources/exclusions/{ref_id}/{genomic_regions}.bed",
    output:
        start="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/{genomic_regions}_start.bed",
        end="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/{genomic_regions}_end.bed",
    log:
        "logs/exclusions/start_end_{bench_id}_{genomic_regions}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/start_end_{bench_id}_{genomic_regions}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    group: "postprocess"
    shell:
        """
        awk '{{FS=OFS="\t"}} {{print $1, $2, $2+1}}' {input.dip} | \
        bedtools intersect -wa -a {input.xregions} -b stdin \
        1> {output.start} 2> {log}

        awk '{{FS=OFS="\t"}} {{print $1, $3-1, $3}}' {input.dip} | \
        bedtools intersect -wa -a {input.xregions} -b stdin \
        1> {output.end} 2>> {log}
        """


# flanks
rule add_flanks:
    input:
        bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_exclusions/flanks.bed",
    log:
        "logs/exclusions/{bench_id}_flanks_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    group: "postprocess"
    shell:
        "bedtools flank -i {input.bed} -g {input.genome} -b 15000 1> {output} 2> {log}"


rule subtract_exclusions:
    input:
        dip_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.bed",
        other_beds=lookup_excluded_region_set,
    output:
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed",
        stats="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded_stats.txt"
    log:
        "logs/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.benchmark"
    conda:
        "../envs/bedtools.yml"
    group: "postprocess"
    shell:
        """
        python scripts/subtract_exclusions.py \
        {input.dip_bed} \
        {output} \
        {input.other_beds} 1> {output.stats} 2> {log}
        """
