import pandas as pd
from common import get_exclusion_inputs


wildcard_constraints:
    ref_id="GRCh38|GRCh37|GRCh38_chr21",
    genomic_region="homopolymers|segdups|tandem-repeats|gaps|self-chains|satellites",


# downloading beds used for exclusion
rule download_bed_gz:
    output:
        "resources/exclusions/{ref_id}/{genomic_region}.bed",
    log:
        "logs/download_bed_gz/{ref_id}-{genomic_region}.log",
    params:
        url=lambda wildcards: config["exclusion_beds"][wildcards.genomic_region],
    shell:
        "curl -L {params.url} 2> {log} | gunzip -c 1> {output} 2>> {log}"


# structural variants - using asm varcalls vcf to identify structural variants for exclusion
rule get_SVs_from_vcf:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_dip_SVs.bed",
    log:
        "logs/exclusions/{bench_id}_div_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    group:
        "exclude_SVs"
    shell:
        """
        gunzip -c {input} | \
            awk 'length($4)>49 || length($5)>49' | \
            awk '{{FS=OFS="\\t"}} {{print $1,$2-1,$2+length($4)}}' \
            1> {output} 2> {log}
        """


rule intersect_SVs_and_homopolymers:
    input:
        sv_bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_dip_SVs_sorted.bed",
        homopoly_bed="resources/exclusions/{ref_id}/homopolymers_sorted.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_svs-and-homopolymers_sorted.bed",
    log:
        "logs/exclusions/{bench_id}_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_SVs_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    group:
        "exclude_SVs"
    shell:
        """
        intersectBed -wa \
                -a {input.homopoly_bed} \
                -b {input.sv_bed} | \
            multiIntersectBed -i stdin {input.sv_bed} | \
            awk '{{FS=OFS="\\t"}} {{print $1,$2-50,$3+50}}' | \
            mergeBed -i stdin -d 1000 |
            sortBed -i stdin -g {input.genome} \
            1> {output} 2>{log} 
        """


## Expanding exclusion regions by 15kb
rule add_slop:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slop.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}.bed",
    params:
        slop=15000,
    conda:
        "../envs/bedtools.yml"
    group:
        "postprocess"
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} |
            bedtools slop -i stdin -g {input.genome} -b {params.slop} \
            1> {output} 2> {log}
        """


## Finding breaks in assemblies for excluded regions
rule intersect_start_and_end:
    input:
        dip=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        xregions="resources/exclusions/{ref_id}/{excluded_region}.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        start="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_{excluded_region}_start_sorted.bed",
        end="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_{excluded_region}_end_sorted.bed",
    log:
        "logs/exclusions/start_end_{bench_id}_{excluded_region}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/start_end_{bench_id}_{excluded_region}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    group:
        "start_end"
    shell:
        """
        awk '{{FS=OFS="\t"}} {{print $1, $2, $2+1}}' {input.dip} \
            | bedtools intersect -wa -a {input.xregions} -b stdin \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.start} 2> {log}

        awk '{{FS=OFS="\t"}} {{print $1, $3-1, $3}}' {input.dip} \
            | bedtools intersect -wa -a {input.xregions} -b stdin  \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.end} 2>> {log}
        """


# flanks
rule add_flanks:
    input:
        bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        genome="resources/exclusions/{ref_id}.genome",
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}_flanks.bed",
    log:
        "logs/exclusions/{bench_id}_flanks_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    group:
        "postprocess"
    shell:
        "bedtools flank -i {input.bed} -g {input.genome} -b 15000 1> {output} 2> {log}"


## Removing excluded genomic regions from asm varcalls bed file
## for draft benchmark regions
rule subtract_exclusions:
    input:
        dip_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        other_beds=get_exclusion_inputs,
    output:
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed",
        stats="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded_stats.txt",
    log:
        "logs/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.benchmark"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        python scripts/subtract_exclusions.py \
        {input.dip_bed} \
        {output.bed} \
        {output.stats} \
        {input.other_beds}
        """
