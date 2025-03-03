import pandas as pd

genomic_regions = [
    "all-tr-and-homopolymers",
    "segdups",
    "tandem-repeats",
    "gaps",
    "self-chains",
    "satellites",
    "hifi-pacbioDV-XY-discrep",
    "imperfecthomopol-gt30",
    "hifiasm-HPRC-T2Tdiscrep",
    "XYelement-homopolymer-T2T-discrep",
    "XYdipcallmanualbugs",
    "VDJ",
    "consecutive-svs",
    "dipcall-bugs-T2TACE",
    "HG002Q100-pav-discrep-smvar",
    "HG002Q100-pav-discrep-stvar",
    "HG002Q100-pav-inversions",
    "HG002Q100-errors",
    "HG002Q100-mosaic",
    "HG002Q100-delins-errors",
    "TSPY2-segdups",
    "self-discrep",
]


wildcard_constraints:
    ref_id="GRCh38|GRCh37|GRCh38_chr21|CHM13v2.0",
    genomic_region="|".join(genomic_regions),


# downloading beds used for exclusion
rule download_bed_gz:
    output:
        "resources/exclusions/{ref_id}/{genomic_region}.bed",
    log:
        "logs/download_bed_gz/{ref_id}-{genomic_region}.log",
    params:
        url=lambda wildcards: config["references"][wildcards.ref_id]["exclusions"][
            wildcards.genomic_region
        ],
    conda:
        "../envs/download_remotes.yml"
    shell:
        "curl -L {params.url} 2> {log} | gunzip -c 1> {output} 2> {log}"


rule make_gaps_bed:
    input:
        fa="resources/references/{ref_id}.fa",
    output:
        bed="resources/exclusions/{ref_id}/gaps.bed",
    log:
        "logs/make_gap_bed/{ref_id}.log",
    params:
        minNs=50,
    conda:
        "../envs/seqtk.yml"
    shell:
        "seqtk gap -l {params.minNs} {input.fa} 1> {output.bed} 2> {log}"


rule get_sv_widen_coords:
    input:
        vcf=get_processed_vcf,
        genome=get_genome_file,
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_dip_SVs.bed",
        tbl="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_dip_SVs.tsv",
    conda:
        "../envs/sv_widen_coords.yml"
    params:
        script=Path(workflow.basedir) / "scripts/get_sv_widen_coords.py",
    log:
        "logs/exclusions/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_dip_SV_coords.log",
    shell:
        """
        python {params.script} \
            --input {input.vcf} \
            --output {output.bed} \
            --verbose \
            --table \
            --sort-merge \
            --genome {input.genome} \
            --log {log} 2>> {log}
        """


rule intersect_SVs_and_simple_repeats:
    input:
        sv_bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_dip_SVs.bed",
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_svs-and-simple-repeats.bed",
    log:
        "logs/exclusions/{bench_id}_SVs_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_SVs_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -wa \
                -a {input.simple_repeat_bed} \
                -b {input.sv_bed} | \
            multiIntersectBed -i stdin {input.sv_bed} | \
            bedtools slop -i stdin -g {input.genome} -b 50 | \
            mergeBed -i stdin -d 1000 \
            1> {output} 2>{log}
        """


## BED file with consecutive insertions and deletions from assembly-assembly bams
rule get_consecutive_svs:
    input:
        hap1_bam=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.hap1.bam",
        hap1_bai=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.hap1.bam.bai",
        hap2_bam=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.hap2.bam",
        hap2_bai=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.hap2.bam.bai",
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_consecutive-svs.bed",
    log:
        "logs/exclusions/consecutive-svs/{bench_id}_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    params:
        min_bp=10,
    conda:
        "../envs/consecutive_svs.yml"
    shell:
        """
        python scripts/get_dipcall_delins.py \
            --min_bp {params.min_bp} \
            --hap1_bam {input.hap1_bam} \
            --hap2_bam {input.hap2_bam} \
            --output_bed {output.bed} \
            &> {log}
        """


## Excluding self comparison
rule self_discrep_happy:
    input:
        vcf=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
        bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        ref=get_ref_file,
        sdf=get_ref_sdf,
    output:
        multiext(
            "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
            ".runinfo.json",
            ".vcf.gz",
            ".vcf.gz.tbi",
            ".summary.csv",
            ".extended.csv",
        ),
    log:
        "logs/exclusions/self-discrep-happy/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    params:
        ## Fix to not hard code and use output prefix
        prefix="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
        engine="vcfeval",
        gender=get_happy_gender_param,
    resources:
        mem_mb=config["_happy_mem"],
    threads: config["_happy_threads"]
    benchmark:
        "benchmark/self-discrep-happy/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/happy.yml"
    shell:
        """
        hap.py \
            {input.vcf} \
            {input.vcf} \
            -R {input.bed}  \
            -r {input.ref}  \
            -o {params.prefix} \
            {params.gender} \
            --pass-only \
            --no-roc \
            --no-json \
            --engine=vcfeval \
            --engine-vcfeval-template {input.sdf} \
            --threads={threads} \
            &> {log}
        """


rule self_discrep_extract_fpfns:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
        faidx=get_ref_index,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.fpfns.bed",
    log:
        "logs/exclusions/self-discrep-fpfns/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    params:
        max_indel=50,
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        bcftools filter \
            --include 'MAX(ILEN)<={params.max_indel} && MIN(ILEN) >= -{params.max_indel} && (FMT/BD=="FN" || FMT/BD=="FP")' \
            {input.vcf} |
                bcftools query -f "%CHROM\t%POS0\t%END\n" |
                bedtools merge -i - | 
                bedtools sort -faidx {input.faidx} -i - \
            1> {output} 2> {log}
        """


rule self_discrep_intersect_slop:
    input:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.fpfns.bed",
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_self-discrep.bed",
    log:
        "logs/exclusions/self-discrep-intersect/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    params:
        slop=50,
        merge_d=1000,
    benchmark:
        "benchmark/exclusions/{bench_id}_self-discrep_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        ## TODO make slop conditional, don't add for intersected vars
        bedtools intersect -wa \
                -a {input.simple_repeat_bed} \
                -b {input.bed} | \
            bedtools multiinter -i stdin {input.bed} | \
            bedtools slop -b {params.slop} -i stdin -g {input.genome} | \
            mergeBed -i stdin -d {params.merge_d} \
            1> {output} 2>{log}
        """


## Expanding exclusion regions by 15kb
rule add_slop:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slop.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}_slop.log",
    params:
        slop=15000,
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} |
            bedtools slop -i stdin -g {input.genome} -b {params.slop} \
            1> {output} 2> {log}
        """


## Expanding regions by 15kb then merging regions within 10kb
rule add_slop_and_merge:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slopmerge.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}_slopmerge.log",
    params:
        slop=15000,
        dist=10000,
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} |
            bedtools slop -i stdin -g {input.genome} -b {params.slop} |
            bedtools merge -i stdin -d {params.dist} \
            1> {output} 2> {log}
        """


## Finding breaks in assemblies for excluded regions
rule intersect_start_and_end:
    input:
        dip=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        xregions="resources/exclusions/{ref_id}/{excluded_region}.bed",
        genome=get_genome_file,
    output:
        start="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_{excluded_region}_start_sorted.bed",
        end="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_{excluded_region}_end_sorted.bed",
    log:
        "logs/exclusions/start_end_{bench_id}_{excluded_region}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/start_end_{bench_id}_{excluded_region}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk '{{FS=OFS="\t"}} {{print $1, $2, $2+1}}' {input.dip} \
            | bedtools intersect -u -wa -a {input.xregions} -b stdin \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.start} 2> {log}

        awk '{{FS=OFS="\t"}} {{print $1, $3-1, $3}}' {input.dip} \
            | bedtools intersect -u -wa -a {input.xregions} -b stdin  \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.end} 2>> {log}
        """


# Generate bed with 15kb regions around assembly breaks (non-dip coverage)
rule get_flanks:
    input:
        bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_flanks.bed",
    params:
        bases=15000,
    log:
        "logs/exclusions/{bench_id}_flanks_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools complement -i {input.bed} -g {input.genome} |
            bedtools flank -i stdin -g {input.genome} -b {params.bases} \
            1> {output} 2> {log}
        """


## Removing excluded genomic regions from asm varcalls bed file
## for draft benchmark regions
rule subtract_exclusions:
    input:
        dip_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        other_beds=get_exclusion_inputs,
    output:
        rpt=report(
            "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_stats.txt",
            caption="../report/exclusion_stats.rst",
            category="Exclusions",
        ),
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed",
    params:
        script=Path(workflow.basedir) / "scripts/subtract_exclusions.py",
    log:
        "logs/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        python {params.script} \
        {input.dip_bed} \
        {output.bed} \
        {output.rpt} \
        {input.other_beds} \
        &> {log}
        """


rule generate_intersection_summary:
    input:
        dip_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
        exclusions=get_exclusion_inputs,
    output:
        summary_table=report(
            "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_intersection_summary.csv",
            caption="../report/exclusion_intersection.rst",
            category="Exclusions",
        ),
    params:
        script=Path(workflow.basedir) / "scripts/exclusion_intersection_summary.py",
        intersect_dir="results/draft_benchmarksets/{bench_id}/exclusions/intersections/",
    log:
        "logs/exclusion-intersect/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """  
        python {params.script} {input.dip_bed} {output.summary_table} {params.intersect_dir} {input.exclusions} &> {log}  
        """


## Used when no exclusions are applied
rule postprocess_bed:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip_sorted.bed",
    output:
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed",
    log:
        "logs/process_benchmark_bed/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "cp {input} {output} &> {log}"
