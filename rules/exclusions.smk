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
    "dipcall-bugs-T2TACE",
    "HG002Q100-pav-discrep-smvar",
    "HG002Q100-pav-discrep-stvar",
    "HG002Q100-pav-inversions",
    "HG002Q100-errors",
    "HG002Q100-mosaic",
    "HG002Q100-delins-errors",
    "TSPY2-segdups",
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
        vcf=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/annotations/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.norm.trfanno.vcf",
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
        hap2_bam=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.hap1.bam"
    output:
       "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_consecutive-svs.bed",
    log:
        "benchmark/exclusions/consecutive-svs/{bench_id}_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.log"
    conda:
        "../envs/consecutive_svs.yml"
    shell:
        """
        python scripts/get_dipcall_delins.py \
            {input.hap1_bam} \
            {input.hap2_bam} \
            {output.bed} \
            &> {log}
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
