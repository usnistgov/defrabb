# import pandas as pd


# wildcard_constraints:
#     ref_id="GRCh38|GRCh37|GRCh38_chr21|CHM13v2.0",
#     genomic_region="homopolymers|segdups|tandem-repeats|gaps|self-chains|satellites",


# structural variants - using asm varcalls vcf to identify structural variants for exclusion
rule get_SVs_from_vcf:
    input:
        "results/{prefix}.vcf.gz",
    output:
        bed="results/{prefix}_SVs.bed",
        tbl="results/{prefix}_SVs.tsv",
    log:
        "logs/exclusions/vc_SVs_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        ## Generating table with SV information and refwiden coordinates
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\%INFO\REFWIDENED\t%INFO\REPTYPE\n'  - > {output.tbl}

        ## Creating bed with SVs for use in defining excluded regions
        # ---- excluding SVs less than 50 bps and reformatting as 0 base tab sep bed file (CHROM\tSTART\tSTOP)
        awk 'length($3)>49 || length($4)>49' {output.tbl} \
            cut -f 5 \
            | sed 's/[:,-]/\t/g'
            awk '{{FS=OFS="\\t"}} {{print $1,$2-1,$3-1}}' \
            1> {output} 2> {log}
        """


rule intersect_SVs_and_simple_repeats:
    input:
        sv_bed="results/{prefix}_SV_sorted.bed",
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/{prefix}_svs-and-simple-repeats.bed",
    log:
        "logs/exclusions/int_SVs_{prefix}.log",
    benchmark:
        "benchmark/exclusions/{prefix}.tsv"
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


## Expanding exclusion regions by 15kb
rule add_slop:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slop.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}.bed",
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


## Finding breaks in assemblies for excluded regions
rule intersect_start_and_end:
    input:
        dip="results/{prefix}.dip_sorted.bed",
        xregions="resources/exclusions/{ref_id}/{excluded_region}.bed",
        genome=get_genome_file,
    output:
        start="results/{prefix}/exclusions/{excluded_region}_start_sorted.bed",
        end="results/{prefix}/exclusions/{excluded_region}_end_sorted.bed",
    log:
        "logs/exclusions/start_end_{excluded_region}_{prefix}.log",
    benchmark:
        "benchmark/exclusions/start_end_{excluded_region}_{prefix}.tsv"
    conda:
        "../envs/bedtools.yml"
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
        bed="results/{prefix}.dip_sorted.bed",
        genome=get_genome_file,
    output:
        "results/{prefix}_flanks.bed",
    log:
        "logs/exclusions/add_flanks)_{prefix}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools flank -i {input.bed} -g {input.genome} -b 15000 1> {output} 2> {log}"


## Removing excluded genomic regions from asm varcalls bed file
## for draft benchmark regions
rule subtract_exclusions:
    input:
        dip_bed=f"results/draft_benchmarksets/intermediates/exclusions/{excluded_bench_space.wildcard_pattern}.dip.sorted.bed",
        other_beds=get_exclusion_inputs,
    output:
        bed=f"results/draft_benchmarksets/{excluded_bench_space.wildcard_pattern}.excluded.bed",
        stats=f"results/draft_benchmarksets/{excluded_bench_space.wildcard_pattern}.excluded_stats.txt",
    params:
        script=workflow.source_path("../scripts/subtract_exclusions.py"),
    log:
        f"logs/exclusions/draft_benchmark_subtract/{excluded_bench_space.wildcard_pattern}.log",
    benchmark:
        f"benchmark/exclusions/draft_benchmark_subtract/{excluded_bench_space.wildcard_pattern}.benchmark"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        python {params.script} \
            {input.dip_bed} \
            {output.bed} \
            {output.stats} \
            {input.other_beds}
        """
