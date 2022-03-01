## Assembly statistics for individual haploids
rule run_assembly_stats:
    input:
        #Input assembly
        assembly="resources/assemblies/{asm_id}/{haplotype}.fa",
    output:
        #Assembly statistics
        assembly_stats=report(
            "results/report/assemblies/{asm_id}_{haplotype}_stats.txt",
            caption="../report/asm_stats.rst",
            category="Assembly",
        ),
    params:
        # Tab delimited output, with a header, is set as the default. Other options are available:
        #   -l <int>
        #       Minimum length cutoff for each sequence.
        #       Sequences shorter than the cutoff will be ignored [1]
        #   -s
        #       Print 'grep friendly' output
        #   -t
        #       Print tab-delimited output
        #   -u
        #       Print tab-delimited output with no header line
        # If you want to add multiple options just delimit them with a space.
        # Note that you can only pick one output format
        # Check https://github.com/sanger-pathogens/assembly-stats for more details
        extra="-t",
    log:
        "logs/run_assembly_stats/{asm_id}_{haplotype}.assembly-stats.log",
    group:
        "report"
    threads: 1
    wrapper:
        "v0.86.0/bio/assembly-stats"


## Potentially modify to calculate for all output beds
rule get_bed_size:
    input:
        "{genomic_region}.bed",
    output:
        "results/bed_size/{genomic_region}.txt",
    # output: report("results/bed_size/{genomic_region}.txt", caption = "report/bed_size.rst", category = "Exclusion Stats")
    log:
        "logs/get_bed_size/{genomic_region}.log",
    group:
        "report"
    shell:
        """
        cat {input} \
            | '{sum+=$3-$2} END {print sum} \
            1> {output} 2> {log}
        """


## TODO - fix output filename and add to rule all
rule get_number_of_variants:
    input:
        vcf="{prefix}.vcf.gz",
        bed="{genomic_region}.bed",
    output:
        "{prefix}/nvars_{genomic_region}.txt",
    log:
        "logs/{prefix}_{genomic_region}_var_counts.tsv",
    conda:
        "../envs/bedtools.yml"
    group:
        "report"
    shell:
        """
        nvar=$(intersectBed \
         -a {input.vcf} \
         -b {input.bed} \
         | wc -l
        echo {input.bed} $nvar 1> {output}  2> {log}
        """


## Summary stats by chromosome, need to test/ debug getting ref_id from input bed
rule get_bed_stats:
    input:
        bed="{bed_dir}/{ref_id}_{genomic_region}.bed",
        genome="resouces/exclusions/{ref_id}.genome",
    output:
        "results/bed_stats/{bed_dir}_{ref_id}/{genomic_region}.tsv",
    # output: report("results/bed_stats/{bed_dir}_{ref_id}/{genomic_region}.tsv", caption = "report/bed_stats.rst", category = "Exclusion Stats")
    log:
        "logs/get_bed/stats/{bed_dir}_{ref_id}_{genomic_region}.log",
    conda:
        "../envs/bedtools.yml"
    group:
        "report"
    shell:
        "bedtools summary -i {input.bed} -g {input.genome} 1> {output} 2> {log}"


rule get_exclusion_coverage:
    input:
        "results/draft_benchmarksets/{draft_bench}_exclusions/{prefix}.bed",
    # output: report("results/draft_benchmarksets/{draft_bench}_exclusions/{prefix}_stats.tsv", caption = "report/exclusion_stats.rst", category = "Exclusion Stats")
    output:
        "results/draft_benchmarksets/{draft_bench}_exclusions/{prefix}_stats.tsv",
    log:
        "logs/get_exclusion_coverage/{draft_bench}_{prefix}.log",
    conda:
        "../envs/bedtools.yml"
    group:
        "report"
    shell:
        """
        multiIntersectBed -header -i
        """


## Variant Callset Stats
rule get_bcftools_stats:
    input:
        "{prefix}.vcf.gz",
    output:
        report(
            "{prefix}_bcftools_stats.txt",
            caption="../report/bcftools_stats.rst",
            category="VCF Stats",
            subcategory="bcftools stats",
        ),
    log:
        "logs/get_bcftools_stats/{prefix}_stats.txt",
    conda:
        "../envs/bcftools.yml"
    group:
        "report"
    shell:
        " bcftools stats {input} > {output} "


rule get_rtg_vcf_stats:
    input:
        "{prefix}.vcf.gz",
    output:
        report(
            "{prefix}_rtg_stats.txt",
            caption="../report/rtg_stats.rst",
            category="VCF Stats",
            subcategory="rtgtools stats",
        ),
    log:
        "logs/get_rtg_vcf_stats/{prefix}_stats.txt",
    conda:
        "../envs/rtgtools.yml"
    group:
        "report"
    shell:
        " rtg vcfstats {input} 1> {output} 2> {log}"


rule summarize_exclusions:
    input:
        lambda wildcards: f"results/draft_benchmarksets/{{bench_id}}/{{ref_id}}_{vc_tbl.loc[(wildcards.vc_id, 'asm_id')]}_{vc_tbl.loc[(wildcards.vc_id, 'vc_cmd')]}-{vc_tbl.loc[(wildcards.vc_id, 'vc_param_id')]}.excluded_stats.txt",
    output:
        "results/report/{bench_id}/{ref_id}_{vc_id}-{exclusion_set}.html",
    conda:
        "../envs/rmd.yml"
    log:
        "logs/summarise_exclusions/{bench_id}_{ref_id}_{vc_id}-{exclusion_set}.log",
    group:
        "report"
    script:
        "../scripts/reports/exclusions.Rmd"
