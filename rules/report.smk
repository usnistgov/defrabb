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
    shell:
        """
        cat {input} \
            | '{sum+=$3-$2} END {print sum} \
            1> {output} 2> {log}
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
    shell:
        "bedtools summary -i {bed} -g {genome}"


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
    shell:
        """
        multiIntersectBed -header -i
        """


## Genome Coverage
## TODO Modify for benchmark set development framework
# rule make_bench_cov_tbls:
#     input:
#         a=ensembl_dir + "/{ref}_mrg_full_{region}.bed",
#         b=benchdir + "/HG002_{ref}_{benchmarkset}_{benchtype}.bed"
#     output: "workflow/results/bench_cov_tbls/HG002_{ref}_{benchmarkset}_{benchtype}_{region}_cov.tsv"
#     threads: 2
#     wrapper: "0.74.0/bio/bedtools/coveragebed"

## TODO Modify for benchmark set development framework
# rule make_strat_cov_tbls:
#     input:
#         a= ensembl_dir + "/{ref}_mrg_full_{region}.bed",
#         b= "workflow/data/strats/{ref}_{strat}.bed.gz"
#     output: "workflow/results/strat_cov_tbls/{strat}_{ref}_mrg_{region}_cov.tsv"
#     threads: 2
#     wrapper: "0.74.0/bio/bedtools/coveragebed"


## Variant Callset Stats
rule get_vcf_stats:
    input:
        "{prefix}.dip.vcf.gz",
    output:
        report(
            "{prefix}_stats.txt",
            caption="../report/vcf_stats.rst",
            category="bcfstats Stats",
        ),
    log:
        "logs/get_vcf_stats/{prefix}_stats.txt",
    conda:
        "../envs/bcftools.yml"
    shell:
        " bcftools stats {input} > {output} "
