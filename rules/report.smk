import yaml
import os


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


## Summary stats by chromosome
rule get_bed_stats:
    input:
        bed="{prefix}.bed",
        genome=get_genome_file,
    output:
        report(
            "{prefix}_bed-summary.tsv",
            caption="../report/bed_stats.rst",
            category="Exclusions",
        ),
    log:
        "logs/get_bed/stats/{prefix}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools summary -i {input.bed} -g {input.genome} 1> {output} 2> {log}"


## Not currently being used
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
    shell:
        " bcftools stats {input} 1> {output} 2> {log}"


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
    shell:
        " rtg vcfstats --allele-lengths {input} 1> {output} 2> {log}"


## Not currently being used
rule summarize_exclusions:
    input:
        lambda wildcards: f"results/draft_benchmarksets/{{bench_id}}/{{ref_id}}_{vc_tbl.loc[(wildcards.vc_id, 'asm_id')]}_{vc_tbl.loc[(wildcards.vc_id, 'vc_cmd')]}-{vc_tbl.loc[(wildcards.vc_id, 'vc_param_id')]}.excluded_stats.txt",
    output:
        "results/report/{bench_id}/{ref_id}_{vc_id}-{exclusion_set}.html",
    conda:
        "../envs/rmd.yml"
    log:
        "logs/summarise_exclusions/{bench_id}_{ref_id}_{vc_id}-{exclusion_set}.log",
    script:
        "../scripts/reports/exclusions.Rmd"


rule write_report_params:
    input:
        asm_var_rtg=expand(
            "results/asm_varcalls/{vc_id}/{ref}_{asm_id}_{vc_cmd}-{vc_param_id}.dip_rtg_stats.txt",
            zip,
            vc_id=dipcall_tbl.index.tolist(),
            ref=dipcall_tbl["ref"].tolist(),
            asm_id=dipcall_tbl["asm_id"].tolist(),
            vc_cmd=dipcall_tbl["vc_cmd"].tolist(),
            vc_param_id=dipcall_tbl["vc_param_id"].tolist(),
        ),
        bench_var_rtg=expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_bench-vars_rtg_stats.txt",
            zip,
            bench_id=bench_tbl.index.tolist(),
            ref=bench_tbl["ref"].tolist(),
            asm_id=bench_tbl["asm_id"].tolist(),
            bench_type=bench_tbl["bench_type"].tolist(),
            vc_cmd=bench_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_tbl["vc_param_id"].tolist(),
        ),
        exclusion_summary=expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_stats.txt",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        exclusion_intersection_summary=expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_intersection_summary.csv",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        bench_cov=expand(
            "results/draft_benchmarksets/{bench_id}/{ref}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark_bed-summary.tsv",
            zip,
            bench_id=bench_excluded_tbl.index.tolist(),
            ref=bench_excluded_tbl["ref"].tolist(),
            asm_id=bench_excluded_tbl["asm_id"].tolist(),
            bench_type=bench_excluded_tbl["bench_type"].tolist(),
            vc_cmd=bench_excluded_tbl["vc_cmd"].tolist(),
            vc_param_id=bench_excluded_tbl["vc_param_id"].tolist(),
        ),
        happy_summary=expand(
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
        happy_extended=expand(
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
        truvari_summary=expand(
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
        truvari_refine_summary=expand(
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
    output:
        report_params="results/analysis_params.yml",
    params:
        param1="value1",
        param2="value2",
    log:
         "logs/write_report_params.log",
    run:
        config_data = config
        analysis_table = f"../{config['analyses']}"

        inputs_converted = {}  
        for key, value in input.items():  
            # Convert each input item to a string and ensure it's a list  
            if isinstance(value, list):  
                inputs_converted[key] = [str(v) for v in value]  
            else:  
                inputs_converted[key] = [str(value)]  

        report_params_dict = {
            "inputs": inputs_converted,
            "config_data": config_data,
            "analysis_table": analysis_table,
            "variables": dict(params),
        }

        with open(output.report_params, "w") as f:
            yaml.safe_dump(report_params_dict, f, default_flow_style=False, sort_keys=False)

rule render_report:
    input:
        report_params=Path(workflow.basedir) / "results/analysis_params.yml",
    output:
        report_html=report(
            "results/analysis.html",
            caption="../report/analysis.rst",
            category="Analysis Report",
        ),
    params:
        qmd="scripts/reports/analysis.qmd",
        results_qmd = "results/analysis.qmd",
        rundir=Path(workflow.basedir),
    log:
         "logs/render_report.log",
    conda:
        "../envs/quarto.yml"
    shell: """
            cp {params.qmd} {params.results_qmd}
            quarto render {params.results_qmd} \
                --execute-dir {params.rundir} \
                -P yaml_path:{input.report_params} \
                --log {log} --log-level info --debug
    """