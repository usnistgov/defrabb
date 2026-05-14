rule exclude_pav_inversions:
    ## TODO fix to work with dipcall as well
    input:
        vcf=get_standardized_vcf,
        vcfidx=get_standardized_vcfidx,
        genome=get_genome_file,
        segdups=get_segdups,
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_pav-inv.bed",
    log:
        "logs/exclusions/exclude_pav_inversions/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    params:
        slop=config["_exclusion_params"]["pav_inv_slop"],
        merge_d=config["_exclusion_params"]["pav_inv_merge_dist"],
    benchmark:
        "benchmark/exclusions/{bench_id}_pav-inv_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        bcftools filter -i 'ALT="<INV>"' {input.vcf} \
            | bcftools query -f '%CHROM\t%POS0\t%END\n' \
            | bedtools slop -b {params.slop} -i stdin -g {input.genome} \
            | bedtools sort -g {input.genome} -i - \
            | bedtools multiinter -i stdin {input.segdups} \
            | mergeBed -i stdin -d {params.merge_d} \
            1> {output.bed} 2>{log}
        """


## Expanding exclusion regions by configured slop (bp or pct of interval length)
rule add_slop:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slop.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}_slop.log",
    params:
        slop=get_slop_value,
        slop_flags=get_slop_flags,
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} |
            bedtools slop -i stdin -g {input.genome} -b {params.slop} {params.slop_flags} \
            1> {output} 2> {log}
        """


## Expanding exclusion regions then merging, with configurable slop and merge distance
rule add_slop_and_merge:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{genomic_region}_slopmerge.bed",
    log:
        "logs/exclusions/{ref_id}_{genomic_region}_slopmerge.log",
    params:
        slop=get_slop_value,
        slop_flags=get_slop_flags,
        dist=get_merge_dist,
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} \
            | bedtools slop -i stdin -g {input.genome} -b {params.slop} {params.slop_flags} \
            | bedtools merge -i stdin -d {params.dist} \
            1> {output} 2> {log}
        """


## Finding breaks in assemblies for excluded regions
rule intersect_start_and_end:
    input:
        baseline_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.baseline.bed",
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
        awk '{{FS=OFS="\t"}} {{print $1, $2, $2+1}}' {input.baseline_bed} \
            | bedtools intersect -u -wa -a {input.xregions} -b stdin \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.start} 2> {log}

        awk '{{FS=OFS="\t"}} {{print $1, $3-1, $3}}' {input.baseline_bed} \
            | bedtools intersect -u -wa -a {input.xregions} -b stdin  \
            | bedtools sort -g {input.genome} -i stdin \
            1> {output.end} 2>> {log}
        """


# Generate bed with 15kb regions around assembly breaks (non-diploid coverage)
rule get_flanks:
    input:
        baseline_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.baseline.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_flanks.bed",
    params:
        bases=config["_exclusion_params"]["flank_bases"],
    log:
        "logs/exclusions/{bench_id}_flanks_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools complement -i {input.baseline_bed} -g {input.genome} \
            | bedtools flank -i stdin -g {input.genome} -b {params.bases} \
            1> {output} 2> {log}
        """


## Removing excluded genomic regions from asm varcalls bed file
## for draft benchmark regions
rule subtract_exclusions:
    input:
        baseline_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.baseline.bed",
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
            {input.baseline_bed} \
            {output.bed} \
            {output.rpt} \
            {input.other_beds} \
            &> {log}
        """


rule generate_intersection_summary:
    input:
        baseline_bed=lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[(wildcards.bench_id, 'vc_id')]}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.baseline.bed",
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
        python {params.script} {input.baseline_bed} {output.summary_table} {params.intersect_dir} {input.exclusions} &> {log}
        """


rule write_exclusion_provenance:
    input:
        intersection_summary="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_intersection_summary.csv",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_provenance.yml",
    params:
        exclusion_set_id=get_bench_exclusion_set_id,
    log:
        "logs/exclusions/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_provenance.log",
    conda:
        "../envs/bedtools.yml"
    script:
        "../scripts/write_exclusion_provenance.py"


## Used when no exclusions are applied
rule postprocess_bed:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.baseline.bed",
    output:
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed",
    log:
        "logs/process_benchmark_bed/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "cp {input} {output} &> {log}"
