rule exclude_pav_inversions:
    ## PAV inversions (symbolic <INV>) plus overlapping segdups, for dipcall or
    ## PAV benchmarks; calls come from the PAV run for this ref + assembly.
    ## Mirrors HG002Q100-pav-inversions (pav_HG002_INV_segdupexpanded_slop50).
    input:
        unpack(get_pav_discrep_inputs),
        genome=get_genome_file,
        segdups=get_segdups,
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_pav-inv.bed",
    log:
        "logs/exclusions/exclude_pav_inversions/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_pav-inv_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bcftools_and_bedtools.yml"
    params:
        slop=config["_exclusion_params"]["pav_inv_slop"],
        merge_d=config["_exclusion_params"]["pav_inv_merge_dist"],
        inv_bed=lambda wildcards, output: f"{output.bed}.inv.tmp",
    shell:
        """
        bcftools filter -i 'ALT="<INV>"' {input.pav_vcf} 2> {log} \
            | bcftools query -f '%CHROM\t%POS0\t%END\n' 2>> {log} \
            | bedtools slop -b {params.slop} -i stdin -g {input.genome} 2>> {log} \
            | bedtools sort -g {input.genome} -i stdin 1> {params.inv_bed} 2>> {log}
        bedtools intersect -wa -a {input.segdups} -b {params.inv_bed} 2>> {log} \
            | cut -f1-3 \
            | cat - {params.inv_bed} \
            | bedtools sort -g {input.genome} -i stdin \
            | mergeBed -i stdin -d {params.merge_d} \
            1> {output.bed} 2>> {log}
        rm {params.inv_bed}
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
    conda:
        "../envs/bedtools.yml"
    params:
        slop=get_slop_value,
        slop_flags=get_slop_flags,
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
    conda:
        "../envs/bedtools.yml"
    params:
        slop=get_slop_value,
        slop_flags=get_slop_flags,
        dist=get_merge_dist,
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} \
            | bedtools slop -i stdin -g {input.genome} -b {params.slop} {params.slop_flags} \
            | bedtools merge -i stdin -d {params.dist} \
            1> {output} 2> {log}
        """


## Profile-aware slop rules for agnostic exclusions (gaps, satellites, etc.)
## These produce profile-scoped paths so different slop values get distinct files.
rule add_slop_with_profile:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{exclusion_profile}/{genomic_region}_slop.bed",
    log:
        "logs/exclusions/{ref_id}_{exclusion_profile}_{genomic_region}_slop.log",
    conda:
        "../envs/bedtools.yml"
    params:
        slop=lambda wildcards: (
            config.get("_exclusion_profiles", {})
            .get(wildcards.exclusion_profile, config["_exclusion_params"])
            .get("slop", config["_exclusion_params"]["slop"])
        ),
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} |
            bedtools slop -i stdin -g {input.genome} -b {params.slop} \
            1> {output} 2> {log}
        """


rule add_slop_and_merge_with_profile:
    input:
        bed="resources/exclusions/{ref_id}/{genomic_region}.bed",
        genome=get_genome_file,
    output:
        "resources/exclusions/{ref_id}/{exclusion_profile}/{genomic_region}_slopmerge.bed",
    log:
        "logs/exclusions/{ref_id}_{exclusion_profile}_{genomic_region}_slopmerge.log",
    conda:
        "../envs/bedtools.yml"
    params:
        slop=lambda wildcards: (
            config.get("_exclusion_profiles", {})
            .get(wildcards.exclusion_profile, config["_exclusion_params"])
            .get("slop", config["_exclusion_params"]["slop"])
        ),
        dist=lambda wildcards: (
            config.get("_exclusion_profiles", {})
            .get(wildcards.exclusion_profile, config["_exclusion_params"])
            .get("slopmerge_dist", config["_exclusion_params"]["slopmerge_dist"])
        ),
    shell:
        """
        bedtools sort -i {input.bed} -g {input.genome} \
            | bedtools slop -i stdin -g {input.genome} -b {params.slop} \
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
    log:
        "logs/exclusions/{bench_id}_flanks_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    params:
        bases=lambda wildcards: _bench_profile_param(wildcards, "flank_bases"),
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
    log:
        "logs/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_subtract_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        python scripts/subtract_exclusions.py \
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
    log:
        "logs/exclusion-intersect/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bedtools.yml"
    params:
        intersect_dir="results/draft_benchmarksets/{bench_id}/exclusions/intersections/",
    shell:
        """
        python scripts/exclusion_intersection_summary.py {input.baseline_bed} {output.summary_table} {params.intersect_dir} {input.exclusions} &> {log}
        """


rule write_exclusion_provenance:
    input:
        intersection_summary="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_intersection_summary.csv",
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.exclusion_provenance.yml",
    log:
        "logs/exclusions/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_provenance.log",
    conda:
        "../envs/bedtools.yml"
    params:
        exclusion_set_id=get_bench_exclusion_set_id,
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
