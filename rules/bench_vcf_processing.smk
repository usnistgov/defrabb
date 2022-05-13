# process T2TXY_v2.7.dip.vcf to match hifiDV GT using JZ sed command `


rule fix_XY_genotype:
    input:
        "results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.fix_XY_genotype.vcf.gz",
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/fix_XY_genotype/{bench_id}_{prefix}.log",
    shell:
        """
        gunzip -c {input} \
            | sed 's/\.|1/1\/1/;s/1|\./1\/1/' \
            | bgzip -c \
            > {output}
        """


## Use when evaluating assembly accuracy with established benchmark
##  This rule is an artifact from previous assembly benchmarking pipeline
##  Not current used - keeping here for potential future use
rule dip_gap2homvarbutfiltered:
    input:
        "results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.gap2homvarbutfiltered.vcf.gz",
    # bgzip is part of samtools, which is part of the dipcall env
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/dip_gap2homvarbutfiltered/{bench_id}_{prefix}.log",
    shell:
        """
        gunzip -c {input} |\
        sed 's/1|\./1|1/' |\
        grep -v 'HET\|GAP1\|DIP' |\
        bgzip -c > {output}
        """


## Primarily for SVs
rule split_multiallelic_sites:
    input:
        "results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.vcf.gz",
    output:
        vcf="results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.split_multi.vcf.gz",
        vcf_tbi="results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.split_multi.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/split_multiallelic_sites/{bench_id}_{prefix}.log",
    shell:
        """
        bcftools norm -m - {input} -Oz -o {output.vcf} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


# Split multi-allelic variants, left-align/normalize, remove duplicates, and
# filter all lines that have REF or ALT > 20bp and no '*' characters. If this
# isn't done, svwiden will choke on commas and star characters
rule normalize_for_svwiden:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
        ref="resources/references/{ref_id}.fa",
    output:
        "results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.gt19_norm.vcf.gz",
    resources:
        mem_mb=8000
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/gt19_norm/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    shell:
        """
<<<<<<< HEAD
        bcftools norm -m- -Ou {input.vcf} \
            | bcftools norm -d exact -Ou \
            | bcftools norm -cs -f {input.ref} -Ov\
            | awk '($4!="*" && $5!="*" && (length($4)>20 || length($5)>20)) || $1~/^#/' \
            | bcftools sort -m{resources.mem_mb}m -Oz > {output} 2> {log}
=======
        bcftools norm -m- {input.vcf} \
            | bcftools norm -f {input.ref} \
            | bcftools norm -d none \
            | awk '($4!="*" && $5!="*" && (length($4)>20 || length($5)>20)) || $1~/^#/' \
            | bcftools view -o {output} -Oz - > {log}
>>>>>>> 06cb5f1... runnable rules
        """


rule run_svwiden:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.gt19_norm.vcf.gz",
        ref="resources/references/{ref_id}.fa",
    output:
<<<<<<< HEAD
        vcf="results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.svwiden.vcf.gz",
    log:
        "logs/svwiden/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/svanalyzer.yml"
    shadow:
        "minimal"
    params:
        prefix="results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.svwiden",
=======
        vcf="results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.svwiden.vcf.gz",
    log:
        "logs/svwiden/{bench_id}_{prefix}.log",
    conda:
        "../envs/svanalyzer.yml"
    shadow:
        "minimal"
    params:
        prefix="results/draft_benchmarksets/{bench_id}/intermediates/{prefix}.svwiden",
>>>>>>> 06cb5f1... runnable rules
    shell:
        """
        svanalyzer widen \
        --variants {input.vcf} \
        --ref {input.ref} \
        --prefix {params.prefix} &> {log} 

<<<<<<< HEAD
        # Removing ".;" at beginning of INFO field introduced by SVwiden
        sed 's/\.;REPTYPE/REPTYPE/' {params.prefix}.vcf \
            | bgzip -c > {params.prefix}.vcf.gz 2>> {log}
=======
        bgzip {params.prefix}.vcf
>>>>>>> 06cb5f1... runnable rules
        """


rule move_asm_vcf_to_draft_bench:
    input:
        lambda wildcards: f"results/asm_varcalls/{bench_tbl.loc[wildcards.bench_id, 'vc_id']}/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}.dip.vcf.gz",
    output:
        "results/draft_benchmarksets/{bench_id}/intermediates/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/process_benchmark_vcf/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    shell:
        "cp {input} {output} &> {log}"


rule move_processed_draft_bench_vcf:
    input:
        get_processed_vcf,
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/move_processed_draft_bench_vcf/{bench_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    shell:
        "cp {input} {output}"
