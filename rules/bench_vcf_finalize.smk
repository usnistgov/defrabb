rule remove_inv:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.no_inv.vcf.gz",
    log:
        "logs/remove_inv/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view -e 'ALT="<INV>"' -Oz -o {output.vcf} {input.vcf} &> {log}
        """


rule copy_std_asm_vcf_to_annotations:
    input:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    output:
        "results/asm_varcalls/{vc_id}/annotations/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/copy_asm_vcf/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "cp {input} {output} &> {log}"


rule move_processed_draft_bench_vcf:
    input:
        get_processed_vcf,
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/move_processed_draft_bench_vcf/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "cp {input} {output} &> {log}"


## Filtering vcf to only include variants in benchmark regions
## - assumes sv benchmark vcfs are annotated with truvari svinfo
rule get_variants_in_benchmark_regions:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
        vcfidx="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
        bed="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.benchmark.bed",
    output:
        vcf="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_bench-vars.vcf.gz",
        vcfidx="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_bench-vars.vcf.gz.tbi",
    log:
        "logs/get_vars_in_bench_regions/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bcftools.yml"
    params:
        filt=lambda wildcards: (
            f"-i 'INFO/SVLEN > 49'" if wildcards.bench_type == "stvar" else ""
        ),
    shell:
        """
    bcftools view -Oz -o {output.vcf} \
        -R {input.bed} --regions-overlap 2 \
        {params.filt} {input.vcf} >> {log}
        bcftools index -t {output.vcf} >> {log}
    """
