rule index_ref:
    input:
        "resources/references/{ref_id}.fa",
    output:
        "resources/references/{ref_id}.fa.fai",
    log:
        "logs/index_ref/{ref_id}.log",
    resources:
        mem_mb=16000,
    wrapper:
        "v3.13.3/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/references/{ref_id}.fa",
    output:
        idx=multiext(
            "resources/references/{ref_id}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"
        ),
    log:
        "logs/bwa_index/{ref_id}.log",
    params:
        prefix=lambda w, input: input[0],
    wrapper:
        "v3.3.3/bio/bwa/index"


rule index_ref_mmi:
    input:
        "resources/references/{ref_id}.fa",
    output:
        "resources/references/{ref_id}.mmi",
    log:
        "logs/index_ref_mmi/{ref_id}.log",
    threads: 4
    resources:
        mem_mb=2400,
    conda:
        "envs/dipcall.yml"
    shell:
        "minimap2 -x asm5 -d {output} {input} &> {log}"


rule index_ref_sdf:
    input:
        "resources/references/{ref_id}.fa",
    output:
        directory("resources/references/{ref_id}.sdf"),
    log:
        "logs/index_ref_sdf/{ref_id}.log",
    conda:
        "envs/rtgtools.yml"
    shell:
        "rtg format -o {output} {input} &>{log}"

## General indexing rule for vcfs
rule tabix:
    input:
        "{filename}.vcf.gz",
    output:
        "{filename}.vcf.gz.tbi",
    params:
        extra="-t",
    log:
        "logs/tabix/{filename}.log",
    wrapper:
        "v3.13.3/bio/bcftools/index"

rule sort_bed:
    input:
        in_file="{prefix}.bed",
        genome=get_genome_file,
    output:
        "{prefix}_sorted.bed",
    log:
        "logs/sort_bed/{prefix}.log",
    wrapper:
        "0.74.0/bio/bedtools/sort"

rule index_dip_bam:
    input:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam",
    output:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam.bai",
    log:
        "logs/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.{hap}.bam.bai.log",
    wrapper:
        "v3.13.3/bio/samtools/index"

## General rule for compressing vcfs
## - bcftools sort based rule to ensure vcf is sorted
rule compress_vcf:
    input:
        "{prefix}.vcf",
    output:
        "{prefix}.vcf.gz",
    log:
        "logs/bcftools/sort/{prefix}.log",
    params:
        extras="-Oz",
    resources:
        mem_mb=8000,
    wrapper:
        "v2.6.0/bio/bcftools/sort"