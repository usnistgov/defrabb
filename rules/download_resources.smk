################################################################################
################################################################################
##
## Downloading Resource Data Files
##
################################################################################
################################################################################


# Get and prepare assemblies
rule get_assemblies:
    output:
        "resources/assemblies/{asm_id}/{haplotype}.fa",
    log:
        "logs/get_assemblies/{asm_id}_{haplotype}.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: asm_config[wildcards.asm_id][wildcards.haplotype],
    shell:
        """
        bash scripts/fetch_resource.sh {params.url} {output} &> {log}
        """


# Get and prepare reference
rule get_ref:
    output:
        "resources/references/{ref_id}.fa",
    log:
        "logs/get_ref/{ref_id}.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: ref_config[wildcards.ref_id]["ref_url"],
    shell:
        """
        bash scripts/fetch_resource.sh {params.url} {output} &> {log}
        """


################################################################################
# Get stratifications


rule get_strats:
    output:
        "resources/strats/{ref_id}/{strat_id}.tar.gz",
    log:
        "logs/get_strats/{ref_id}_{strat_id}.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: f"{config['references'][wildcards.ref_id]['stratifications']['url']}",
    shell:
        "bash scripts/fetch_resource.sh {params.url} {output} &> {log}"


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations


rule get_comparison_vcf:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_vcf.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "vcf_url"
        ],
    shell:
        "bash scripts/fetch_resource.sh {params.url} {output} &> {log}"


rule get_comparison_bed:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_bed.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "bed_url"
        ],
    shell:
        "bash scripts/fetch_resource.sh {params.url} {output} &> {log}"
