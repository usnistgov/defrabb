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
    params:
        url=lambda wildcards: asm_config[wildcards.asm_id][wildcards.haplotype],
    log:
        "logs/get_assemblies/{asm_id}_{haplotype}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        """
        curl -f -L {params.url} 2> {log} | gunzip -c 1> {output} 2>> {log};
        """


# Get and prepare reference
rule get_ref:
    output:
        "resources/references/{ref_id}.fa",
    params:
        url=lambda wildcards: ref_config[wildcards.ref_id]["ref_url"],
    log:
        "logs/get_ref/{ref_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        """
        curl -f --connect-timeout 120 -L {params.url} 2> {log} \
            | gunzip -c 1> {output} 2>> {log}
        """


################################################################################
# Get stratifications


rule get_strats:
    output:
        "resources/strats/{ref_id}/{strat_id}.tar.gz",
    params:
        url=lambda wildcards: f"{config['references'][wildcards.ref_id]['stratifications']['url']}",
    log:
        "logs/get_strats/{ref_id}_{strat_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


################################################################################
# Get vcf and bed files used in draft benchmark set evaluations


rule get_comparison_vcf:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "vcf_url"
        ],
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_vcf.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "curl -f -L -o {output} {params.url} &> {log}"


rule get_comparison_bed:
    output:
        "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",
    params:
        url=lambda wildcards: comp_config[wildcards.ref_id][wildcards.comp_id][
            "bed_url"
        ],
    log:
        "logs/get_comparisons/{ref_id}_{comp_id}_bed.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        """
        if [[ "{params.url}" == *gz ]]; then
            curl -f -L {params.url} |
                gunzip - 1> {output} 2> {log}
        else
            curl -f -L -o {output} {params.url} &> {log}
        fi
    """
