# process T2TXY_v2.7.dip.vcf to match hifiDV GT using JZ sed command
rule fix_XY_genotype:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        par_bed=get_par_bed,
        genome=get_genome_file,
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.fix_XY_genotype.vcf.gz",
        tbi="results/asm_varcalls/{vc_id}/annotations/{prefix}.fix_XY_genotype.vcf.gz.tbi",
    log:
        "logs/fix_XY_genotype/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        bash scripts/fix_xy_gt.sh \
            -i {input.vcf} -o {output.vcf} \
            -p {input.par_bed} \
            -g {input.genome} \
            &> {log}
        """


## Use when evaluating assembly accuracy with established benchmark
##  This rule is an artifact from previous assembly benchmarking pipeline
##  Not current used - keeping here for potential future use
rule dip_gap2homvarbutfiltered:
    input:
        "results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
    output:
        "results/asm_varcalls/{vc_id}/annotations/{prefix}.gap2homvarbutfiltered.vcf.gz",
    log:
        "logs/dip_gap2homvarbutfiltered/{vc_id}_{prefix}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        """
        gunzip -c {input} |\
        sed 's/1|\./1|1/' |\
        grep -v 'HET\|GAP1\|DIP' |\
        bgzip -c 1> {output} 2> {log}
        """


## Primarily for SVs
rule split_multiallelic_sites:
    input:
        "results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.split_multi.vcf.gz",
    log:
        "logs/split_multiallelic_sites/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools norm -m-any {input} -Oz -o {output.vcf} &> {log}
        """


# Split multi-allelic variants, left-align/normalize, remove duplicates, and
# filter all lines that have REF or ALT > 20bp and no '*' characters. If this
# isn't done, svwiden will choke on commas and star characters
rule filter_lt19_and_norm:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        ref=get_ref_file,
    output:
        "results/asm_varcalls/{vc_id}/annotations/{prefix}.gt19_norm.vcf.gz",
    log:
        "logs/gt19_norm/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    resources:
        mem_mb=8000,
    shell:
        """
        bcftools norm -m-any -Ou {input.vcf} 2> {log} \
            | bcftools norm -d exact -Ou 2>> {log} \
            | bcftools norm -cs -f {input.ref} -Ov 2>> {log} \
            | awk '($4!="*" && $5!="*" && (length($4)>20 || length($5)>20)) || $1~/^#/' \
            | bcftools sort -m{resources.mem_mb}m -Oz > {output} 2>> {log}
        """


# Split multi-allelic variants, left-align/normalize, and remove duplicates.
rule normalize_vars:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        ref=get_ref_file,
    output:
        "results/asm_varcalls/{vc_id}/annotations/{prefix}.norm.vcf.gz",
    log:
        "logs/normalize_vars/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    resources:
        mem_mb=8000,
    shell:
        """
        bcftools norm -m-any -Ou {input.vcf} 2> {log} \
            | bcftools norm -d exact -Ou 2>> {log} \
            | bcftools norm -cs -f {input.ref} -Ov 2>> {log} \
            | bcftools sort -m{resources.mem_mb}m -Oz > {output} 2>> {log}
        """


## Adding END INFO tag missing from PAV callsets
rule add_end_info_header:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.end_info.vcf.gz",
    log:
        "logs/add_end_info/{vc_id}_{prefix}.log",
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t END &> {log}
        """
