# process T2TXY_v2.7.dip.vcf to match hifiDV GT using JZ sed command `
rule fix_XY_genotype:
    input:
        "results/{prefix}.vcf.gz",
    output:
        "results/{prefix}.fix_XY_genotype.vcf.gz",
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/fix_XY_genotype/{prefix}.log",
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
        "results/{prefix}.vcf.gz",
    output:
        "results/{prefix}.gap2homvarbutfiltered.vcf.gz",
    # bgzip is part of samtools, which is part of the dipcall env
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/dip_gap2homvarbutfiltered/{prefix}.log",
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
        "results/{prefix}.vcf.gz",
    output:
        vcf="results/{prefix}.split_multi.vcf.gz",
        vcf_tbi="results/{prefix}.split_multi.vcf.gz.tbi",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/split_multiallelic_sites/{prefix}.log",
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
        vcf="results/{component_dir}/intermediates/{prefix}.vcf.gz",
        ref=get_ref_file,
    output:
        "results/{component_dir}/intermediates/{prefix}.gt19_norm.vcf.gz",
    resources:
        mem_mb=8000,
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/gt19_norm/{component_dir}/{prefix}.log",
    shell:
        """
        bcftools norm -m- -Ou {input.vcf} \
            | bcftools norm -d exact -Ou \
            | bcftools norm -cs -f {input.ref} -Ov\
            | awk '($4!="*" && $5!="*" && (length($4)>20 || length($5)>20)) || $1~/^#/' \
            | bcftools sort -m{resources.mem_mb}m -Oz > {output} 2> {log}
        """


rule run_svwiden:
    input:
        vcf="results/{prefix}.gt19_norm.vcf.gz",
        ref=get_ref_file,
    output:
        vcf="results/{prefix}.svwiden.vcf.gz",
    log:
        "logs/svwiden/{prefix}.log",
    conda:
        "../envs/svanalyzer.yml"
    shadow:
        "minimal"
    params:
        prefix="results/{prefix}.svwiden",
    shell:
        """
        svanalyzer widen \
            --variants {input.vcf} \
            --ref {input.ref} \
            --prefix {params.prefix} &> {log} 

        # Removing ".;" at beginning of INFO field introduced by SVwiden
        sed 's/\.;REPTYPE/REPTYPE/' {params.prefix}.vcf \
            | bgzip -c > {params.prefix}.vcf.gz 2>> {log}
        """
