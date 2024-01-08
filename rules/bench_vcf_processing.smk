# process T2TXY_v2.7.dip.vcf to match hifiDV GT using JZ sed command `
rule fix_XY_genotype:
    input:
        "results/asm_varcalls/{vc_id}/{prefix}.vcf.gz",
    output:
        "results/asm_varcalls/{vc_id}/{prefix}.fix_XY_genotype.vcf.gz",
    conda:
        "../envs/dipcall.yml"
    log:
        "logs/fix_XY_genotype/{vc_id}_{prefix}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        """
        gunzip -c {input} \
            | sed 's/\.|1:0,1/1:1/;s/1|\.:0,1/1:1/'  \
            | bgzip -c \
            1> {output} 2> {log}
        """


## Use when evaluating assembly accuracy with established benchmark
##  This rule is an artifact from previous assembly benchmarking pipeline
##  Not current used - keeping here for potential future use
rule dip_gap2homvarbutfiltered:
    input:
        "results/asm_varcalls/{vc_id}/{prefix}.vcf.gz",
    output:
        "results/asm_varcalls/{vc_id}/{prefix}.gap2homvarbutfiltered.vcf.gz",
    # bgzip is part of samtools, which is part of the dipcall env
    conda:
        "../envs/download_remotes.yml"
    log:
        "logs/dip_gap2homvarbutfiltered/{vc_id}_{prefix}.log",
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
        "results/asm_varcalls/{vc_id}/{prefix}.vcf.gz",
    output:
        vcf="results/asm_varcalls/{vc_id}/{prefix}.split_multi.vcf.gz",
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/split_multiallelic_sites/{vc_id}_{prefix}.log",
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
    resources:
        mem_mb=8000,
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/gt19_norm/{vc_id}_{prefix}.log",
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
    resources:
        mem_mb=8000,
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/normalize_vars/{vc_id}_{prefix}.log",
    shell:
        """
        bcftools norm -m-any -Ou {input.vcf} 2> {log} \
            | bcftools norm -d exact -Ou 2>> {log} \
            | bcftools norm -cs -f {input.ref} -Ov 2>> {log} \
            | bcftools sort -m{resources.mem_mb}m -Oz > {output} 2>> {log}
        """


## Using Adotto as tr catalogue for SV annotations - currently only for GRCh38
## https://github.com/ACEnglish/adotto
rule get_adotto_tr_anno_db:
    output:
        adotto_db="resources/references/{ref_id}_adotto_db.bed.gz",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=get_addoto_tr_anno_db_url,
    log:
        "logs/get_addoto_tr_anno_db/{ref_id}.log",
    shell:
        """
        curl -L {params.url} 1> {output.adotto_db} 2> {log}
        """


rule make_db_for_truvari_anno_trf:
    input:
        adotto_db="resources/references/{ref_id}_adotto_db.bed.gz",
        genome=get_genome_file,
    output:
        trfdb="resources/references/{ref_id}_adotto_trf.bed.gz",
        trfdbtbi="resources/references/{ref_id}_adotto_trf.bed.gz.tbi",
    log:
        "logs/make_db_for_truvari_anno_trf/{ref_id}.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        echo "Getting number of columns in {input.adotto_db}" >{log}
        ## Using python one-liner as using `zcat {input.adotto_db} | head -n 1 | awk -v FS='\\t' '{{print NF}}'
        ## - causes a pipe error (captured using `trap '' PIPE`), using `set +e` allowed the rule to run
        ## - python one-liner avoid this error and avoids having to use `set +e`
        last_col=$(awk -v FS='\\t' 'NR==1 {{print NF; exit}}' <(gzip -dc {input.adotto_db}))
        echo "Number of columns $last_col" >> {log}
        zcat {input.adotto_db} \
            | cut -f1-3,${{last_col}} \
            | bedtools sort -i stdin -g {input.genome} \
            | bgzip 1> {output.trfdb} 2>>{log}
        tabix {output.trfdb} 2>>{log}
    """


rule run_truvari_anno_trf:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
        ref=get_ref_file,
        trdb=get_ref_trdb,
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.trfanno.vcf",
    log:
        "logs/truvari_anno_trf/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari.yml"
    params:
        min_length=20,
    threads: 5
    shell:
        """
        truvari anno trf \
            -i {input.vcf} \
            -o {output.vcf} \
            -r {input.trdb} \
            -f {input.ref} \
            -t {threads} \
            --min-length {params.min_length} \
            -e trf &> {log}
        """


rule run_truvari_anno_svinfo:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.svinfo.vcf",
    log:
        "logs/truvari_anno_svinfo/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari.yml"
    params:
        minsize=20,
    shell:
        """
        truvari anno svinfo \
            -o {output.vcf} \
            --minsize {params.minsize} \
            {input.vcf} \
            &> {log}
        """


rule run_truvari_anno_repmask:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.repmask.vcf",
    log:
        "logs/truvari_anno_repmask/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari.yml"
    threads: 5
    params:
        min_length=20,
    shell:
        """
        truvari anno repmask \
            -i {input.vcf} \
            -o {output.vcf} \
            -e RepeatMasker \
            --min-length {params.min_length} \
            -T {threads} \
            --debug \
            2>&1 | tee -a {log}
        """


rule run_truvari_anno_remap:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
        ref=get_ref_file,
        refidx=get_ref_index,
        refbwaidx=get_ref_bwaindex,
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.remap.vcf",
    log:
        "logs/truvari_anno_remap/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari_remap.yml"
    params:
        min_length=20,
    shell:
        """
        truvari anno remap \
            -r {input.ref} \
            -o {output.vcf} \
            --minlength {params.min_length} \
            --debug \
            {input.vcf} \
            &> {log}
        """


rule run_truvari_anno_lcr:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.lcr.vcf",
    log:
        "logs/truvari_anno_lcr/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari.yml"
    shell:
        """
        truvari anno lcr \
            -o {output.vcf} \
            {input.vcf} \
            &> {log}
        """


rule copy_asm_vcf_to_annotations:
    input:
        "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
    output:
        "results/asm_varcalls/{vc_id}/annotations/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/copy_asm_vcf/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/download_remotes.yml"
    shell:
        "cp {input} {output} &> {log}"


rule rename_and_move_processed_draft_bench_vcf:
    input:
        get_processed_vcf,
    output:
        "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
    log:
        "logs/rename_and_move_processed_draft_bench_vcf/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bcftools.yml"
    params:
        get_sample_id,
    shell:
        """
        echo "syndip {params}\n" | \
            bcftools reheader \
            -s -\
            -o {output} \
            {input} \
            &> {log}
        """
