## Using Adotto as tr catalogue for SV annotations
rule get_adotto_tr_anno_db:
    output:
        adotto_db="resources/references/{ref_id}_adotto_db.bed.gz",
    log:
        "logs/get_addoto_tr_anno_db/{ref_id}.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=get_addoto_tr_anno_db_url,
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
        "../envs/truvari_trf.yml"
    threads: 5
    params:
        min_length=20,
    shell:
        """
    OPENSSL_CONF=/dev/null \
        truvari anno trf \
            -i {input.vcf} \
            -o {output.vcf} \
            -r {input.trdb} \
            -f {input.ref} \
            -t {threads} \
            --min-length {params.min_length} \
            --trf-params "3 7 7 80 5 40 500 -h -ngs -l 1" \
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
        "../envs/truvari_core.yml"
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


rule install_dfam_hmm:
    output:
        touch("installed_dfam_hmm"),
    conda:
        "../envs/truvari_repmask.yml"
    run:
        rmdir = f"{env.CONDA_PREFIX}/share/RepeatMasker"
        # fetch the curated‐only HMM
        shell(
            f"""
          cd {rmdir}/Libaries/famdb/
          wget -qO- https://www.dfam.org/releases/current/families/FamDB/dfam39_full.5.h5.gz
          gunzip dfam39_full.5.h5.gz
          cd {rmdir}
          perl ./configure
        """
        )


rule run_truvari_anno_repmask:
    input:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/annotations/{prefix}.vcf.gz.tbi",
    output:
        vcf="results/asm_varcalls/{vc_id}/annotations/{prefix}.repmask.vcf",
    log:
        "logs/truvari_anno_repmask/{vc_id}_{prefix}.log",
    conda:
        "../envs/truvari_repmask.yml"
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
        "../envs/truvari_core.yml"
    shell:
        """
        truvari anno lcr \
            -o {output.vcf} \
            {input.vcf} \
            &> {log}
        """
