################################################################################
################################################################################
##
## Component: Assembly Based Variant Calling
##
################################################################################
################################################################################

## Run Dipcall


rule run_dipcall:
    input:
        h1="resources/assemblies/{asm_id}/paternal.fa",
        h2="resources/assemblies/{asm_id}/maternal.fa",
        ref="resources/references/{ref_id}.fa",
        ref_idx="resources/references/{ref_id}.fa.fai",
        ref_mmi="resources/references/{ref_id}.mmi",
        par=get_par_bed,
    output:
        make="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.mak",
        vcf="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
        bed="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.bed",
        bam1="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.hap1.bam",
        bam2="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.hap2.bam",
    conda:
        "../envs/dipcall.yml"
    params:
        prefix="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
        male_bed=get_dipcall_par_param,
        ts=config["_dipcall_threads"],
        make_jobs=config["_dipcall_jobs"],
        extra=lambda wildcards: (
            ""
            if vc_tbl.loc[wildcards.vc_id]["vc_params"] == "default"
            else vc_tbl.loc[wildcards.vc_id]["vc_params"]
        ),
    log:
        multiext(
            "results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}",
            ".hap1.paf.gz.log",
            ".hap2.paf.gz.log",
            ".hap1.sam.gz.log",
            ".hap2.sam.gz.log",
        ),
        rulelog="logs/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/asm_varcalls/{vc_id}_{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    resources:
        ## GB per make job run in parallel
        mem_mb=config["_dipcall_jobs"] * config["_dipcall_mem"],
    threads: config["_dipcall_threads"] * config["_dipcall_jobs"]
    shell:
        """
        echo "Writing Makefile defining dipcall pipeline" > {log.rulelog}
        run-dipcall \
            -t {params.ts} \
            -d {input.ref_mmi} \
            {params.extra} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            1> {output.make} 2>> {log.rulelog}

        echo "Running dipcall pipeline" >> {log.rulelog}
        make -j{params.ts} -f {output.make} &>>{log.rulelog}
        """


## TODO run_pav rule
rule run_pav:
    input:
        ref=get_ref_file,
        refidx=get_ref_index,
        hap1=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/paternal.fa",
        hap2=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/maternal.fa",
    output:
        vcf="results/asm_varcalls/{vc_id}/results/{asm_id}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/results/{asm_id}.vcf.gz.tbi",
    params:
        outdir="results/asm_varcalls/{vc_id}",
        pav_config=lambda wildcards: config["_pav_config"][
            vc_tbl.loc[wildcards.vc_id]["vc_param_id"]
        ],
        name=lambda wildcards: vc_tbl.loc[wildcards.vc_id]["asm_id"],
    container:
        "docker://becklab/pav:latest"
    threads: 8
    script: "../scripts/run_pav.py"
