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
    conda:
        "../envs/dipcall.yml"
    threads: config["_dipcall_threads"] * config["_dipcall_jobs"]
    resources:
        ## GB per make job run in parallel
        mem_mb=config["_dipcall_jobs"] * config["_dipcall_mem"],
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
        ## -j must match the jobs knob the mem_mb reservation is based on
        ## (mem_mb = _dipcall_jobs * _dipcall_mem); driving it off the thread
        ## count over-parallelizes memory-heavy minimap2 and OOMs (see
        ## docs/issues/run_pav_run_dipcall_failures.md).
        make -j{params.make_jobs} -f {output.make} &>>{log.rulelog}
        """


rule rename_dipcall_vcf_sample:
    input:
        vcf="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.vcf.gz",
    output:
        vcf="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip.rename.vcf.gz",
    log:
        "logs/rename_dipcall/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bcftools.yml"
    params:
        get_sample_id,
    shell:
        """
        echo "syndip {params}\n" \
            | bcftools reheader -s - -o {output.vcf} {input.vcf} &> {log}
        """


## Generate PAV input files (assemblies.tsv + config.json).
##
## Runs on the HOST (no container) so Snakemake's script: machinery uses the
## host Snakemake. It must NOT run inside the becklab/pav container: that
## container ships its own older Snakemake, and injecting the host-generated
## script preamble there crashes with a version skew
## (`No module named 'snakemake.io.container'`). See
## docs/issues/run_pav_run_dipcall_failures.md.
rule pav_config:
    input:
        ref=get_ref_file,
        refidx=get_ref_index,
        hap1=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/paternal.fa",
        hap2=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/maternal.fa",
    output:
        assemblies="results/asm_varcalls/{vc_id}/assemblies.tsv",
        config="results/asm_varcalls/{vc_id}/config.json",
    log:
        "logs/asm_varcalls/{vc_id}_pav_config.log",
    params:
        outdir="results/asm_varcalls/{vc_id}",
        pav_config=lambda wildcards: config["_pav_config"][
            vc_tbl.loc[wildcards.vc_id]["vc_param_id"]
        ],
        name=lambda wildcards: asm_config[vc_tbl.loc[wildcards.vc_id, "asm_id"]][
            "sample_id"
        ],
    script:
        "../scripts/setup_pav.py"


## Running the PAV assembly variant caller.
##
## A bare shell: directive (NOT a script:) — the nested PAV Snakemake must run
## with the container's own Snakemake, and Snakemake does not inject its
## script-machinery preamble for shell: rules, so there is no host/container
## version skew. Input config files are produced on the host by pav_config.
rule run_pav:
    input:
        assemblies="results/asm_varcalls/{vc_id}/assemblies.tsv",
        config="results/asm_varcalls/{vc_id}/config.json",
        ref=get_ref_file,
        refidx=get_ref_index,
        hap1=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/paternal.fa",
        hap2=lambda wildcards: f"resources/assemblies/{vc_tbl.loc[wildcards.vc_id]['asm_id']}/maternal.fa",
    output:
        vcf="results/asm_varcalls/{vc_id}/{sample_id}.vcf.gz",
        vcfidx="results/asm_varcalls/{vc_id}/{sample_id}.vcf.gz.tbi",
        h1_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h1_500.bed.gz",
        h2_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h2_500.bed.gz",
    log:
        "logs/asm_varcalls/{vc_id}_{sample_id}_pav.log",
    benchmark:
        "benchmark/run_pav/{vc_id}_{sample_id}.tsv"
    container:
        config.get("_pav_container", "docker://becklab/pav:latest")
    threads: config["_pav_threads"]
    resources:
        mem_mb=config["_pav_mem"],
    params:
        outdir="results/asm_varcalls/{vc_id}",
    shell:
        # Capture the absolute log path before cd; PAV runs in outdir.
        # On FIPS hosts, run_defrabb auto-binds /proc/sys/crypto/fips_enabled -> 0
        # inside apptainer containers to prevent the PAV pysam FIPS self-test crash.
        # See docs/issues/run_pav_fips_selftest.md.
        # The nested PAV output is captured to the rule log; on failure PAV's
        # own .snakemake logs are surfaced too so the cause is debuggable.
        """
        LOG="$(pwd)/{log}"
        cd {params.outdir}
        if ! snakemake -s /opt/pav/Snakefile --ri -k -w 20 \
            --rerun-triggers mtime -c {threads} --config ignore_env_file=True \
            > "$LOG" 2>&1; then
            echo '--- inner PAV .snakemake/log ---' >> "$LOG"
            cat .snakemake/log/*.log >> "$LOG" 2>/dev/null || true
            exit 1
        fi
        """


rule intersect_pav_callable_regions:
    input:
        h1_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h1_500.bed.gz",
        h2_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h2_500.bed.gz",
        faidx=get_ref_index,
    output:
        diploid_regions="results/asm_varcalls/{vc_id}/{sample_id}.diploid_regions.bed",
    log:
        "logs/asm_varcalls/{vc_id}/{sample_id}_intersect.log",
    conda:
        "../envs/bedtools.yml"
    threads: config.get("intersect_threads", 1)
    resources:
        mem_mb=config.get("intersect_mem_mb", 1024),
    params:
        intersect_opts=config.get("intersect_opts", ""),
    shell:
        """
        bedtools intersect {params.intersect_opts} \
            -a {input.h1_bed} -b {input.h2_bed} \
            | bedtools sort -faidx {input.faidx} -i stdin \
            1> {output.diploid_regions} 2> {log}
        """


rule standardize_vcasm_output:
    input:
        unpack(branch(is_pav, then=get_pav_outputs, otherwise=get_dipcall_outputs)),
    output:
        standardized_vcf="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
        standardized_vcfidx="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
        standardized_bed="results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.baseline.bed",
    shell:
        """
        cp {input.vcf} {output.standardized_vcf}
        cp {input.vcfidx} {output.standardized_vcfidx}
        cp {input.bed} {output.standardized_bed}
        """
