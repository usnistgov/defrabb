################################################################################
################################################################################
##
## Component: Evaluating Draft Benchmarksets
##
################################################################################
################################################################################

## Run happy


rule run_happy:
    input:
        unpack(partial(get_happy_inputs, analyses, config)),
    output:
        multiext(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
            ".runinfo.json",
            ".vcf.gz",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
        report(
            "results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.summary.csv",
            caption="report/happy_summary.rst",
            category="Happy",
        ),
    params:
        prefix="results/evaluations/happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
        strat_tsv=lambda wildcards: f"{config['references'][wildcards.ref_id]['stratifications']['tsv']}",
        threads=config["_happy_threads"],
        engine="vcfeval",
        engine_extra=lambda wildcards: f"--engine-vcfeval-template resources/references/{wildcards.ref_id}.sdf",
        gender_param=get_happy_gender_param
    resources:
        mem_mb=config["_happy_mem"],
    threads: config["_happy_threads"]
    log:
        "logs/run_happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/happy.log",
    benchmark:
        "benchmark/run_happy/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/happy.yml"
    script:
        "../scripts/run_happy.py"


################################################################################
## Run Truvari


rule run_truvari:
    input:
        unpack(partial(get_truvari_inputs, analyses, config)),
    output:
        report(
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/summary.json",
            caption="report/truvari_summary.rst",
            category="Truvari",
        ),
        multiext(
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/",
            "fn.vcf.gz",
            "fn.vcf.gz.tbi",
            "fp.vcf.gz",
            "fp.vcf.gz.tbi",
            "tp-base.vcf.gz",
            "tp-base.vcf.gz.tbi",
            "tp-comp.vcf.gz",
            "tp-comp.vcf.gz.tbi",
            "candidate.refine.bed",
        ),
    log:
        "logs/run_travari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/truvari.log",
    benchmark:
        "benchmark/run_truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    params:
        dir=lambda wildcards, output: Path(output[0]).parent,
        tmpdir=lambda wildcards: expand("truvari_{eval_id}", eval_id=wildcards.eval_id),
    conda:
        "../envs/truvari.yml"
    shell:
        """
        ## Removing temp directory before starting run
        rm -rf {params.tmpdir}

        ## Running truvari
        truvari bench \
            -b {input.truth} \
            -c {input.query} \
            -o {params.tmpdir} \
            --pick ac \
            --passonly \
            -r 2000 \
            -C 5000 \
            -f {input.genome} \
            --includebed {input.truth_regions} \
        &> {log}

        mv {params.tmpdir}/* {params.dir}
        rm -r {params.tmpdir}
        """


rule truvari_refine:
    input:
        unpack(partial(get_truvari_inputs, analyses, config)),
        refine_bed="results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/candidate.refine.bed",
        json="results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/summary.json",
        ref="resources/references/{ref_id}.fa",
        refidx="resources/references/{ref_id}.fa.fai",
    output:
        report(
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/refine.variant_summary.json",
            caption="report/truvari_refine_summary.rst",
            category="Truvari",
        ),
        multiext(
            "results/evaluations/truvari/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/",
            "phab.output.vcf.gz",
            "phab.output.vcf.gz.tbi",
            "refine.log.txt",
            "refine.region_summary.json",
            "refine.regions.txt",
        ),
    log:
        "logs/run_travari_refine/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}/truvari_refine.log",
    params:
        bench_output=lambda w, output: os.path.dirname(output[0]),
    benchmark:
        "benchmark/run_truvari_refine/{eval_id}_{bench_id}/{ref_id}_{comp_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/truvari.yml"
    threads: config["_truvari_refine_threads"]
    script:
        "../scripts/run_truvari_refine.py"
