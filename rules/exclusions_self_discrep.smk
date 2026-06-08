## Excluding self comparison
rule self_discrep_filter_symbolic:
    # hap.py / vcfeval cannot parse symbolic or breakend ALT alleles (e.g. PAV
    # inversions represented as <INV>), which crashes the self-comparison with
    # "Invalid allele ... <INV>" (#192). Drop those records before hap.py. They
    # cannot participate in vcfeval haplotype comparison anyway; SV-aware
    # self-discrepancy handling is tracked separately in #193.
    input:
        vcf=get_standardized_vcf,
        vcfidx=get_standardized_vcfidx,
    output:
        # Index is produced by the generic `tabix` rule to avoid a rule clash.
        vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.no-symbolic.vcf.gz",
    log:
        "logs/exclusions/self-discrep-filter-symbolic/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        bcftools view -e 'ALT~"<" || ALT~"\\[" || ALT~"\\]"' -Oz -o {output.vcf} {input.vcf} 2> {log}
        """


rule self_discrep_happy:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.no-symbolic.vcf.gz",
        vcfidx="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.no-symbolic.vcf.gz.tbi",
        bed=get_standardized_bed,
        ref=get_ref_file,
        sdf=get_ref_sdf,
    output:
        multiext(
            "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
            ".runinfo.json",
            ".vcf.gz",
            ".vcf.gz.tbi",
            ".summary.csv",
            ".extended.csv",
        ),
    log:
        "logs/exclusions/self-discrep-happy/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/self-discrep-happy/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/happy.yml"
    threads: config["_happy_threads"]
    resources:
        mem_mb=config["_happy_mem"],
    params:
        ## Fix to not hard code and use output prefix
        prefix="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}",
        engine="vcfeval",
        gender=get_happy_gender_param,
    shell:
        """
        hap.py \
            {input.vcf} \
            {input.vcf} \
            -R {input.bed}  \
            -r {input.ref}  \
            -o {params.prefix} \
            {params.gender} \
            --pass-only \
            --no-roc \
            --no-json \
            --engine=vcfeval \
            --engine-vcfeval-template {input.sdf} \
            --threads={threads} \
            &> {log}
        """


rule self_discrep_extract_fpfns:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
        faidx=get_ref_index,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.fpfns.bed",
    log:
        "logs/exclusions/self-discrep-fpfns/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    params:
        max_indel=config["_exclusion_params"]["self_discrep_max_indel"],
    shell:
        """
        echo "Filtering VCF for indels <={params.max_indel} and extracting FP/FN regions" >> {log}
        bcftools filter \
            --include 'ABS(ILEN)<={params.max_indel} && (FMT/BD=="FN" || FMT/BD=="FP")' {input.vcf} 2>> {log} \
                | bcftools query -f "%CHROM\t%POS0\t%END\n" 2>> {log} \
                | bedtools merge -i - 2>> {log} \
                | bedtools sort -faidx {input.faidx} -i - 1> {output} 2>> {log}
        echo "Completed successfully" >> {log}
        """


rule self_discrep_intersect_slop:
    input:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.fpfns.bed",
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_self-discrep.bed",
    log:
        "logs/exclusions/self-discrep-intersect/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_self-discrep_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    params:
        slop=config["_exclusion_params"]["self_discrep_slop"],
        merge_d=config["_exclusion_params"]["self_discrep_merge_dist"],
    shell:
        """
        ## TODO make slop conditional, don't add for intersected vars
        bedtools intersect -wa \
                -a {input.simple_repeat_bed} \
                -b {input.bed} \
            | bedtools multiinter -i stdin {input.bed} \
            | bedtools slop -b {params.slop} -i stdin -g {input.genome} \
            | mergeBed -i stdin -d {params.merge_d} \
            1> {output} 2> {log}
        """
