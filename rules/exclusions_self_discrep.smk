## Excluding self comparison


def get_self_discrep_fpfns_bed(wildcards):
    """Route self-discrepancy FP/FN BED to the correct intermediate file.

    stvar benchmarks use the truvari bench path (SV-aware, #193);
    smvar benchmarks use the hap.py/vcfeval path.
    """
    base = (
        f"results/draft_benchmarksets/{wildcards.bench_id}/exclusions/self-discrep/"
        f"{wildcards.ref_id}_{wildcards.asm_id}_{wildcards.bench_type}"
        f"_{wildcards.vc_cmd}-{wildcards.vc_param_id}"
    )
    if wildcards.bench_type == "stvar":
        return f"{base}_truvari.fpfns.bed"
    return f"{base}.fpfns.bed"


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
        bed=get_self_discrep_fpfns_bed,
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


## SV-aware self-discrepancy via truvari bench (#193)
## hap.py/vcfeval is fundamentally wrong for SVs: it drops symbolic/BND alleles
## before running and cannot evaluate large insertion/deletion concordance.
## For stvar benchmarks, truvari bench (VCF vs itself) detects discordant SV
## records in the same window that cannot be mutually matched, producing FP/FN
## VCFs that drive the same intersect+slop exclusion pipeline.
##
## Parameter rationale:
##   --passonly    : restrict to PASS variants (same as smvar path)
##   --sizemin 50  : SV threshold — sub-50bp variants are handled by smvar path
##   -p 0          : disable pctseq (symbolic alleles have no sequence to compare)
##   -P 0.7        : require 70% size overlap for a match
##   -B -1         : disable BND distance matching (no BNDs in SV benchmark,
##                   per stvar_v5 profile reasoning in resources.yml #194)
##   -r 2000       : reference window (consistent with stvar_v5 / default profile)


rule self_discrep_truvari:
    input:
        vcf=get_standardized_vcf,
        vcfidx=get_standardized_vcfidx,
        bed=get_standardized_bed,
        genome=get_ref_file,
        genomeidx=get_ref_index,
    output:
        multiext(
            "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_truvari/",
            "fn.vcf.gz",
            "fn.vcf.gz.tbi",
            "fp.vcf.gz",
            "fp.vcf.gz.tbi",
            "tp-base.vcf.gz",
            "tp-base.vcf.gz.tbi",
            "tp-comp.vcf.gz",
            "tp-comp.vcf.gz.tbi",
            "summary.json",
            "candidate.refine.bed",
        ),
    log:
        "logs/exclusions/self-discrep-truvari/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_self-discrep-truvari_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    wildcard_constraints:
        bench_type="stvar",
    conda:
        "../envs/truvari_core.yml"
    params:
        dir=lambda wildcards, output: Path(output[0]).parent,
        tmpdir=lambda wildcards: (
            f"truvari_sd_{wildcards.bench_id}_{wildcards.ref_id}"
            f"_{wildcards.asm_id}_{wildcards.bench_type}"
            f"_{wildcards.vc_cmd}-{wildcards.vc_param_id}"
        ),
    shell:
        """
        rm -rf {params.tmpdir}
        truvari bench \
            -b {input.vcf} \
            -c {input.vcf} \
            -o {params.tmpdir} \
            -f {input.genome} \
            --includebed {input.bed} \
            --passonly \
            --sizemin 50 \
            -p 0 \
            -P 0.7 \
            -B -1 \
            -r 2000 \
        &> {log}
        mv {params.tmpdir}/* {params.dir}
        rm -r {params.tmpdir}
        """


rule self_discrep_truvari_extract_fpfns:
    input:
        fn_vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_truvari/fn.vcf.gz",
        fp_vcf="results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_truvari/fp.vcf.gz",
        faidx=get_ref_index,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/self-discrep/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_truvari.fpfns.bed",
    log:
        "logs/exclusions/self-discrep-truvari-fpfns/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    wildcard_constraints:
        bench_type="stvar",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        echo "Extracting FP/FN regions from truvari self-comparison" >> {log}
        (
            bcftools query -f "%CHROM\t%POS0\t%END\n" {input.fn_vcf} 2>> {log}
            bcftools query -f "%CHROM\t%POS0\t%END\n" {input.fp_vcf} 2>> {log}
        ) \
            | bedtools merge -i - 2>> {log} \
            | bedtools sort -faidx {input.faidx} -i - 1> {output} 2>> {log}
        echo "Completed successfully" >> {log}
        """
