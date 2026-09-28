## PAV vs dipcall discrepancy exclusions
## Regenerates, for any assembly, the HG002Q100-pav-discrep-{smvar,stvar} beds
## that were originally curated for HG002 Q100 v1.0/v1.1 (downloaded beds):
##   smvar: hap.py/vcfeval PAV vs dipcall FP/FN (<=50bp), + overlapping repeats,
##          slop 5, merge 1000 (pav10vsbench10corrXY_fnfp_slop5_merge1000_withoverlaprepeats)
##   stvar: truvari PAV vs dipcall FP/FN (>=50bp), + overlapping repeats,
##          slop 50, merge 10000 (pav10vsbench10_truvarimafft_fnfp_slop50_merge10000_withoverlaprepeats)
## The stvar comparison uses truvari bench without the original mafft refine step.
## Both comparisons are restricted to regions covered by both callers.
## Requires dipcall and PAV rows for the ref + assembly in the analyses table.

## Comparison intermediates are shared by all benchmarks of a ref + assembly.
PAV_DISCREP_ID = "{ref_id}_{asm_id}_dipcall-{dipcall_param}_pav-{pav_param}"
PAV_DISCREP_DIR = f"results/pav_discrep/{PAV_DISCREP_ID}/"


rule pav_discrep_prep:
    # hap.py / vcfeval cannot parse symbolic or breakend ALT alleles (#192);
    # truvari handles them, but both comparisons use the same filtered inputs.
    input:
        unpack(get_pav_discrep_inputs),
        genome=get_genome_file,
    output:
        # Indexes are produced by the generic `tabix` rule to avoid a rule clash.
        dip_vcf=f"{PAV_DISCREP_DIR}dipcall.no-symbolic.vcf.gz",
        pav_vcf=f"{PAV_DISCREP_DIR}pav.no-symbolic.vcf.gz",
        bed=f"{PAV_DISCREP_DIR}shared-regions.bed",
    log:
        f"logs/exclusions/pav-discrep-prep/{PAV_DISCREP_ID}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        bcftools view -e 'ALT~"<" || ALT~"\\[" || ALT~"\\]"' -Oz -o {output.dip_vcf} {input.dip_vcf} 2> {log}
        bcftools view -e 'ALT~"<" || ALT~"\\[" || ALT~"\\]"' -Oz -o {output.pav_vcf} {input.pav_vcf} 2>> {log}
        bedtools intersect -a {input.dip_bed} -b {input.pav_bed} \
            | bedtools sort -g {input.genome} -i stdin \
            | bedtools merge -i stdin \
            1> {output.bed} 2>> {log}
        """


rule pav_discrep_happy:
    input:
        truth=f"{PAV_DISCREP_DIR}dipcall.no-symbolic.vcf.gz",
        truthidx=f"{PAV_DISCREP_DIR}dipcall.no-symbolic.vcf.gz.tbi",
        query=f"{PAV_DISCREP_DIR}pav.no-symbolic.vcf.gz",
        queryidx=f"{PAV_DISCREP_DIR}pav.no-symbolic.vcf.gz.tbi",
        bed=f"{PAV_DISCREP_DIR}shared-regions.bed",
        ref=get_ref_file,
        sdf=get_ref_sdf,
    output:
        multiext(
            f"{PAV_DISCREP_DIR}happy",
            ".runinfo.json",
            ".vcf.gz",
            ".vcf.gz.tbi",
            ".summary.csv",
            ".extended.csv",
        ),
    log:
        f"logs/exclusions/pav-discrep-happy/{PAV_DISCREP_ID}.log",
    benchmark:
        f"benchmark/exclusions/pav-discrep-happy_{PAV_DISCREP_ID}.tsv"
    conda:
        "../envs/happy.yml"
    threads: config["_happy_threads"]
    resources:
        mem_mb=config["_happy_mem"],
    params:
        prefix=f"{PAV_DISCREP_DIR}happy",
        gender=get_happy_gender_param,
    shell:
        """
        hap.py \
            {input.truth} \
            {input.query} \
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


rule pav_discrep_happy_extract_fpfns:
    input:
        vcf=f"{PAV_DISCREP_DIR}happy.vcf.gz",
        faidx=get_ref_index,
    output:
        f"{PAV_DISCREP_DIR}happy.fpfns.bed",
    log:
        f"logs/exclusions/pav-discrep-happy-fpfns/{PAV_DISCREP_ID}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    params:
        max_indel=config["_exclusion_params"]["pav_discrep_smvar_max_indel"],
    shell:
        """
        bcftools filter \
            --include 'ABS(ILEN)<={params.max_indel} && (FMT/BD=="FN" || FMT/BD=="FP")' {input.vcf} 2> {log} \
                | bcftools query -f "%CHROM\t%POS0\t%END\n" 2>> {log} \
                | bedtools sort -faidx {input.faidx} -i - 2>> {log} \
                | bedtools merge -i - 1> {output} 2>> {log}
        """


rule pav_discrep_truvari:
    input:
        base=f"{PAV_DISCREP_DIR}dipcall.no-symbolic.vcf.gz",
        baseidx=f"{PAV_DISCREP_DIR}dipcall.no-symbolic.vcf.gz.tbi",
        comp=f"{PAV_DISCREP_DIR}pav.no-symbolic.vcf.gz",
        compidx=f"{PAV_DISCREP_DIR}pav.no-symbolic.vcf.gz.tbi",
        bed=f"{PAV_DISCREP_DIR}shared-regions.bed",
        genome=get_ref_file,
        genomeidx=get_ref_index,
    output:
        fn=f"{PAV_DISCREP_DIR}truvari/fn.vcf.gz",
        fp=f"{PAV_DISCREP_DIR}truvari/fp.vcf.gz",
        summary=f"{PAV_DISCREP_DIR}truvari/summary.json",
    log:
        f"logs/exclusions/pav-discrep-truvari/{PAV_DISCREP_ID}.log",
    benchmark:
        f"benchmark/exclusions/pav-discrep-truvari_{PAV_DISCREP_ID}.tsv"
    conda:
        "../envs/truvari_core.yml"
    params:
        dir=lambda wildcards, output: Path(output.fn).parent,
        tmpdir=lambda wildcards: f"truvari_pd_{PAV_DISCREP_ID.format(** dict(wildcards.items()))}",
    shell:
        """
        rm -rf {params.tmpdir}
        truvari bench \
            -b {input.base} \
            -c {input.comp} \
            -o {params.tmpdir} \
            -f {input.genome} \
            --includebed {input.bed} \
            --passonly \
            --pick ac \
            --sizemin 50 \
            -B -1 \
            -r 2000 \
        &> {log}
        mkdir -p {params.dir}
        mv {params.tmpdir}/* {params.dir}
        rm -r {params.tmpdir}
        """


rule pav_discrep_truvari_extract_fpfns:
    input:
        fn=f"{PAV_DISCREP_DIR}truvari/fn.vcf.gz",
        fp=f"{PAV_DISCREP_DIR}truvari/fp.vcf.gz",
        faidx=get_ref_index,
    output:
        f"{PAV_DISCREP_DIR}truvari.fpfns.bed",
    log:
        f"logs/exclusions/pav-discrep-truvari-fpfns/{PAV_DISCREP_ID}.log",
    conda:
        "../envs/bcftools_and_bedtools.yml"
    shell:
        """
        (
            bcftools query -f "%CHROM\t%POS0\t%END\n" {input.fn} 2>> {log}
            bcftools query -f "%CHROM\t%POS0\t%END\n" {input.fp} 2>> {log}
        ) \
            | bedtools sort -faidx {input.faidx} -i - 2>> {log} \
            | bedtools merge -i - 1> {output} 2>> {log}
        """


rule pav_discrep_intersect_slop:
    input:
        bed=get_pav_discrep_fpfns_bed,
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_pav-discrep-{discrep_type}.bed",
    log:
        "logs/exclusions/pav-discrep-intersect/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_{discrep_type}.log",
    wildcard_constraints:
        discrep_type="smvar|stvar",
    conda:
        "../envs/bedtools.yml"
    params:
        slop=lambda wildcards: config["_exclusion_params"][
            f"pav_discrep_{wildcards.discrep_type}_slop"
        ],
        merge_d=lambda wildcards: config["_exclusion_params"][
            f"pav_discrep_{wildcards.discrep_type}_merge_dist"
        ],
    shell:
        """
        bedtools intersect -wa \
                -a {input.simple_repeat_bed} \
                -b {input.bed} \
            | bedtools multiinter -i stdin {input.bed} \
            | bedtools slop -b {params.slop} -i stdin -g {input.genome} \
            | bedtools sort -g {input.genome} -i stdin \
            | mergeBed -i stdin -d {params.merge_d} \
            1> {output} 2> {log}
        """
