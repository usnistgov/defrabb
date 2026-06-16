import pandas as pd

genomic_regions = [
    "all-tr-and-homopolymers",
    "segdups",
    "tandem-repeats",
    "gaps",
    "self-chains",
    "satellites",
    "hifi-pacbioDV-XY-discrep",
    "imperfecthomopol-gt30",
    "hifiasm-HPRC-T2Tdiscrep",
    "XYelement-homopolymer-T2T-discrep",
    "XYdipcallmanualbugs",
    "VDJ",
    "consecutive-svs",
    "dipcall-bugs-T2TACE",
    "HG002Q100-pav-discrep-smvar",
    "HG002Q100-pav-discrep-stvar",
    "HG002Q100-pav-inversions",
    "HG002Q100-errors",
    "HG002Q100-mosaic",
    "HG002Q100-delins-errors",
    "TSPY2-segdups",
    "self-discrep",
    "pav-inv",
]


wildcard_constraints:
    ref_id="GRCh38|GRCh37|GRCh38_chr21|CHM13v2.0",
    genomic_region="|".join(genomic_regions),


# downloading beds used for exclusion
rule download_bed_gz:
    output:
        "resources/exclusions/{ref_id}/{genomic_region}.bed",
    log:
        "logs/download_bed_gz/{ref_id}-{genomic_region}.log",
    conda:
        "../envs/download_remotes.yml"
    params:
        url=lambda wildcards: config["references"][wildcards.ref_id]["exclusions"][
            wildcards.genomic_region
        ],
    shell:
        "bash scripts/fetch_resource.sh {params.url} {output} &> {log}"


rule make_gaps_bed:
    input:
        fa="resources/references/{ref_id}.fa",
    output:
        bed="resources/exclusions/{ref_id}/gaps.bed",
    log:
        "logs/make_gap_bed/{ref_id}.log",
    conda:
        "../envs/seqtk.yml"
    params:
        minNs=config["_exclusion_params"]["gap_min_ns"],
    shell:
        "seqtk gap -l {params.minNs} {input.fa} 1> {output.bed} 2> {log}"


rule get_sv_widen_coords:
    input:
        vcf=get_processed_vcf,
        genome=get_genome_file,
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_SVs.bed",
        tbl="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_SVs.tsv",
    log:
        "logs/exclusions/sv_widen_coords/{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_SV_coords.log",
    conda:
        "../envs/sv_widen_coords.yml"
    params:
        verbose="--verbose" if config.get("debug") else "",
    shell:
        """
        python scripts/get_sv_widen_coords.py \
            --input {input.vcf} \
            --output {output.bed} \
            {params.verbose} \
            --table \
            --sort-merge \
            --genome {input.genome} \
            --log {log}
        """


rule intersect_SVs_and_simple_repeats:
    input:
        sv_bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_SVs.bed",
        simple_repeat_bed="resources/exclusions/{ref_id}/all-tr-and-homopolymers_sorted.bed",
        genome=get_genome_file,
    output:
        "results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_svs-and-simple-repeats.bed",
    log:
        "logs/exclusions/SVs/{bench_id}_SVs_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    benchmark:
        "benchmark/exclusions/{bench_id}_SVs_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.tsv"
    conda:
        "../envs/bedtools.yml"
    params:
        slop=config["_exclusion_params"]["sv_repeat_slop"],
        merge_d=config["_exclusion_params"]["sv_repeat_merge_dist"],
    shell:
        """
        intersectBed -wa \
                -a {input.simple_repeat_bed} \
                -b {input.sv_bed} \
            | multiIntersectBed -i stdin {input.sv_bed} \
            | bedtools slop -i stdin -g {input.genome} -b {params.slop} \
            | mergeBed -i stdin -d {params.merge_d} \
            1> {output} 2>{log}
        """


## BED file with consecutive insertions and deletions from assembly-assembly bams
rule get_consecutive_svs:
    input:
        ## Reuses the existing dipcall run for this ref+asm (see
        ## get_consecutive_svs_bams) so PAV benchmarks do not trigger a
        ## duplicate dipcall run inside their own vc_id directory.
        unpack(get_consecutive_svs_bams),
    output:
        bed="results/draft_benchmarksets/{bench_id}/exclusions/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}_consecutive-svs.bed",
    log:
        "logs/exclusions/consecutive-svs/{bench_id}_{ref_id}_{bench_type}_{asm_id}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/consecutive_svs.yml"
    params:
        min_bp=config["_exclusion_params"]["consecutive_sv_min_bp"],
    shell:
        """
        python scripts/get_dipcall_delins.py \
            --min_bp {params.min_bp} \
            --hap1_bam {input.hap1_bam} \
            --hap2_bam {input.hap2_bam} \
            --output_bed {output.bed} \
            &> {log}
        """
