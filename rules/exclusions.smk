import pandas as pd
from pathlib import Path
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.io import apply_wildcards
from functools import partial

# results paths
dip_bed_path = "results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed"
exclusions_dir = Path(dip_bed_path).parent / "exclusions"

# downloading stuff

rule download_bed_gz:
    output: "resources/exclusions/{ref_refix}_{genome_index}.bed",
    params:
        url=lambda wildcards: config["exclusion_beds"][wildcards.index],
    shell:
        "curl -L {params.url} | gunzip -c > {output}"

# TODO hack since this is the only bed file that isn't processed according to
# the output from dipcall
rule link_gaps:
    input: "resources/exclusions/{ref_prefix}_gaps.bed",
    output:
        exclusions_dir / "{ref_prefix}_gaps.bed"
    shell:
        "cp {input} {output}"
    

# get genome.txt (used for bedtools slop and flank)

mysql_login_params = "--user=genome --host=genome-mysql.soe.ucsc.edu -P 3306 -D hg38"


rule get_genome:
    output: "resources/exclusions/{ref_prefix}.genome",
    params:
        login=mysql_login_params,
        extra="-N",
    conda:
        "../envs/mysql.yml"
    shell:
        """
        mysql {params.login} {params.extra} -A -B -e \
        \"SELECT chrom,size FROM chromInfo
        ORDER BY 'chrom';\" \
        > {output} 
        """


# get satellites


rule get_satellites:
    output: "resources/exclusions/{ref_prefix}_satellites.bed",
    params:
        login=mysql_login_params,
        extra="-N",
    conda:
        "../envs/mysql.yml"
    shell:
        """
        mysql {params.login} {params.extra} -A -B -e \
        \"SELECT genoName,genoStart,genoEnd,repClass FROM rmsk
        WHERE repClass = 'Satellite'
        ORDER BY 'genoName';\" | \
        awk 'OFS="\t" {{print $1, $2, $3}}' \
        > {output} 
        """


# structural variants


rule get_SVs_from_vcf:
    input: "results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.vcf.gz",
    output:
        exclusions_dir / "dip_SVs.bed",
    shell:
        """
        gunzip -c {input} | \
            awk 'length($4)>49 || length($5)>49' | \
            awk '{{FS=OFS="\\t"}} {{print $1,$2-1,$2+length($4)}}' \
            > {output}
        """


rule intersect_SVs_and_homopolymers:
    input:
        sv_bed=exclusions_dir / "dip_SVs.bed",
        homopoly_bed="resources/exclusions/{ref_prefix}_homopolymers.bed",
    output:
        exclusions_dir / "structural_variants.bed",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        intersectBed -wa \
        -a {input.homopoly_bed} \
        -b {input.sv_bed} | \
        multiIntersectBed -i stdin {input.sv_bed} | \
        awk '{{FS=OFS="\\t"}} {{print $1,$2-50,$3+50}}' | \
        mergeBed -i stdin -d 1000 \
        > {output}
        """


# segdups/tandemreps/self chains


rule add_slop:
    input:
        bed="resources/exclusions/{ref_prefix}_{genomic_regions}.bed",
        genome="resources/exclusions/{ref_prefix}.genome",
    output: exclusions_dir / "slop" / "{genomics_regions}.bed",
    conda: "../envs/bedtools.yml"
    shell: "bedtools slop -i {input.bed} -g {input.genome} -b 15000 > {output}"


rule intersect_start_and_end:
    input:
        dip="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed",
        xregions=exclusions_dir / "slop" / "{genomic_regions}.bed",
    output:
        start=exclusions_dir / "{genomic_regions}_start.bed",
        end=exclusions_dir / "{genomic_regions}_end.bed",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk '{{FS=OFS="\t"}} {{print $1, $2, $2+1}}' {input.dip} | \
        bedtools intersect -wa -a {input.xregions} -b stdin \
        > {output.start}

        awk '{{FS=OFS="\t"}} {{print $1, $3-1, $3}}' {input.dip} | \
        bedtools intersect -wa -a {input.xregions} -b stdin \
        > {output.end}
        """


# flanks


rule add_flanks:
    input:
        bed="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed",
        genome="resources/references/{ref_prefix}.fa",
    output:
        exclusions_dir / "flanks.bed",
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools flank -i {input.bed} -g {input.genome} -b 15000 > {output}"

rule subtract_exclusions:
    input:
        dip_bed="results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed",
        other_beds = lookup_excluded_region_set
    output:
        exclusions_dir / "{bench_prefix}" / "excluded.bed"
    log: "logs/subtract_exclusions/{bench_prefix}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}_exclusions.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        python scripts/subtract_exclusions.py \
        {input.dip_bed} \
        {output} \
        {input.other_beds} > {log}
        """
