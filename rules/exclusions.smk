import pandas as pd
from pathlib import Path
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.io import apply_wildcards
from functools import partial

# resource paths

downloads_dir = resource_dir / "exclusions"

# dip_vcf_path = rules.run_dipcall.vcf
# dip_bed_path = rules.run_dipcall.bed

# results paths

exclusions_dir = Path(dip_bed_path).parent / "exclusions"
slop_dir = exclusions_dir / "slop"

# excluded_dip_file = exclusions_dir / "excluded.bed"
# excluded_log_file = exclusions_dir / "excluded.log"
dip_dir = exclusions_dir / "{bench_prefix}"
excluded_dip_file = dip_dir / "excluded.bed"
excluded_log_file = dip_dir / "excluded.log"

# downloading stuff

rule download_bed_gz:
    output:
        downloads_dir / "{index}.bed",
    params:
        url=lambda wildcards: config["exclusion_beds"][wildcards.index],
    shell:
        "curl -L {params.url} | gunzip -c > {output}"

# TODO hack since this is the only bed file that isn't processed according to
# the output from dipcall
rule link_gaps:
    input:
        downloads_dir / "gaps.bed",
    output:
        exclusions_dir / "gaps.bed"
    shell:
        "cp {input} {output}"
    

# get genome.txt (used for bedtools slop and flank)

mysql_login_params = "--user=genome --host=genome-mysql.soe.ucsc.edu -P 3306 -D hg38"


rule get_genome:
    output:
        downloads_dir / "genome.txt",
    params:
        login=mysql_login_params,
        extra="-N",
    conda:
        "envs/mysql.yml"
    shell:
        """
        mysql {params.login} {params.extra} -A -B -e \
        \"SELECT chrom,size FROM chromInfo
        ORDER BY 'chrom';\" \
        > {output} 
        """


# get satellites


rule get_satellites:
    output:
        downloads_dir / "satellites.bed",
    params:
        login=mysql_login_params,
        extra="-N",
    conda:
        "envs/mysql.yml"
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
    input:
        dip_vcf_path,
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
        sv_bed=rules.get_SVs_from_vcf.output,
        homopoly_bed=downloads_dir / "homopolymers.bed",
    output:
        exclusions_dir / "structural_variants.bed",
    conda:
        "envs/bedtools.yml"
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
        bed=downloads_dir / "{index}.bed",
        genome=rules.get_genome.output,
    output:
        slop_dir / "{index}.bed",
    conda:
        "envs/bedtools.yml"
    shell:
        "bedtools slop -i {input.bed} -g {input.genome} -b 15000 > {output}"


rule intersect_start_and_end:
    input:
        dip=dip_bed_path,
        xregions=rules.add_slop.output,
    output:
        start=exclusions_dir / "{index}_start.bed",
        end=exclusions_dir / "{index}_end.bed",
    conda:
        "envs/bedtools.yml"
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
        bed=dip_bed_path,
        genome=rules.get_genome.output,
    output:
        exclusions_dir / "flanks.bed",
    conda:
        "envs/bedtools.yml"
    shell:
        "bedtools flank -i {input.bed} -g {input.genome} -b 15000 > {output}"

def lookup_exclusions(wildcards):
    xset = analyses.loc[(wildcards.bench_prefix, "exclusion_set")]
    return [exclusions_dir / p for p in config["exclusion_sets"][xset]]

rule subtract_exclusions:
    input:
        dip_bed=dip_bed_path,
        # exclude these bed files in this order:
        # other_beds=[
        #     exclusions_dir / "segdups_start.bed",
        #     exclusions_dir / "segdups_end.bed",
        #     exclusions_dir / "tandemreps_start.bed",
        #     exclusions_dir / "tandemreps_end.bed",
        #     exclusions_dir / "satellites_start.bed",
        #     exclusions_dir / "satellites_end.bed",
        #     downloads_dir / "gaps.bed",
        #     exclusions_dir / "self_chains_start.bed",
        #     exclusions_dir / "self_chains_end.bed",
        #     rules.intersect_SVs_and_homopolymers.output,
        #     rules.add_flanks.output,
        # ],
        other_beds = lookup_exclusions
    output:
        excluded_dip_file
    log:
        excluded_log_file
    conda:
        "envs/bedtools.yml"
    shell:
        """
        python scripts/subtract_exclusions.py \
        {input.dip_bed} \
        {output} \
        {input.other_beds} > {log}
        """
