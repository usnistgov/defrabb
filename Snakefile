# import pandas as pd
from snakemake.utils import min_version
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from os.path import join, basename

include: "rules/common.smk"

## File download
FTP = FTPRemoteProvider()

min_version("6.0")

configfile: "config/config.yaml"

## Read table of samples and set wildcard prefix and constraints
# asm = pd.read_table(config["assemblies"]).set_index(["prefix"], drop = False)
# ASM_prefix = list(set(asm["prefix"]))
ASM_prefix = config["ASM_prefix"]

wildcard_constraints:
    prefix="|".join(ASM_prefix)

## Set reference to be used. Pipeline only uses GRCh38 or GRCh37.
ref_id = config["reference"]

ref_dependent_data = config["tool_data"][config["reference"]]

ref_url = ref_dependent_data["ref_url"]
par_ref = join(config["par_bed_root"], ref_dependent_data["par_bed"])
# strat_url = ref_dependent_data["strat_url"]
# strat_tsv = ref_dependent_data["strat_tsv"]
# strat_id = ref_dependent_data["strat_id"]

## Define basename for benchmark files
# base = os.path.basename(config["benchmark_vcfgz"])
# benchmark_name = base.split(".vcf.gz")[0]


rule all:
    input:
        "resources/references/{}.fa".format(ref_id),
        "resources/references/{}.fa.fai".format(ref_id),
        # "resources/benchmark/{}.vcf.gz".format(benchmark_name),
        # "resources/benchmark/{}.vcf.gz.tbi".format(benchmark_name),
        # "resources/benchmark/{}.bed".format(benchmark_name),
        expand("results/dipcall/{prefix}.mak", prefix = ASM_prefix),
        expand("results/dipcall/{prefix}.dip.vcf.gz", prefix = ASM_prefix),
        expand("results/dipcall/{prefix}.dip.bed", prefix = ASM_prefix),
        expand("results/dipcall/{prefix}.hap1.bam", prefix = ASM_prefix),
        expand("results/dipcall/{prefix}.hap2.bam", prefix = ASM_prefix),
        # expand("results/dipcall/{prefix}.dip.vcf.gz.tbi", prefix = ASM_prefix),
        # expand("results/dipcall/{prefix}.hap1.bam.bai", prefix = ASM_prefix),
        # expand("results/dipcall/{prefix}.hap2.bam.bai", prefix = ASM_prefix),
        expand("results/dipcall/{prefix}.dip.gap2homvarbutfiltered.vcf.gz", prefix = ASM_prefix),
        expand("resources/assemblies/{hap}.fa", hap = ["maternal", "paternal"])
        # "resources/stratifications/",
        # "resources/stratifications/{}/{}".format(config["reference"], strat_tsv),
        # expand("results/happy/nontargeted/{prefix}_benchmark-v" + config["benchmark_version"] + ".extended.csv", prefix = ASM_prefix),
        # expand("results/happy/targeted/{prefix}_benchmark-v" + config["benchmark_version"] + ".extended.csv", prefix = ASM_prefix)

################################################################################
## Get and prepare assemblies
################################################################################

rule get_assemblies:
    output:
        "resources/assemblies/{haplotype}.fa"
    params:
        url = lambda wildcards: config["assemblies"][wildcards["haplotype"]]
    shell:
        "curl -L {params.url} | gunzip -c > {output}"

################################################################################
## Get and prepare reference
################################################################################

## TODO these FTP calls sometimes cause timeout errors depending on how long we
## wait between rule calls
rule get_ref:
    input: FTP.remote(ref_url, keep_local=True)
    output: "resources/references/{}.fa".format(ref_id)
    shell: "gunzip -c {input} > {output}"

rule index_ref:
    input: "resources/references/{}.fa".format(ref_id)
    output: "resources/references/{}.fa.fai".format(ref_id)
    wrapper: "0.61.0/bio/samtools/faidx"

################################################################################
## Run Dipcall
################################################################################

rule dipcall_makefile:
    input:
        h1="resources/assemblies/paternal.fa",
        h2="resources/assemblies/maternal.fa",
        ref= rules.get_ref.output,
        ref_idx= rules.index_ref.output,
    output:
        "results/dipcall/{prefix}.mak"
    conda:
        "envs/dipcall.yml"
    params:
        prefix = "results/dipcall/{prefix}",
        male_bed = "-x " + par_ref if config["male"] else ""
    shell: """
    run-dipcall \
    {params.male_bed} \
    {params.prefix} \
    {input.ref} \
    {input.h1} \
    {input.h2} \
    > {output}
    """

## conda env not needed here since the .mak file uses absolute paths to the
## executables
rule run_dipcall:
    input:
        h1="resources/assemblies/paternal.fa",
        h2="resources/assemblies/maternal.fa",
        # h1=get_hap1,
        # h2=get_hap2,
        make= "results/dipcall/{prefix}.mak"
    output:
        vcf="results/dipcall/{prefix}.dip.vcf.gz",
        bed="results/dipcall/{prefix}.dip.bed",
        bam1="results/dipcall/{prefix}.hap1.bam",
        bam2="results/dipcall/{prefix}.hap2.bam"
    params:
        ts = config["dipcall_threads"]
    log:
        "results/dipcall/{prefix}_dipcall.log"
    threads:
        config["dipcall_threads"]
    shell: """
    date
    make -j{params.ts} -f {input.make}
    date
    """

rule dip_gap2homvarbutfiltered:
    input:
        vcf="results/dipcall/{prefix}.dip.vcf.gz"
    output:
        "results/dipcall/{prefix}.dip.gap2homvarbutfiltered.vcf.gz"
    shell: """
    gunzip -c {input.vcf} |\
    sed 's/1|\./1|1/' |\
    grep -v 'HET\|GAP1\|DIP' |\
    bgzip -c > {output}
    """
