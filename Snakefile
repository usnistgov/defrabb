import pandas as pd
from snakemake.utils import min_version
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from os.path import join, basename

include: "rules/common.smk"

## File download
FTP = FTPRemoteProvider()

min_version("6.0")

configfile: "config/config.yaml"

## Read table of samples and set wildcard prefix and constraints
benchmark_params = pd.read_table(config["benchmark_params"]).set_index(["prefix"], drop = False)
benchmark_prefix = list(set(benchmark_params["prefix"]))
# ASM_prefix = list(set(asm["prefix"]))
ASM_prefix = config["ASM_prefix"]

wildcard_constraints:
    prefix="|".join(ASM_prefix)

## Set reference to be used. Pipeline only uses GRCh38 or GRCh37.
ref_id = config["reference"]

ref_dependent_data = config["tool_data"][config["reference"]]

ref_url = ref_dependent_data["ref_url"]
par_ref = join(config["par_bed_root"], ref_dependent_data["par_bed"])
strat_url = ref_dependent_data["strat_url"]
strat_tsv = ref_dependent_data["strat_tsv"]
strat_id = ref_dependent_data["strat_id"]

rule all:
    input:
        # stratifications
        "resources/stratifications/",
        "resources/stratifications/{}/{}".format(config["reference"], strat_tsv),

        # reference
        "resources/references/{}.fa".format(ref_id),
        "resources/references/{}.fa.fai".format(ref_id),

        # assemblies
        expand("resources/assemblies/{hap}.fa", hap = ["maternal", "paternal"]),

        # dipcall
        expand("results/dipcall/{v}.mak", v = ASM_prefix),
        expand("results/dipcall/{v}.dip.vcf.gz", v = ASM_prefix),
        expand("results/dipcall/{v}.dip.bed", v = ASM_prefix),
        expand("results/dipcall/{v}.hap1.bam", v = ASM_prefix),
        expand("results/dipcall/{v}.hap2.bam", v = ASM_prefix),
        expand("results/dipcall/{v}.dip.gap2homvarbutfiltered.vcf.gz", v = ASM_prefix),

        # benchmark truthsets
        expand("resources/benchmark/{b}.vcf.gz", b = benchmark_prefix),
        expand("resources/benchmark/{b}.bed", b = benchmark_prefix),

        # benchmark output
        expand("results/happy/{v}/{b}/happy_out.extended.csv", v = ASM_prefix, b = benchmark_prefix),

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
    output: "resources/references/{}.fa".format(ref_id)
    params:
        url = "http://{}".format(ref_dependent_data["ref_url"])
    shell:
        "curl --connect-timeout 120 -L {params.url} | gunzip -c > {output}"

rule index_ref:
    input: "resources/references/{}.fa".format(ref_id)
    output: "resources/references/{}.fa.fai".format(ref_id)
    wrapper: "0.61.0/bio/samtools/faidx"

################################################################################
## Get benchmark vcf.gz and .bed
################################################################################

# TODO wet code
rule get_benchmark_vcf:
    output: "resources/benchmark/{bm_prefix}.vcf.gz"
    params:
        url = lambda wildcards: benchmark_params.loc[wildcards.bm_prefix, "truth_vcf_url"]
    shell: "curl -L -o {output} {params.url}"

rule get_benchmark_bed:
    output: "resources/benchmark/{bm_prefix}.bed"
    params:
        url = lambda wildcards: benchmark_params.loc[wildcards.bm_prefix, "truth_bed_url"]
    shell: "curl -L -o {output} {params.url}"

# TODO add rule to get the tbi file as well when we need it
# (analogous to these rules)

################################################################################
## Get v2.0 stratifications
################################################################################

rule get_strats:
    output:
        dir = directory("resources/stratifications/"),
        tsv = "resources/stratifications/{}/{}".format(ref_id, strat_tsv)
    params: strats = {strat_url}
    shell: "wget -r {params.strats} -nH --cut-dirs=5 -P {output.dir}"

################################################################################
## Run Dipcall
################################################################################

rule run_dipcall:
    input:
        h1 = "resources/assemblies/paternal.fa",
        h2 = "resources/assemblies/maternal.fa",
        ref = rules.get_ref.output,
        ref_idx = rules.index_ref.output
    output:
        make = "results/dipcall/{vc_prefix}.mak",
        vcf = "results/dipcall/{vc_prefix}.dip.vcf.gz",
        bed = "results/dipcall/{vc_prefix}.dip.bed",
        bam1 = "results/dipcall/{vc_prefix}.hap1.bam",
        bam2 = "results/dipcall/{vc_prefix}.hap2.bam"
    conda: "envs/dipcall.yml"
    params:
        prefix = "results/dipcall/{vc_prefix}",
        male_bed = "-x " + par_ref if config["male"] else "",
        ts = config["dipcall_threads"],
        zdrop = "-z " + config["dipcall_zdrop"] if config["dipcall_zdrop"] else ""
    log: "results/dipcall/{vc_prefix}_dipcall.log"
    resources: mem_mb = config["dipcall_threads"] * 2 * 32000 ## GB per thread
    threads: config["dipcall_threads"] * 2 ## For diploid
    shell: """
        echo "Writing Makefile defining dipcall pipeline"
        run-dipcall \
            {params.zdrop} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            > {output.make}
        
        echo "Running dipcall pipeline"
        make -j{params.ts} -f {output.make}
    """

rule dip_gap2homvarbutfiltered:
    input: rules.run_dipcall.output.vcf
    output: "results/dipcall/{vc_prefix}.dip.gap2homvarbutfiltered.vcf.gz"
    # bgzip is part of samtools, which is part of the diptcall env
    conda: "envs/dipcall.yml"
    shell: """
        gunzip -c {input} |\
        sed 's/1|\./1|1/' |\
        grep -v 'HET\|GAP1\|DIP' |\
        bgzip -c > {output}
    """

rule run_happy:
    input:
        query = rules.dip_gap2homvarbutfiltered.output,
        truth = rules.get_benchmark_vcf.output,
        truth_regions = rules.get_benchmark_bed.output,
        strats = rules.get_strats.output.tsv,
        genome = rules.get_ref.output
    # NOTE many files will be produced as output but this will be used to
    # signify completion
    output: "results/happy/{vc_prefix}/{bm_prefix}/happy_out.extended.csv",
    priority: 1
    params:
        prefix = lambda _, output: output [0][:-13],
        threads = 6,
        engine = "vcfeval",
        extra = "--target-regions " + rules.run_dipcall.output.bed if (lambda wildcards: benchmark_params.loc[wildcards.bm_prefix, "targeted"]) else ""
    log: "results/happy/{vc_prefix}/{bm_prefix}/happy.log"
    wrapper: "0.78.0/bio/hap.py/hap.py"
