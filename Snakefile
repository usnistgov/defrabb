import pandas as pd
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from os.path import join, basename

include: "rules/common.smk"

## File download
FTP = FTPRemoteProvider()

min_version("6.0")

configfile: "config/static.yml"
configfile: "config/dynamic.yml"
validate(config, "config/schema.yml")

# TODO add validations

asm_config = config["assemblies"]
bmk_config = config["benchmarks"]
str_config = config["stratifications"]
vcr_config = config["variant_caller_runs"]
hpy_config = config["happy_runs"]
ref_config = config["_references"]

hpy_combinations = unzip(
    (hkey, vkey, bkey, s["bed"], s["tsv"])
    for k, v in hpy_config.items()
    for hkey, vkey, bkey, s in product(
            [k],
            v["variant_caller_runs"],
            v["benchmarks"],
            v["strats"]
    )
)
hpy_keys, vcr_keys, bmk_keys, bed_keys, tsv_keys = tuple(
    map(list, hpy_combinations)
)

ref_keys = list(set([v["reference"] for v in str_config.values()]))

wildcard_constraints:
    asm_prefix="|".join(list(asm_config)),
    bmk_prefix="|".join(bmk_keys),
    bed_prefix="|".join(bed_keys),
    tsv_prefix="|".join(tsv_keys),
    vcr_prefix="|".join(vcr_keys),
    hpy_prefix="|".join(hpy_keys),
    ref_prefix="|".join(ref_keys)

# path initialization

resource_dir = "resources"
output_dir = "results"

asm_full_path = join(resource_dir, "assemblies", "{asm_prefix}")
ref_full_prefix = join(resource_dir, "references", "{ref_prefix}")
benchmark_full_prefix = join(resource_dir, "benchmarks", "{bmk_prefix}")
strats_full_path = join(resource_dir, "stratifications", "{bed_prefix}")
tsv_full_path = join(strats_full_path, "{tsv_prefix}.tsv")

vcr_full_prefix = join(
    output_dir,
    "dipcall",
    "{ref_prefix}",
    "{asm_prefix}",
    "{vcr_prefix}",
    "dipcall"
)

hpy_full_path = join(
    output_dir,
    "happy",
    "{vcr_prefix}",
    "{bmk_prefix}",
    "{bed_prefix}-{tsv_prefix}",
    "{hpy_prefix}",
)

rule all:
    input:
        expand(
            join(hpy_full_path, "happy_out.extended.csv"),
            zip,
            hpy_prefix = hpy_keys,
            vcr_prefix = vcr_keys,
            bmk_prefix = bmk_keys,
            bed_prefix = bed_keys,
            tsv_prefix = tsv_keys,
        )

################################################################################
## Get and prepare assemblies
################################################################################

rule get_assemblies:
    output: join(asm_full_path, "{haplotype}.fa")
    params:
        url = lambda wildcards: asm_config[wildcards.asm_prefix][wildcards.haplotype]
    shell:
        "curl -f -L {params.url} | gunzip -c > {output}"

################################################################################
## Get and prepare reference
################################################################################

## TODO these FTP calls sometimes cause timeout errors depending on how long we
## wait between rule calls
rule get_ref:
    output: "{}.fa".format(ref_full_prefix)
    params:
        # TODO not sure if this underscore thing is necessary here
        url = lambda wildcards: ref_config[wildcards.ref_prefix]["_ref_url"]
    shell:
        "curl -f --connect-timeout 120 -L {params.url} | gunzip -c > {output}"

rule index_ref:
    input: rules.get_ref.output
    output: "{}.fai".format(ref_full_prefix)
    wrapper: "0.61.0/bio/samtools/faidx"

################################################################################
## Get benchmark vcf.gz and .bed
################################################################################

# TODO wet code
rule get_benchmark_vcf:
    output: "{}.vcf.gz".format(benchmark_full_prefix)
    params:
        url = lambda wildcards: bmk_config[wildcards["bmk_prefix"]]["vcf_url"]
    shell: "curl -f -L -o {output} {params.url}"

rule get_benchmark_bed:
    output: "{}.bed".format(benchmark_full_prefix)
    params:
        url = lambda wildcards: bmk_config[wildcards["bmk_prefix"]]["bed_url"]
    shell: "curl -f -L -o {output} {params.url}"

# TODO add rule to get the tbi file as well when we need it
# (analogous to these rules)

################################################################################
## Get v2.0 stratifications
################################################################################

def read_tsv_bed_files (wildcards):
    b = wildcards.bed_prefix
    t = wildcards.tsv_prefix
    o = checkpoints.get_strat_tsv.get(bed_prefix=b, tsv_prefix=t).output[0]
    with o.open() as f:
        beds = pd.read_table(f, header=None)[1].tolist()
        return [join(strats_full_path, b) for b in beds]

checkpoint get_strat_tsv:
    output: tsv_full_path
    params:
        root = lambda wildcards: str_config[wildcards.bed_prefix]["root"],
        tsv = lambda wildcards: str_config[wildcards.bed_prefix]["tsv"][wildcards.tsv_prefix],
    shell: "curl -f -L -o {output} {params.root}/{params.tsv}"

rule get_strat_beds:
    output: join(strats_full_path, "{bed}")
    wildcard_constraints:
        # NOTE: not all files in the tsv are "*.bed.gz"; a few just have "*.gz"
        bed = ".*\.gz"
    params:
        root = lambda wildcards: str_config[wildcards["bed_prefix"]]["root"],
    shell: "curl -f -L -o {output} {params.root}/{wildcards.bed}"

################################################################################
## Run Dipcall
################################################################################

# TODO this seems brittle
def get_zdrop(wildcards):
    v = wildcards["vcr_prefix"]
    z = vcr_config[v]["zdrop"]
    return "-z {},{}".format(z[0], z[1]) if z else ""

def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_prefix]["is_male"]
    root = config["_par_bed_root"]
    par_path = join(root, ref_config[wildcards.ref_prefix]["_par_bed"])
    return "-x " + par_path if is_male else ""

rule run_dipcall:
    input:
        h1 = join(asm_full_path, "paternal.fa"),
        h2 = join(asm_full_path, "maternal.fa"),
        ref = rules.get_ref.output,
        ref_idx = rules.index_ref.output
    output:
        make = "{}.mak".format(vcr_full_prefix),
        vcf = "{}.dip.vcf.gz".format(vcr_full_prefix),
        bed = "{}.dip.bed".format(vcr_full_prefix),
        bam1 = "{}.hap1.bam".format(vcr_full_prefix),
        bam2 = "{}.hap2.bam".format(vcr_full_prefix)
    conda: "envs/dipcall.yml"
    params:
        prefix = vcr_full_prefix,
        male_bed = get_male_bed,
        ts = config["_dipcall_threads"],
        zdrop = get_zdrop
    log: "{}.log".format(vcr_full_prefix)
    resources:
        mem_mb = config["_dipcall_threads"] * 2 * 32000 ## GB per thread
    threads: config["_dipcall_threads"] * 2 ## For diploid
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
    output: "{}.dip.gap2homvarbutfiltered.vcf.gz".format(vcr_full_prefix)
    # bgzip is part of samtools, which is part of the diptcall env
    conda: "envs/dipcall.yml"
    shell: """
        gunzip -c {input} |\
        sed 's/1|\./1|1/' |\
        grep -v 'HET\|GAP1\|DIP' |\
        bgzip -c > {output}
    """

def get_targeted (wildcards):
    bed = "--target-regions " + rules.run_dipcall.output.bed
    ws = {
        "ref_prefix": str_config[wildcards.bed_prefix]["reference"],
        # TODO this is a list and this will make the expand return not a singleton list below
        "asm_prefix": vcr_config[wildcards.vcr_prefix]["assemblies"],
        "vcr_prefix": wildcards.vcr_prefix
        }
    # TODO this is ridiculous
    return expand(bed, **ws, )[0] if hpy_config[wildcards.hpy_prefix]["use_targeted"] else ""

rule run_happy:
    input:
        query = lambda wildcards: expand(
            rules.dip_gap2homvarbutfiltered.output,
            ref_prefix = str_config[wildcards.bed_prefix]["reference"],
            asm_prefix = vcr_config[wildcards.vcr_prefix]["assemblies"],
            allow_missing = True,
        ),
        truth = rules.get_benchmark_vcf.output,
        truth_regions = rules.get_benchmark_bed.output,
        strats = rules.get_strat_tsv.output,
        genome = lambda wildcards: expand(
            rules.get_ref.output,
            ref_prefix = str_config[wildcards.bed_prefix]["reference"],
        ),
        # NOTE this isn't actually required by the wrapper, but hap.py itself
        # will be quite sad without them
        strat_beds = read_tsv_bed_files
    # NOTE many files will be produced as output but this will be used to
    # signify completion
    output: join(hpy_full_path, "happy_out.extended.csv")
    priority: 1
    params:
        prefix = lambda _, output: output [0][:-13],
        threads = 6,
        engine = "vcfeval",
        extra = get_targeted
    log: join(hpy_full_path, "happy.log")
    wrapper: "0.78.0/bio/hap.py/hap.py"
