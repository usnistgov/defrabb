import pandas as pd
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.io import apply_wildcards
from os.path import join, basename, splitext
from functools import partial

include: "rules/common.smk"

min_version("6.0")

################################################################################
# init resources

configfile: "config/resources.yml"
validate(config, "config/resources-schema.yml")

asm_config = config["assemblies"]
bmk_config = config["benchmarks"]
ref_config = config["references"]

################################################################################
# init analyses

ANALYSES_TSV = "config/analyses.tsv"
_analyses = pd.read_table(ANALYSES_TSV)
validate(_analyses, "config/analyses-schema.yml")

try:
    analyses = _analyses.set_index("happy_id", verify_integrity=True)
except ValueError:
    print("All keys in column 'happy_id' must by unique")

################################################################################
# init paths

resource_dir = "resources"
output_dir = "results"
manual_dir = "manual"

manual_target_regions_path = join(manual_dir, "target_regions")

asm_full_path = join(resource_dir, "assemblies", "{asm_prefix}")
ref_full_prefix = join(resource_dir, "references", "{ref_prefix}")
benchmark_full_prefix = join(resource_dir, "benchmarks", "{bmk_prefix}")
strats_base_path = join(
    resource_dir,
    "stratifications",
)
strats_full_path = join(
    strats_base_path,
    "{str_ver}",
    "{ref_prefix}",
)
tsv_full_path = join(strats_full_path, "{str_tsv}")

vcr_full_prefix = join(
    output_dir,
    "dipcall",
    "{ref_prefix}",
    "{asm_prefix}",
    "{vcr_cmd}_{vcr_params}",
    "dipcall"
)

hpy_full_path = join(output_dir, "happy", "{hpy_prefix}")

################################################################################
# init wildcard constraints

def format_constraint (xs):
    return "|".join(set(xs))

# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_prefix=format_constraint(asm_config),
    bmk_prefix=format_constraint(bmk_config),
    ref_prefix=format_constraint(ref_config),

################################################################################
# main rule
#
# Define what files we want hap.py to make, and these paths will contain the
# definitions for the assemblies, variant caller, etc to use in upstream rules.

rule all:
    input:
        expand(
            join(hpy_full_path, "happy_out.extended.csv"),
            hpy_prefix = analyses.index.tolist()
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

rule get_ref:
    output: "{}.fa".format(ref_full_prefix)
    params:
        url = lambda wildcards: ref_config[wildcards.ref_prefix]["ref_url"]
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

# def read_tsv_bed_files (wildcards):
#     b = wildcards.bed_prefix
#     t = wildcards.tsv_prefix
#     o = checkpoints.get_strat_tsv.get(bed_prefix=b, tsv_prefix=t).output[0]
#     with o.open() as f:
#         beds = pd.read_table(f, header=None)[1].tolist()
#         return [join(strats_full_path, b) for b in beds]

# checkpoint get_strat_tsv:
#     output: tsv_full_path
#     params:
#         root = lambda wildcards: str_config[wildcards.bed_prefix]["root"],
#         tsv = lambda wildcards: str_config[wildcards.bed_prefix]["tsv"][wildcards.tsv_prefix],
#     shell: "curl -f -L -o {output} {params.root}/{params.tsv}"

# rule get_strat_beds:
#     output: join(strats_full_path, "{bed}")
#     wildcard_constraints:
#         # NOTE: not all files in the tsv are "*.bed.gz"; a few just have "*.gz"
#         bed = ".*\.gz"
#     params:
#         root = lambda wildcards: str_config[wildcards["bed_prefix"]]["root"],
#     shell: "curl -f -L -o {output} {params.root}/{wildcards.bed}"

rule get_strats:
    output: tsv_full_path
    params:
        url_root = config["_strats_root"],
        prefix = strats_base_path
    shell: """
    wget -r {params.url_root}/{wildcards.str_ver}/{wildcards.ref_prefix} \
    -nH --cut-dirs=4 \
    -P {params.prefix}
    """

################################################################################
## Run Dipcall
################################################################################


def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_prefix]["is_male"]
    root = config["_par_bed_root"]
    par_path = join(root, ref_config[wildcards.ref_prefix]["par_bed"])
    return "-x " + par_path if is_male else ""

def get_extra (wildcards):
    print(wildcards.vcr_params)
    # TODO this seems brittle
    return "" if "nan" == wildcards.vcr_params else wildcards.vcr_params

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
        extra = get_extra
        # extra = lambda wildcards: "" if pd.isna(wildcards.vcr_params) else wildcards.vcr_params
    log: "{}.log".format(vcr_full_prefix)
    resources:
        mem_mb = config["_dipcall_threads"] * 2 * 32000 ## GB per thread
    threads: config["_dipcall_threads"] * 2 ## For diploid
    shell: """
        echo "Writing Makefile defining dipcall pipeline"
        run-dipcall \
            {params.extra} \
            {params.male_bed} \
            {params.prefix} \
            {input.ref} \
            {input.h1} \
            {input.h2} \
            > {output.make}
        
        echo "Running dipcall pipeline"
        make -j{params.ts} -f {output.make}
    """

# rule dip_gap2homvarbutfiltered:
#     input: rules.run_dipcall.output.vcf
#     output: "{}.dip.gap2homvarbutfiltered.vcf.gz".format(vcr_full_prefix)
#     # bgzip is part of samtools, which is part of the diptcall env
#     conda: "envs/dipcall.yml"
#     shell: """
#         gunzip -c {input} |\
#         sed 's/1|\./1|1/' |\
#         grep -v 'HET\|GAP1\|DIP' |\
#         bgzip -c > {output}
#     """

def apply_analyses_wildcards(s, keyvals, wildcards):
    print(type(keyvals))
    ws = {k: analyses.loc[(wildcards.hpy_prefix, v)] for k, v in keyvals.items()}
    return expand(s, **ws)

def get_targeted(wildcards):
    # ASSUME: target_regions is either True, False, or a string
    h = wildcards.hpy_prefix
    trs = analyses.loc[(h, "target_regions")]
    if trs == False:
        return ""
    else:
        if trs == True:
            # TODO not dry
            bed = apply_wildcards(
                rules.run_dipcall.output.bed,
                {
                    "ref_prefix": analyses.loc[(h, "ref")],
                    "asm_prefix": analyses.loc[(h, "asm_id")],
                    "vcr_cmd": analyses.loc[(h, "varcaller")],
                    "vcr_params": analyses.loc[(h, "vc_params")]
                }
            )
        else:
            bed = join(manual_target_regions_path, trs)
        return "--target-regions {}".format(bed)

rule run_happy:
    input:
        query = partial(
            apply_analyses_wildcards,
            rules.run_dipcall.output.vcf,
            {
                "ref_prefix": "ref",
                "asm_prefix": "asm_id",
                "vcr_cmd": "varcaller",
                "vcr_params": "vc_params",
            }
        ),
        # TODO not dry
        truth = partial(
            apply_analyses_wildcards,
            rules.get_benchmark_vcf.output,
            {"bmk_prefix": "truth_var_id"}
        ),
        truth_regions = partial(
            apply_analyses_wildcards,
            rules.get_benchmark_bed.output,
            {"bmk_prefix": "truth_var_id"}
        ),
        strats = partial(
            apply_analyses_wildcards,
            rules.get_strats.output,
            {
                "str_ver": "strat_version",
                "str_tsv": "strat_list",
                "ref_prefix": "ref",
            }
        ),
        genome = partial(
            apply_analyses_wildcards,
            rules.get_ref.output,
            {"ref_prefix": "ref"}
        )
    output: join(hpy_full_path, "happy_out.extended.csv")
    priority: 1
    params:
        prefix = lambda _, output: output [0][:-13],
        threads = 6,
        engine = "vcfeval",
        extra = get_targeted
    log: join(hpy_full_path, "happy.log")
    wrapper: "0.78.0/bio/hap.py/hap.py"
