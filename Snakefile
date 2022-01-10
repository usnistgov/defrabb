import pandas as pd
from pathlib import Path
from itertools import product
from more_itertools import unzip, flatten
from snakemake.utils import min_version, validate
from snakemake.io import apply_wildcards
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


def get_analyses(path):
    # target_regions must be a string even though it might only contain
    # 'boolean' values
    analyses = pd.read_table(path, dtype={"target_regions": str})
    validate(analyses, "config/analyses-schema.yml")

    try:
        return analyses.set_index("bench_id", verify_integrity=True)
    except ValueError:
        print("All keys in column 'bench_id' must by unique")


analyses = get_analyses("config/analyses.tsv")

################################################################################
# init paths

resource_dir = Path("resources")
output_dir = Path("results")

manual_target_regions_path = resource_dir / "manual" / "target_regions"

asm_full_path = resource_dir / "assemblies" / "{asm_prefix}"
ref_full_prefix = resource_dir / "references" / "{ref_prefix}"
benchmark_full_prefix = resource_dir / "benchmarks" / "{bmk_prefix}"
strats_base_path = resource_dir / "stratifications"

strats_full_path = strats_base_path / "v3.0"
tsv_full_path = (
    strats_full_path / "{ref_prefix}" / "v3.0-{ref_prefix}-all-stratifications.tsv"
)

vcr_full_prefix = (
    output_dir
    / "dipcall"
    / "{ref_prefix}"
    / "{asm_prefix}"
    / "{vcr_cmd}_{vcr_params}"
    / "dipcall"
)
bench_full_path = output_dir / "bench" / "{bench_prefix}"

hpy_full_path = bench_full_path / "happy"

tvi_full_path = bench_full_path / "truvari"

################################################################################
# init wildcard constraints


def format_constraint(xs):
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


def expand_bench_output(path, cmd):
    bps = analyses[analyses["bench_cmd"] == cmd].index.tolist()
    return expand(path, bench_prefix=bps)

## Rules to run locally
localrules: get_ref, get_assemblies

rule all:
    input:
        expand_bench_output(hpy_full_path / "happy_out.extended.csv", "happy"),
        expand_bench_output(tvi_full_path / "out" / "summary.txt", "truvari"),


################################################################################
# Get and prepare assemblies


rule get_assemblies:
    output:
        asm_full_path / "{haplotype}.fa",
    params:
        url=lambda wildcards: asm_config[wildcards.asm_prefix][wildcards.haplotype],
    shell:
        "curl -f -L {params.url} | gunzip -c > {output}"


################################################################################
# Get and prepare reference


rule get_ref:
    output:
        ref_full_prefix.with_suffix(".fa"),
    params:
        url=lambda wildcards: ref_config[wildcards.ref_prefix]["ref_url"],
    shell:
        "curl -f --connect-timeout 120 -L {params.url} | gunzip -c > {output}"


rule index_ref:
    input: "resources/references/{ref}.fa"
    output: "resources/references/{ref}.fa.fai"
    wrapper: "0.79.0/bio/samtools/faidx"

################################################################################
# Get stratifications


rule get_strats:
    output:
        tsv_full_path,
    params:
        root=config["_strats_root"],
        target=strats_full_path,
    shell:
        """
        curl -L \
            {params.root}/v3.0/v3.0-stratifications-{wildcards.ref_prefix}.tar.gz | \
            gunzip -c | \
            tar x -C {params.target}
        """


################################################################################
# Run Dipcall


def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_prefix]["is_male"]
    root = config["_par_bed_root"]
    par_path = Path(root) / ref_config[wildcards.ref_prefix]["par_bed"]
    return "-x " + str(par_path) if is_male else ""


def get_extra(wildcards):
    # TODO this seems brittle
    return "" if "nan" == wildcards.vcr_params else wildcards.vcr_params


rule run_dipcall:
    input:
        h1=asm_full_path / "paternal.fa",
        h2=asm_full_path / "maternal.fa",
        ref=rules.get_ref.output,
        ref_idx=rules.index_ref.output,
    output:
        make=vcr_full_prefix.with_suffix(".mak"),
        vcf=vcr_full_prefix.with_suffix(".dip.vcf.gz"),
        bed=vcr_full_prefix.with_suffix(".dip.bed"),
        bam1=vcr_full_prefix.with_suffix(".hap1.bam"),
        bam2=vcr_full_prefix.with_suffix(".hap2.bam"),
    conda:
        "envs/dipcall.yml"
    params:
        prefix=str(vcr_full_prefix),
        male_bed=get_male_bed,
        ts=config["_dipcall_threads"],
        extra=get_extra,
    log:
        vcr_full_prefix.with_suffix(".log"),
    resources:
        mem_mb=config["_dipcall_threads"] * 2 * 128000,  ## GB per thread
    threads: config["_dipcall_threads"] * 2  ## For diploid
    shell:
        """
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


################################################################################
# Postprocess variant caller output

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


rule split_multiallelic_sites:
    input:
        rules.run_dipcall.output.vcf,
    output:
        vcf=vcr_full_prefix.with_suffix(".dip.split_multi.vcf.gz"),
        vcf_tbi=vcr_full_prefix.with_suffix(".dip.split_multi.vcf.gz.tbi"),
    conda:
        "envs/bcftools.yml"
    shell:
        """
        bcftools norm -m - {input} -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


################################################################################
## Run happy


def apply_analyses_wildcards(s, keyvals, wildcards):
    ws = {k: analyses.loc[(wildcards.bench_prefix, v)] for k, v in keyvals.items()}
    return expand(s, **ws)


# TODO happy will break in non-obvious ways if this file doesn't exist
def get_targeted(wildcards):
    # ASSUME: target_regions is either "true," "false," or a filename (all
    # strings); the schema itself defines either a string or boolean type for
    # this field, but the dataframe when parsed will contain all strings for
    # this column
    h = wildcards.bench_prefix
    trs = analyses.loc[(h, "target_regions")]
    if trs == "false":
        return ""
    else:
        if trs == "true":
            # TODO not dry
            bed = apply_wildcards(
                rules.run_dipcall.output.bed,
                {
                    "ref_prefix": analyses.loc[(h, "ref")],
                    "asm_prefix": analyses.loc[(h, "asm_id")],
                    "vcr_cmd": analyses.loc[(h, "varcaller")],
                    "vcr_params": analyses.loc[(h, "vc_params")],
                },
            )
        else:
            bed = manual_target_regions_path / trs
        return "--target-regions {}".format(bed)


rule run_happy:
    input:
        query=partial(
            apply_analyses_wildcards,
            rules.run_dipcall.output.vcf,
            {
                "ref_prefix": "ref",
                "asm_prefix": "asm_id",
                "vcr_cmd": "varcaller",
                "vcr_params": "vc_params",
            },
        ),
        # TODO not dry
        truth=partial(
            apply_analyses_wildcards,
            rules.get_benchmark_vcf.output,
            {"bmk_prefix": "truth_var_id"},
        ),
        truth_regions=partial(
            apply_analyses_wildcards,
            rules.get_benchmark_bed.output,
            {"bmk_prefix": "truth_var_id"},
        ),
        strats=partial(
            apply_analyses_wildcards,
            rules.get_strats.output,
            {
                "ref_prefix": "ref",
            },
        ),
        genome=partial(
            apply_analyses_wildcards, rules.get_ref.output, {"ref_prefix": "ref"}
        ),
    # TODO not dry
    output:
        hpy_full_path / "happy_out.extended.csv",
    params:
        prefix=str(hpy_full_path / "happy_out"),
        threads=6,
        engine="vcfeval",
        extra=get_targeted,
    log:
        hpy_full_path / "happy.log",
    wrapper:
        "0.78.0/bio/hap.py/hap.py"


################################################################################
## Run Truvari


rule run_truvari:
    input:
        query=partial(
            apply_analyses_wildcards,
            rules.split_multiallelic_sites.output.vcf,
            {
                "ref_prefix": "ref",
                "asm_prefix": "asm_id",
                "vcr_cmd": "varcaller",
                "vcr_params": "vc_params",
            },
        ),
        # TODO not dry
        truth=partial(
            apply_analyses_wildcards,
            rules.get_benchmark_vcf.output,
            {"bmk_prefix": "truth_var_id"},
        ),
        truth_regions=partial(
            apply_analyses_wildcards,
            rules.get_benchmark_bed.output,
            {"bmk_prefix": "truth_var_id"},
        ),
        # NOTE this isn't actually fed to the command but is still necessary
        truth_tbi=partial(
            apply_analyses_wildcards,
            rules.get_benchmark_tbi.output,
            {"bmk_prefix": "truth_var_id"},
        ),
        genome=partial(
            apply_analyses_wildcards, rules.get_ref.output, {"ref_prefix": "ref"}
        ),
    output:
        tvi_full_path / "out" / "summary.txt",
    log:
        tvi_full_path / "truvari.log",
    params:
        extra=lambda wildcards: analyses.loc[(wildcards.bench_prefix, "bench_params")],
        prefix=str(tvi_full_path / "out"),
        tmpdir="/tmp/truvari",
    conda:
        "envs/truvari.yml"
    # TODO this tmp thing is a workaround for the fact that snakemake
    # over-zealously makes output directories when tools like truvari expect
    # them to not exist. Also, /tmp is only a thing on Linux (if that matters)
    shell:
        """
        truvari bench \
            -b {input.truth} \
            -c {input.query} \
            -o {params.tmpdir} \
            -f {input.genome} \
            --includebed {input.truth_regions} \
            {params.extra}
        mv {params.tmpdir}/* {params.prefix}
        rm -r {params.tmpdir}
        """
