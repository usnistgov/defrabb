import sys
from pathlib import Path
from snakemake.utils import validate
import pandas as pd


################################################################################
## Config processing functions


def load_df(path, schema):
    df = pd.read_table(path, dtype={"eval_target_regions": bool})
    validate(df, schema)
    return df


def load_analyses(path, schema):
    return load_df(path, schema).astype(dtype={"eval_target_regions": bool})


def _filter_subtable(df, filter_re, id_cols, new_index):
    params = df.filter(regex=filter_re).drop_duplicates()
    ids = df[[new_index] + id_cols].drop_duplicates()
    tbl = pd.merge(ids, params, how="inner", on=new_index).set_index(new_index)
    return (params, tbl)


def analyses_to_vc_tbl(analyses):
    return _filter_subtable(analyses, "vc_", ["asm_id", "ref"], "vc_id")


def analyses_to_bench_tbls(analyses):
    id_cols = [
        "asm_id",
        "vc_id",
        "vc_cmd",
        "vc_param_id",
        "ref",
        "exclusion_set",
    ]
    params, tbl = _filter_subtable(analyses, "bench_", id_cols, "bench_id")
    excluded_tbl = tbl[tbl.exclusion_set != "none"]
    return (params, tbl, excluded_tbl)


################################################################################
# load config


configfile: workflow.source_path("../config/resources.yml")


validate(config, "../schema/resources-schema.yml")

asm_config = config["assemblies"]
comp_config = config["comparisons"]
ref_config = config["references"]

################################################################################
# init analyses

## Loading analysis table with run information
analyses = load_analyses(
    workflow.source_path(f"../{config['analyses']}"), "../schema/analyses-schema.yml"
)


## Generating seperate tables for individual framework components
## asm variant calls
vc_params, vc_tbl = analyses_to_vc_tbl(analyses)

## draft benchmark set generation
bench_params, bench_tbl, bench_excluded_tbl = analyses_to_bench_tbls(analyses)

## Setting index for analysis run lookup
analyses = analyses.set_index("eval_id")

# Cross-reference validation: catches asm_id/ref/comp_id/exclusion_set
# typos in analyses.tsv before any rule expansion.
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from validate_configs import validate_cross_references

validate_cross_references(config, analyses)

################################################################################
# init wildcard constraints

## Wildcard variables and ids

## Variables for assembly based variant calling
VCIDS = set(vc_tbl.index.tolist())
REFIDS = set(vc_tbl["ref"].tolist())
ASMIDS = set(vc_tbl["asm_id"].tolist())
VCCMDS = set(vc_tbl["vc_cmd"].tolist())
VCPARAMIDS = set(vc_tbl["vc_param_id"].tolist())
SAMPLEIDS = set([asm_config[asm]["sample_id"] for asm in ASMIDS])

## Draft benchmark set generation variables
BENCHIDS = set(bench_tbl.index.tolist())
BENCHTYPS = set(bench_tbl["bench_type"].tolist())


## Evaluations
EVALIDS = set(analyses.index.tolist())
EVALCOMPIDS = set(analyses["eval_comp_id"].tolist())


# Only constrain the wildcards to match what is in the resources file. Anything
# else that can be defined on the command line or in the analyses.tsv can is
# unconstrained (for now).
wildcard_constraints:
    asm_id="|".join(ASMIDS),
    comp_id="|".join(EVALCOMPIDS),
    ref_id="|".join(REFIDS),
    bench_id="|".join(BENCHIDS),
    bench_type="|".join(BENCHTYPS),
    eval_id="|".join(EVALIDS),
    vc_id="|".join(VCIDS),
    vc_cmd="|".join(VCCMDS),
    vc_param_id="|".join(VCPARAMIDS),
    sample_id="|".join(SAMPLEIDS),


## Using zip in rule all to get config sets by config table rows

# defining variables for cleaner rule all
happy_analyses = analyses[analyses["eval_cmd"] == "happy"]
truvari_analyses = analyses[analyses["eval_cmd"] == "truvari"]
truvari_refine_analyses = analyses[analyses["eval_cmd"] == "truvari_refine"]
dipcall_tbl = vc_tbl[vc_tbl["vc_cmd"] == "dipcall"]


################################################################################
# Helper modules (must be included AFTER module-level globals above so the
# helper functions can reference vc_tbl, bench_tbl, ref_config, asm_config,
# REFIDS, etc.)


include: "helpers_ref.smk"
include: "helpers_varcall.smk"
include: "helpers_eval.smk"
include: "helpers_bench.smk"
