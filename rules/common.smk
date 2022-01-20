################################################################################
## Utility Functions

## Loading and validating analysis tables
def get_analyses(path):
    # target_regions must be a string even though it might only contain
    # 'boolean' values
    analyses = pd.read_table(path, dtype={"target_regions": str})
    validate(analyses, "config/analyses-schema.yml")

    try:
        return analyses.set_index("bench_id", verify_integrity=True)
    except ValueError:
        print("All keys in column 'bench_id' must by unique")


## Generating concatenated string for wildcard constraints
def format_constraint(xs):
    return "|".join(set(xs))


################################################################################
## Functions for Setting up variable paths
################################################################################

# Define what files we want hap.py to make, and these paths will contain the
# definitions for the assemblies, variant caller, etc to use in upstream rules.


def expand_bench_output(path, cmd):
    bps = analyses[analyses["bench_cmd"] == cmd].index.tolist()
    return expand(path, bench_prefix=bps)

################################################################################
## Rule parameters

def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_prefix]["is_male"]
    root = config["_par_bed_root"]
    par_path = Path(root) / ref_config[wildcards.ref_prefix]["par_bed"]
    return f"-x {str(par_path)}" if is_male else ""