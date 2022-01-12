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

