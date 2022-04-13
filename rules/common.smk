################################################################################
## Utility Functions

## Loading and validating analysis tables
def get_analyses(path):
    # target_regions must be a string even though it might only contain
    # 'boolean' values
    analyses = pd.read_table(path, dtype={"eval_target_regions": str})
    validate(analyses, "schema/analyses-schema.yml")

    try:
        return analyses.set_index("eval_id", verify_integrity=True)
    except ValueError:
        print("All keys in column 'eval_id' must by unique")


## Generating concatenated string for wildcard constraints
def format_constraint(xs):
    return "|".join(set(xs))


################################################################################
## Rule parameters
def get_genome_file(wildcards):
    ## get ref_id from prefix
    for id in REFIDS:
        if id in wildcards.prefix:
            return f"resources/references/{id}.genome"
    print("ref_id not found in bed file prefix")

def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    root = config["_par_bed_root"]
    par_path = Path(root) / ref_config[wildcards.ref_id]["par_bed"]
    return f"-x {str(par_path)}" if is_male else ""


## Happy Inputs and Parameters
def get_happy_inputs(wildcards):
    ## Creating empty dictionary for storing inputs
    inputs = {}

    ## Reference genome and stratifications
    ref_id = analyses.loc[wildcards.eval_id, "ref"]
    inputs["genome"] = f"resources/references/{ref_id}.fa"
    inputs["genome_index"] = f"resources/references/{ref_id}.fa.fai"
    strat_tb = config["references"][wildcards.ref_id]["stratifications"]["tarball"]
    inputs["strat_tb"] = f"resources/strats/{ref_id}/{strat_tb}"

    ## draft benchmark variant calls
    draft_bench_vcf = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz"

    ## draft benchmark regions
    if analyses.loc[(wildcards.eval_id, "exclusion_set")] == "none":
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed"
    else:
        draft_bench_bed = "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.excluded.bed"

    ## comparison variant call paths
    comp_vcf = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz"
    comp_vcfidx = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi"
    comp_bed = "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed"

    ## Determining which callsets and regions are used as truth
    if analyses.loc[wildcards.eval_id, "eval_comp_id_is_truth"] == True:
        query = "draft_bench"
    else:
        query = "comp"

    ## Defining truth calls and regions along with query calls
    if query == "draft_bench":
        inputs["query"] = draft_bench_vcf
        inputs["truth"] = comp_vcf
        inputs["truth_vcfidx"] = comp_vcfidx
        inputs["truth_regions"] = comp_bed
    else:
        inputs["query"] = comp_vcf
        inputs["query_vcfidx"] = comp_vcfidx
        inputs["truth"] = draft_bench_vcf
        inputs["truth_regions"] = draft_bench_bed

    ## Determining Target regions
    trs = eval_tbl.loc[wildcards.eval_id, "eval_target_regions"]
    if trs.lower() != "false":
        if trs.lower() == "true":
            if query == "draft_bench":
                inputs["target_regions"] = draft_bench_bed
            else:
                inputs["target_regions"] = comp_bed
        else:
            inputs["target_regions"] = f"resources/manual/target_regions/{trs}"

    ## Returning happy inputs
    return inputs


## Exclusions
def get_exclusion_inputs(wildcards):

    ## Getting list of excluded regions
    exclusion_set_id = bench_tbl.loc[wildcards.bench_id, "exclusion_set"]
    exclusion_set = config["exclusion_set"][exclusion_set_id]

    ## Initiating empty list for storing paths for beds to excluded from
    ## diploid assembled regions
    exc_paths = []
    for exclusion in exclusion_set:

        ## Determining path for asm specific exclusions and asm agnostic exclusions
        if exclusion in config["exclusion_asm_agnostic"]:
            exc_path = f"resources/exclusions/{{ref_id}}/{exclusion}"
        else:
            exc_path = f"results/draft_benchmarksets/{{bench_id}}/exclusions/{{ref_id}}_{{asm_id}}_{{vc_cmd}}-{{vc_param_id}}_{exclusion}"

        ## Adding slop - currently a 15kb hard coded buffer around excluded repeat regions
        if exclusion in config["exclusion_slop_regions"]:
            exc_path = f"{exc_path}_slop"
        elif  exclusion in config["exclusion_asm_intersect"]:
            ## Ensuring bed files are sorted before intersect
            exc_path = f"{exc_path}_sorted"

        ## Defining which regions are excluded based on diploid assembly breaks
        if exclusion in config["exclusion_asm_intersect"]:
            exc_paths += [f"{exc_path}_start", f"{exc_path}_end"]
        else:
            exc_paths = exc_paths + [exc_path]

    ## Adding to exc_paths list and ensuring all beds are sorted
    ## prior to exclusion from dip assembled regions
    exclusion_paths = [f"{exc}_sorted.bed" for exc in exc_paths]

    ## Returning list of bed paths for exclusion
    return exclusion_paths
