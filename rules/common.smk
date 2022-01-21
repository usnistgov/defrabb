################################################################################
## Utility Functions

## Loading and validating analysis tables
def get_analyses(path):
    # target_regions must be a string even though it might only contain
    # 'boolean' values
    analyses = pd.read_table(path, dtype={"target_regions": str})
    validate(analyses, "schema/analyses-schema.yml")

    try:
        return analyses.set_index("bench_id", verify_integrity=True)
    except ValueError:
        print("All keys in column 'bench_id' must by unique")


## Generating concatenated string for wildcard constraints
def format_constraint(xs):
    return "|".join(set(xs))

################################################################################
## Rule parameters

def get_male_bed(wildcards):
    is_male = asm_config[wildcards.asm_prefix]["is_male"]
    root = config["_par_bed_root"]
    par_path = Path(root) / ref_config[wildcards.ref_prefix]["par_bed"]
    return f"-x {str(par_path)}" if is_male else ""

## Happy Inputs and Parameters
def get_happy_inputs(wildcards):
    ## Creating empty dictionary for storing inputs
    inputs = {}

    ## asm variant call output - TODO link to post processing output
    asm_vcf = "results/dipcall-{bench_id}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.vcf.gz"
    #asm_vcfidx = "results/dipcall/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.vcf.gz.tbi"

    ## Asm regions
    if analyses.loc[(wildcards.bench_id, "exclusion_set")] == "none":
        asm_bed = "results/dipcall-{bench_id}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/dipcall.dip.bed"
    else: 
        asm_bed = "results/dipcall-{bench_id}/{ref_prefix}_{asm_prefix}_{varcaller}-{vc_param_id}/exclusions/excluded.bed"

    ## comparison variant call paths
    comp_id = analyses.loc[wildcards.bench_id, "compare_var_id"]
    comp_vcf = f"resources/benchmarks/{comp_id}.vcf.gz"
    comp_vcfidx = f"resources/benchmarks/{comp_id}.vcf.gz.tbi"
    comp_bed = f"resources/benchmarks/{comp_id}.bed"

    ## Determining which callsets and regions are used as truth
    if analyses.loc[wildcards.bench_id, "vcr_is_query"] == "true":
        query = "asm"
    else:
        query = "comp"
    
    ## Defining truth calls and regions along with query calls
    if query == "asm":
        inputs["query"] = asm_vcf
        #inputs["query_vcfidx"] = asm_vcfidx
        inputs["truth"] = comp_vcf
        inputs["truth_vcfidx"] = comp_vcfidx
        inputs["truth_regions"] = comp_bed
    else:
        inputs["query"] = comp_vcf
        inputs["query_vcfidx"] = comp_vcfidx
        inputs["truth"] = asm_vcf
        #inputs["truth_vcfidx"] = asm_vcfidx
        inputs["truth_regions"] = asm_bed

    ## Determining Target regions
    trs = analyses.loc[(wildcards.bench_id, "target_regions")]
    if trs != "false":
        if trs == "true":
            if query == "asm":
                inputs["target_regions"] = asm_bed
            else:
                inputs["target_regions"] = comp_bed
        else:
            inputs["target_regions"]=f"resources/manual/target_regions/{trs}"

    ## Returning happy inputs
    return inputs

## Exclusions
def lookup_excluded_region_set(wildcards):
    xset = analyses.loc[(wildcards.bench_id, "exclusion_set")]
    return [f"results/dipcall-{{bench_id}}/{{ref_prefix}}_{{asm_prefix}}_{{varcaller}}-{{vc_param_id}}/exclusions/{p}" for p in config["exclusion_set"][xset]]