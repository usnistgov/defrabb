## Reference / sample / PAR / segdup helpers
def get_ref_id(wildcards):
    ref_id = ""
    if wildcards.get("ref", ""):
        return wildcards.get("ref", "")
    if wildcards.get("ref_id", ""):
        return wildcards.get("ref_id", "")
    if wildcards.get("vc_id", ""):
        return vc_tbl.loc[wildcards.vc_id]["ref"]
    if wildcards.get("prefix", ""):
        prefix = wildcards.get("prefix", "")
        for id in REFIDS:
            if id in prefix:
                return id
    try:
        ref_id != ""
    except:
        return f"Ref ID could not be determined from {wildcards}"


def get_ref_file(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.fa"


def get_ref_index(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.fa.fai"


def get_ref_bwaindex(wildcards):
    ref_id = get_ref_id(wildcards)
    return [
        f"resources/references/{ref_id}.fa.{ext}"
        for ext in ["amb", "ann", "bwt", "pac", "sa"]
    ]


def get_genome_file(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.genome"


def get_ref_sdf(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}.sdf"


def get_ref_trdb(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/references/{ref_id}_adotto_trf.bed.gz"


def get_addoto_tr_anno_db_url(wildcards):
    return ref_config[wildcards.ref_id]["adotto_db"]


def get_trf_db_url(wildcards):
    return ref_config[wildcards.ref_id]["trf_db"]


def get_sample_id(wildcards):
    return asm_config[wildcards.asm_id]["sample_id"]


def get_par_bed(wildcards):
    root = config["_par_bed_root"]
    ref_id = get_ref_id(wildcards)
    filename = ref_config[ref_id]["par_bed"]
    return Path(root) / filename


def get_dipcall_par_param(wildcards):
    is_male = asm_config[wildcards.asm_id]["is_male"]
    par_path = get_par_bed(wildcards)
    return f"-x {par_path}" if is_male else ""


def get_segdups(wildcards):
    ref_id = get_ref_id(wildcards)
    return f"resources/exclusions/{ref_id}/segdups_slopmerge_sorted.bed"
