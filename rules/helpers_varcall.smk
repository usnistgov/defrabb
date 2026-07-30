## PAV / dipcall variant-call helpers
def get_pav_basename(wildcards):
    vc_id = wildcards.vc_id
    asm_id = vc_tbl.loc[vc_id, "asm_id"]
    sample_id = asm_config[asm_id]["sample_id"]
    base_name = f"results/asm_varcalls/{vc_id}/{sample_id}"
    return base_name


def get_pav_hap1_bed(wildcards):
    base_name = get_pav_basename(wildcards.vc_id)
    return f"{base_name}/callable_regions_h1_500.bed.gz"


def get_pav_hap2_bed(wildcards):
    base_name = get_pav_basename(wildcards.vc_id)
    return f"{base_name}/callable_regions_h2_500.bed.gz"


def get_pav_outputs(wildcards):
    base_name = get_pav_basename(wildcards)
    outdict = {
        "vcf": f"{base_name}.vcf.gz",
        "vcfidx": f"{base_name}.vcf.gz.tbi",
        "bed": f"{base_name}.diploid_regions.bed",
    }
    return outdict


def get_dipcall_basename(wildcards):
    vc_id = wildcards.vc_id
    ref_id = wildcards.ref_id
    asm_id = wildcards.asm_id
    vc_cmd = wildcards.vc_cmd
    vc_param_id = wildcards.vc_param_id
    return f"results/asm_varcalls/{vc_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.dip"


def get_dipcall_outputs(wildcards):
    base_name = get_dipcall_basename(wildcards)
    return {
        "vcf": f"{base_name}.rename.vcf.gz",
        "vcfidx": f"{base_name}.rename.vcf.gz.tbi",
        "bed": f"{base_name}_sorted.bed",
    }


def is_pav(wildcards):
    return wildcards.vc_cmd == "pav"


from varcall_lookup import resolve_asm_varcall_run


def get_asm_varcall_run(ref, asm_id, vc_cmd, prefer_vc_param_id=None):
    """Resolve the existing asm_varcall run for a (ref, asm_id, vc_cmd).

    Used by exclusions that need variant calls from *another* assembly variant
    caller (e.g. ``consecutive-svs`` needs dipcall hap BAMs). Rather than
    triggering a duplicate caller run inside the current benchmark's vc_id
    directory, this reuses the results from the matching run already declared in
    the analyses table.

    If multiple runs match (parameter sweeps), prefer ``prefer_vc_param_id``.

    Returns ``(vc_id, vc_param_id)`` for the match.
    """
    return resolve_asm_varcall_run(vc_tbl, ref, asm_id, vc_cmd, prefer_vc_param_id)
