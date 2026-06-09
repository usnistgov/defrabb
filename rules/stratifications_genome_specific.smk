################################################################################
##
## Component: Genome-specific stratifications (#59 / #173)
##
## Ports the complex/overlapping-variant stratifications from the GIAB
## genome-stratifications "GenomeSpecific" pipeline. These are derived from the
## draft benchmark VCF (genome-specific) and can be merged into the hap.py
## stratification TSV to break evaluation metrics down by hard-to-call complex
## variant contexts.
##
## OPT-IN: nothing here runs unless its outputs are requested. hap.py only
## consumes these strata when config["genome_specific_strats"] is true (see
## rules/helpers_eval.smk + scripts/run_happy.py), so default/pinned evaluations
## are byte-for-byte unchanged.
##
## NOTE: the SV/CNV stratum (GS rules 7-10: complement of the intermediate
## small-variant benchmark bed, unioned with complex variants and the GIAB
## all-difficult regions) is documented in docs/issues/stratification-59-173-
## design.md but not yet wired here; it needs the defrabb-specific intermediate
## bed and the version-specific GIAB all-difficult strat path.
################################################################################

GENOME_SPECIFIC_STRATA = [
    "comphetsnp10bp",
    "comphetindel10bp",
    "complexindel10bp",
    "snpswithin10bp",
    "othercomplexwithin10bp",
]

_gs_strat_dir = "results/draft_benchmarksets/{bench_id}/genome_specific_strats"
_gs_strat_base = (
    _gs_strat_dir + "/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}."
)


## Combine variants within 10bp into haplotype records (GS rule 1) so that
## compound hets / nearby complex variants can be classified together.
rule genome_specific_geno2haplo:
    input:
        vcf="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
        ref="resources/references/{ref_id}.fa",
        refidx="resources/references/{ref_id}.fa.fai",
    output:
        _gs_strat_base + "geno2haplo10.vcf.gz",
    log:
        "logs/genome_specific_strats/geno2haplo_{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/genome_strats.yml"
    shell:
        """
        ( tmp=$(mktemp --suffix .vcf)
          zcat {input.vcf} > "$tmp"
          vcfgeno2haplo -w 10 -r {input.ref} "$tmp" | bgzip -c > {output}
          rm -f "$tmp" ) &> {log}
        """


## Classify variants into the five complex/overlapping-variant strata
## (GS rules 2-6). Logic is in scripts/genome_specific_strats.py, validated
## byte-identical to the upstream awk/bedtools pipeline on HG002 chr21.
rule genome_specific_complex_strats:
    input:
        geno2haplo=_gs_strat_base + "geno2haplo10.vcf.gz",
        raw="results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",
    output:
        multiext(
            _gs_strat_base,
            "comphetsnp10bp_slop50.bed",
            "comphetindel10bp_slop50.bed",
            "complexindel10bp_slop50.bed",
            "snpswithin10bp_slop50.bed",
            "othercomplexwithin10bp_slop50.bed",
        ),
    log:
        "logs/genome_specific_strats/complex_{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/genome_strats.yml"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        prefix=lambda wildcards: f"{wildcards.ref_id}_{wildcards.asm_id}_{wildcards.bench_type}_{wildcards.vc_cmd}-{wildcards.vc_param_id}.",
    shell:
        """
        python scripts/genome_specific_strats.py \
            --geno2haplo {input.geno2haplo} \
            --raw {input.raw} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            &> {log}
        """


## Assemble the genome-specific stratification TSV (name<TAB>bed, paths relative
## to the TSV) that hap.py consumes alongside the GIAB strat TSV.
rule genome_specific_strat_tsv:
    input:
        beds=multiext(
            _gs_strat_base,
            "comphetsnp10bp_slop50.bed",
            "comphetindel10bp_slop50.bed",
            "complexindel10bp_slop50.bed",
            "snpswithin10bp_slop50.bed",
            "othercomplexwithin10bp_slop50.bed",
        ),
    output:
        _gs_strat_base + "genome_specific_strats.tsv",
    log:
        "logs/genome_specific_strats/tsv_{bench_id}_{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.log",
    conda:
        "../envs/genome_strats.yml"
    params:
        strat_args=lambda wildcards, input: " ".join(
            f"--strat {stratum}={bed}"
            for stratum, bed in zip(GENOME_SPECIFIC_STRATA, input.beds)
        ),
    shell:
        """
        python scripts/build_stratification_tsv.py \
            {params.strat_args} \
            -o {output} \
            &> {log}
        """
