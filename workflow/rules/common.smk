################################################################################
## Functions to get assembly haplotype paths from tsv
################################################################################


def _get_hap(key, wildcards):
    path = asm.loc[(wildcards.prefix), key]
    return path


def get_hap1(wildcards):
    return _get_hap("paternal_haplotype", wildcards)


def get_hap2(wildcards):
    return _get_hap("maternal_haplotype", wildcards)
