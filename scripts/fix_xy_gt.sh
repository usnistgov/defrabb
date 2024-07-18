#!/bin/bash

########### Fix XY genotype ####################################################
# Fix genotypes to match VCF specifications
# Convert ./1, 1|., and .|1 variant genotypes in non-PAR XY to 1
#
# Approach
# 1. Extract non-PAR XY variants and all others
# 2. Convert passing variant genotypes in non-PAR XY
# 3. Concatenate other variants vcf and non-PAR XY gt fixed VCF, sort and index
#
# QC: XY GT counts before and after fix
#
# Dependencies
## Tools
# - bcftools
# - Others: sed, awk, ...
#
## Input
# - VCF file to fix
# - PAR XY bed
# - genome file

# Usage
# ./fix_xy_gt.sh -i <input.vcf.gz> -o <output.vcf.gz> -p <non-par-xy.bed> -g <genome>

# Author: ND Olson
# Date: 7/10/2024

# Parse arguments
while getopts ":i:o:p:g:" opt; do
    case $opt in
    i)
        input_vcf="$OPTARG" # input vcf
        ;;
    o)
        output_vcf="$OPTARG" # output vcf
        ;;
    p)
        parbed="$OPTARG" # non-PAR XY bed
        ;;
    g)
        genome="$OPTARG" # ref genome file
        ;;
    \?)
        echo "Invalid option -$OPTARG" >&2
        ;;
    esac
done

printf "\n\nProvenance Info...\n" >&2
echo "Input VCF: $(md5sum ${input_vcf})" >&2
echo "PAR XY bed: $(md5sum ${parbed})" >&2
echo "Ref Genome File: $(md5sum ${genome})" >&2
echo "Output VCF: $output_vcf" >&2

echo "Tool versions: " >&2
bcftools --version >&2

## Temp directory for intermediates
TEMPDIR=$(mktemp -d)
echo "Temp dir: $TEMPDIR" >&2
xy=${TEMPDIR}/xy.bed
nonparxy=${TEMPDIR}/nonparxy.bed
nonparxy_comp=${TEMPDIR}/nonparxy_comp.bed
nonparxyvcf=${TEMPDIR}/non-par-xy.vcf.gz
nonparxycompvcf=${TEMPDIR}/non-par-xy_comp.vcf.gz
nonparxyvcffixed=${TEMPDIR}/non-par-xy_fixed.vcf.gz

# Check inputs
[ -f ${input_vcf} ] || {
    echo "Input VCF file not found"
    exit 1
}
[ -f ${parbed} ] || {
    echo "PAR XY bed file not found"
    exit 1
}
[ -f ${genome} ] || {
    echo "Ref genome file not found"
    exit 1
}

# Creating non-par xy bed and compliment
grep -P "[XY]\t" ${genome} \
    | awk '{OFS="\t"} {print $1,0,$2}' \
    > ${xy}

bedtools subtract \
    -a ${xy} \
    -b ${parbed} \
    > ${nonparxy}

bedtools complement \
    -i ${nonparxy} \
    -g ${genome} \
    > ${nonparxy_comp}



# Extract non-PAR XY variants
printf "\n\nSpliting input vcf..." >&2
if [ ! -f ${input_vcf}.tbi ]; then
    bcftools index --tbi ${input_vcf}
fi

bcftools view \
    -R "$nonparxy" \
    ${input_vcf} \
    -Oz -o ${nonparxyvcf}
bcftools index -f ${nonparxyvcf}

# Extract all other variants
bcftools view \
    -R ${nonparxy_comp} \
    ${input_vcf} \
    -Oz -o ${nonparxycompvcf}
bcftools index -f ${nonparxycompvcf}

## Sanity check
printf "\n\nSanity check..." >&2
## Comparing total XY input variants to variants in the vcf subsets
echo "Input XY variants:" >&2
bcftools index -s -a ${input_vcf} | grep -P "^chr[XY]"$'\t'
echo "non-PAR XYvariants:" >&2
bcftools index -s -a ${nonparxyvcf}
echo "Others:" >&2
bcftools index -s -a ${nonparxycompvcf} | grep "^chrX"

# Fix genotypes in non-PAR XY
printf "\n\nFixing genotypes in non-PAR XY...\n" >&2
bcftools +setGT -Oz -o ${nonparxyvcffixed} \
    ${nonparxyvcf} \
    -- -n c:'M' -t ./x
bcftools index ${nonparxyvcffixed}

# Concatenate other variants and non-PAR XY fixed
printf "\n\nConcatenating other variants and non-PAR XY fixed...\n" >&2
bcftools concat -a \
    ${nonparxycompvcf} \
    ${nonparxyvcffixed} |
    bcftools sort -Oz -o "$output_vcf"

bcftools index --tbi "$output_vcf"

# QC - number of par and non-par XY variants in input and output VCF
#      by GT and FILTER
printf "\n\nQC GT Counts Pre and Post Fix:..." >&2
for i in $input_vcf $output_vcf; do
    vcfname=$(basename ${i} .vcf.gz)

    # PAR and non-PAR XY genotypes
    bcftools query -R ${parbed} \
        -f '%CHROM\_par\t%FILTER[\t%GT]\n' \
        ${i} >${TEMPDIR}/${vcfname}_xy_gt.tsv

    bcftools query -R ${nonparxy} \
        -f '%CHROM\_nonpar\t%FILTER[\t%GT]\n' \
        ${i} >>${TEMPDIR}/${vcfname}_xy_gt.tsv

    printf "\n\nchrom	filter	gt	${vcfname}\n" >&2
    sort ${TEMPDIR}/${vcfname}_xy_gt.tsv |
        awk 'BEGIN {OFS="\t"}; {!seen[$0]++}END{for (i in seen) print i, seen[i]}' |
        sort >&2
done

# Clean up
printf "\n\nCleaning up...\n" >&2
rm -f $TEMPDIR/*

echo "Done." >&2
