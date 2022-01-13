## Make Test Stratification Set
for strat in GRCh38/ancestry/GRCh38_ancestry_Neanderthal.bed.gz \
	GRCh38/union/GRCh38_alldifficultregions.bed.gz \
	GRCh38/union/GRCh38_notinalldifficultregions.bed.gz \
	GRCh38/OtherDifficult/GRCh38_population_CNV_FP_regions.bed.gz ; 
	do
		bedtools intersect -a ${strat} -b test/resources/chr21.bed \
			> test/resources/stratifications/${strat}
	done

## Commands used to generate ch21 reference
## Region coordinates manually extracted from GRCh38_chr21.fa.fai
echo "chr21:1-46709983" > GRCh38_chr21.txt
samtools faidx -o GRCh38_chr21.fa -r GRCh38_chr21.txt GRCh38.fa