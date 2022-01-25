## TODO Test
rule run_assembly_stats:
    input:
        #Input assembly
        assembly="resources/assemblies/{asm_id}/{haplotype}.fa",
    output:
        #Assembly statistics
        assembly_stats="results/report/assemblies/{asm_id}_{haplotype}_stats.txt",
    params:
        # Tab delimited output, with a header, is set as the default. Other options are available:
        #   -l <int>
        #       Minimum length cutoff for each sequence.
        #       Sequences shorter than the cutoff will be ignored [1]
        #   -s
        #       Print 'grep friendly' output
        #   -t
        #       Print tab-delimited output
        #   -u
        #       Print tab-delimited output with no header line
        # If you want to add multiple options just delimit them with a space.
        # Note that you can only pick one output format
        # Check https://github.com/sanger-pathogens/assembly-stats for more details
        extra="-t",
    log:
        "logs/run_assembly_stats/{asm_id}_{haplotype}.assembly-stats.log",
    threads: 1
    wrapper:
        "v0.86.0/bio/assembly-stats"

## Potentially modify to calculate for all output beds
rule get_bed_size:
	input: "{genomic_region}.bed"
	output: report("results/bed_size/{genomic_region}.txt", caption = "report/bed_size.rst", category = "Exclusion Stats")
	log: "logs/get_bed_size/{genomic_region}.log"
	shell: """
		cat {input} \
			| '{sum+=$3-$2} END {print sum} \
			1> {output} 2> {log}
		"""
## Genome Coverage
## TODO Modify for benchmark set development framework
# rule make_bench_cov_tbls:
#     input:
#         a=ensembl_dir + "/{ref}_mrg_full_{region}.bed",
#         b=benchdir + "/HG002_{ref}_{benchmarkset}_{benchtype}.bed"
#     output: "workflow/results/bench_cov_tbls/HG002_{ref}_{benchmarkset}_{benchtype}_{region}_cov.tsv"
#     threads: 2
#     wrapper: "0.74.0/bio/bedtools/coveragebed"

## TODO Modify for benchmark set development framework
# rule make_strat_cov_tbls:
#     input:
#         a= ensembl_dir + "/{ref}_mrg_full_{region}.bed",
#         b= "workflow/data/strats/{ref}_{strat}.bed.gz"
#     output: "workflow/results/strat_cov_tbls/{strat}_{ref}_mrg_{region}_cov.tsv"
#     threads: 2
#     wrapper: "0.74.0/bio/bedtools/coveragebed"

## Variant Callset Stats
rule get_vcf_stats:
	input: "{variant_callset}.vcf.gz"
	output: 
		stats="results/report/{variant_callset}_stats.txt"
	log: "logs/get_vcf_stats/{variant_callset}.log"
	conda: "envs/bcftools.yml"
	shell: """
		bcftools stats {input} > {output.txt}
	"""

### Small Var Table
## TODO Modify for benchmark set development framework
# rule make_smallvar_tbls:
#     input: 
#         anno_vcf="workflow/data/anno_vcf/HG002_{ref}_{benchmarkset}_smallvar_anno.vcf"
#     output: "workflow/results/anno_vcf_tbls/HG002_{ref}_{benchmarkset}_smallvar_anno.tsv"
#     conda: "envs/bcftools.yml"
#     shell: """
#         bcftools query \
#             -f '%CHROM\\t%POS\\t[ %GT]\\t%TYPE\\t%INFO/GENE\\t%INFO/EXON\\t%REF\\t%ALT\\n' \
#             {input.anno_vcf} \
#             > {output}
#     """

## Structural Var Table
## TODO Modify for benchmark set development framework
# rule make_sv_tbls:
#     input: 
#         anno_vcf="workflow/data/anno_vcf/HG002_{ref}_{benchmarkset}_SV_anno.vcf"
#     output: "workflow/results/anno_vcf_tbls/HG002_{ref}_{benchmarkset}_SV_anno.tsv"
#     conda: "envs/bcftools.yml"
#     shell: """
#         bcftools query \
#             -f '%CHROM\\t%POS\\t[ %GT]\\t%TYPE\\t%REPTYPE\\t%BREAKSIMLENGTH\\t%INFO/GENE\\t%INFO/EXON\\t%REF\\t%ALT\\n' \
#             {input.anno_vcf} \
#             > {output}
#     """