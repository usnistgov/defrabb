import boto3
import os

# Initialize the S3 client
s3 = boto3.client('s3')
bucket_name = 'giab-data'
prefix = 'defrabb_runs/20240523_v0.016_HG002Q100/'

# Function to determine the analysis type
def determine_analysis(key):
    subdirs = key.split('/')
    if "results" in key:
        return subdirs[3]
    if len(subdirs) > 3:
        return subdirs[2]
    return 'NA'

# Function to determine the file type
def determine_file_type(key):
    if key.endswith('tar'):
        return 'snakemake_archive'
    if 'README' in key and key.endswith('.md'):
        return 'README'
    if 'config/analyses_' in key and key.endswith('.tsv'):
        return 'analysis table'
    if 'config/resources.yml' in key:
        return 'config'
    if 'run.log' in key:
        return 'run.log'
    if 'environment.yml' in key:
        return 'mamba env'
    if 'resources/comparison' in key:
        return 'comparison vcfs'
    if 'resources/exclusions' in key:
        return "exclusion bed"
    if 'results/asm_varcalls' in key and key.endswith('dip.bed'):
        return "dip.bed"
    if 'results/asm_varcalls' in key and key.endswith('hap1.bam'):
        return "hap1 bam"
    if 'results/asm_varcalls' in key and key.endswith('hap1.bam.bai'):
        return "hap1 bam idx"
    if 'results/asm_varcalls' in key and key.endswith('hap2.bam'):
        return "hap2 bam"
    if 'results/asm_varcalls' in key and key.endswith('hap2.bam.bai'):
        return "hap2 bam idx"
    if 'draft_benchmarksets' in key and 'bench-vars' not in key and key.endswith('vcf.gz'):
        if "smvar" in key:
            return "smvar benchmark vcf"
        if "stvar" in key:
            return "stvar benchmark vcf"
    if 'draft_benchmarksets' in key and 'bench-vars' not in key and key.endswith('vcf.gz.tbi'):
        if "smvar" in key:
            return "smvar benchmark vcf idx"
        if "stvar" in key:
            return "stvar benchmark vcf idx"
    if 'draft_benchmarksets' in key and key.endswith('benchmark.bed'):
        if "smvar" in key:
            return "smvar benchmark bed"
        if "stvar" in key:
            return "stvar benchmark bed"
    else:
        return "other files"

# Function to determine ref_id
def determine_ref_id(key):
    for ref in ['GRCh37', 'GRCh38', 'GRCh38-GIABv3', 'v1.0.1mat', 'v1.0.1pat', 'v1.0.1dip']:
        if ref in key:
            return ref.replace('v1.0.1', 'HG002Q100v1.0.1')
    return 'NA'

# Open a file to write the TSV data
with open('data_manifest.tsv', 'w') as file:
    file.write("analysis\tfile_type\tref_id\tsize_Gb\ts3_uri\turl\n")

    # Initialize variables for pagination
    continuation_token = None
    is_truncated = True

    while is_truncated:
        # Fetch the list of files from S3 with pagination
        if continuation_token:
            response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix, ContinuationToken=continuation_token)
        else:
            response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

        for item in response.get('Contents', []):
            key = item['Key']
            size_gb = item['Size'] / (1024 ** 3)

            analysis = determine_analysis(key)
            file_type = determine_file_type(key)
            ref_id = determine_ref_id(key)
            s3_uri = f's3://{bucket_name}/{key}'
            url = f'https://{bucket_name}.s3.amazonaws.com/{key}'

            file.write(f"{analysis}\t{file_type}\t{ref_id}\t{size_gb:.2f}\t{s3_uri}\t{url}\n")

        # Check if the response is truncated and set the continuation token
        is_truncated = response.get('IsTruncated', False)
        continuation_token = response.get('NextContinuationToken', None)
