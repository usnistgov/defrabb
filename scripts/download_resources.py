from email.message import Message
from logging import warning
import os
import pathlib
import tempfile
from re import sub
from warnings import WarningMessage
from snakemake.common import get_file_hash
from snakemake.utils import makedirs
from snakemake.shell import shell
from snakemake.remote import AUTO
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

S3 = S3RemoteProvider(keep_local=True)


log = snakemake.log_fmt_shell(stdout=False, stderr=True)

source_uri = snakemake.params.source_uri
source_hash = snakemake.params.source_hash
hash_algorithm = snakemake.params.hash_algo
outfile = snakemake.output
outfmt = snakemake.params.outfmt  # alternatively use decompress

## Downloading File
in_path = pathlib.Path(source_uri)
infile = in_path.name
infile_ext = in_path.suffix

tmpremote = tempfile.mktemp(suffix=infile_ext)
### Determining download protocol
if "s3://" in source_uri:
    ## Issue with aws s3 command on workstation, using curl only for now
    shell("aws s3 cp {source_uri} {tmpremote}")
elif "s3.amazonaws.com" in source_uri:
    # shell("aws s3 cp {source_uri} {tmpremote}")
    ## Issue with aws s3 command on workstation, using curl only for now
    shell("curl --insecure -L {source_uri} > {tmpremote}")
elif "https://" in source_uri or "ftp://" in source_uri:
    shell("curl --insecure -L {source_uri} > {tmpremote}")
else:
    print(
        "Assuming remote from S3, ftp, or https.",
        "URI does not contain expected protocol.",
        f"Downloading using {source_uri} curl.",
    )

## MD5 check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## - getting md5 of remote file
local_hash = get_file_hash(tmpremote, algorithm=hash_algorithm)

## - Verifying checksum
if local_hash != source_hash:
    print(
        f"Checksum for {infile} not match the expected hash, please verify {hash_algorithm} was used to compute digest\n"
    )


## Decompressing or compressing file as appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~
if infile_ext == ".gz" and (outfmt == "bgzip" or outfmt == "gzip"):
    print("moving compressed file to output")
    ## Keeping file compressed and moving file to output path
    shell("cp {tmpremote} {outfile}")
elif infile_ext != ".gz" and outfmt == "decompress":
    print("moving uncompressed file to output")
    ## Keeping file uncompressed
    shell("mv {tmpremote} {outfile}")
elif infile_ext == ".gz" and outfmt == "decompressed":
    print("uncompressing file and moving to output")
    ## decompressing file and moving to output path
    shell("gunzip -cf {tmpremote} > {outfile}")
elif infile_ext != ".gz" and outfmt == "bgzip":
    print("compressing file and moving to output")
    ## bgzip compressing file and moving to specified output directory
    shell("bgzip -fc {tmpremote} > {outfile}")
elif infile_ext == ".bed" and outfmt == "decompressed":
    print("moving uncompressed file to output")
    ## bgzip compressing file and moving to specified output directory
    shell("cp {tmpremote} {outfile}")
else:
    print(f"no rename step applied, bug in code for downloading {outfile}")
