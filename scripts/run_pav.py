import os
from snakemake.script import snakemake
from snakemake import shell
import pandas as pd
import json

## Defining directory paths
outdir = snakemake._params_store.outdir

## Generating assemblies table
dat = {
    "NAME": [snakemake._params_store.name],
    "HAP1": [os.path.relpath(snakemake.input.hap1, outdir)],
    "HAP2": [os.path.relpath(snakemake.input.hap2, outdir)],
}
df = pd.DataFrame(dat)

# Write assemblies table
df.to_csv(f"{outdir}/assemblies.tsv", sep="\t", index=False)

## Generating pav config.json ------------------------------------------
config_dict = snakemake._params_store.pav_config
## Adding reference file
# Ref

config_dict["reference"] = os.path.relpath(snakemake.input.ref, outdir)

## write config to json
with open(f"{outdir}/config.json", "w") as file:
    json_string = json.dumps(
        config_dict, default=lambda o: o.__dict__, sort_keys=True, indent=2
    )
    file.write(json_string)

## calling pav ----------------------------------------------------------
# OPENSSL_CONF=/dev/null disables OpenSSL FIPS self-test for the inner
# snakemake invocation. NIST FIPS hosts trigger "FATAL FIPS SELFTEST FAILURE"
# inside the becklab/pav container otherwise. Same pattern is used for
# `truvari anno trf` in rules/bench_vcf_processing.smk.
#
# The nested PAV snakemake output is captured to the rule log (it was
# previously lost to the parent's stdout/stderr, making PAV failures
# undebuggable). On failure we also surface PAV's own inner .snakemake logs.
log_path = os.path.abspath(snakemake.log[0])
shell(
    "cd {snakemake._params_store.outdir}; "
    "OPENSSL_CONF=/dev/null snakemake -s /opt/pav/Snakefile --ri -k -w 20 "
    "--rerun-triggers mtime -c {snakemake.threads} --config ignore_env_file=True "
    "> {log_path} 2>&1 "
    "|| {{ echo '--- inner PAV .snakemake/log ---' >> {log_path}; "
    "cat .snakemake/log/*.log >> {log_path} 2>/dev/null; exit 1; }}"
)
