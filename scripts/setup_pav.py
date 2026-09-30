"""Generate PAV3 input files (assemblies.tsv + pav.json) for a run_pav job.

Split out of the former ``run_pav.py``: this runs on the *host* (no container),
so Snakemake's ``script:`` machinery uses the host Snakemake. The container step
(the nested PAV Snakemake) is now a plain ``shell:`` directive in ``run_pav``.

Previously the whole thing was a single ``script:`` rule with a
``container:``. Snakemake injects a host-generated preamble
(``from snakemake.script import snakemake``) that then executes *inside* the
becklab/pav container, which ships its own older Snakemake. The version skew
crashed at import with::

    ModuleNotFoundError: No module named 'snakemake.io.container';
    'snakemake.io' is not a package

You cannot reliably run Snakemake's ``script:`` machinery inside a container
whose Snakemake differs from the host's. Generating the config on the host and
invoking PAV with a bare ``shell:`` avoids the injection entirely. See
``docs/issues/run_pav_run_dipcall_failures.md``.

PAV3 (``pav3 batch``) input format: assembly table columns are ``name`` plus
one ``hap_<NAME>`` column per haplotype; column order sets the VCF GT allele
order. Paternal is ``hap_h1`` and maternal ``hap_h2`` (GT = pat|mat, same as the
PAV2 HAP1/HAP2 layout). ``pav.json`` holds ``reference`` plus any PAV3 config
parameters from the selected ``_pav_config`` profile.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import os
import json
import pandas as pd

# PAV runs with its working directory == outdir; all paths in the generated
# files are therefore relative to outdir.
outdir = snakemake.params.outdir

## Assemblies table -----------------------------------------------------------
df = pd.DataFrame(
    {
        "name": [snakemake.params.name],
        "hap_h1": [os.path.relpath(snakemake.input.hap1, outdir)],
        "hap_h2": [os.path.relpath(snakemake.input.hap2, outdir)],
    }
)
df.to_csv(os.path.join(outdir, "assemblies.tsv"), sep="\t", index=False)

## PAV pav.json (+ reference) ------------------------------------------------
config_dict = dict(snakemake.params.pav_config or {})
config_dict["reference"] = os.path.relpath(snakemake.input.ref, outdir)
with open(os.path.join(outdir, "pav.json"), "w") as fh:
    json.dump(config_dict, fh, default=lambda o: o.__dict__, sort_keys=True, indent=2)
