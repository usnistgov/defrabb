#!/usr/bin/env python3
"""Generate report parameters YAML from snakemake inputs and config."""

from pathlib import Path
import yaml
from snakemake.script import snakemake

# Resolve analysis table path
config_data = snakemake.config
analysis_path = Path(snakemake.config["analyses"]).expanduser()
analysis_table = (
    str(analysis_path)
    if analysis_path.is_absolute()
    else f"../{snakemake.config['analyses']}"
)

# Convert input items to string lists
inputs_converted = {}
for key, value in snakemake.input.items():
    if isinstance(value, list):
        inputs_converted[key] = [str(v) for v in value]
    else:
        inputs_converted[key] = [str(value)]

# Build report params dict
report_params_dict = {
    "inputs": inputs_converted,
    "config_data": config_data,
    "analysis_table": analysis_table,
    "variables": dict(snakemake.params),
}

# Write to YAML
with open(snakemake.output.report_params, "w") as f:
    yaml.safe_dump(report_params_dict, f, default_flow_style=False, sort_keys=False)
