import yaml
import shutil
from pathlib import Path

print("starting quarto_params.py script")
config_data = snakemake.config
analysis_table = f"../{snakemake.config['analyses']}"

inputs_converted = {}
for key, value in snakemake.input.items():
    # Convert each input item to a string and ensure it's a list
    if isinstance(value, list):
        inputs_converted[key] = [str(v) for v in value]
    else:
        inputs_converted[key] = [str(value)]

report_params_dict = {
    "inputs": inputs_converted,
    "config_data": config_data,
    "analysis_table": analysis_table,
    "variables": dict(snakemake.params),
}

with open(snakemake.output.report_params, "w") as f:
    yaml.safe_dump(
        report_params_dict, f, default_flow_style=False, sort_keys=False
    )

# Define the source and destination paths
src_path = Path(snakemake.params.qmd)
dest_path = Path(snakemake.output.analysis_qmd)

# Copy the file
shutil.copy(src_path, dest_path)

