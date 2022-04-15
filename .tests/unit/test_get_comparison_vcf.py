import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_get_comparison_vcf():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/get_comparison_vcf/data")
        expected_path = PurePosixPath(".tests/unit/get_comparison_vcf/expected")
        config_path = PurePosixPath("config")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")

        # dbg
        print("resources/comparison_variant_callsets/GRCh38_chr21_v4.2.1chr21.vcf.gz", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "resources/comparison_variant_callsets/GRCh38_chr21_v4.2.1chr21.vcf.gz",
            "-f", 
            "-j1",
            "--keep-target-files",
    
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
