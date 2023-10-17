import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_run_assembly_stats():
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/run_assembly_stats/data")
        expected_path = PurePosixPath(".tests/unit/run_assembly_stats/expected")
        config_path = PurePosixPath("config")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")

        # dbg
        print(
            "results/report/assemblies/asm17aChr21_maternal_stats.txt", file=sys.stderr
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/report/assemblies/asm17aChr21_maternal_stats.txt",
                "-f",
                "-j1",
                "--keep-target-files",
                "--use-conda",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
