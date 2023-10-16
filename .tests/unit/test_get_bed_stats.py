import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_get_bed_stats():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/get_bed_stats/data")
        expected_path = PurePosixPath(".tests/unit/get_bed_stats/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "results/draft_benchmarksets/testB/GRCh38_chr21_asm17aChr21_dipcall-default.benchmark_bed-summary.tsv",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/draft_benchmarksets/testB/GRCh38_chr21_asm17aChr21_dipcall-default.benchmark_bed-summary.tsv",
                "-f",
                "-j1",
                "--keep-target-files",
                "--touch",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
