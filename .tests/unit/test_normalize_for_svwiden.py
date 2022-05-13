import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_run_normalize_for_svwiden():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/normalize_for_svwiden/data")
        expected_path = PurePosixPath(".tests/unit/normalize_for_svwiden/expected")
        config_path = PurePosixPath("config")
        references_path = PurePosixPath(".tests/integration/resources/references")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")
        shutil.copytree(references_path, workdir / "resources" / "references")

        # dbg
        output = "results/draft_benchmarksets/testC/intermediates/GRCh38_chr21_asm17aChr21_dipcall-default.gt19_norm.vcf.gz"
        print(
            output,
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                output,
                "-f",
                "-j1",
                "--keep-target-files",
                "--use-conda",
<<<<<<< HEAD
                "--touch",
=======
>>>>>>> c5d7be7... adding unit tests
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(
            data_path, expected_path, workdir, ignore_unexpected=True
        ).check()
