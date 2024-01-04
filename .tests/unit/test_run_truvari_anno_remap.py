import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_run_truvari_anno_remap():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/run_truvari_anno_remap/data")
        expected_path = PurePosixPath(".tests/unit/run_truvari_anno_remap/expected")
        resource_path = PurePosixPath(".tests/integration/resources")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(resource_path, workdir / "resources" )

        # dbg
        output="results/draft_benchmarksets/testI/intermediates/GRCh38_chr21_asm17aChr21_stvar_dipcall-default.remap.vcf"
        print(output, file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            output,
            "-f", 
            "-j1",
            "--keep-target-files",
            "--touch",
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
