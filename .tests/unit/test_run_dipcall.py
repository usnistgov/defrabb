import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_run_dipcall():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/run_dipcall/data")
        expected_path = PurePosixPath(".tests/unit/run_dipcall/expected")
        config_path = PurePosixPath("config")
        references_path = PurePosixPath(".tests/integration/resources/references")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)
        shutil.copytree(config_path, workdir / "config")
        shutil.copytree(references_path, workdir / "resources" / "references")

        # dbg
        print(
            "results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.mak results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.dip.vcf.gz results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.dip.bed results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.hap1.bam results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.hap2.bam",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/asm_varcalls/vc1/GRCh38_chr21_asm17aChr21_dipcall-default.mak",
                "-f",
                "-j1",
                "--keep-target-files",
                "--verbose",
                "--use-conda",
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
