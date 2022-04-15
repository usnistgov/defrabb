"""
Common code for unit testing of rules generated with Snakemake 7.3.0.
"""

from pathlib import Path
import subprocess as sp
import os
import hashlib
import gzip

class OutputChecker:
    def __init__(self, data_path, expected_path, workdir, ignore_unexpected = False):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir
        self.ignore_unexpected = ignore_unexpected

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if (
                    str(f).startswith(".snakemake")
                    or str(f).startswith("config")
                    or str(f).startswith("resources")
                    or str(f).startswith("logs")
                    or str(f).startswith("benchmark")
                    or str(f).endswith("log")
                    or str(f).endswith("json")
                    or str(f).endswith("json.gz")
                ):
                    continue
                elif f in expected_files:
                    if str(f).endswith("mak"):
                        print(f"Not comparing makefile: {f}")
                    elif str(f).endswith("vcf.gz"):
                        compare_vcfs(self.workdir / f, self.expected_path / f)
                    elif str(f).endswith("bam"):
                        compare_bams(self.workdir / f, self.expected_path / f)
                    else:
                        self.compare_files(self.workdir / f, self.expected_path / f)
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
        if unexpected_files:
            if self.ignore_unexpected: 
                "Ignoring unexpected files:\n{}".format(
                    "\n".join(sorted(map(str, unexpected_files)))
                )
            else:
                raise ValueError(
                    "Unexpected files:\n{}".format(
                        "\n".join(sorted(map(str, unexpected_files)))
                    )
                )

    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])


## Comparison rules from https://github.com/nikostr/dna-seq-deepvariant-glnexus-variant-calling
def calc_md5sum(file, comment_char):
    with gzip.open(file, mode = "rt", encoding="latin_1") as f:
        file_string="".join([l for l in f if not l.startswith(comment_char)]).encode('utf-8')
        return hashlib.md5(
           
        ).hexdigest()


def compare_vcfs(generated_file, expected_file):
    assert calc_md5sum(generated_file, comment_char = "##") == calc_md5sum(
        expected_file, comment_char = "##"
    ), "md5sum of vcfs do not match"

def compare_bams(generated_file, expected_file):
    assert calc_md5sum(generated_file, comment_char = "@") == calc_md5sum(
        expected_file, comment_char = "@"
    ), "md5sum of vcfs do not match"
