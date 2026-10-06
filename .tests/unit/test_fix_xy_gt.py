"""fix_XY_genotype must not silently emit an empty VCF from a stale input index.

2026-10-05 v1.2 fulltest: normalize_vars regenerated `.norm.vcf.gz` but the
untracked `.tbi` from the previous run was reused; fix_xy_gt.sh region queries
hit BGZF errors, the script exited 0, and every downstream PAV annotation step
ran on an empty VCF.
"""

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "fix_xy_gt.sh"
RULES = REPO_ROOT / "rules" / "bench_vcf_normalize.smk"


def test_script_fails_fast():
    assert "set -euo pipefail" in SCRIPT.read_text()


def test_script_reindexes_stale_index():
    text = SCRIPT.read_text()
    assert '-nt "${input_vcf}"' in text
    assert "bcftools index -f --tbi" in text


def test_rule_declares_input_index():
    block = RULES.read_text().split("rule fix_XY_genotype:", 1)[1]
    block = block.split("\nrule ", 1)[0]
    assert re.search(r'tbi="[^"]*\{prefix\}\.vcf\.gz\.tbi"', block)


TOOLS = ("bcftools", "bedtools", "bgzip")


@pytest.mark.skipif(
    any(shutil.which(t) is None for t in TOOLS), reason="bcftools/bedtools needed"
)
def test_stale_index_end_to_end(tmp_path: Path):
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=1000>\n##contig=<ID=X,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )

    def write_vcf(n: int) -> None:
        recs = "".join(f"1\t{i}\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n" for i in range(1, n))
        recs += "X\t500\t.\tA\tC\t.\tPASS\t.\tGT\t.|1\n"
        plain = tmp_path / "in.vcf"
        plain.write_text(header + recs)
        subprocess.run(["bgzip", "-f", str(plain)], check=True)

    vcf = tmp_path / "in.vcf.gz"
    write_vcf(400)
    subprocess.run(["bcftools", "index", "--tbi", str(vcf)], check=True)
    time.sleep(1)
    write_vcf(50)  # rewrite data; index now stale
    os.utime(vcf)

    (tmp_path / "genome").write_text("1\t1000\nX\t1000\n")
    (tmp_path / "par.bed").write_text("X\t0\t10\n")
    out = tmp_path / "out.vcf.gz"
    subprocess.run(
        ["bash", str(SCRIPT), "-i", str(vcf), "-o", str(out)]
        + ["-p", str(tmp_path / "par.bed"), "-g", str(tmp_path / "genome")],
        check=True,
        cwd=tmp_path,
        capture_output=True,
    )
    view = subprocess.run(
        ["bcftools", "view", "-H", str(out)], check=True, capture_output=True, text=True
    )
    lines = view.stdout.splitlines()
    assert len(lines) == 50
    assert lines[-1].split("\t")[-1] == "1"
