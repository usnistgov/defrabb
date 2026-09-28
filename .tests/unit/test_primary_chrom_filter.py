"""Primary-chromosome filter applied when standardizing assembly variant calls.

dipcall baseline beds and PAV VCF/beds include alt/random/Un/decoy contigs;
standardize_vcasm_output restricts both to the reference's primary chromosomes.
See docs/issues/truvari-refine-primary-chr-filter.md.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HELPERS_REF = REPO_ROOT / "rules" / "helpers_ref.smk"
ASM_VARCALL = REPO_ROOT / "rules" / "asm-varcall.smk"


def load_get_primary_chromosomes(ref_config: dict):
    """Exec get_primary_chromosomes from helpers_ref.smk with a stub ref_config."""
    src = HELPERS_REF.read_text()
    match = re.search(
        r"^def get_primary_chromosomes\(.*?(?=^def |\Z)", src, re.MULTILINE | re.DOTALL
    )
    assert match, "get_primary_chromosomes not found in helpers_ref.smk"
    namespace = {"ref_config": ref_config}
    exec(match.group(0), namespace)
    return namespace["get_primary_chromosomes"]


def test_grch38_chm13_chr_prefixed():
    fn = load_get_primary_chromosomes({"GRCh38": {}, "CHM13v2.0": {}})
    for ref_id in ["GRCh38", "CHM13v2.0"]:
        chroms = fn(ref_id)
        assert chroms[0] == "chr1" and chroms[-2:] == ["chrX", "chrY"]
        assert len(chroms) == 24
        assert "chrM" not in chroms


def test_grch37_unprefixed():
    chroms = load_get_primary_chromosomes({"GRCh37": {}})("GRCh37")
    assert chroms[0] == "1" and chroms[-2:] == ["X", "Y"]
    assert len(chroms) == 24
    assert "MT" not in chroms


def test_config_override():
    fn = load_get_primary_chromosomes(
        {"GRCh38_chr21": {"primary_chromosomes": ["chr21"]}}
    )
    assert fn("GRCh38_chr21") == ["chr21"]


def test_standardize_filters_vcf_and_bed():
    block = ASM_VARCALL.read_text().split("rule standardize_vcasm_output:")[1]
    block = re.split(r"^rule ", block, flags=re.MULTILINE)[0]
    assert "get_primary_chrom_param" in block
    assert re.search(r"bcftools view -t \{params\.chroms\}", block)
    assert "{input.bed}" in block and "keep[" in block
    assert not re.search(r"^\s*cp \{input", block, re.MULTILINE)
