"""Tests for scripts/diploidize_gt.py (#59 / #173).

Guards the haploid-GT diploidization that keeps vcfgeno2haplo from segfaulting
on hemizygous chrX / chrY calls. See
docs/issues/genome-specific-geno2haplo-haploid-segfault.md.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from diploidize_gt import (  # noqa: E402
    diploidize_gt_field,
    diploidize_line,
)


# --- GT field ----------------------------------------------------------------


def test_haploid_alt_becomes_homozygous():
    assert diploidize_gt_field("1") == "1|1"


def test_haploid_ref_becomes_homozygous():
    assert diploidize_gt_field("0") == "0|0"


def test_haploid_multiallelic_index_preserved():
    # dipcall emits allele indices >1 on chrY; keep the index, just diploidize.
    assert diploidize_gt_field("2") == "2|2"
    assert diploidize_gt_field("3") == "3|3"


def test_haploid_missing_becomes_missing_pair():
    assert diploidize_gt_field(".") == ".|."


def test_phased_diploid_untouched():
    assert diploidize_gt_field("0|1") == "0|1"
    assert diploidize_gt_field("1|2") == "1|2"


def test_unphased_diploid_untouched():
    assert diploidize_gt_field("0/1") == "0/1"
    assert diploidize_gt_field("1/2") == "1/2"


def test_missing_diploid_untouched():
    assert diploidize_gt_field(".|.") == ".|."
    assert diploidize_gt_field("./.") == "./."


# --- whole line --------------------------------------------------------------


def test_header_passes_through():
    line = "##fileformat=VCFv4.2\n"
    assert diploidize_line(line) == line
    col = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
    assert diploidize_line(col) == col


def test_haploid_record_only_gt_token_rewritten():
    # GT is the first FORMAT token; AD and the rest must be preserved verbatim.
    line = "chrY\t11423451\t.\tA\tG\t30\tHET1;DIPY\t.\tGT:AD\t1:1,5\n"
    expected = "chrY\t11423451\t.\tA\tG\t30\tHET1;DIPY\t.\tGT:AD\t1|1:1,5\n"
    assert diploidize_line(line) == expected


def test_diploid_record_untouched():
    line = "chr21\t5000\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0|1:3,4\n"
    assert diploidize_line(line) == line


def test_record_without_sample_column_untouched():
    line = "chr1\t100\t.\tA\tG\t.\tPASS\t.\n"
    assert diploidize_line(line) == line
