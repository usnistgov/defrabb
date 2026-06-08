"""Tests for scripts/dryrun_summary.py — pure-parse units, no snakemake invocation.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from dryrun_summary import (
    parse_job_stats,
    summarize_analyses,
    unique_ordered,
)  # noqa: E402


SAMPLE_OUTPUT = """\
Building DAG of jobs...
Job stats:
job                count
---------------  -------
all                    1
run_dipcall            3
tabix                 79
total                 83

Reasons:
    (check individual jobs above for details)
This was a dry-run (flag -n).
"""


def test_parse_job_stats_basic():
    counts, total = parse_job_stats(SAMPLE_OUTPUT)
    assert counts == {"all": 1, "run_dipcall": 3, "tabix": 79}
    assert total == 83


def test_parse_job_stats_no_table():
    counts, total = parse_job_stats("nothing here\n")
    assert counts == {}
    assert total == 0


def test_unique_ordered_preserves_first_occurrence():
    assert unique_ordered(["a", "b", "a", "c", "", "b"]) == ["a", "b", "c"]


def test_summarize_analyses_collects_unique_columns():
    rows = [
        {
            "ref": "GRCh38",
            "asm_id": "HG002",
            "vc_cmd": "dipcall",
            "bench_type": "smvar",
            "exclusion_set": "ex1",
            "eval_comp_id": "v4.2.1",
            "eval_cmd": "happy",
        },
        {
            "ref": "GRCh38",
            "asm_id": "HG002",
            "vc_cmd": "dipcall",
            "bench_type": "stvar",
            "exclusion_set": "ex2",
            "eval_comp_id": "v0.6",
            "eval_cmd": "truvari",
        },
        {
            "ref": "CHM13",
            "asm_id": "HG002",
            "vc_cmd": "dipcall",
            "bench_type": "smvar",
            "exclusion_set": "ex1",
            "eval_comp_id": "v4.2.1",
            "eval_cmd": "happy",
        },
    ]
    s = summarize_analyses(rows)
    assert s["ref"] == ["GRCh38", "CHM13"]
    assert s["asm_id"] == ["HG002"]
    assert s["bench_type"] == ["smvar", "stvar"]
    assert s["exclusion_set"] == ["ex1", "ex2"]
    assert s["eval_comp_id"] == ["v4.2.1", "v0.6"]
    assert s["eval_cmd"] == ["happy", "truvari"]
