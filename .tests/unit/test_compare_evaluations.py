"""Tests for scripts/compare_evaluations.py — pure-parse units, no I/O beyond tmp.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from compare_evaluations import (  # noqa: E402
    add_deltas,
    discover_eval_files,
    load_happy_summary,
    load_truvari_summary,
    parse_eval_path,
    resolve_baseline,
    sort_records,
)

HAPPY_CSV = (
    "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,"
    "FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score\n"
    "INDEL,ALL,500,490,10,520,8,5,2,1,0.98,0.984,0.01,0.982\n"
    "INDEL,PASS,500,489,11,515,6,4,1,1,0.978,0.988,0.01,0.983\n"
    "SNP,ALL,1000,995,5,1010,4,3,1,0,0.995,0.996,0.01,0.9955\n"
    "SNP,PASS,1000,994,6,1008,3,2,1,0,0.994,0.997,0.01,0.9955\n"
)

TRUVARI_JSON = (
    '{"TP-base":900,"TP-comp":905,"FP":20,"FN":100,"precision":0.978,'
    '"recall":0.9,"f1":0.937,"gt_concordance":0.95,"weighted":{"x":1},'
    '"gt_matrix":{"a":1}}'
)


def _write_happy(tmp_path, group, callset):
    d = tmp_path / "evaluations" / "happy" / group
    d.mkdir(parents=True)
    p = d / f"{callset}.summary.csv"
    p.write_text(HAPPY_CSV)
    return p


def _write_truvari(tmp_path, group, callset):
    d = tmp_path / "evaluations" / "truvari" / group / callset
    d.mkdir(parents=True)
    p = d / "summary.json"
    p.write_text(TRUVARI_JSON)
    return p


def test_parse_eval_path_happy():
    p = Path(
        "results/evaluations/happy/evalA_benchA/"
        "GRCh38_HG002_asm1_smvar_dipcall-default.summary.csv"
    )
    meta = parse_eval_path(p)
    assert meta["eval_tool"] == "happy"
    assert meta["bench_type"] == "smvar"
    assert meta["analysis_id"] == "evalA_benchA/GRCh38_HG002_asm1_smvar_dipcall-default"


def test_parse_eval_path_truvari():
    p = Path(
        "results/evaluations/truvari/evalS_benchS/"
        "GRCh38_compS_asm1_stvar_dipcall-default/summary.json"
    )
    meta = parse_eval_path(p)
    assert meta["eval_tool"] == "truvari"
    assert meta["bench_type"] == "stvar"
    assert meta["analysis_id"].endswith("stvar_dipcall-default")


def test_load_happy_summary_keeps_pass_snp_indel(tmp_path):
    p = _write_happy(
        tmp_path, "evalA_benchA", "GRCh38_HG002_asm1_smvar_dipcall-default"
    )
    records = load_happy_summary(p)
    # Only the two PASS rows (SNP, INDEL); ALL rows dropped.
    assert {r["variant_type"] for r in records} == {"SNP", "INDEL"}
    assert len(records) == 2
    indel = next(r for r in records if r["variant_type"] == "INDEL")
    assert indel["precision"] == 0.988
    assert indel["recall"] == 0.978
    assert indel["f1"] == 0.983
    assert indel["tp"] == 489.0  # TRUTH.TP
    assert indel["fp"] == 6.0  # QUERY.FP
    assert indel["fn"] == 11.0  # TRUTH.FN
    assert indel["gt_concordance"] is None


def test_load_truvari_summary(tmp_path):
    p = _write_truvari(
        tmp_path, "evalS_benchS", "GRCh38_compS_asm1_stvar_dipcall-default"
    )
    records = load_truvari_summary(p)
    assert len(records) == 1
    r = records[0]
    assert r["variant_type"] == "SV"
    assert r["precision"] == 0.978
    assert r["f1"] == 0.937
    assert r["tp"] == 905.0  # TP-comp preferred
    assert r["fp"] == 20.0
    assert r["fn"] == 100.0
    assert r["gt_concordance"] == 0.95


def test_discover_eval_files(tmp_path):
    _write_happy(tmp_path, "evalA_benchA", "GRCh38_HG002_asm1_smvar_dipcall-default")
    _write_truvari(tmp_path, "evalS_benchS", "GRCh38_compS_asm1_stvar_dipcall-default")
    files = discover_eval_files(str(tmp_path))
    names = sorted(p.name for p in files)
    assert "summary.json" in names
    assert any(n.endswith(".summary.csv") for n in names)


def _records():
    return [
        {
            "analysis_id": "base",
            "variant_type": "SNP",
            "precision": 0.99,
            "recall": 0.99,
            "f1": 0.99,
        },
        {
            "analysis_id": "cand",
            "variant_type": "SNP",
            "precision": 0.995,
            "recall": 0.992,
            "f1": 0.993,
        },
        {
            "analysis_id": "worse",
            "variant_type": "SNP",
            "precision": 0.98,
            "recall": 0.98,
            "f1": 0.98,
        },
    ]


def test_add_deltas_flags_improved_regressed_baseline():
    recs = add_deltas(_records(), "base", "f1", threshold=0.001)
    status = {r["analysis_id"]: r["status"] for r in recs}
    assert status["base"] == "baseline"
    assert status["cand"] == "improved"
    assert status["worse"] == "regressed"
    cand = next(r for r in recs if r["analysis_id"] == "cand")
    assert cand["delta_f1"] == pytest.approx(0.003)


def test_add_deltas_threshold_marks_same():
    # cand delta +0.003 is below the 0.005 threshold -> "same";
    # worse delta -0.01 exceeds it -> "regressed".
    recs = add_deltas(_records(), "base", "f1", threshold=0.005)
    status = {r["analysis_id"]: r["status"] for r in recs}
    assert status["cand"] == "same"
    assert status["worse"] == "regressed"


def test_resolve_baseline_substring_and_errors():
    recs = _records()
    assert resolve_baseline(recs, "base") == "base"
    assert resolve_baseline(recs, "wor") == "worse"
    with pytest.raises(SystemExit):
        resolve_baseline(recs, "nomatch")


def test_resolve_baseline_ambiguous():
    recs = [
        {"analysis_id": "run_a", "variant_type": "SNP"},
        {"analysis_id": "run_b", "variant_type": "SNP"},
    ]
    with pytest.raises(SystemExit):
        resolve_baseline(recs, "run_")


def test_sort_records_by_type_then_metric_desc():
    recs = sort_records(_records(), "f1")
    # all SNP; highest f1 first
    assert [r["analysis_id"] for r in recs] == ["cand", "base", "worse"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
