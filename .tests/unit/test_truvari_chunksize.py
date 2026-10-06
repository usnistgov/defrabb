"""Hardcoded `truvari bench` calls must set --chunksize >= -r/--refdist.

Truvari v5 aborts with "--chunksize must be >= --refdist" (default chunksize
1000); pav_discrep_truvari failed this way in the 2026-10-06 v1.2 fulltest.
"""

import re
from pathlib import Path

RULES_DIR = Path(__file__).resolve().parents[2] / "rules"


def _bench_calls():
    for smk in RULES_DIR.glob("*.smk"):
        for m in re.finditer(r"truvari bench(.*?)&>", smk.read_text(), re.S):
            yield smk.name, m.group(1)


def test_refdist_within_chunksize():
    checked = 0
    for name, call in _bench_calls():
        refdist = re.search(r"(?:-r|--refdist) (\d+)", call)
        if not refdist:
            continue
        chunk = re.search(r"(?:-C|--chunksize) (\d+)", call)
        assert chunk, f"{name}: truvari bench sets -r without --chunksize"
        assert int(chunk.group(1)) >= int(refdist.group(1)), name
        checked += 1
    assert checked >= 2
