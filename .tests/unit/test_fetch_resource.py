"""Tests for scripts/fetch_resource.sh — local-path / compression handling (#186).

Exercises the local (non-network) branches: copy, decompress-by-dest-extension,
file:// URIs, and error on missing source.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""

import gzip
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FETCH = REPO_ROOT / "scripts" / "fetch_resource.sh"


def run(src, dest):
    return subprocess.run(
        ["bash", str(FETCH), str(src), str(dest)],
        capture_output=True,
        text=True,
    )


def test_local_plain_to_plain_copies(tmp_path):
    src = tmp_path / "src.bed"
    src.write_text("chr1\t0\t100\n")
    dest = tmp_path / "out.bed"
    assert run(src, dest).returncode == 0
    assert dest.read_text() == "chr1\t0\t100\n"


def test_local_gzipped_to_plain_decompresses(tmp_path):
    src = tmp_path / "src.bed.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("chr1\t0\t100\n")
    dest = tmp_path / "out.bed"  # plain dest -> decompress
    assert run(src, dest).returncode == 0
    assert dest.read_text() == "chr1\t0\t100\n"


def test_local_gzipped_to_gz_keeps_bytes(tmp_path):
    src = tmp_path / "src.bed.gz"
    with gzip.open(src, "wt") as fh:
        fh.write("chr1\t0\t100\n")
    dest = tmp_path / "out.bed.gz"  # .gz dest -> keep compressed
    assert run(src, dest).returncode == 0
    with gzip.open(dest, "rt") as fh:
        assert fh.read() == "chr1\t0\t100\n"


def test_file_uri_is_treated_as_local(tmp_path):
    src = tmp_path / "ref.fa.gz"
    with gzip.open(src, "wt") as fh:
        fh.write(">c\nACGT\n")
    dest = tmp_path / "out.fa"
    assert run(f"file://{src}", dest).returncode == 0
    assert dest.read_text() == ">c\nACGT\n"


def test_missing_local_source_errors(tmp_path):
    result = run(tmp_path / "nope.bed", tmp_path / "out.bed")
    assert result.returncode != 0
    assert "not found" in (result.stdout + result.stderr)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
