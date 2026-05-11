"""Tests for the ``run_defrabb`` wrapper script.

Scope (TODO #18, run-id slice): exercise the run-id format contract enforced
by ``validate_and_set_defaults``. S3 release-rule and manifest coverage live
in separate test modules.
"""

import importlib.machinery
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace


REPO_ROOT = Path(__file__).resolve().parents[2]


def load_run_defrabb_module():
    if "boto3" not in sys.modules:
        boto3_stub = types.ModuleType("boto3")
        boto3_stub.client = None
        sys.modules["boto3"] = boto3_stub

    loader = importlib.machinery.SourceFileLoader(
        "run_defrabb_module", str(REPO_ROOT / "run_defrabb")
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def make_args(runid, analyses, outdir):
    return SimpleNamespace(
        runid=runid,
        analyses=str(analyses) if analyses is not None else None,
        outdir=str(outdir),
        jobs=1,
        steps="all",
    )


class RunIDFormatTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)
        self.analyses = self.tmpdir / "analyses.tsv"
        self.analyses.write_text("placeholder\n")

    def tearDown(self):
        self._tmp.cleanup()

    def test_valid_runid_passes(self):
        args = make_args("20260507_v0.022_smoke", self.analyses, self.tmpdir)
        result = self.module.validate_and_set_defaults(args)
        self.assertEqual(result.runid, "20260507_v0.022_smoke")

    def test_valid_runid_with_complex_suffix(self):
        args = make_args(
            "20260101_v1.999_HG002-v1.1-stratify", self.analyses, self.tmpdir
        )
        result = self.module.validate_and_set_defaults(args)
        self.assertEqual(result.runid, "20260101_v1.999_HG002-v1.1-stratify")

    def test_missing_runid_raises(self):
        args = make_args(None, self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError) as ctx:
            self.module.validate_and_set_defaults(args)
        self.assertIn("--runid is required", str(ctx.exception))

    def test_empty_runid_raises(self):
        args = make_args("", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_short_date_prefix_rejected(self):
        args = make_args("2026507_v0.022_smoke", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError) as ctx:
            self.module.validate_and_set_defaults(args)
        self.assertIn("Invalid RUN ID format", str(ctx.exception))

    def test_missing_version_segment_rejected(self):
        args = make_args("20260507_smoke", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_two_digit_major_version_rejected(self):
        # Pattern pins the major version to a single digit.
        args = make_args("20260507_v10.022_smoke", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_wrong_minor_width_rejected(self):
        # Minor version must be exactly 3 digits.
        args = make_args("20260507_v0.22_smoke", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_missing_trailing_underscore_rejected(self):
        args = make_args("20260507_v0.022", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_alpha_date_prefix_rejected(self):
        args = make_args("YYYYMMDD_v0.022_smoke", self.analyses, self.tmpdir)
        with self.assertRaises(self.module.InvalidRunIDError):
            self.module.validate_and_set_defaults(args)

    def test_default_analyses_path_uses_runid(self):
        # When --analyses is not supplied, the wrapper derives it from runid.
        # Create the expected default file inside a fake repo layout.
        fake_repo = self.tmpdir / "repo"
        (fake_repo / "config").mkdir(parents=True)
        runid = "20260507_v0.022_smoke"
        expected = fake_repo / "config" / f"analyses_{runid}.tsv"
        expected.write_text("placeholder\n")

        cwd = Path.cwd()
        try:
            import os

            os.chdir(fake_repo)
            args = make_args(runid, None, self.tmpdir)
            result = self.module.validate_and_set_defaults(args)
            self.assertEqual(result.analyses, f"config/analyses_{runid}.tsv")
        finally:
            import os

            os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
