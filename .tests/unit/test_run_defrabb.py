"""Tests for the ``run_defrabb`` wrapper script.

Scope (TODO #18): exercise run-id format, S3 release rules, and manifest generation.
"""

import importlib.machinery
import importlib.util
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch


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


class ReleasePatternTests(unittest.TestCase):
    """Tests for S3 upload include/exclude pattern matching."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def test_matches_release_patterns_exact_match(self):
        patterns = ["run.log", "environment.yml"]
        self.assertTrue(
            self.module.matches_release_patterns(
                "run.log", "20260507_v0.022_smoke", patterns
            )
        )

    def test_matches_release_patterns_wildcard(self):
        patterns = ["results/**", "config/**"]
        self.assertTrue(
            self.module.matches_release_patterns(
                "results/draft_benchmarksets/test.vcf.gz",
                "20260507_v0.022_smoke",
                patterns,
            )
        )

    def test_matches_release_patterns_runid_substitution(self):
        patterns = ["config/analyses_{RUNID}.tsv"]
        runid = "20260507_v0.022_smoke"
        self.assertTrue(
            self.module.matches_release_patterns(
                f"config/analyses_{runid}.tsv", runid, patterns
            )
        )

    def test_matches_release_patterns_no_match(self):
        patterns = ["run.log", "environment.yml"]
        self.assertFalse(
            self.module.matches_release_patterns(
                "random_file.txt", "20260507_v0.022_smoke", patterns
            )
        )

    def test_should_release_file_included(self):
        release_rules = {
            "include": ["results/**", "config/**", "run.log"],
            "exclude": ["*"],
        }
        runid = "20260507_v0.022_smoke"
        self.assertTrue(
            self.module.should_release_file("results/test.vcf.gz", runid, release_rules)
        )
        self.assertTrue(
            self.module.should_release_file("config/analyses.tsv", runid, release_rules)
        )
        self.assertTrue(
            self.module.should_release_file("run.log", runid, release_rules)
        )

    def test_should_release_file_excluded(self):
        release_rules = {
            "include": ["results/**", "config/**"],
            "exclude": ["*"],
        }
        runid = "20260507_v0.022_smoke"
        # File not matching any include pattern should be excluded
        self.assertFalse(
            self.module.should_release_file("random_file.txt", runid, release_rules)
        )

    def test_should_release_file_with_runid_substitution(self):
        release_rules = {
            "include": ["config/analyses_{RUNID}.tsv"],
            "exclude": ["*"],
        }
        runid = "20260507_v0.022_smoke"
        self.assertTrue(
            self.module.should_release_file(
                f"config/analyses_{runid}.tsv", runid, release_rules
            )
        )
        # Different runid should not match
        self.assertFalse(
            self.module.should_release_file(
                "config/analyses_20260101_v0.001_test.tsv", runid, release_rules
            )
        )


class ReleaseRulesExpansionTests(unittest.TestCase):
    """Tests for release rules validation and expansion."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def test_same_as_local_expansion(self):
        release_config = {
            "release_rules": {
                "local": {"include": ["*.log"], "exclude": ["*"]},
                "s3": "same as local",
            }
        }
        expanded = self.module.validate_and_expand_release_rules(release_config, "both")
        self.assertEqual(expanded["s3"]["include"], ["*.log"])
        self.assertEqual(expanded["s3"]["exclude"], ["*"])

    def test_same_as_s3_expansion(self):
        release_config = {
            "release_rules": {
                "s3": {"include": ["*.log"], "exclude": ["*"]},
                "local": "same as s3",
            }
        }
        expanded = self.module.validate_and_expand_release_rules(release_config, "both")
        self.assertEqual(expanded["local"]["include"], ["*.log"])
        self.assertEqual(expanded["local"]["exclude"], ["*"])

    def test_missing_include_raises(self):
        release_config = {
            "release_rules": {
                "local": {"exclude": ["*"]},
            }
        }
        with self.assertRaises(ValueError) as ctx:
            self.module.validate_and_expand_release_rules(release_config, "local")
        self.assertIn(
            "Include and exclude patterns must be defined", str(ctx.exception)
        )

    def test_missing_exclude_raises(self):
        release_config = {
            "release_rules": {
                "s3": {"include": ["*.log"]},
            }
        }
        with self.assertRaises(ValueError) as ctx:
            self.module.validate_and_expand_release_rules(release_config, "s3")
        self.assertIn(
            "Include and exclude patterns must be defined", str(ctx.exception)
        )


class ProfileLoadingTests(unittest.TestCase):
    """Tests for profile loading and application (TODO #14)."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def test_load_profile_config_success(self):
        """Test loading a valid profile config."""
        with tempfile.TemporaryDirectory() as tmpdir:
            profile_dir = Path(tmpdir) / "profiles" / "test"
            profile_dir.mkdir(parents=True)
            config_path = profile_dir / "config.json"
            config_path.write_text('{"defaults": {"outdir": "/test/path"}}')

            # Temporarily change to tmpdir so relative path works
            import os

            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                config = self.module.load_profile_config("test")
                self.assertEqual(config["defaults"]["outdir"], "/test/path")
            finally:
                os.chdir(original_cwd)

    def test_load_profile_config_missing(self):
        """Test loading a non-existent profile returns None."""
        with tempfile.TemporaryDirectory() as tmpdir:
            import os

            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                config = self.module.load_profile_config("nonexistent")
                self.assertIsNone(config)
            finally:
                os.chdir(original_cwd)

    def test_apply_profile_defaults_only_when_not_set(self):
        """Test profile defaults are only applied when args not explicitly set."""
        import argparse

        # Create mock args with parser defaults (not set by user) and explicit values
        args = argparse.Namespace(
            outdir=".",  # Parser default (updated from "./defrabb_runs/"), should use profile default
            archive_dir="/explicit/path",  # Explicitly set, should NOT be overridden
            release_config="config/release.json",  # Parser default, should use profile default
            release_type="s3",  # Parser default (matches profile, but demonstrates logic)
        )

        profile_config = {
            "defaults": {
                "outdir": "/profile/outdir",
                "archive_dir": "/profile/archive",
                "release_config": "profiles/test/release.json",
                "release_type": "local",
            }
        }

        self.module.apply_profile_defaults(args, profile_config)

        # Check that parser defaults were replaced with profile defaults
        self.assertEqual(args.outdir, "/profile/outdir")
        self.assertEqual(args.release_config, "profiles/test/release.json")
        self.assertEqual(args.release_type, "local")

        # Check that explicit value was NOT overridden
        self.assertEqual(args.archive_dir, "/explicit/path")


class SubcommandParsingTests(unittest.TestCase):
    """Tests for new subcommand-based CLI (TODO #6)."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def test_subcommand_run_maps_to_pipe_step(self):
        """'run' subcommand should map to 'pipe' step internally."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "run", "-r", "20260515_v0.022_test"]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.command, "run")
            self.assertEqual(args.steps, "pipe")
            self.assertFalse(args.validate_only)
        finally:
            sys.argv = original_argv

    def test_subcommand_validate_sets_validate_only(self):
        """'validate' subcommand should set validate_only=True."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "validate", "-a", "config/analyses.tsv"]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.command, "validate")
            self.assertTrue(args.validate_only)
            self.assertIsNone(args.steps)
        finally:
            sys.argv = original_argv

    def test_legacy_steps_syntax_still_works(self):
        """Old -s/--steps syntax should still work for backward compatibility."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = [
                "run_defrabb",
                "-s",
                "pipe",
                "-r",
                "20260515_v0.022_test",
            ]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.steps, "pipe")
            self.assertFalse(args.validate_only)
        finally:
            sys.argv = original_argv

    def test_legacy_validate_only_still_works(self):
        """Old --validate-only flag should still work."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "--validate-only", "-a", "config/analyses.tsv"]
            args, _ = self.module.parse_arguments()
            self.assertTrue(args.validate_only)
        finally:
            sys.argv = original_argv

    def test_subcommand_report_maps_correctly(self):
        """'report' subcommand should map to 'report' step."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "report", "-r", "20260515_v0.022_test"]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.command, "report")
            self.assertEqual(args.steps, "report")
        finally:
            sys.argv = original_argv

    def test_subcommand_archive_maps_correctly(self):
        """'archive' subcommand should map to 'archive' step."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "archive", "-r", "20260515_v0.022_test"]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.command, "archive")
            self.assertEqual(args.steps, "archive")
        finally:
            sys.argv = original_argv

    def test_subcommand_release_maps_correctly(self):
        """'release' subcommand should map to 'release' step."""
        import sys

        original_argv = sys.argv.copy()
        try:
            sys.argv = ["run_defrabb", "release", "-r", "20260515_v0.022_test"]
            args, _ = self.module.parse_arguments()
            self.assertEqual(args.command, "release")
            self.assertEqual(args.steps, "release")
        finally:
            sys.argv = original_argv


class ManifestGenerationTests(unittest.TestCase):
    """Tests for data manifest generation."""

    @classmethod
    def setUpClass(cls):
        cls.module = load_run_defrabb_module()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    @patch("boto3.client")
    def test_create_data_manifest_basic_structure(self, mock_boto3_client):
        """Test that manifest has correct headers and structure."""
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3
        mock_s3.list_objects_v2.return_value = {
            "Contents": [
                {"Key": "defrabb_runs/20260507_v0.022_test/run.log", "Size": 1024}
            ],
            "IsTruncated": False,
        }

        manifest_path = self.tmpdir / "manifest.tsv"
        self.module.create_data_manifest(
            "20260507_v0.022_test", "test-bucket", "defrabb_runs", str(manifest_path)
        )

        # Verify manifest exists and has correct header
        self.assertTrue(manifest_path.exists())
        content = manifest_path.read_text()
        lines = content.strip().split("\n")
        self.assertEqual(lines[0], "analysis\tfile_type\tref_id\tsize_Gb\ts3_uri\turl")

    @patch("boto3.client")
    def test_create_data_manifest_file_type_detection(self, mock_boto3_client):
        """Test that different file types are correctly identified."""
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        test_files = [
            ("defrabb_runs/20260507_v0.022_test/run.log", "run.log"),
            ("defrabb_runs/20260507_v0.022_test/environment.yml", "mamba env"),
            ("defrabb_runs/20260507_v0.022_test/archive.tar.gz", "snakemake_archive"),
            (
                "defrabb_runs/20260507_v0.022_test/config/analyses_20260507_v0.022_test.tsv",
                "analysis table",
            ),
            (
                "defrabb_runs/20260507_v0.022_test/results/draft_benchmarksets/bench1/GRCh38_HG002_smvar.vcf.gz",
                "smvar benchmark vcf",
            ),
            (
                "defrabb_runs/20260507_v0.022_test/results/draft_benchmarksets/bench1/GRCh38_HG002_stvar.vcf.gz",
                "stvar benchmark vcf",
            ),
        ]

        mock_s3.list_objects_v2.return_value = {
            "Contents": [{"Key": key, "Size": 1024} for key, _ in test_files],
            "IsTruncated": False,
        }

        manifest_path = self.tmpdir / "manifest.tsv"
        self.module.create_data_manifest(
            "20260507_v0.022_test", "test-bucket", "defrabb_runs", str(manifest_path)
        )

        content = manifest_path.read_text()
        lines = content.strip().split("\n")[1:]  # Skip header

        for i, (_, expected_type) in enumerate(test_files):
            cols = lines[i].split("\t")
            actual_type = cols[1]
            self.assertEqual(
                actual_type, expected_type, f"File type mismatch for {test_files[i][0]}"
            )

    @patch("boto3.client")
    def test_create_data_manifest_ref_id_detection(self, mock_boto3_client):
        """Test that reference genome IDs are correctly extracted."""
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        test_files = [
            ("defrabb_runs/test/results/GRCh38_HG002_smvar.vcf.gz", "GRCh38"),
            ("defrabb_runs/test/results/GRCh37_HG002_smvar.vcf.gz", "GRCh37"),
            ("defrabb_runs/test/results/CHM13_HG002_smvar.vcf.gz", "CHM13"),
            ("defrabb_runs/test/run.log", "NA"),
        ]

        mock_s3.list_objects_v2.return_value = {
            "Contents": [{"Key": key, "Size": 1024} for key, _ in test_files],
            "IsTruncated": False,
        }

        manifest_path = self.tmpdir / "manifest.tsv"
        self.module.create_data_manifest(
            "test", "test-bucket", "defrabb_runs", str(manifest_path)
        )

        content = manifest_path.read_text()
        lines = content.strip().split("\n")[1:]  # Skip header

        for i, (_, expected_ref) in enumerate(test_files):
            cols = lines[i].split("\t")
            actual_ref = cols[2]
            self.assertEqual(
                actual_ref, expected_ref, f"Ref ID mismatch for {test_files[i][0]}"
            )

    @patch("boto3.client")
    def test_create_data_manifest_pagination(self, mock_boto3_client):
        """Test that manifest generation handles S3 pagination."""
        mock_s3 = MagicMock()
        mock_boto3_client.return_value = mock_s3

        # Simulate paginated response
        mock_s3.list_objects_v2.side_effect = [
            {
                "Contents": [{"Key": "defrabb_runs/test/file1.log", "Size": 1024}],
                "IsTruncated": True,
                "NextContinuationToken": "token1",
            },
            {
                "Contents": [{"Key": "defrabb_runs/test/file2.log", "Size": 2048}],
                "IsTruncated": False,
            },
        ]

        manifest_path = self.tmpdir / "manifest.tsv"
        self.module.create_data_manifest(
            "test", "test-bucket", "defrabb_runs", str(manifest_path)
        )

        content = manifest_path.read_text()
        lines = content.strip().split("\n")[1:]  # Skip header

        # Should have 2 entries (one from each page)
        self.assertEqual(len(lines), 2)

        # Verify both pagination calls were made
        self.assertEqual(mock_s3.list_objects_v2.call_count, 2)


if __name__ == "__main__":
    unittest.main()
