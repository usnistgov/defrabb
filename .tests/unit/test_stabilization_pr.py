import importlib.machinery
import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


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


class FakeS3Client:
    def __init__(self):
        self.uploads = []

    def upload_file(self, file_path, bucket_name, key, ExtraArgs=None):
        self.uploads.append(
            {
                "file_path": str(file_path),
                "bucket_name": bucket_name,
                "key": key,
                "extra_args": ExtraArgs,
            }
        )


class StabilizationPRTests(unittest.TestCase):
    def test_intersect_pav_callable_regions_uses_distinct_haplotype_inputs(self):
        asm_varcall_rules = (REPO_ROOT / "rules" / "asm-varcall.smk").read_text()
        intersect_rule = asm_varcall_rules.split("rule intersect_pav_callable_regions:")[
            1
        ]
        intersect_rule = intersect_rule.split("rule standardize_vcasm_output:")[0]

        self.assertIn(
            'h1_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h1_500.bed.gz"',
            intersect_rule,
        )
        self.assertIn(
            'h2_bed="results/asm_varcalls/{vc_id}/results/{sample_id}/callable/callable_regions_h2_500.bed.gz"',
            intersect_rule,
        )

    def test_render_report_uses_root_analysis_qmd(self):
        report_rules = (REPO_ROOT / "rules" / "report.smk").read_text()

        self.assertIn('qmd=Path(workflow.basedir) / "analysis.qmd"', report_rules)
        self.assertTrue((REPO_ROOT / "analysis.qmd").exists())
        self.assertFalse((REPO_ROOT / "scripts" / "reports" / "analysis.qmd").exists())

    def test_release_config_includes_data_manifest(self):
        release_config = json.loads((REPO_ROOT / "config" / "release.json").read_text())

        self.assertIn(
            "data_manifest.tsv", release_config["release_rules"]["local"]["include"]
        )

    def test_upload_to_s3_keeps_included_files_with_default_exclude(self):
        module = load_run_defrabb_module()
        fake_s3_client = FakeS3Client()

        with tempfile.TemporaryDirectory() as tmpdir:
            run_id = "20260317_v0.001_test"
            run_dir = Path(tmpdir) / run_id
            (run_dir / "results").mkdir(parents=True)
            (run_dir / "config").mkdir(parents=True)
            (run_dir / "results" / "kept.txt").write_text("keep me")
            (run_dir / "config" / f"analyses_{run_id}.tsv").write_text("analysis_id\n")
            (run_dir / "data_manifest.tsv").write_text("analysis\tfile_type\n")
            (run_dir / "skip.txt").write_text("skip me")

            with mock.patch.object(
                module.boto3, "client", return_value=fake_s3_client
            ):
                module.upload_to_s3(
                    run_id,
                    str(run_dir),
                    "example-bucket",
                    "defrabb_runs",
                    {
                        "include": [
                            "results/**",
                            "config/analyses_{RUNID}.tsv",
                            "data_manifest.tsv",
                        ],
                        "exclude": ["*"],
                    },
                )

        self.assertEqual(
            {upload["key"] for upload in fake_s3_client.uploads},
            {
                f"defrabb_runs/{run_id}/config/analyses_{run_id}.tsv",
                f"defrabb_runs/{run_id}/data_manifest.tsv",
                f"defrabb_runs/{run_id}/results/kept.txt",
            },
        )

    def test_release_run_uploads_manifest_after_generating_it(self):
        module = load_run_defrabb_module()
        fake_s3_client = FakeS3Client()

        with tempfile.TemporaryDirectory() as tmpdir:
            run_id = "20260317_v0.001_test"
            run_dir = Path(tmpdir) / run_id
            run_dir.mkdir()
            events = []

            def fake_upload_to_s3(*args, **kwargs):
                events.append("upload")

            def fake_create_data_manifest(
                manifest_run_id, bucket_name, s3_path, manifest_path
            ):
                events.append("manifest")
                Path(manifest_path).write_text("analysis\tfile_type\n")

            with mock.patch.object(
                module.boto3, "client", return_value=fake_s3_client
            ), mock.patch.object(
                module, "upload_to_s3", side_effect=fake_upload_to_s3
            ), mock.patch.object(
                module,
                "create_data_manifest",
                side_effect=fake_create_data_manifest,
            ):
                module.release_run(
                    SimpleNamespace(runid=run_id),
                    str(run_dir),
                    str(Path(tmpdir) / "archive"),
                    "example-bucket",
                    "defrabb_runs",
                    "s3",
                    {
                        "local": {"include": [], "exclude": ["*"]},
                        "s3": {"include": ["data_manifest.tsv"], "exclude": ["*"]},
                    },
                )

        self.assertEqual(events, ["upload", "manifest"])
        self.assertEqual(
            fake_s3_client.uploads,
            [
                {
                    "file_path": str(run_dir / "data_manifest.tsv"),
                    "bucket_name": "example-bucket",
                    "key": f"defrabb_runs/{run_id}/data_manifest.tsv",
                    "extra_args": {"ACL": "public-read"},
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
