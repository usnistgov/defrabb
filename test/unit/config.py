import unittest as ut
import pandas as pd
from snakemake.io import load_configfile
from snakemake.utils import validate
from rules.common import get_happy_inputs_inner, load_analyses


RESOURCES_PATH = "test/unit/config/resources.yml"
ANALYSES_PATH = "test/unit/config/analyses.tsv"

R_SCHEMA_PATH = "schema/resources-schema.yml"
A_SCHEMA_PATH = "schema/analyses-schema.yml"


def import_yaml(path, schema):
    config = load_configfile(path)
    validate(config, schema)
    return config


def assert_happy(self, truth, eval_id, ref_id):
    # "I assert you shall be happy with truth (or else...purple burglar alarm)"
    keys = [
        "genome",
        "genome_index",
        "strat_tb",
        "comp_idx",
        "query",
        "truth",
        "truth_regions",
        "target_regions",
    ]
    config = import_yaml(RESOURCES_PATH, R_SCHEMA_PATH)
    analyses = load_analyses(ANALYSES_PATH, A_SCHEMA_PATH)
    actual = get_happy_inputs_inner(ref_id, eval_id, analyses, config)
    for k in keys:
        assert k in actual, "Key '{}' not present".format(k)
        assert truth[k] == actual[k], "Key '{}' mismatch".format(k)


class TestConfig(ut.TestCase):
    def test_happy_config(self):
        """Ensure happy runs are properly configured"""
        truth = {
            "genome": "resources/references/GRCh38_chr21.fa",
            "genome_index": "resources/references/GRCh38_chr21.fa.fai",
            "strat_tb": "resources/strats/GRCh38_chr21/v3.0-stratifications-GRCh38_chr21.tar.gz",
            "comp_idx": "resources/comparison_variant_callsets/{comp_id}.vcf.gz.tbi",
            "query": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
            "truth": "resources/comparison_variant_callsets/{comp_id}.vcf.gz",
            "truth_regions": "resources/comparison_variant_callsets/{comp_id}.bed",
            "target_regions": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed",
        }
        assert_happy(self, truth, "eval1", "GRCh38_chr21")


if __name__ == "__main__":
    ut.main()
