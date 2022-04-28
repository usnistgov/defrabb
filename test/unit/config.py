import unittest as ut
import pandas as pd
from snakemake.io import load_configfile
from snakemake.utils import validate
from rules.common import get_happy_inputs_inner, get_truvari_inputs_inner, load_analyses


RESOURCES_PATH = "config/resources.yml"
ANALYSES_PATH = "config/analyses.tsv"

R_SCHEMA_PATH = "schema/resources-schema.yml"
A_SCHEMA_PATH = "schema/analyses-schema.yml"


def import_yaml(path, schema):
    config = load_configfile(path)
    validate(config, schema)
    return config


def assert_eval_config(get_fun, keys, self, truth, eval_id, ref_id):
    config = import_yaml(RESOURCES_PATH, R_SCHEMA_PATH)
    analyses = load_analyses(ANALYSES_PATH, A_SCHEMA_PATH)
    analyses = analyses.set_index("eval_id")
    actual = get_fun(ref_id, eval_id, analyses, config)
    for k in keys:
        assert k in actual, "Key '{}' not present".format(k)
        assert (
            truth[k] == actual[k]
        ), f"Key {k} mismatch, expect {truth[k]} got {actual[k]}"


def assert_happy(self, truth, eval_id, ref_id):
    # "I assert you shall be happy with truth (or else...purple burglar alarm)"
    keys = [
        "genome",
        "genome_index",
        "strat_tb",
        "query",
        "query_vcfidx",
        "truth",
        "truth_vcfidx",
        "truth_regions",
        "target_regions",
    ]
    assert_eval_config(
        get_happy_inputs_inner,
        keys,
        self,
        truth,
        eval_id,
        ref_id,
    )
    # config = import_yaml(RESOURCES_PATH, R_SCHEMA_PATH)
    # analyses = load_analyses(ANALYSES_PATH, A_SCHEMA_PATH)
    # analyses = analyses.set_index("eval_id")
    # actual = get_happy_inputs_inner(ref_id, eval_id, analyses, config)
    # for k in keys:
    #     assert k in actual, "Key '{}' not present".format(k)
    #     assert (
    #         truth[k] == actual[k]
    #     ), f"Key {k} mismatch, expect {truth[k]} got {actual[k]}"


def assert_truvari(self, truth, eval_id, ref_id):
    keys = [
        "genome",
        "genome_index",
        "query",
        "query_vcfidx",
        "truth",
        "truth_vcfidx",
        "truth_regions",
    ]
    assert_eval_config(
        get_truvari_inputs_inner,
        keys,
        self,
        truth,
        eval_id,
        ref_id,
    )


class TestConfig(ut.TestCase):
    def test_happy_config(self):
        """Ensure happy runs are properly configured"""
        truth = {
            "genome": "resources/references/GRCh38_chr21.fa",
            "genome_index": "resources/references/GRCh38_chr21.fa.fai",
            "strat_tb": "resources/strats/GRCh38_chr21/v3.0-stratifications-GRCh38_chr21.tar.gz",
            "query": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
            "query_vcfidx": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
            "truth": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",
            "truth_vcfidx": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi",
            "truth_regions": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",
            "target_regions": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.bed",
        }
        assert_happy(self, truth, "eval1", "GRCh38_chr21")

    def test_truvari_config(self):
        """Ensure truvar runs are properly configured"""
        truth = {
            "genome": "resources/references/GRCh38_chr21.fa",
            "genome_index": "resources/references/GRCh38_chr21.fa.fai",
            "query": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz",
            "query_vcfidx": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",
            "truth": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",
            "truth_vcfidx": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi",
            "truth_regions": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",
        }
        assert_truvari(self, truth, "eval5", "GRCh38_chr21")


if __name__ == "__main__":
    ut.main()
