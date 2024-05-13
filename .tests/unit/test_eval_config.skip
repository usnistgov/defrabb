import pytest
from snakemake.io import load_configfile  
from snakemake.utils import validate 
import importlib.machinery
import importlib.util

## Importing smk files as module
loader = importlib.machinery.SourceFileLoader('common', 'rules/common.smk')
spec = importlib.util.spec_from_loader('common', loader)
common = importlib.util.module_from_spec(spec)
loader.exec_module(common)
  
RESOURCES_PATH = "config/resources.yml"  
ANALYSES_PATH = "config/analyses.tsv"  
  
R_SCHEMA_PATH = "schema/resources-schema.yml"  
A_SCHEMA_PATH = "schema/analyses-schema.yml"  
  
@pytest.fixture  
def yaml_configs():  
    config = load_configfile(RESOURCES_PATH)  
    validate(config, R_SCHEMA_PATH)  
    analyses = common.load_analyses(ANALYSES_PATH, A_SCHEMA_PATH)  
    analyses = analyses.set_index("eval_id")  
    return config, analyses  
  
def assert_eval_config(get_fun, keys, truth, eval_id, ref_id, config, analyses):  
    actual = get_fun(ref_id, eval_id, analyses, config)  
    for k in keys:  
        assert k in actual, f"Key '{k}' not present"  
        assert truth[k] == actual[k], f"Key {k} mismatch, expect {truth[k]} got {actual[k]}"  
  
def test_happy_config(yaml_configs):  
    config, analyses = yaml_configs  
    truth = {  
        "genome": "resources/references/GRCh38_chr21.fa",  
        "genome_index": "resources/references/GRCh38_chr21.fa.fai",  
        "strat_tb": "resources/strats/GRCh38_chr21/genome-stratifications-GRCh38@all.tar.gz",  
        "query": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",  
        "query_vcfidx": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",  
        "truth": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",  
        "truth_vcfidx": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi",  
        "truth_regions": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",  
        "target_regions": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.bed",  
    }  
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
    assert_eval_config(common.get_happy_inputs_inner, keys, truth, "eval1", "GRCh38_chr21", config, analyses)  
  
def test_truvari_config(yaml_configs):  
    config, analyses = yaml_configs  
    truth = {  
        "genome": "resources/references/GRCh38_chr21.fa",  
        "genome_index": "resources/references/GRCh38_chr21.fa.fai",  
        "query": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz",  
        "query_vcfidx": "results/draft_benchmarksets/{bench_id}/{ref_id}_{asm_id}_{bench_type}_{vc_cmd}-{vc_param_id}.vcf.gz.tbi",  
        "truth": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz",  
        "truth_vcfidx": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.vcf.gz.tbi",  
        "truth_regions": "resources/comparison_variant_callsets/{ref_id}_{comp_id}.bed",  
    }  
    keys = [  
        "genome",  
        "genome_index",  
        "query",  
        "query_vcfidx",  
        "truth",  
        "truth_vcfidx",  
        "truth_regions",  
    ]  
    assert_eval_config(common.get_truvari_inputs_inner, keys, truth, "eval6", "GRCh38_chr21", config, analyses)  
