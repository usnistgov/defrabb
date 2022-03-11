import unittest as ut
import pandas as pd
from snakemake.io import load_configfile
from snakemake.utils import validate


def import_yaml(path, schema):
    config = load_configfile(path)
    validate(config, schema)
    return config


def import_df(path, schema):
    df = pd.read_csv(path, header=0, sep="\t")
    validate(df, schema)
    return df


class TestConfig(ut.TestCase):
    def test_load(self):
        """Dummy test to make sure pytest is set up correctly"""
        x = import_df("test/unit/config/analyses.tsv", "schema/analyses-schema.yml")
        assert x.shape[0] > 0, "learn to program you dummy"

        y = import_yaml("test/unit/config/resources.yml", "schema/resources-schema.yml")
        assert len(y) > 0, "learn to program you dummy"


if __name__ == "__main__":
    ut.main()
