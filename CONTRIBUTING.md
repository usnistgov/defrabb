# Contributing to DeFrABB

DeFrABB is developed by the NIST Genome in a Bottle team, primarily for
internal benchmark development. The code is public for transparency and
reproducibility, and support is best-effort.

## Questions and issues

Please open an issue on the [GitHub repository](https://github.com/usnistgov/defrabb/issues)
for questions about the pipeline, how a benchmark was generated, or suspected
bugs. Include the DeFrABB version (git tag or commit), the analyses table, and
relevant log output.

Questions about the benchmark sets themselves (for example, a suspected error in
the HG002 v5.0q benchmark) can also be sent to the GIAB team through the
contact information in the benchmark release README.

## Pull requests

Pull requests are welcome but are reviewed as time allows. Development happens
on NIST-internal GitLab, and accepted changes are applied there and appear on
GitHub at the next release. Before submitting:

```sh
pytest .tests          # unit tests
snakefmt --check .     # Snakemake formatting
black --check scripts/ # Python formatting
```

Describe the change's effect on workflow outputs, and note any changes to
`config/resources.yml` or the analyses-table schema.
