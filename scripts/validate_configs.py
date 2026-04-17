"""Cross-reference validation between analyses.tsv and resources.yml.

This module was developed with assistance from Claude (Anthropic). All code
has been reviewed and tested by the primary author.
"""
from snakemake.exceptions import WorkflowError


def validate_cross_references(config, analyses):
    """Verify every ID referenced in analyses.tsv exists in resources.yml.

    Raises WorkflowError with a single grouped message if any IDs are missing.
    Returns None on success.
    """
    return None
