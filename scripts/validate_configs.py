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
    errors = {
        "assemblies": {},
        "references": {},
        "exclusion_sets": {},
        "comparisons": {},
    }

    asm_ids = set(config["assemblies"])
    ref_ids = set(config["references"])
    excl_set_ids = set(config["exclusion_set"])

    for eval_id, row in analyses.iterrows():
        if row["asm_id"] not in asm_ids:
            errors["assemblies"].setdefault(row["asm_id"], []).append(eval_id)
        if row["ref"] not in ref_ids:
            errors["references"].setdefault(row["ref"], []).append(eval_id)
        if row["exclusion_set"] != "none" and row["exclusion_set"] not in excl_set_ids:
            errors["exclusion_sets"].setdefault(row["exclusion_set"], []).append(
                eval_id
            )

    if any(errors.values()):
        raise WorkflowError(_format_grouped_errors(errors))
    return None


def _format_grouped_errors(errors):
    sections = []
    if errors["assemblies"]:
        sections.append(
            _format_section(
                "Missing assemblies (resources.yml:assemblies)",
                errors["assemblies"],
            )
        )
    if errors["references"]:
        sections.append(
            _format_section(
                "Missing references (resources.yml:references)",
                errors["references"],
            )
        )
    if errors["exclusion_sets"]:
        sections.append(
            _format_section(
                "Missing exclusion sets (resources.yml:exclusion_set)",
                errors["exclusion_sets"],
            )
        )
    return (
        "Config validation failed: missing cross-references in resources.yml.\n\n"
        + "\n\n".join(sections)
    )


def _format_section(title, missing):
    lines = [f"{title}:"]
    for missing_id, eval_ids in missing.items():
        lines.append(f"  - {missing_id} (used by eval_ids: {', '.join(eval_ids)})")
    return "\n".join(lines)
