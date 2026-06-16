"""Pure lookup helpers for resolving assembly variant-call runs.

Used by the Snakemake helpers so that exclusions which depend on *another*
assembly variant caller's output reuse an existing run (same reference +
assembly + caller) instead of triggering a duplicate caller run. Kept free of
Snakemake globals so it can be unit-tested directly.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""


def resolve_asm_varcall_run(vc_tbl, ref, asm_id, vc_cmd):
    """Resolve the unique asm_varcall run matching (ref, asm_id, vc_cmd).

    ``vc_tbl`` is the variant-call table indexed by ``vc_id`` with ``ref``,
    ``asm_id``, ``vc_cmd``, and ``vc_param_id`` columns.

    Returns ``(vc_id, vc_param_id)``. Raises ``ValueError`` when no run or more
    than one run matches, so the "appropriate" parameters are unambiguous.
    """
    matches = vc_tbl[
        (vc_tbl["ref"] == ref)
        & (vc_tbl["asm_id"] == asm_id)
        & (vc_tbl["vc_cmd"] == vc_cmd)
    ]
    if matches.empty:
        raise ValueError(
            f"No {vc_cmd} run found for ref={ref} asm_id={asm_id}; an exclusion "
            f"needs its variant calls. Add a {vc_cmd} variant-call row for this "
            f"reference + assembly to the analyses table."
        )
    if len(matches) > 1:
        raise ValueError(
            f"Multiple {vc_cmd} runs for ref={ref} asm_id={asm_id}: "
            f"{matches.index.tolist()}. Cannot pick the appropriate one for a "
            f"cross-caller exclusion; disambiguate the analyses table."
        )
    return matches.index[0], matches.iloc[0]["vc_param_id"]
