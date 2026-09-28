"""Pure lookup helpers for resolving assembly variant-call runs.

Used by the Snakemake helpers so that exclusions which depend on *another*
assembly variant caller's output reuse an existing run (same reference +
assembly + caller) instead of triggering a duplicate caller run. Kept free of
Snakemake globals so it can be unit-tested directly.

Developed with assistance from Claude (Anthropic); reviewed by the primary author.
"""


def resolve_asm_varcall_run(vc_tbl, ref, asm_id, vc_cmd, prefer_vc_param_id=None):
    """Resolve the unique asm_varcall run matching (ref, asm_id, vc_cmd).

    ``vc_tbl`` is the variant-call table indexed by ``vc_id`` with ``ref``,
    ``asm_id``, ``vc_cmd``, and ``vc_param_id`` columns.

    If multiple runs match (parameter sweeps), use the one whose
    ``vc_param_id`` equals ``prefer_vc_param_id``.

    Returns ``(vc_id, vc_param_id)``. Raises ``ValueError`` when no run matches,
    or when several match and none can be selected via ``prefer_vc_param_id`` —
    guessing would silently tie exclusions to the wrong caller parameters.
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
        preferred = matches[matches["vc_param_id"] == prefer_vc_param_id]
        if len(preferred) != 1:
            raise ValueError(
                f"Multiple {vc_cmd} runs found for ref={ref} asm_id={asm_id} "
                f"({', '.join(matches.index)}) and none uniquely matches "
                f"vc_param_id={prefer_vc_param_id!r}; cannot choose which run's "
                f"calls the exclusion should use."
            )
        matches = preferred
    return matches.index[0], matches.iloc[0]["vc_param_id"]
