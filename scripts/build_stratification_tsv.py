"""Build a hap.py stratification TSV from named BED files.

hap.py's ``--stratification`` argument takes a 2-column TSV of
``<name>\\t<bed_path>`` where each bed_path is resolved relative to the TSV's own
directory. This utility assembles such a TSV from named BED inputs, writing each
path relative to the output TSV location so hap.py can find them.

It is the shared building block for:
  - #173: stratify a hap.py run by the per-benchmark exclusion BEDs to quantify
    the impact of each exclusion.
  - #59: feed genome-specific stratifications derived from dipcall output.

Standalone and side-effect-free apart from writing the output TSV, so it can be
unit tested and called from a Snakemake rule. NOTE: not yet wired into the
run_happy rule — that wiring (and, for #59, the exact dipcall-derived strats)
needs hap.py end-to-end validation and JZ input; see the issues.

This module was developed with assistance from Claude (Anthropic). All code has
been reviewed and tested by the primary author.
"""

import argparse
import os
import sys
from typing import List, Sequence, Tuple

Entry = Tuple[str, str]


def parse_strat_args(values: Sequence[str]) -> List[Entry]:
    """Parse ``name=bed_path`` strings into (name, path) entries."""
    entries: List[Entry] = []
    for v in values:
        if "=" not in v:
            raise ValueError(f"--strat expects name=path, got: {v!r}")
        name, path = v.split("=", 1)
        if not name or not path:
            raise ValueError(f"--strat expects non-empty name and path, got: {v!r}")
        entries.append((name, path))
    return entries


def format_strat_rows(
    entries: Sequence[Entry], out_dir: str, relative: bool = True
) -> List[str]:
    """Return ``name<TAB>path`` rows, paths relative to out_dir when requested.

    Duplicate names raise — hap.py would otherwise silently use the last one.
    """
    seen = set()
    rows = []
    for name, path in entries:
        if name in seen:
            raise ValueError(f"duplicate stratification name: {name}")
        seen.add(name)
        out_path = os.path.relpath(path, out_dir) if relative else path
        rows.append(f"{name}\t{out_path}")
    return rows


def write_strat_tsv(
    entries: Sequence[Entry], out_path: str, relative: bool = True
) -> None:
    out_dir = os.path.dirname(os.path.abspath(out_path))
    rows = format_strat_rows(entries, out_dir, relative)
    with open(out_path, "w") as fh:
        fh.write("\n".join(rows) + ("\n" if rows else ""))


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--strat",
        action="append",
        default=[],
        metavar="NAME=BED",
        help="A stratification as name=bed_path (repeatable)",
    )
    p.add_argument("-o", "--out", required=True, help="Output stratification TSV path")
    p.add_argument(
        "--no-relative",
        action="store_true",
        help="Write bed paths verbatim instead of relative to the output TSV",
    )
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    entries = parse_strat_args(args.strat)
    if not entries:
        sys.exit("Error: at least one --strat name=bed is required")
    write_strat_tsv(entries, args.out, relative=not args.no_relative)
    print(f"Wrote {len(entries)} stratifications to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
