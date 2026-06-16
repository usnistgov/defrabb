# `run_dipcall` and `run_pav` failures (full-pipeline test 20260615_v0.022)

Diagnosis and remediation plan from the whole-genome regression run
`/defrabb_runs/runs/20260615_v0.022_fulltest-HG002v1.1` (HG002 T2T-Q100 v1.1,
GRCh37 / GRCh38 / CHM13v2.0; one PAV benchmark on GRCh38).

## TL;DR

Three independent root causes, two of which interact:

| # | Rule | Symptom | Root cause |
|---|------|---------|------------|
| A | `run_dipcall` | `make` exits 2, `samtools sort: failed to read header from "-"`, `Killed` in log | **OOM.** `make -j` over-parallelizes memory-heavy minimap2, and the per-rule `mem_mb` reservation is never enforced because the pipeline is launched with `--cores` only (no global `mem_mb` resource), so multiple dipcall jobs run concurrently. |
| B | `run_pav` + `run_dipcall` | PAV output dir contains dipcall files; nested PAV snakemake fails immediately (exit 1) | **Duplicate caller run into a shared directory.** A PAV benchmark whose exclusion set contains `consecutive-svs` forces a *second, redundant* `run_dipcall` to write into the *same* `results/asm_varcalls/{vc_id}/` directory that `run_pav` runs its nested snakemake in. The two collide. |
| C | `run_pav` | Actual failure cause not recoverable from logs | **No logging + fragile nested-snakemake design.** `run_pav` invokes `snakemake -s /opt/pav/Snakefile` via a `script:` with no `log:` directive and no output redirection. The exit-1 cause was never captured. |

## Evidence

### Bug A — dipcall OOM

`logs/asm_varcalls/CHM13_..._dipcall-z2k.log` (the run that ran a full ~31 min then failed):

```
... minimap2 -a -xasm5 ... (4 mapping jobs launched in parallel) ...
Killed
Killed
... samtools sort -m4G --threads 4 -o ...hap1.bam -
[W::hts_set_opt] Cannot change block size for this format
samtools sort: failed to read header from "-"
make: *** [...hap1.bam] Error 1
```

`Killed` = SIGKILL from the Linux OOM killer. The killed minimap2 left a
truncated `hap1.sam.gz`, so the downstream `samtools sort` reading that stream
got an empty/garbage header.

Per-haplotype minimap2 logs show **Peak RSS ~25–28 GB each** (minimap2 2.31,
whole-genome asm5). Contributing factors:

- `rules/asm-varcall.smk` runs `make -j{params.ts}` where `params.ts =
  _dipcall_threads = 5`. dipcall's generated makefile has 4 independent
  minimap2 targets (2 PAF + 2 SAM), so all 4 run at once → ~4 × 28 GB ≈ 112 GB
  for a single dipcall job.
- The rule reserves `mem_mb = _dipcall_jobs * _dipcall_mem = 4 * 32000 =
  128 GB`, but `run_defrabb` (and the documented `snakemake … --cores N`
  invocation) **never passes a global `--resources mem_mb=…`**, so Snakemake
  ignores the reservation and schedules **4 `run_dipcall` jobs concurrently**
  (4 × 20 threads = 80 ≤ 84 cores). Aggregate demand ≈ 4 × 112 GB → OOM.
- The `make -j` count uses the *thread* param (`params.ts`), not the *jobs*
  param (`params.make_jobs` / `_dipcall_jobs`) the memory math is based on — the
  two knobs are wired inconsistently.

### Bug B — shared output directory between PAV and dipcall

The PAV analysis row uses `vc_id = GRCh38_HG002-T2TQ100v1.1-pav` and exclusion
set `HG002Q100stvarv0.022`, which includes `consecutive-svs` and
`dipcall-bugs-T2TACE`.

- `rules/exclusions_download.smk::get_consecutive_svs` requires
  `results/asm_varcalls/{vc_id}/{ref}_{asm}_{vc_cmd}-{vc_param}.hap1.bam` /
  `.hap2.bam` — produced **only by `run_dipcall`**.
- For `vc_cmd=pav` that path lives **inside the PAV `vc_id` directory**, so
  Snakemake schedules `run_dipcall` (jobid 288, `vc_cmd=pav`) into
  `results/asm_varcalls/GRCh38_HG002-T2TQ100v1.1-pav/` at the same time as
  `run_pav` (jobid 95).
- `run_pav` does `cd results/asm_varcalls/{vc_id}` and runs a nested
  `snakemake -s /opt/pav/Snakefile` there. Confirmed collision: the PAV dir
  contains dipcall artifacts (`*_pav-giab.mak`, `*.hap2.sam.gz`,
  `*.hap2.paf.gz`, `*.hap2.var.gz`, per-hap `.log`).

So even ignoring the OOM, a PAV benchmark that needs the assembly-BAM-derived
exclusions runs dipcall and PAV in the same directory — fragile by design.

### Bug C — `run_pav` is undebuggable and fragile

`rules/asm-varcall.smk::run_pav` (script `scripts/run_pav.py`):

- No `log:` directive; the inner
  `OPENSSL_CONF=/dev/null snakemake -s /opt/pav/Snakefile … -c {threads}` writes
  only to stdout/stderr. The exit-1 cause is not in any log file.
- Runs a **nested Snakemake** (PAV's own workflow) inside the container with
  `--rerun-triggers mtime -k -w 20`, cwd set to the (shared) output dir. This is
  the "wonky" usage noted in the task: PAV is a standalone pipeline, not a
  module/subworkflow, so there is no DAG integration, no per-rule resource
  control, no isolation, and no log capture.

## Remediation plan

### Fix A — bound dipcall memory (highest priority; this is what failed)

1. **Wire `make -j` to the jobs knob, not the thread knob.** In
   `rules/asm-varcall.smk` change `make -j{params.ts}` → `make
   -j{params.make_jobs}` so the parallelism that drives memory matches the
   `mem_mb` reservation math.
2. **Reduce default per-job parallelism / raise per-job mem.** With minimap2
   2.31 at ~28 GB/mapping, `_dipcall_jobs: 2` (one per haplotype-ish) with
   `_dipcall_mem: 32000` keeps a single dipcall job ≲ 64 GB. Re-tune in
   `config/resources.yml`.
3. **Enforce a global memory budget.** Make `run_defrabb` (and the documented
   invocation) pass `--resources mem_mb=<host_total_mb>` so the per-rule
   `mem_mb` reservations actually serialize concurrent `run_dipcall` jobs.
   Source the value from `resources.yml` (e.g. a new `_host_mem_mb`) or detect
   it. Without this, no amount of per-rule tuning prevents N-way concurrent OOM.

### Fix B — never run a duplicate caller into the PAV directory (decided)

Two coordinated changes (per maintainer direction):

- **B1 — drop dipcall-specific exclusions from PAV exclusion sets.**
  `consecutive-svs` and `dipcall-bugs-T2TACE` are dipcall artifacts and do not
  belong in a PAV benchmark's exclusion set. Add a PAV-specific exclusion set
  (e.g. `HG002Q100stvarv0.022pav`) that omits them — mirroring the existing
  `HG008stvarv0.021pav` convention — and point PAV analysis rows at it.
- **B2 — cross-caller exclusions reuse the existing run, never regenerate.**
  Any exclusion that needs the *other* caller's variant calls (today only
  `consecutive-svs`, which needs dipcall hap BAMs) must resolve its input to the
  **existing `results/asm_varcalls/` run for the same reference + assembly +
  appropriate caller params**, instead of the current benchmark's own `vc_id`.
  This prevents a duplicate caller run (and the shared-directory collision)
  entirely. Implemented via a `vc_tbl` lookup helper
  (`get_asm_varcall_run` in `rules/helpers_varcall.smk`) that finds the unique
  matching run and asserts uniqueness; `get_consecutive_svs` now consumes that
  run's hap BAMs. For a dipcall benchmark the lookup returns its own run, so
  behavior is unchanged.

### Fix C — make `run_pav` robust and debuggable

1. Add a `log:` directive and redirect the nested snakemake's stdout/stderr to
   it (`&> {log}` / `> {log} 2>&1`).
2. Run PAV in an isolated working directory (ties in with B1).
3. On failure, surface the inner PAV `.snakemake/log/*.log` into the rule log
   (e.g. `cat` the latest PAV log on non-zero exit) so the real cause is always
   captured.
4. Longer term, evaluate running PAV via a clean CLI entrypoint /
   per-run temp dir wrapper rather than a `script:`-driven nested snakemake.

## Status

Fix A3 is implemented: `run_defrabb` now passes `--resources
mem_mb=<budget>` to the pipeline, defaulting to **80% of system memory**
(`find_mem_limit`), overridable with `--mem_mb`. This enforces the per-rule
`mem_mb` reservations so concurrent `run_dipcall` jobs serialize instead of
collectively OOMing.

Remaining longer-term item: PAV is still a fragile nested-snakemake (run via a
`script:`); consider running it in an isolated working dir / via a clean
entrypoint.

## Changes applied (worktree `6d6f8ea` dev line)

- `rules/asm-varcall.smk`: `make -j{params.make_jobs}` (Fix A1); `run_pav` gains
  a `log:` directive (Fix C).
- `scripts/run_pav.py`: nested PAV snakemake output captured to the rule log;
  inner PAV `.snakemake/log/*` surfaced on failure (Fix C).
- `scripts/varcall_lookup.py` (new) + `rules/helpers_varcall.smk`
  `get_asm_varcall_run`: resolve the existing run for a (ref, asm, caller),
  asserting uniqueness (Fix B2).
- `rules/helpers_bench.smk` `get_consecutive_svs_bams` +
  `rules/exclusions_download.smk` `get_consecutive_svs`: consume the resolved
  dipcall run's hap BAMs instead of the benchmark's own `{vc_cmd}` path (Fix B2).
- `config/resources.yml`: new `HG002Q100stvarv0.022pav` exclusion set omitting
  `consecutive-svs` and `dipcall-bugs-T2TACE` (Fix B1).
- `run_defrabb`: `find_mem_limit` + `--mem_mb` arg; pipeline command now passes
  `--resources mem_mb=<80% of system memory>` (Fix A3).

### Follow-up on the run/test branch (not in this worktree)

The failed run used analyses table
`analyses_20260615_v0.022_fulltest-HG002v1.1.tsv` (on the test clone), whose PAV
row points at `HG002Q100stvarv0.022-selfdiscrep`. When regenerating that table,
point the PAV row at a PAV-specific set (e.g. `HG002Q100stvarv0.022pav`, plus the
`self-discrep` member if desired) so the dipcall-only exclusions are dropped.

## Tests added

- `.tests/unit/test_varcall_lookup.py` — behavior of the cross-caller resolver:
  a PAV benchmark resolves to the matching dipcall run; a dipcall benchmark
  resolves to itself; missing/ambiguous runs raise.
- `.tests/unit/test_asm_varcall_rules.py` — structural invariants: `run_pav`
  has a `log:` and the script captures PAV output (Bug C); `run_dipcall`'s
  `make -j` uses the jobs knob (Bug A); PAV exclusion sets omit dipcall-only
  exclusions (Mod 1); `get_consecutive_svs` reuses the resolved dipcall run and
  never references the benchmark's own `{vc_cmd}` hap BAMs (Mod 2 / Bug B).
