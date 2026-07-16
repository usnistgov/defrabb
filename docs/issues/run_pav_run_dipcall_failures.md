# `run_dipcall`, `run_pav`, and `truvari anno trf` failures

Diagnosis and remediation from full-pipeline regression tests:
- 20260615_v0.022_fulltest-HG002v1.1 (`run_dipcall` OOM, `run_pav` failures)
- 20260617_v0.022_fulltest-HG002v1.1 (`truvari anno trf` 4-day stall on PAV)

## TL;DR

Three independent root causes, two of which interact:

| # | Rule | Symptom | Root cause |
|---|------|---------|------------|
| A | `run_dipcall` | `make` exits 2, `samtools sort: failed to read header from "-"`, `Killed` in log | **OOM.** `make -j` over-parallelizes memory-heavy minimap2, and the per-rule `mem_mb` reservation is never enforced because the pipeline is launched with `--cores` only (no global `mem_mb` resource), so multiple dipcall jobs run concurrently. |
| B | `run_pav` + `run_dipcall` | PAV output dir contains dipcall files; nested PAV snakemake fails immediately (exit 1) | **Duplicate caller run into a shared directory.** A PAV benchmark whose exclusion set contains `consecutive-svs` forces a *second, redundant* `run_dipcall` to write into the *same* `results/asm_varcalls/{vc_id}/` directory that `run_pav` runs its nested snakemake in. The two collide. |
| C | `run_pav` | Actual failure cause not recoverable from logs | **No logging + fragile nested-snakemake design.** `run_pav` invokes `snakemake -s /opt/pav/Snakefile` via a `script:` with no `log:` directive and no output redirection. The exit-1 cause was never captured. |
| D | `run_pav` | `ModuleNotFoundError: No module named 'snakemake.io.container'; 'snakemake.io' is not a package` | **Snakemake version skew: `script:`-in-container.** `run_pav` was a `script:` rule *with* a `container:`. Snakemake injects a host-generated preamble (`from snakemake.script import snakemake`) that then runs **inside** the becklab/pav container, which ships its own older Snakemake where `snakemake.io` is a module, not a package. The preamble (host Snakemake 8.x) and the container's Snakemake disagree → import crash before the script body runs. |

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

### Bug D — `script:`-in-container Snakemake version skew

The follow-up full-pipeline run (after Fixes A–C) reached `run_pav` and crashed
during the script preamble import:

```
ModuleNotFoundError: No module named 'snakemake.io.container';
'snakemake.io' is not a package
```

`run_pav` was declared as a `script:` rule (`scripts/run_pav.py`) with a
`container: docker://becklab/pav:latest`. For a `script:` rule, Snakemake
generates a preamble on the **host** (`from snakemake.script import snakemake`,
plus a pickled `snakemake` object) and executes it **inside** the container. The
host runs Snakemake 8.x (Python 3.13) where `snakemake.io` is a *package* with a
`container` submodule; the becklab/pav container ships its own older Snakemake
where `snakemake.io` is a single *module*. The host-generated preamble therefore
fails to import against the container's Snakemake — before any of the script body
(config generation, PAV invocation) runs. The host site-packages are bind-mounted
in, but the container's own Snakemake still wins on `sys.path`, so the two are
mismatched either way.

This is the concrete failure mode of the "wonky" nested-snakemake design flagged
in the task: you cannot reliably run Snakemake's `script:` machinery inside a
container whose Snakemake differs from the host's.

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

### Fix D — split host config-generation from the in-container PAV run (decided)

Root cause D means `run_pav` must not be a `script:` rule while it carries a
`container:`. Split the single rule in two:

- **`pav_config`** (no `container:`) — runs `scripts/setup_pav.py` via `script:`
  on the **host**, so Snakemake's script machinery uses the host Snakemake.
  Generates `assemblies.tsv` + `config.json` (declared as outputs).
- **`run_pav`** (`container:`, **`shell:`** not `script:`) — runs the nested
  `snakemake -s /opt/pav/Snakefile` with a plain `shell:` directive. Snakemake
  does **not** inject its script preamble for `shell:` rules, so the nested PAV
  workflow runs against the container's own Snakemake (which is what PAV needs)
  with no host/container version skew. Keeps Fix C's log capture (`> "$LOG"
  2>&1`, absolute log path captured before `cd`, inner PAV `.snakemake/log/*`
  surfaced on failure).

This removes the `script:`-in-container failure mode entirely and is the
"clean entrypoint" follow-up from Fix C item 4.

## Status

Fixes A3, B1, B2, C, and D are implemented.

- **A3** — `run_defrabb` passes `--resources mem_mb=<budget>` (default **80% of
  system memory**, `find_mem_limit`, overridable with `--mem_mb`), so per-rule
  `mem_mb` reservations serialize concurrent `run_dipcall` jobs instead of
  collectively OOMing.
- **D** — `run_pav` no longer uses `script:`-in-container; config generation is a
  host-side `pav_config` rule and PAV runs via a bare `shell:`. Dry-run confirms
  the `pav_config → run_pav` edge resolves with no duplicate dipcall run in the
  PAV directory.

## Changes applied (worktree `6d6f8ea` dev line)

- `rules/asm-varcall.smk`: `make -j{params.make_jobs}` (Fix A1); `run_pav` split
  into `pav_config` (host `script:`) + `run_pav` (in-container `shell:`) with a
  `log:` directive (Fix C + D).
- `scripts/setup_pav.py` (new, replaces `scripts/run_pav.py`): host-side
  generation of `assemblies.tsv` + `config.json` (Fix D). The nested PAV
  snakemake is now invoked directly from `run_pav`'s `shell:`, capturing output
  to the rule log and surfacing inner PAV `.snakemake/log/*` on failure (Fix C).
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
  has a `log:` and captures PAV output (Bug C); `run_pav` uses `shell:` (not
  `script:`) inside the container and `pav_config` generates inputs on the host
  (Bug D); `run_dipcall`'s
  `make -j` uses the jobs knob (Bug A); PAV exclusion sets omit dipcall-only
  exclusions (Mod 1); `get_consecutive_svs` reuses the resolved dipcall run and
  never references the benchmark's own `{vc_cmd}` hap BAMs (Mod 2 / Bug B).
- `.tests/unit/test_trfanno_merge.py` — VCF header merging and BCF roundtrip for
  `merge_trfanno_vcfs.py`.

## E — truvari anno trf 4-day stall on PAV (20260617 fulltest)

**Root cause:** `truvari anno trf` (v4.3.0) stalled 4+ days on PAV callsets due
to O(n²) `edlib.align()` in Python estimation path (`TRFAnno.ins_estimate_anno`)
when processing 420 kb insertions. PAV emits full-length insertions with no size
cap; the largest 15 are ≥100 kb (max 420 kb). Truvari parallelizes by reference
chunks (`-C 5MB` default), so a single chunk holding a pathological variant
stalls the entire job (classic tail latency).

**Solution (20260708, commit 633196d):**
- Size-cap insertions routed to TRF (100 kb default, configurable via
  `truvari_anno_max_ins_length` in `_vcf_processing_params`).
- Oversized insertions bypass TRF annotation but remain in output (un-annotated).
- Three-rule workflow:
  1. `split_vcf_by_ins_size`: bcftools split on `STRLEN(ALT) > cap`
  2. `run_truvari_anno_trf`: input `trf_insize.vcf.gz`, output
     `trf_insize_annotated.vcf`; threads now configurable
     (`_truvari_anno_threads: 8`, was hardcoded 5); added `-C 1` (1 MB chunks for
     load balance) and `-l 0.5` (TRF memory cap)
  3. `merge_trfanno`: Python/pysam merge with proper BCF header handling (create
     new records in output header context instead of `entry.translate()` to avoid
     BCF INFO tag ID corruption when headers have different field counts); added
     END to header (pysam auto-generates for symbolic alleles with SVLEN)

**Verified:** 20260617 fulltest HG002v1.1 PAV: 6,627,937 variants (6,627,915
annotated + 22 oversize), pipeline completes successfully. TRF annotation now
bounds at ~5 min (was 4+ days).

## F — run_pav `call_cigar` concurrency failure inside the nested PAV Snakemake (20260708 fulltest)

**Symptom:** `run_pav` for `GRCh38_HG002-T2TQ100v1.1-pav` failed ~39% into the
nested PAV workflow. Three `call_cigar` jobs (h2 batches 1, 4, 8) reported
`Error in rule call_cigar:` with **no message and no traceback**, all at the
exact same second, ~24 s after launch — while their sibling h2 batches ran a
further ~4 min and succeeded, and every h1 batch (launched ~11 s later)
succeeded. h1's downstream ran to completion; h2 stalled because its
`call_cigar_merge` needs all 10 batches. With `-k` the run drained all h1 work,
then exited failed. Every subsequent rerun then hit a `LockException` because
the killed job left a stale `.snakemake/locks/` entry (and a 0-byte incomplete
`aligned_tig_h2.sam.gz`).

**Root cause — the failure is concurrency-conditional; most likely a transient
memory spike (OOM).** Evidence:

- An empty error message with **no Python traceback** is the signature of a
  signal kill (SIGKILL), not a Python exception (which PAV's `run:`-directive
  jobs surface with a full traceback — confirmed separately). Three independent
  jobs dying at the *same second* points to a shared external event; the OOM
  killer reaping the top-RSS processes fits.
- **Reproduced the conditionality:** re-running the failing batch (h2 batch 4)
  *in isolation* (single `call_cigar`, no concurrency, fresh alignment)
  **succeeds**. So there is no data/code defect in the batch — the failure only
  occurs under the full ~20-way `call_cigar` + `call_lg_discover` concurrency of
  the real run.

In the container PAV (`/opt/pav/rules/call.snakefile`), `call_cigar` declares
**no `threads:` and no `resources:`**, so the nested `snakemake -c {_pav_threads}`
(24) schedules up to 24 at once. Each accumulates **one `pandas.Series` per
called variant** in a Python list (`pavlib/cigarcall.py::make_insdel_snv_calls`,
`df_snv_list.append(pd.Series(...))`) and concatenates at the end, so both
steady-state and the terminal concat spike scale with the batch's variant count
(millions for T2T-Q100 vs GRCh38). The highest-variant batches (the 3 that died)
peak highest; ~20 batches started together would spike together. The batch count
is a hardcoded constant (`pavlib.cigarcall.CALL_CIGAR_BATCH_COUNT = 10`), not
config-tunable, and the outer `--resources mem_mb` budget governs only *outer*
jobs — it does not propagate into the nested PAV Snakemake, so nothing bounds the
inner concurrency.

**Not fully confirmed:** kernel OOM logs (dmesg/journalctl) were not accessible
on this host, and observed *steady-state* per-batch RSS was modest (hundreds of
MB to ~2.7 GB); the transient concat spike was not captured. So "OOM" is the
most probable mechanism given the signature, not a log-verified fact. What *is*
verified is that the failure is triggered by concurrency and disappears without
it — which is exactly what the fix removes.

**Fix (config-driven concurrency cap):** `run_pav` now passes
`--set-threads call_cigar={_pav_cigar_threads}` to the nested PAV Snakemake.
Since `call_cigar` has no declared threads, `--set-threads` is the only lever
that bounds its concurrency (PAV exposes no batch-count or per-rule memory
config). With `_pav_threads: 24` and `_pav_cigar_threads: 8`, at most
`floor(24/8) = 3` `call_cigar` jobs run at once. This directly removes the
proven trigger (high `call_cigar` concurrency) and bounds any aggregate memory
spike. Lower `_pav_cigar_threads` for more concurrency (faster), raise it to
serialize harder on memory-tight hosts.

**Changes applied:**
- `config/resources.yml`: new `_pav_cigar_threads: 8`.
- `schema/resources-schema.yml`: `_pav_cigar_threads` added (required + typed).
- `rules/asm-varcall.smk` `run_pav`: `--set-threads call_cigar={params.cigar_threads}`.

**Operational note (stale lock / incomplete temp on restart):** a killed PAV run
leaves `.snakemake/locks/` populated (next run fails with `LockException`) and
can leave 0-byte incomplete temp outputs. Recovery: remove the lock files (or run
the nested PAV with `--unlock`) and delete any 0-byte `temp/.../*.sam.gz` before
rerunning. Note that PAV's `temp()` intermediates (e.g. `contigs_{hap}.fa.gz`)
are cleaned on success, so a `--ri` resume re-runs alignment from scratch — an
interrupted PAV run effectively restarts the whole ~8 h job.

**Validation:**
- Failing h2 batch 4 runs to completion in isolation (single batch, fresh
  alignment) — confirms the failure is concurrency-conditional, not a data bug.
- `--set-threads call_cigar={n}` parses in the container PAV Snakemake and a
  dry-run shows `call_cigar` assigned the capped thread count.
- Structural regression tests added (`.tests/unit/test_asm_varcall_rules.py`).
- **Confirmed end-to-end (20260715 rerun, `_pav_cigar_threads: 8`):** `run_pav`
  cleared the `call_cigar` stage with 0 errors and `call_cigar_merge` ran for
  *both* haplotypes (previously h2 could never merge — batches 1/4/8 died). PAV
  ran to 550/550 steps and produced `HG002.vcf.gz` (6,744,518 records). The full
  pipeline completed cleanly (outer `all` rule, 34/34 steps, no errors under
  `--keep-going`) in ~11h20m; all 7 evaluations produced (3 hap.py smvar + 4
  Truvari stvar, incl. the PAV benchmark: precision 0.844 / recall 0.934 /
  F1 0.887 vs v5.0q).
