#!/usr/bin/env bash
###############################################################################
# run_full_pipeline_test.sh
#
# Stand up a *whole-genome* DeFrABB regression run that exercises the larger
# features CI does not cover. CI only runs the bundled chr21 data plus
# formatting/lint and a dry-run, so changes to the real variant-calling,
# exclusion, stratification, and evaluation paths can pass CI and still break
# on a full genome. Run this by hand AFTER a large change to the codebase.
#
# WARNING: this configures a compute- and time-intensive whole-genome run
# (3 references x dipcall smvar+stvar, plus a PAV stvar row). It is NOT for
# routine use. By default it only *sets up and validates* the run and prints
# the command to launch it (see --go to launch immediately).
#
# What it exercises (this dev round):
#   #194      stvar_v5 truvari bench profile (bnddist -1) on the stvar rows
#   #59/#173  genome_specific_strats (opt-in) on the dipcall smvar hap.py rows
#   #188      HG002Q100stvarv0.022 exclusion set (mosaic dropped) on stvar rows
#   #192      self-discrepancy symbolic-allele (<INV>) filter, via a PAV stvar
#             row whose exclusion set includes self-discrep
#
# How it runs (regular clone-into-runid workflow, per the README):
#   1. Clones the current repo (committed HEAD) into <DEST>/<RUNID>.
#   2. Inside that clone, generates the dated analyses table, enables the
#      genome_specific_strats config flag, adds the self-discrep exclusion set,
#      and COMMITS those changes on a run branch.
#   3. Validates + dry-runs in the clone, then prints the run command.
#   All run artifacts (results/, logs/, the commit) live in the clone, keeping
#   the source checkout clean.
#
# IMPORTANT: the clone copies *committed* HEAD only. Commit the codebase
# changes you want to test BEFORE running this, or they will not be exercised.
# The script warns if the source tree has uncommitted tracked changes.
#
# Usage:
#   conda activate snakemake   # or ~/miniforge3/envs/snakemake
#   scripts/run_full_pipeline_test.sh [RUNID] [--dest-dir DIR] [--go] [--cores N]
#
#   RUNID        run id, format YYYYMMDD_v#.###_brief-id
#                (default: <today>_v0.022_fulltest-HG002v1.1)
#   --dest-dir   parent directory for the clone (default: repo parent dir, ..)
#   --go         after a clean validate + dry-run, launch the real run
#   --cores      cores/jobs for the run + dry-run (default: 1)
###############################################################################
set -euo pipefail

TEMPLATE_REL="config/analyses_20260323_v0.020_HG002v1.1-CI.tsv"
NEW_EXCLUSION_SET="HG002Q100stvarv0.022-selfdiscrep"

GO=0
CORES=1
RUNID=""
DEST_DIR=".."

die() { echo "ERROR: $*" >&2; exit 1; }
info() { echo ">> $*"; }
usage() { sed -n '2,51p' "$0" | sed 's/^# \{0,1\}//'; }

# --- parse args ------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --go)    GO=1; shift ;;
        --cores|-j) CORES="${2:?--cores needs a value}"; shift 2 ;;
        --dest-dir) DEST_DIR="${2:?--dest-dir needs a value}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1 (try --help)" ;;
        *) [[ -z "$RUNID" ]] || die "unexpected extra argument: $1"; RUNID="$1"; shift ;;
    esac
done
RUNID="${RUNID:-$(date +%Y%m%d)_v0.022_fulltest-HG002v1.1}"

[[ "$RUNID" =~ ^[0-9]{8}_v[0-9]+\.[0-9]{3}_[A-Za-z0-9.-]+$ ]] \
    || die "RUNID '$RUNID' does not match YYYYMMDD_v#.###_brief-id"

# --- preflight -------------------------------------------------------------
[[ -f Snakefile && -x run_defrabb ]] \
    || die "run this from the repository root (Snakefile / run_defrabb not found)"
[[ -f "$TEMPLATE_REL" ]] || die "template not found: $TEMPLATE_REL"
command -v snakemake >/dev/null 2>&1 \
    || die "snakemake not on PATH; activate the env first: conda activate snakemake (or ~/miniforge3/envs/snakemake)"
command -v git >/dev/null 2>&1 || die "git not found on PATH"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "not inside a git work tree"

REPO_ROOT="$(pwd)"

# The clone takes committed HEAD; warn if tracked changes are still uncommitted.
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
    echo "WARNING: the source repo has uncommitted tracked changes. The clone"
    echo "         copies committed HEAD only, so these will NOT be tested:"
    git status --short --untracked-files=no | sed 's/^/         /'
    echo "         Commit them first if you want this run to exercise them."
    echo
fi

CLONE_DIR="${DEST_DIR%/}/${RUNID}"
[[ -e "$CLONE_DIR" ]] && die "clone target already exists: $CLONE_DIR (remove it or pick a new RUNID)"

# --- 1. clone the repo into the run directory ------------------------------
info "cloning $REPO_ROOT -> $CLONE_DIR"
git clone --quiet "$REPO_ROOT" "$CLONE_DIR"
CLONE_ABS="$(cd "$CLONE_DIR" && pwd)"
cd "$CLONE_ABS"
git checkout --quiet -b "run/${RUNID}"

OUTFILE="config/analyses_${RUNID}.tsv"
RESOURCES="config/resources.yml"

# --- 2. generate the analyses table (in the clone) ------------------------
info "generating $OUTFILE from $TEMPLATE_REL"
# Keep smvar rows as-is (genome_specific_strats is driven by the config flag,
# not the table). For stvar rows: switch to the stvar_v5 truvari profile (#194)
# and the v0.022 exclusion set (#188).
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
    if ($10 == "stvar") { $4 = "stvar_v5"; $13 = "HG002Q100stvarv0.022" }
    print
}' "$TEMPLATE_REL" > "$OUTFILE"

# Append one PAV stvar row (GRCh38) to exercise the PAV caller and the #192
# symbolic-allele self-discrepancy filter (its exclusion set includes
# self-discrep; PAV emits <INV> records). The assembly HG2-T2TQ100-V1.1 is
# caller-agnostic; only vc_cmd/vc_id/vc_param differ from the dipcall rows.
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "GRCh38_HG002_T~T2TQ100v1.1_Q~v5.0q-stvar-pav_TR~v5.0q" \
    "GRCh38_HG002-T2TQ100v1.1-pav_stvar-excluded" \
    "truvari" "stvar_v5" "v5.0q-stvar" "FALSE" "TRUE" "TRUE" \
    "GRCh38_HG002-T2TQ100v1.1-pav" "stvar" "norm_xy_trf_sv_lcr" "exclude" \
    "$NEW_EXCLUSION_SET" "HG2-T2TQ100-V1.1" "GRCh38" "pav" "giab" "default" \
    >> "$OUTFILE"

info "wrote $(($(wc -l < "$OUTFILE") - 1)) analysis rows"

# --- 3. patch resources.yml (in the clone) --------------------------------
# 3a. add the new self-discrep-bearing exclusion set after the v0.022 block.
grep -qE '^  HG002Q100stvarv0\.022:' "$RESOURCES" \
    || die "could not find HG002Q100stvarv0.022 in $RESOURCES to anchor the new set"
BLOCK_FILE="$(mktemp)"
cat > "$BLOCK_FILE" <<EOF
  # full-pipeline regression harness (#192): stvarv0.022 + self-discrep so the
  # symbolic-allele self-discrepancy filter runs on PAV <INV> records.
  ${NEW_EXCLUSION_SET}:
    - segdups
    - tandem-repeats
    - satellites
    - gaps
    - flanks
    - VDJ
    - consecutive-svs
    - dipcall-bugs-T2TACE
    - HG002Q100-pav-discrep-smvar
    - HG002Q100-pav-discrep-stvar
    - HG002Q100-pav-inversions
    - HG002Q100-errors
    - TSPY2-segdups
    - self-discrep
EOF
TMP="$(mktemp)"
awk -v bf="$BLOCK_FILE" '
    function emit(   l){ while ((getline l < bf) > 0) print l; close(bf) }
    /^  HG002Q100stvarv0\.022:/ { seen=1 }
    seen && !ins && /^  [^ ]/ && !/^  HG002Q100stvarv0\.022:/ { emit(); ins=1 }
    { print }
    END { if (!ins) emit() }
' "$RESOURCES" > "$TMP"
mv "$TMP" "$RESOURCES"
rm -f "$BLOCK_FILE"
info "added exclusion set $NEW_EXCLUSION_SET to $RESOURCES"

# 3b. enable the genome_specific_strats config flag (#59/#173), opt-in on
# master, turned on here so the smvar hap.py evals build the genome-specific
# strata. Confined to this run clone.
sed -i -E '/^debug:/a\
\
# Enabled by run_full_pipeline_test.sh to exercise #59/#173 on smvar hap.py evals.\
genome_specific_strats: true' "$RESOURCES"
info "set genome_specific_strats: true in $RESOURCES"

# --- 4. commit the tested config in the clone -----------------------------
info "committing test config on branch run/${RUNID}"
git add "$OUTFILE" "$RESOURCES"
git commit --quiet -m "test(fullpipe): whole-genome regression run ${RUNID}

Exercises stvar_v5 truvari profile (#194), genome_specific_strats (#59/#173),
HG002Q100stvarv0.022 exclusion set (#188), and the symbolic-allele
self-discrepancy filter via a PAV stvar row (#192).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"

# --- 5. validate config ----------------------------------------------------
info "validating config (run_defrabb validate)"
./run_defrabb validate -a "$OUTFILE"

# --- 6. dry-run to confirm the DAG resolves --------------------------------
info "dry-run (snakemake -n) to confirm the DAG resolves with the new config"
DRYLOG="$(mktemp)"
if snakemake -n --use-conda --use-apptainer --cores "$CORES" \
        --config analyses="$OUTFILE" > "$DRYLOG" 2>&1; then
    tail -n 30 "$DRYLOG"
    info "dry-run OK. Feature paths present in the DAG:"
    for rule in self_discrep_filter_symbolic genome_specific_ run_pav; do
        if grep -q "$rule" "$DRYLOG"; then echo "   [x] $rule"; else echo "   [ ] $rule (not matched)"; fi
    done
    rm -f "$DRYLOG"
else
    echo "----- dry-run output (tail) -----" >&2
    tail -n 60 "$DRYLOG" >&2
    rm -f "$DRYLOG"
    die "dry-run failed; fix the config/code before launching the full run"
fi

# --- 7. write the post-run comparison helper (in the clone) ---------------
COMPARE_HELPER="compare_results_${RUNID}.sh"
PLOT_LINE=""
command -v Rscript >/dev/null 2>&1 && PLOT_LINE=" --plot plots_${RUNID}" || true
cat > "$COMPARE_HELPER" <<EOF
#!/usr/bin/env bash
# Aggregate this run's hap.py + truvari results into a comparison table.
# Run from this clone after the pipeline finishes.
set -euo pipefail
RESULTS_DIR="\${RESULTS_DIR:-results}"
python scripts/compare_evaluations.py \\
    --results-dir "\$RESULTS_DIR" \\
    --out "comparison_${RUNID}.tsv"${PLOT_LINE}
echo "wrote comparison_${RUNID}.tsv"
EOF
chmod +x "$COMPARE_HELPER"
info "wrote $COMPARE_HELPER (run it from the clone after the pipeline completes)"

# --- 8. launch or print the run command ------------------------------------
echo
echo "=========================================================================="
echo "Setup complete for run: $RUNID"
echo "  clone dir      : $CLONE_ABS  (branch run/${RUNID})"
echo "  analyses table : $OUTFILE  (committed)"
echo "  resources.yml  : + $NEW_EXCLUSION_SET, genome_specific_strats: true (committed)"
echo "  cleanup        : rm -rf '$CLONE_ABS'  (when done)"
echo "=========================================================================="
if [[ "$GO" -eq 1 ]]; then
    info "launching the full run (--go)"
    exec ./run_defrabb run -r "$RUNID" -a "$OUTFILE" -j "$CORES"
else
    echo "Setup-only (no --go). To launch the (compute-intensive) whole-genome run:"
    echo
    echo "    cd '$CLONE_ABS' && ./run_defrabb run -r ${RUNID} -j ${CORES}"
    echo
    echo "Then aggregate results (from the clone) with:  ./${COMPARE_HELPER}"
fi
