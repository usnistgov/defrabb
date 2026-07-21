#!/usr/bin/env bash
###############################################################################
# release_test.sh
#
# Pre-release validation: clone repo, run full-genome pipeline test in tmux.
# This script prepares for a release by running the comprehensive whole-genome
# regression test that CI cannot cover (time/compute constraints).
#
# Features:
#   - Uses today's date for RUNID (e.g., 20260721_v0.022_fulltest-HG002v1.1)
#   - Clones repo into dated directory
#   - Launches run_defrabb in a detachable tmux session
#   - Monitors progress until completion or error
#   - Generates comparison results automatically
#
# IMPORTANT: Commit all changes to master BEFORE running this script.
# The test runs against committed HEAD.
#
# Usage:
#   conda activate snakemake   # or ~/miniforge3/envs/snakemake
#   scripts/release_test.sh [--cores N] [--dest-dir DIR] [--version v0.022]
#
#   --cores      cores/jobs for the run (default: 88 = all available)
#   --dest-dir   parent directory for the clone (default: ../runs)
#   --version    version string for RUNID (default: v0.022)
#   --detach     create tmux session but don't attach (returns immediately)
#   --help       show this help
#
# The script will:
#   1. Clone repo to DEST_DIR/<RUNID>
#   2. Generate analyses table with fulltest config
#   3. Validate and dry-run
#   4. Launch in tmux session named "defrabb-<RUNID>"
#   5. Wait for completion (unless --detach)
#   6. Generate comparison results
#
# To check status: tmux attach -t defrabb-<RUNID>
# To kill run: tmux kill-session -t defrabb-<RUNID>
###############################################################################
set -euo pipefail

# Defaults
VERSION="v0.022"
DEST_DIR="../runs"
CORES=$(nproc 2>/dev/null || echo 88)
DETACH=0
TEMPLATE_REL="config/analyses_20260323_v0.020_HG002v1.1-CI.tsv"
NEW_EXCLUSION_SET="HG002Q100stvarv0.022pav"

die() { echo "ERROR: $*" >&2; exit 1; }
info() { echo ">> $*"; }
usage() { sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; }

# --- parse args ------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores|-j) CORES="${2:?--cores needs a value}"; shift 2 ;;
        --dest-dir) DEST_DIR="${2:?--dest-dir needs a value}"; shift 2 ;;
        --version)  VERSION="${2:?--version needs a value}"; shift 2 ;;
        --detach)   DETACH=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        -*) die "unknown option: $1 (try --help)" ;;
        *) die "unexpected argument: $1 (try --help)" ;;
    esac
done

RUNID="$(date +%Y%m%d)_${VERSION}_fulltest-HG002v1.1"
SESSION_NAME="defrabb-${RUNID}"

[[ "$RUNID" =~ ^[0-9]{8}_v[0-9]+\.[0-9]{3}_[A-Za-z0-9.-]+$ ]] \
    || die "RUNID '$RUNID' does not match YYYYMMDD_v#.###_brief-id"

# --- preflight -------------------------------------------------------------
info "Pre-flight checks"
[[ -f Snakefile && -x run_defrabb ]] \
    || die "run this from the repository root (Snakefile / run_defrabb not found)"
[[ -f "$TEMPLATE_REL" ]] || die "template not found: $TEMPLATE_REL"
command -v snakemake >/dev/null 2>&1 \
    || die "snakemake not on PATH; activate the env first: conda activate snakemake"
command -v git >/dev/null 2>&1 || die "git not found on PATH"
command -v tmux >/dev/null 2>&1 || die "tmux not found on PATH (required for detachable session)"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "not inside a git work tree"

REPO_ROOT="$(pwd)"

# Warn about uncommitted changes
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
    echo "WARNING: uncommitted tracked changes detected. Test runs against committed HEAD."
    git status --short --untracked-files=no | sed 's/^/         /'
    echo "         Commit changes first to include them in the test."
    read -p "Continue anyway? [y/N] " -n 1 -r
    echo
    [[ $REPLY =~ ^[Yy]$ ]] || die "aborted by user"
fi

mkdir -p "$DEST_DIR"
CLONE_DIR="${DEST_DIR%/}/${RUNID}"
[[ -e "$CLONE_DIR" ]] && die "clone target already exists: $CLONE_DIR (remove it or wait until tomorrow)"

# Check if tmux session already exists
if tmux has-session -t "$SESSION_NAME" 2>/dev/null; then
    die "tmux session '$SESSION_NAME' already exists. Attach with: tmux attach -t $SESSION_NAME"
fi

# --- 1. clone the repo -----------------------------------------------------
info "Cloning $REPO_ROOT -> $CLONE_DIR"
git clone --quiet "$REPO_ROOT" "$CLONE_DIR"
CLONE_ABS="$(cd "$CLONE_DIR" && pwd)"
cd "$CLONE_ABS"
git checkout --quiet -b "run/${RUNID}"

OUTFILE="config/analyses_${RUNID}.tsv"
RESOURCES="config/resources.yml"

# --- 2. generate analyses table --------------------------------------------
info "Generating $OUTFILE from $TEMPLATE_REL"
awk 'BEGIN{FS=OFS="\t"}
NR==1 {print; next}
{
    if ($10 == "stvar") { $4 = "stvar_v5"; $13 = "HG002Q100stvarv0.022" }
    print
}' "$TEMPLATE_REL" > "$OUTFILE"

# Append PAV stvar row
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "GRCh38_HG002_T~T2TQ100v1.1_Q~v5.0q-stvar-pav_TR~v5.0q" \
    "GRCh38_HG002-T2TQ100v1.1-pav_stvar-excluded" \
    "truvari" "stvar_v5" "v5.0q-stvar" "FALSE" "TRUE" "TRUE" \
    "GRCh38_HG002-T2TQ100v1.1-pav" "stvar" "norm_xy_trf_sv_lcr" "exclude" \
    "$NEW_EXCLUSION_SET" "HG2-T2TQ100-V1.1" "GRCh38" "pav" "giab" "default" \
    >> "$OUTFILE"

info "Wrote $(($(wc -l < "$OUTFILE") - 1)) analysis rows"

# --- 3. enable genome_specific_strats --------------------------------------
grep -qE "^  ${NEW_EXCLUSION_SET}:" "$RESOURCES" \
    || die "PAV exclusion set ${NEW_EXCLUSION_SET} not found in $RESOURCES"

sed -i -E '/^debug:/a\
\
# Enabled by release_test.sh to exercise #59/#173 on smvar hap.py evals.\
genome_specific_strats: true' "$RESOURCES"
info "Set genome_specific_strats: true in $RESOURCES"

# --- 4. commit test config -------------------------------------------------
info "Committing test config on branch run/${RUNID}"
git add "$OUTFILE" "$RESOURCES"
git commit --quiet -m "test(release): whole-genome validation run ${RUNID}

Exercises stvar_v5 truvari profile (#194), genome_specific_strats (#59/#173),
HG002Q100stvarv0.022 exclusion set (#188), PAV-specific exclusion,
TRF contig filter (#197), and SV-aware self-discrepancy (#193).

Release validation for ${VERSION}.

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"

# --- 5. validate -----------------------------------------------------------
info "Validating config"
./run_defrabb validate -a "$OUTFILE" || die "config validation failed"

# --- 6. dry-run ------------------------------------------------------------
info "Dry-run to confirm DAG resolves"
DRYLOG="$(mktemp)"
if snakemake -n --use-conda --use-apptainer --cores "$CORES" \
        --config analyses="$OUTFILE" > "$DRYLOG" 2>&1; then
    tail -n 20 "$DRYLOG"
    info "Dry-run OK. Feature coverage:"
    for rule in self_discrep_filter_symbolic sv_self_discrep_truvari genome_specific_ run_pav merge_trfanno; do
        if grep -q "$rule" "$DRYLOG"; then
            echo "   ✓ $rule"
        else
            echo "   ✗ $rule (not in DAG - may be optional)"
        fi
    done
    rm -f "$DRYLOG"
else
    tail -n 60 "$DRYLOG" >&2
    rm -f "$DRYLOG"
    die "dry-run failed"
fi

# --- 7. create comparison helper -------------------------------------------
COMPARE_HELPER="compare_results_${RUNID}.sh"
cat > "$COMPARE_HELPER" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
RESULTS_DIR="${RESULTS_DIR:-results}"
python scripts/compare_evaluations.py \
    --results-dir "$RESULTS_DIR" \
    --regions \
    --out "comparison_${RUNID}.tsv"
echo "Wrote comparison_${RUNID}.tsv"
if command -v Rscript >/dev/null 2>&1; then
    mkdir -p plots_${RUNID}
    python scripts/compare_evaluations.py \
        --results-dir "$RESULTS_DIR" \
        --regions \
        --out "comparison_${RUNID}.tsv" \
        --plot plots_${RUNID}
    echo "Wrote plots to plots_${RUNID}/"
fi
EOF
sed -i "s/\${RUNID}/$RUNID/g" "$COMPARE_HELPER"
chmod +x "$COMPARE_HELPER"

# --- 8. create tmux session ------------------------------------------------
info "Creating tmux session: $SESSION_NAME"

# Create tmux session with run command
tmux new-session -d -s "$SESSION_NAME" -c "$CLONE_ABS"

# Set up panes: main pane runs pipeline, split pane monitors logs
tmux send-keys -t "$SESSION_NAME" "# DeFrABB Release Test: $RUNID" C-m
tmux send-keys -t "$SESSION_NAME" "# Clone: $CLONE_ABS" C-m
tmux send-keys -t "$SESSION_NAME" "# Start: $(date)" C-m
tmux send-keys -t "$SESSION_NAME" "" C-m

# Activate conda env (try common paths)
if [[ -n "${CONDA_PREFIX:-}" ]]; then
    tmux send-keys -t "$SESSION_NAME" "conda activate snakemake || source ~/miniforge3/bin/activate snakemake" C-m
fi

# Launch pipeline
RUN_CMD="./run_defrabb run -r ${RUNID} -a ${OUTFILE} -j ${CORES}"
tmux send-keys -t "$SESSION_NAME" "echo 'Launching: $RUN_CMD'" C-m
tmux send-keys -t "$SESSION_NAME" "$RUN_CMD" C-m

# After pipeline completes, run comparison
tmux send-keys -t "$SESSION_NAME" "echo ''" C-m
tmux send-keys -t "$SESSION_NAME" "echo 'Pipeline complete. Generating comparison results...'" C-m
tmux send-keys -t "$SESSION_NAME" "./$COMPARE_HELPER" C-m
tmux send-keys -t "$SESSION_NAME" "echo ''" C-m
tmux send-keys -t "$SESSION_NAME" "echo 'Release test complete: $(date)'" C-m
tmux send-keys -t "$SESSION_NAME" "echo 'Review results in: $CLONE_ABS'" C-m

# Split window for log monitoring
tmux split-window -v -t "$SESSION_NAME" -c "$CLONE_ABS"
tmux send-keys -t "$SESSION_NAME:0.1" "tail -f run.log 2>/dev/null || echo 'Waiting for run.log...'; sleep 5; tail -f run.log" C-m

# Resize panes (main 70%, log 30%)
tmux resize-pane -t "$SESSION_NAME:0.0" -y 70%

echo
echo "=========================================================================="
echo "Release test setup complete"
echo "=========================================================================="
echo "  Run ID        : $RUNID"
echo "  Clone dir     : $CLONE_ABS"
echo "  Tmux session  : $SESSION_NAME"
echo "  Cores         : $CORES"
echo "  Analyses      : $(($(wc -l < "$OUTFILE") - 1)) rows"
echo "=========================================================================="
echo
echo "Commands:"
echo "  Attach        : tmux attach -t $SESSION_NAME"
echo "  Detach        : Ctrl-b d"
echo "  Kill session  : tmux kill-session -t $SESSION_NAME"
echo "  List sessions : tmux ls"
echo
echo "After completion, results will be in:"
echo "  $CLONE_ABS/results/"
echo "  $CLONE_ABS/comparison_${RUNID}.tsv"
echo "=========================================================================="

if [[ $DETACH -eq 0 ]]; then
    info "Attaching to tmux session (Ctrl-b d to detach)"
    sleep 2
    tmux attach -t "$SESSION_NAME"
else
    info "Session running in background (--detach)"
    echo "Attach with: tmux attach -t $SESSION_NAME"
fi
