#!/usr/bin/env bash
# =============================================================================
#  run_all.sh - the whole pipeline, in order
# =============================================================================
#  Every stage is idempotent: a finished stage prints that it is skipping and
#  returns. So a run that died at stage 6 is resumed by running this again, and
#  --from restarts from a named stage after a code change.
#
#  Order and what each stage costs on the machine this was developed on, which
#  is 16 threads with 8 GB visible to the kernel. The measured figures are in
#  logs/*.resources.tsv and the README quotes them from there rather than from
#  this comment.
#
#      00_configure       seconds
#      01_install         three conda solves plus a 2.1 GB GNINA download
#      02_fetch_benchmarks  a 54 MB archive, then seconds
#      03_build_crossdock   around 700 RCSB API calls
#      04_prepare         643 receptors and 786 ligand preparations
#      05_define_boxes    fpocket on 393 receptors
#      06_dock            the long one. Vina measured 20 s and smina 25 s per
#                         complex on one core; GNINA with its CNN on CPU is the
#                         slow one. Multiply by the arms, methods and datasets in
#                         project.conf, then divide by JOBS.
#      07_score_poses     PoseBusters at roughly 0.5 s per pose
#      08_analyse         seconds
#      09_figures         seconds
#      10_report          one quarto render
#
#  Usage:
#      bash run_all.sh                             # everything, resuming
#      bash run_all.sh --from 06_dock              # restart at docking
#      bash run_all.sh --arm A2_genconf_refbox     # one arm only
#      bash run_all.sh --method vina --dataset astex
#      bash run_all.sh --jobs 12                   # override JOBS for this run
#      bash run_all.sh --smoke                     # 8 complexes, one method
#      bash run_all.sh --help
#
#  Options:
#      --from STAGE    start at this stage name (00_configure .. 10_report)
#      --only STAGE    run just this stage
#      --arm NAME      restrict docking and scoring to one arm
#      --method NAME   restrict to one method
#      --dataset NAME  restrict to one dataset
#      --jobs N        override JOBS from project.conf for this run
#      --smoke         8 complexes, one arm, one method: proves the wiring
#      --force         pass --force to every stage it applies to
#
#  Environment:
#      VGB_RUN_CONVERGENCE=1   also run the exhaustiveness grid in stage 6. It
#                              is the evidence behind EXHAUSTIVENESS in
#                              project.conf and costs about as much as one arm,
#                              so it is opt-in rather than part of every run.
#      -h, --help      this text
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="${ROOT}/scripts"

STAGES=(00_configure 01_install 02_fetch_benchmarks 03_build_crossdock 04_prepare
        05_define_boxes 06_dock 07_score_poses 08_analyse 09_figures 10_report)

FROM=""; ONLY=""; ARM=""; METHOD=""; DATASET=""; JOBS_OVERRIDE=""
SMOKE=0; FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)    FROM="$2"; shift 2 ;;
        --only)    ONLY="$2"; shift 2 ;;
        --arm)     ARM="$2"; shift 2 ;;
        --method)  METHOD="$2"; shift 2 ;;
        --dataset) DATASET="$2"; shift 2 ;;
        --jobs)    JOBS_OVERRIDE="$2"; shift 2 ;;
        --smoke)   SMOKE=1; shift ;;
        --force)   FORCE=1; shift ;;
        -h|--help) sed -n '2,45p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

stage_index() {
    local want="$1" i
    for i in "${!STAGES[@]}"; do
        [[ "${STAGES[$i]}" == "$want" ]] && { echo "$i"; return; }
    done
    echo "-1"
}
START=0
if [[ -n "$FROM" ]]; then
    START="$(stage_index "$FROM")"
    (( START >= 0 )) || { echo "[error] unknown stage: ${FROM}" >&2; exit 1; }
fi
should_run() {
    local name="$1"
    if [[ -n "$ONLY" ]]; then [[ "$ONLY" == "$name" ]] && return 0 || return 1; fi
    local i; i="$(stage_index "$name")"
    (( i >= START ))
}

banner() {
    echo
    echo "=============================================================="
    echo " $1"
    echo "=============================================================="
}

FORCE_FLAG=(); (( FORCE == 1 )) && FORCE_FLAG=(--force)

# -----------------------------------------------------------------------------
# 00 configure. Only run when project.conf is absent, because re-running it
# would silently overwrite a deliberately edited thread count or data directory
# on a resumed run. Pass --only 00_configure to force it.
# -----------------------------------------------------------------------------
if should_run 00_configure; then
    if [[ ! -f "${ROOT}/project.conf" || "$ONLY" == "00_configure" ]]; then
        banner "00_configure"
        bash "${SCRIPTS}/00_configure.sh" --yes
    else
        echo "[00_configure] project.conf exists; skipping (--only 00_configure to redo)."
    fi
fi

[[ -f "${ROOT}/project.conf" ]] || {
    echo "[error] project.conf missing. Run: bash scripts/00_configure.sh --yes" >&2; exit 1; }
# shellcheck source=/dev/null
source "${ROOT}/project.conf"
[[ -n "$JOBS_OVERRIDE" ]] && JOBS="$JOBS_OVERRIDE"

# Activate the analysis environment once for the python stages. The docking
# programs live in their own environments and 06_dock.sh calls them by absolute
# path, so this activation does not reach them.
activate_py() {
    set +u
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    conda activate "$CONDA_ENV_NAME"
    set -u
}

if should_run 01_install; then
    banner "01_install"
    bash "${SCRIPTS}/01_install.sh" "${FORCE_FLAG[@]}"
fi
if should_run 02_fetch_benchmarks; then
    banner "02_fetch_benchmarks"
    bash "${SCRIPTS}/02_fetch_benchmarks.sh" "${FORCE_FLAG[@]}"
fi

activate_py

SMOKE_ARGS=()
if (( SMOKE == 1 )); then
    SMOKE_ARGS=(--limit 8)
    [[ -n "$ARM" ]]    || ARM="A2_genconf_refbox"
    [[ -n "$METHOD" ]] || METHOD="vina"
    [[ -n "$DATASET" ]] || DATASET="astex"
    echo "[run_all] smoke run: ${ARM} / ${METHOD} / ${DATASET}, 8 complexes"
fi

if should_run 03_build_crossdock; then
    banner "03_build_crossdock_set"
    if [[ -s "${CONFIG_DIR}/crossdock_pairs.tsv" && $FORCE -eq 0 ]]; then
        echo "[03_build_crossdock] config/crossdock_pairs.tsv exists; skipping."
    else
        python "${SCRIPTS}/03_build_crossdock_set.py" --config "${ROOT}/project.conf" \
            "${FORCE_FLAG[@]}"
    fi
fi

if should_run 04_prepare; then
    banner "04_prepare"
    ARGS=(--config "${ROOT}/project.conf" --jobs "$JOBS" "${FORCE_FLAG[@]}")
    [[ -n "$DATASET" && "$DATASET" != "crossdock" ]] && ARGS+=(--dataset "$DATASET")
    (( SMOKE == 1 )) && ARGS+=("${SMOKE_ARGS[@]}")
    python "${SCRIPTS}/04_prepare.py" "${ARGS[@]}"
fi

if should_run 05_define_boxes; then
    banner "05_define_boxes"
    ARGS=(--config "${ROOT}/project.conf" --jobs "$JOBS" "${FORCE_FLAG[@]}")
    [[ -n "$DATASET" && "$DATASET" != "crossdock" ]] && ARGS+=(--dataset "$DATASET")
    (( SMOKE == 1 )) && ARGS+=("${SMOKE_ARGS[@]}")
    python "${SCRIPTS}/05_define_boxes.py" "${ARGS[@]}"
fi

# -----------------------------------------------------------------------------
# 06 docking. The cross product of arms, methods and datasets from project.conf,
# narrowed by whatever was passed on the command line. A0 is run once per
# dataset rather than once per method, because it involves no docking program and
# running it three times would put three identical rows in the table under
# different method names.
# -----------------------------------------------------------------------------
if should_run 06_dock; then
    banner "06_dock"
    IFS=',' read -r -a ARM_LIST    <<< "${ARM:-$ARMS}"
    IFS=',' read -r -a METHOD_LIST <<< "${METHOD:-$METHODS}"
    IFS=',' read -r -a DS_LIST     <<< "${DATASET:-$DATASETS}"
    for arm in "${ARM_LIST[@]}"; do
        [[ -n "$arm" ]] || continue
        if [[ "$arm" == "A4_crossdock" ]]; then
            for method in "${METHOD_LIST[@]}"; do
                bash "${SCRIPTS}/06_dock.sh" --arm "$arm" --method "$method" \
                    --jobs "$JOBS" "${FORCE_FLAG[@]}" "${SMOKE_ARGS[@]}"
            done
            continue
        fi
        for ds in "${DS_LIST[@]}"; do
            [[ -n "$ds" ]] || continue
            if [[ "$arm" == "A0_null" ]]; then
                bash "${SCRIPTS}/06_dock.sh" --arm "$arm" --method "${METHOD_LIST[0]}" \
                    --dataset "$ds" --jobs "$JOBS" "${FORCE_FLAG[@]}" "${SMOKE_ARGS[@]}"
                continue
            fi
            for method in "${METHOD_LIST[@]}"; do
                bash "${SCRIPTS}/06_dock.sh" --arm "$arm" --method "$method" \
                    --dataset "$ds" --jobs "$JOBS" "${FORCE_FLAG[@]}" "${SMOKE_ARGS[@]}"
            done
        done
    done

    # The convergence grid and the seed-variance experiment. Both run after the
    # arms so a failure in either cannot cost the main result. The convergence
    # grid is off unless asked for, because it is the evidence for a protocol
    # choice that has already been made and costs about as much as one arm.
    if [[ "${VGB_RUN_CONVERGENCE:-0}" == "1" ]]; then
        banner "06_dock --convergence"
        for method in "${METHOD_LIST[@]}"; do
            bash "${SCRIPTS}/06_dock.sh" --convergence --method "$method"                 --jobs "$JOBS" "${FORCE_FLAG[@]}" ||                 echo "[run_all] convergence grid for ${method} failed; continuing"
        done
    fi

    # The seed-variance experiment: one complex, five seeds, one method. Run
    # after the arms so that a failure here cannot cost the main result.
    banner "06_dock --seed-variance"
    for method in "${METHOD_LIST[@]}"; do
        bash "${SCRIPTS}/06_dock.sh" --seed-variance --method "$method" \
            --jobs 1 "${FORCE_FLAG[@]}" || \
            echo "[run_all] seed variance for ${method} failed; continuing"
    done
fi

if should_run 07_score_poses; then
    banner "07_score_poses"
    ARGS=(--config "${ROOT}/project.conf" --jobs "$JOBS" "${FORCE_FLAG[@]}")
    [[ -n "$ARM" ]] && ARGS+=(--arm "$ARM")
    [[ -n "$METHOD" ]] && ARGS+=(--method "$METHOD")
    python "${SCRIPTS}/07_score_poses.py" "${ARGS[@]}"
fi

if should_run 08_analyse; then
    banner "08_analyse"
    python "${SCRIPTS}/08_analyse.py" --config "${ROOT}/project.conf"
fi

if should_run 09_figures; then
    banner "09_figures"
    python "${SCRIPTS}/09_figures.py" --config "${ROOT}/project.conf"
fi

if should_run 10_report; then
    banner "10_report"
    if command -v quarto >/dev/null 2>&1; then
        mkdir -p "${RESULTS_DIR}/report"
        # Rendered from the repository root. The report finds project.conf by
        # walking up from the working directory rather than taking it as a
        # quarto parameter, because quarto renders from the qmd's own directory
        # in some versions and from the project root in others, and a report
        # that cannot find its own config depending on how it was invoked is not
        # reproducible. VGB_CONFIG overrides if the config ever moves.
        ( cd "$ROOT" && VGB_CONFIG="${ROOT}/project.conf" \
            quarto render "scripts/10_report.qmd" \
                --to html --output-dir "${RESULTS_DIR}/report" ) \
            || echo "[10_report] quarto render failed; the tables in results/ are unaffected"
    else
        echo "[10_report] quarto not found on PATH; skipping the report."
        echo "            Every number it would show is already in results/."
    fi
fi

banner "done"
echo "Tables:  ${RESULTS_DIR}"
echo "Figures: ${FIGURES_DIR}"
echo "Logs:    ${LOG_DIR}"
if [[ -s "${LOG_DIR}/summary.tsv" ]] || ls "${LOG_DIR}"/*.resources.tsv >/dev/null 2>&1; then
    echo
    echo "Per-stage elapsed time, peak resident memory and data growth:"
    { printf 'stage\telapsed_s\tpeak_rss_mb\tdata_growth_mb\tdata_total_mb\tjobs\tnote\n'
      for f in "${LOG_DIR}"/*.resources.tsv; do
          [[ -f "$f" ]] && tail -n +2 "$f"
      done
    } > "${LOG_DIR}/summary.tsv"
    column -t -s "$(printf '\t')" "${LOG_DIR}/summary.tsv" | sed 's/^/  /'
fi
