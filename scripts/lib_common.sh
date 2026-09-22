#!/usr/bin/env bash
# =============================================================================
#  lib_common.sh - sourced by every stage. Not executable on its own.
# =============================================================================
#  Four jobs:
#    1. find and source project.conf, failing loudly if 00_configure.sh has
#       not run, because a stage that silently uses shell defaults for THREADS
#       and SEED produces results nobody can reproduce
#    2. activate the conda environment without needing a login shell
#    3. record elapsed time, peak resident memory and peak disk per stage into
#       logs/, so the README can quote measurements
#    4. stage stamps, so re-running a finished stage skips and says so
# =============================================================================

# shellcheck disable=SC2148
# (no shebang effect: this file is sourced, and the shebang above is for
# editors and shellcheck rather than for execution)

VGB_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VGB_REPO_DIR="$(cd "${VGB_LIB_DIR}/.." && pwd)"

# ---- 1. project.conf --------------------------------------------------------
vgb_load_conf() {
    local conf="${VGB_REPO_DIR}/project.conf"
    if [[ ! -r "$conf" ]]; then
        echo "[error] project.conf not found. Run scripts/00_configure.sh first." >&2
        echo "        Every stage reads its thread count, seed, box size and" >&2
        echo "        paths from there; nothing here has a silent default." >&2
        return 3
    fi
    # shellcheck source=/dev/null
    source "$conf"
    local v
    for v in THREADS RAM_GB DISK_GB JOBS SEED DATA_DIR RESULTS_DIR LOG_DIR \
             BOX_SIZE EXHAUSTIVENESS TOP_N_POSES RMSD_PASS RMSD_STRICT; do
        [[ -n "${!v:-}" ]] || { echo "[error] ${v} missing from project.conf; re-run 00_configure.sh" >&2; return 3; }
    done
    mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$FIGURES_DIR" "$LOG_DIR" "$CONFIG_DIR"
}

# ---- 2. conda --------------------------------------------------------------
# Activating inside a non-interactive script needs conda.sh sourced by hand.
# `conda run` is the documented alternative and it was tried first; it buffers
# stdout until the child exits, which makes a four-hour docking stage look
# hung, so the environment is activated in-process instead.
vgb_activate() {
    [[ -n "${CONDA_SH:-}" && -r "$CONDA_SH" ]] || {
        echo "[error] CONDA_SH='${CONDA_SH:-}' is not readable; re-run 00_configure.sh" >&2; return 3; }
    # conda.sh trips -u on unset variables of its own
    set +u
    # shellcheck source=/dev/null
    source "$CONDA_SH"
    conda activate "$CONDA_ENV_NAME" || { set -u; echo "[error] conda env '${CONDA_ENV_NAME}' not found. Run scripts/01_install.sh." >&2; return 3; }
    set -u
}

# vgb_tool_path <env_name> <binary>
# Absolute path to a binary in a sibling conda environment. 06_dock.sh needs
# this because Vina and smina live in different environments (libboost 1.82
# against 1.86) and activating one would deactivate the other mid-loop.
vgb_tool_path() {
    local env="$1" bin="$2" root
    root="$(dirname "$(dirname "$CONDA_SH")")"      # <root>/etc/profile.d -> <root>
    root="$(dirname "$root")"                        # -> conda root
    local p="${root}/envs/${env}/bin/${bin}"
    [[ -x "$p" ]] && { echo "$p"; return 0; }
    # Fall back to PATH, so a system install of the tool still works.
    command -v "$bin" 2>/dev/null && return 0
    echo ""; return 1
}

# ---- 3. measurement --------------------------------------------------------
# Peak RSS comes from GNU time's %M, which reports kilobytes. /usr/bin/time is
# a separate package from the shell builtin `time` and the builtin cannot
# report memory at all, so the absence of the binary is reported rather than
# quietly skipped: a resources table with NA in it is honest, one with a
# guessed number is not.
vgb_time_bin() {
    local t
    for t in /usr/bin/time /bin/time "$(command -v gtime 2>/dev/null || true)"; do
        [[ -n "$t" && -x "$t" ]] && { echo "$t"; return; }
    done
    echo ""
}

# vgb_dir_kb <dir> - apparent size in kilobytes, 0 if absent
vgb_dir_kb() {
    [[ -d "$1" ]] || { echo 0; return; }
    du -sk "$1" 2>/dev/null | awk '{print $1}' || echo 0
}

# vgb_stage_start <stage>
# Records the wall clock and the data directory size so vgb_stage_end can
# subtract them.
vgb_stage_start() {
    VGB_STAGE="$1"
    VGB_T0="$(date +%s)"
    VGB_DISK0="$(vgb_dir_kb "$DATA_DIR")"
    VGB_PEAK_RSS_KB=0
    echo "[${VGB_STAGE}] started $(date -Iseconds)"
}

# vgb_run <label> <command...>
# Runs a command under GNU time, appends its elapsed seconds and peak RSS to
# the stage's per-command log, and keeps the largest RSS seen so the stage's
# peak is the peak of its children rather than of the shell.
vgb_run() {
    local label="$1"; shift
    local tb; tb="$(vgb_time_bin)"
    local tf; tf="$(mktemp)"
    local rc=0
    if [[ -n "$tb" ]]; then
        "$tb" -f '%e %M' -o "$tf" "$@" || rc=$?
        local elapsed peak
        read -r elapsed peak < "$tf" || { elapsed=NA; peak=NA; }
        printf '%s\t%s\t%s\t%s\n' "$VGB_STAGE" "$label" "$elapsed" "$peak" \
            >> "${LOG_DIR}/${VGB_STAGE}.commands.tsv"
        if [[ "$peak" =~ ^[0-9]+$ ]] && (( peak > VGB_PEAK_RSS_KB )); then
            VGB_PEAK_RSS_KB="$peak"
        fi
    else
        local t0 t1; t0="$(date +%s)"
        "$@" || rc=$?
        t1="$(date +%s)"
        printf '%s\t%s\t%s\tNA\n' "$VGB_STAGE" "$label" "$(( t1 - t0 ))" \
            >> "${LOG_DIR}/${VGB_STAGE}.commands.tsv"
    fi
    rm -f "$tf"
    return $rc
}

# vgb_stage_end [note]
vgb_stage_end() {
    local note="${1:-ok}"
    local t1 disk1 elapsed disk_delta peak
    t1="$(date +%s)"
    disk1="$(vgb_dir_kb "$DATA_DIR")"
    elapsed=$(( t1 - VGB_T0 ))
    disk_delta=$(( disk1 - VGB_DISK0 ))
    (( disk_delta < 0 )) && disk_delta=0
    peak="$VGB_PEAK_RSS_KB"
    [[ "$peak" == "0" ]] && peak=NA
    local out="${LOG_DIR}/${VGB_STAGE}.resources.tsv"
    printf 'stage\telapsed_s\tpeak_rss_mb\tdata_dir_growth_mb\tdata_dir_total_mb\tjobs\tnote\n' > "$out"
    printf '%s\t%d\t%s\t%d\t%d\t%s\t%s\n' \
        "$VGB_STAGE" "$elapsed" \
        "$( [[ "$peak" == NA ]] && echo NA || echo $(( peak / 1024 )) )" \
        "$(( disk_delta / 1024 ))" "$(( disk1 / 1024 ))" "${JOBS}" "$note" >> "$out"
    echo "[${VGB_STAGE}] finished in ${elapsed} s, data dir grew $(( disk_delta / 1024 )) MB to $(( disk1 / 1024 )) MB"
}

# ---- 4. idempotency --------------------------------------------------------
# A stamp file per stage. Stages check it and skip, which makes run_all.sh
# restartable after a failure without redoing hours of docking. --force on any
# stage removes its own stamp.
vgb_stamp() { echo "${LOG_DIR}/${1}.done"; }
vgb_is_done() { [[ -f "$(vgb_stamp "$1")" ]]; }
vgb_mark_done() { date -Iseconds > "$(vgb_stamp "$1")"; }
vgb_skip_if_done() {  # vgb_skip_if_done <stage> <force 0|1>
    local stage="$1" force="${2:-0}"
    if [[ "$force" == "1" ]]; then rm -f "$(vgb_stamp "$stage")"; return 1; fi
    if vgb_is_done "$stage"; then
        echo "[${stage}] already completed on $(cat "$(vgb_stamp "$stage")"); skipping (--force to redo)."
        return 0
    fi
    return 1
}

# ---- misc ------------------------------------------------------------------
vgb_need() {
    command -v "$1" >/dev/null 2>&1 || { echo "[error] ${1} not found on PATH" >&2; return 3; }
}
# Comma-separated list to newline-separated, dropping blanks.
vgb_split() { tr ',' '\n' <<<"$1" | sed '/^$/d'; }
# Is <needle> in the comma-separated <haystack>?
vgb_in_list() { case ",${2}," in *",${1},"*) return 0 ;; *) return 1 ;; esac; }
