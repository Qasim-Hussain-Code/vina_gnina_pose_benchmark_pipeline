#!/usr/bin/env bash
# =============================================================================
#  06_dock.sh - run one arm with one method
# =============================================================================
#  One invocation covers one arm, one method, one dataset. run_all.sh loops over
#  the combinations named in project.conf. Keeping the loop outside means a
#  failed combination can be re-run on its own without redoing the rest.
#
#  Arms
#      A0_null            no search at all. The generated conformer is dropped
#                         into the box at a random orientation about the box
#                         centre. This is the floor every success rate should be
#                         read against, and without it a reader cannot tell
#                         whether the method did anything.
#      A1_refconf_refbox  crystal ligand coordinates in, box on the crystal
#                         ligand. The standard protocol and the flattering one.
#      A2_genconf_refbox  conformer generated from SMILES, box still on the
#                         crystal ligand. The primary result.
#      A3_genconf_detbox  generated conformer, box from fpocket on the receptor
#                         alone. No crystal information anywhere.
#      A4_crossdock       generated conformer docked into a different structure
#                         of the same protein, box on that structure's own
#                         ligand.
#
#  Methods
#      vina      AutoDock Vina, Vina scoring
#      vinardo   smina with --scoring vinardo. smina is a fork of Vina 1.1.2, so
#                this differs from the vina arm in both scoring function and
#                search. The README says so rather than letting a reader assume
#                the search is held constant.
#      gnina     GNINA with the CNN applied as configured in project.conf.
#                GNINA_CNN_SCORING is recorded next to every number because
#                rescore and refinement are not the same experiment.
#
#  Disk. After the first ten complexes the stage measures what they actually
#  cost, projects the total for the remaining complexes in this arm plus the
#  arms still to come, and refuses to continue if that will not fit inside
#  DISK_GB. A run that dies two thirds of the way through has wasted more than
#  it saved. Intermediate PDBQT poses are deleted as soon as the pose and the
#  score have been written into the per-run record, so the footprint is the
#  compressed pose archive and nothing else.
#
#  Timing. JOBS from project.conf defaults to 1. Serial is the default because
#  the per-complex wall clock in the README is a measurement, and N docking
#  processes competing for memory bandwidth inflate it by an amount that depends
#  on the machine. Raise --jobs for a production run; the timing figures quoted
#  in the README come from the serial subset, which is marked in the results
#  table by its jobs column.
#
#  Usage:
#      bash scripts/06_dock.sh --arm A2_genconf_refbox --method vina --dataset posebusters
#      bash scripts/06_dock.sh --arm A4_crossdock --method gnina --jobs 8
#      bash scripts/06_dock.sh --seed-variance --method vina
#
#  Options:
#      --arm NAME        one of the five arms above
#      --method NAME     vina, vinardo or gnina
#      --dataset NAME    posebusters or astex; ignored for A4_crossdock
#      --jobs N          override JOBS from project.conf
#      --limit N         stop after N complexes (smoke test)
#      --seed-variance   re-dock one complex once per seed in SEED_REPLICATES
#                        and write the spread. Ignores --arm.
#      --force           ignore the stage stamp and redo
#      -h, --help        this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
vgb_load_conf

ARM=""; METHOD=""; DATASET=""; LIMIT=0; FORCE=0; SEED_VAR=0
JOBS_OVERRIDE=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --arm)            ARM="$2"; shift 2 ;;
        --method)         METHOD="$2"; shift 2 ;;
        --dataset)        DATASET="$2"; shift 2 ;;
        --jobs)           JOBS_OVERRIDE="$2"; shift 2 ;;
        --limit)          LIMIT="$2"; shift 2 ;;
        --seed-variance)  SEED_VAR=1; shift ;;
        --force)          FORCE=1; shift ;;
        -h|--help)        sed -n '2,66p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$JOBS_OVERRIDE" ]] && JOBS="$JOBS_OVERRIDE"
[[ -n "$METHOD" ]] || { echo "[error] --method is required" >&2; exit 1; }
if [[ $SEED_VAR -eq 0 ]]; then
    [[ -n "$ARM" ]] || { echo "[error] --arm is required" >&2; exit 1; }
    if [[ "$ARM" == "A4_crossdock" ]]; then
        DATASET="crossdock"
    else
        [[ -n "$DATASET" ]] || { echo "[error] --dataset is required for ${ARM}" >&2; exit 1; }
    fi
    STAGE="06_dock_${ARM}_${METHOD}_${DATASET}"
else
    STAGE="06_dock_seedvariance_${METHOD}"
fi

vgb_skip_if_done "$STAGE" "$FORCE" && exit 0

vgb_activate
PY="$(command -v python)"

# Each docking program lives in its own conda environment because Vina 1.2.7
# and smina cannot share one (libboost 1.86 against 1.82). They are called by
# absolute path rather than by activating and deactivating inside the loop.
VINA_BIN="$(vgb_tool_path "$CONDA_ENV_VINA" vina || true)"
SMINA_BIN="$(vgb_tool_path "$CONDA_ENV_SMINA" smina || true)"
# GNINA_BIN is set in project.conf, not next to SMINA_BIN above.
# shellcheck disable=SC2153
case "$METHOD" in
    vina)    [[ -n "$VINA_BIN"  ]] || { echo "[error] vina not installed; run 01_install.sh" >&2; exit 3; } ;;
    vinardo) [[ -n "$SMINA_BIN" ]] || { echo "[error] smina not installed; run 01_install.sh" >&2; exit 3; } ;;
    gnina)   [[ -x "$GNINA_BIN" ]] || { echo "[error] ${GNINA_BIN} not present; run 01_install.sh" >&2; exit 3; } ;;
    *) echo "[error] unknown method: ${METHOD}" >&2; exit 1 ;;
esac

vgb_stage_start "$STAGE"

# A killed run leaves per-complex temporary directories behind. They are under
# the arm's own work directory so the trap can clear them without touching
# anything another arm is using.
WORK="${DATA_DIR}/work/${STAGE}"
mkdir -p "$WORK"
cleanup() {
    local rc=$?
    rm -rf "$WORK"
    (( rc != 0 )) && echo "[${STAGE}] exited with status ${rc}; no stamp written" >&2
    return 0
}
trap cleanup EXIT INT TERM

ARGS=(--config "${REPO_DIR}/project.conf" --method "$METHOD" --jobs "$JOBS"
      --work-dir "$WORK" --stage "$STAGE")
[[ -n "$VINA_BIN"  ]] && ARGS+=(--vina-bin "$VINA_BIN")
[[ -n "$SMINA_BIN" ]] && ARGS+=(--smina-bin "$SMINA_BIN")
[[ -x "$GNINA_BIN" ]] && ARGS+=(--gnina-bin "$GNINA_BIN")
(( LIMIT > 0 )) && ARGS+=(--limit "$LIMIT")
(( FORCE == 1 )) && ARGS+=(--force)
if (( SEED_VAR == 1 )); then
    ARGS+=(--seed-variance)
else
    ARGS+=(--arm "$ARM" --dataset "$DATASET")
fi

# The python side does the per-complex work: it owns the disk gate, the process
# pool, the pose extraction and the deletion of intermediates. Bash owns the
# orchestration, the stage stamp and the resource log, which is the split used
# throughout this repository.
# SCRIPTS_DIR is set in project.conf, which vgb_load_conf sourced above.
# shellcheck disable=SC2153
vgb_run "dock_${METHOD}" "$PY" "${SCRIPTS_DIR}/lib_dock.py" "${ARGS[@]}"

vgb_stage_end "arm=${ARM:-seedvariance} method=${METHOD} dataset=${DATASET:-na} jobs=${JOBS}"
vgb_mark_done "$STAGE"
