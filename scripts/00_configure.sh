#!/usr/bin/env bash
# =============================================================================
#  00_configure.sh - measure the machine, take the operator's answers, project
#  the disk footprint of the requested arms, write project.conf
# =============================================================================
#  Every later script sources project.conf rather than hardcoding a thread
#  count, a seed, a box size or a path. Re-run this script to change any of
#  them; it overwrites project.conf and nothing else.
#
#  The disk projection here is coarse and deliberately so. It multiplies a
#  per-complex figure by the number of complex-runs the requested arms imply.
#  Those per-complex figures are defaults until something has been measured on
#  this machine. The real gate is in 06_dock.sh, which measures the first ten
#  complexes of an arm and refuses to continue if the projection from those ten
#  does not fit. A projection from zero measurements cannot be trusted and this
#  script says so rather than pretending otherwise.
#
#  Usage:
#      bash scripts/00_configure.sh --threads 16 --ram 16 --disk 40 --yes
#      bash scripts/00_configure.sh --yes                  # detected defaults
#      bash scripts/00_configure.sh --data-dir /scratch/bench --yes
#
#  Options:
#      --threads N     CPU threads the machine may use (default: detected)
#      --ram GB        RAM in gigabytes (default: detected)
#      --disk GB       disk budget for data/ in gigabytes (default: detected
#                      free space on the data filesystem, minus 5 GB margin)
#      --exhaustiveness N
#                      Vina, smina and GNINA search effort (default: 32). This
#                      is NOT the programs' default of 8, and the reason is in
#                      results/convergence.tsv rather than in an opinion: on the
#                      25 Angstrom box this pipeline uses, exhaustiveness 8
#                      leaves the search unconverged badly enough that the
#                      random seed decides the verdict. One complex measured
#                      4.57 Angstrom on one seed and 0.33 on another with
#                      everything else fixed. Run
#                      scripts/06_dock.sh --convergence to reproduce the grid
#                      this was chosen from.
#      --jobs N        concurrent docking processes (default: 1). Serial is the
#                      default because the timing distribution in the README is
#                      a measurement, and sixteen docking processes competing
#                      for memory bandwidth on one socket inflate per-complex
#                      wall clock by an amount that depends on the machine.
#                      Memory is not the reason: one Vina process on a 25 A box
#                      peaks near 200 MB. Raise --jobs for a production run and
#                      read timings off the serial subset instead.
#      --data-dir DIR  where structures, prepared inputs and poses go. Default
#                      is <repo>/data. Point it at a larger filesystem if the
#                      repository sits on a small one, which is the case on the
#                      machine this was developed on.
#      --arms LIST     comma-separated subset of A0_null,A1_refconf_refbox,
#                      A2_genconf_refbox,A3_genconf_detbox,A4_crossdock
#      --methods LIST  comma-separated subset of vina,vinardo,gnina
#      --datasets LIST comma-separated subset of posebusters,astex
#      --yes, -y       accept every default without asking
#      -h, --help      this text
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
CONF="${REPO_DIR}/project.conf"
CONF_TMP="${CONF}.tmp.$$"

# A killed run must not leave a half-written project.conf behind, because the
# next script would source it and inherit a truncated variable list.
cleanup() { rm -f "$CONF_TMP"; }
trap cleanup EXIT INT TERM

THREADS=""; RAM_GB=""; DISK_GB=""; JOBS=""; DATA_DIR=""; ASSUME_YES=0
EXHAUSTIVENESS=32
ARMS="A0_null,A1_refconf_refbox,A2_genconf_refbox,A3_genconf_detbox,A4_crossdock"
METHODS="vina,vinardo,gnina"
DATASETS="posebusters,astex"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)  THREADS="$2";  shift 2 ;;
        --ram)      RAM_GB="$2";   shift 2 ;;
        --disk)     DISK_GB="$2";  shift 2 ;;
        --jobs)     JOBS="$2";     shift 2 ;;
        --exhaustiveness) EXHAUSTIVENESS="$2"; shift 2 ;;
        --data-dir) DATA_DIR="$2"; shift 2 ;;
        --arms)     ARMS="$2";     shift 2 ;;
        --methods)  METHODS="$2";  shift 2 ;;
        --datasets) DATASETS="$2"; shift 2 ;;
        --yes|-y)   ASSUME_YES=1;  shift ;;
        -h|--help)  sed -n '2,46p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

# -----------------------------------------------------------------------------
# 1. Measure the machine. Three platforms because this was developed under WSL2
#    on Windows and is meant to run on a Linux workstation without edits.
# -----------------------------------------------------------------------------
detect_threads() {
    nproc 2>/dev/null && return
    getconf _NPROCESSORS_ONLN 2>/dev/null && return
    sysctl -n hw.ncpu 2>/dev/null && return
    echo 2
}
detect_ram_gb() {
    # MemTotal, not MemAvailable. Under WSL2 the kernel is handed a fraction of
    # the host's RAM, so this reports what the analysis can actually use rather
    # than what the laptop has on the motherboard, and the two differ here.
    if [[ -r /proc/meminfo ]]; then
        awk '/^MemTotal/ {printf "%d", ($2 + 524288) / 1048576}' /proc/meminfo; return
    fi
    if sysctl -n hw.memsize >/dev/null 2>&1; then
        echo $(( $(sysctl -n hw.memsize) / 1073741824 )); return
    fi
    echo 4
}
detect_free_gb() {  # detect_free_gb <path>
    df -Pk "$1" 2>/dev/null | awk 'NR==2 {printf "%d", $4 / 1048576}' || echo 0
}

DET_THREADS="$(detect_threads)"
DET_RAM_GB="$(detect_ram_gb)"

[[ -n "$DATA_DIR" ]] || DATA_DIR="${REPO_DIR}/data"
mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"
DET_FREE_GB="$(detect_free_gb "$DATA_DIR")"
REPO_FREE_GB="$(detect_free_gb "$REPO_DIR")"

echo "=============================================================="
echo " vina/vinardo/gnina pose benchmark - machine configuration"
echo "=============================================================="
echo "  CPU threads visible     : ${DET_THREADS}"
echo "  RAM visible to kernel   : ${DET_RAM_GB} GB"
echo "  free disk at data dir   : ${DET_FREE_GB} GB  (${DATA_DIR})"
echo "  free disk at repo       : ${REPO_FREE_GB} GB  (${REPO_DIR})"
echo

ask() {  # ask <prompt> <default> <varname>
    local prompt="$1" default="$2" __var="$3" reply
    if [[ $ASSUME_YES -eq 1 || ! -t 0 ]]; then
        printf -v "$__var" '%s' "$default"; return
    fi
    read -r -p "  ${prompt} [${default}]: " reply || reply=""
    printf -v "$__var" '%s' "${reply:-$default}"
}

DEFAULT_DISK=$(( DET_FREE_GB - 5 )); (( DEFAULT_DISK < 1 )) && DEFAULT_DISK=1
[[ -n "$THREADS" ]] || ask "CPU threads to use"         "$DET_THREADS"  THREADS
[[ -n "$RAM_GB"  ]] || ask "RAM to use (GB)"            "$DET_RAM_GB"   RAM_GB
[[ -n "$DISK_GB" ]] || ask "disk budget for data/ (GB)" "$DEFAULT_DISK" DISK_GB
[[ -n "$JOBS"    ]] || ask "concurrent docking jobs"    "1"             JOBS

# -----------------------------------------------------------------------------
# 2. Validate. A wrong number here costs hours of docking later.
# -----------------------------------------------------------------------------
for pair in "THREADS:$THREADS" "RAM_GB:$RAM_GB" "DISK_GB:$DISK_GB" "JOBS:$JOBS"             "EXHAUSTIVENESS:$EXHAUSTIVENESS"; do
    name="${pair%%:*}"; val="${pair#*:}"
    [[ "$val" =~ ^[0-9]+$ ]] || { echo "[error] ${name} must be a non-negative integer, got '${val}'" >&2; exit 1; }
done
(( THREADS >= 1 )) || { echo "[error] need at least 1 thread" >&2; exit 1; }
(( JOBS >= 1 ))    || JOBS=1
if (( JOBS > THREADS )); then
    echo "[warn] --jobs ${JOBS} exceeds ${THREADS} threads; using ${THREADS}."
    JOBS=$THREADS
fi
# One docking process peaks near 200 MB for Vina and smina and above 1 GB for
# GNINA running its CNN on CPU. 2 GB per job is the budget assumed here; the
# measured peaks per stage land in logs/*.resources.tsv.
MAX_JOBS_BY_RAM=$(( RAM_GB / 2 )); (( MAX_JOBS_BY_RAM < 1 )) && MAX_JOBS_BY_RAM=1
if (( JOBS > MAX_JOBS_BY_RAM )); then
    echo "[warn] ${JOBS} jobs at 2 GB each would exceed ${RAM_GB} GB; using ${MAX_JOBS_BY_RAM}."
    JOBS=$MAX_JOBS_BY_RAM
fi
if (( RAM_GB < 6 )); then
    echo "[warn] ${RAM_GB} GB visible. GNINA's CNN on CPU has been seen above 1 GB;"
    echo "       if the gnina method is in --methods, expect it to be the tight one."
fi

# -----------------------------------------------------------------------------
# 3. Project the disk footprint of the requested arms.
#    Per-complex-run figures, in kilobytes so the arithmetic stays integer.
#    These are the defaults used before anything has been measured here;
#    06_dock.sh replaces them with its own measurement after ten complexes and
#    re-runs the same arithmetic.
#      prepared receptor PDBQT   ~1200 kB, deleted once the arm finishes
#      prepared ligand PDBQT     ~8 kB
#      top-N poses, gzipped SDF  ~25 kB per run
#      per-run log and score     ~4 kB per run
#    Structures are downloaded once and shared across arms and methods.
# -----------------------------------------------------------------------------
KB_PER_RUN_DEFAULT=29        # poses plus per-run log and score
KB_PER_COMPLEX_PREP=1210     # receptor and ligand PDBQT, transient
MB_STRUCTURES_DEFAULT=1400   # benchmark archive unpacked plus cross-docking entries
# The GNINA release is a single statically linked CUDA binary and it is large:
# 2.1 GB for v1.3.3. It lands under DATA_DIR/tools and dominates the footprint
# of a small run, which is why it is counted separately rather than folded into
# the per-complex figure.
MB_GNINA_BINARY=2100

n_items() { tr ',' '\n' <<<"$1" | grep -c . ; }
N_METHODS="$(n_items "$METHODS")"
N_ARMS="$(n_items "$ARMS")"
N_DATASETS="$(n_items "$DATASETS")"

# 308 PoseBusters plus 85 Astex is 393 complexes if both datasets are asked
# for. The cross-docking pair count is not known until
# 03_build_crossdock_set.py has run, so assume 200 here and let the manifest
# correct it once it exists. Two hundred is a guess, not a measurement.
N_COMPLEXES=0
case ",${DATASETS}," in *,posebusters,*) N_COMPLEXES=$(( N_COMPLEXES + 308 )) ;; esac
case ",${DATASETS}," in *,astex,*)       N_COMPLEXES=$(( N_COMPLEXES + 85 ))  ;; esac
N_CROSS_PAIRS=200
CROSS_MANIFEST="${REPO_DIR}/config/crossdock_pairs.tsv"
if [[ -s "$CROSS_MANIFEST" ]]; then
    N_CROSS_PAIRS="$(grep -cv '^#' "$CROSS_MANIFEST" || true)"
    N_CROSS_PAIRS=$(( N_CROSS_PAIRS - 1 ))   # header row
    (( N_CROSS_PAIRS < 0 )) && N_CROSS_PAIRS=0
fi

# Arms A1 to A3 run on the chosen datasets; A4 runs on the cross-docking pairs.
SELF_ARMS=0
case ",${ARMS}," in *,A1_refconf_refbox,*) SELF_ARMS=$(( SELF_ARMS + 1 )) ;; esac
case ",${ARMS}," in *,A2_genconf_refbox,*) SELF_ARMS=$(( SELF_ARMS + 1 )) ;; esac
case ",${ARMS}," in *,A3_genconf_detbox,*) SELF_ARMS=$(( SELF_ARMS + 1 )) ;; esac
CROSS_ARM=0
case ",${ARMS}," in *,A4_crossdock,*) CROSS_ARM=1 ;; esac
# A0 is the null floor. It places the generated conformer in the box without
# searching, so it costs one pose per complex and no docking time at all.
NULL_ARM=0
case ",${ARMS}," in *,A0_null,*) NULL_ARM=1 ;; esac

N_RUNS=$(( N_COMPLEXES * SELF_ARMS * N_METHODS \
         + N_CROSS_PAIRS * CROSS_ARM * N_METHODS \
         + N_COMPLEXES * NULL_ARM ))
MB_TOOLS=0
case ",${METHODS}," in *,gnina,*) MB_TOOLS=$(( MB_TOOLS + MB_GNINA_BINARY )) ;; esac
PROJ_MB=$(( N_RUNS * KB_PER_RUN_DEFAULT / 1024 \
          + N_COMPLEXES * KB_PER_COMPLEX_PREP / 1024 \
          + MB_STRUCTURES_DEFAULT + MB_TOOLS ))
PROJ_GB=$(( PROJ_MB / 1024 + 1 ))

echo "  arms requested           : ${N_ARMS} (${ARMS})"
echo "  methods requested        : ${N_METHODS} (${METHODS})"
echo "  datasets requested       : ${N_DATASETS} (${DATASETS})"
echo "  complexes                : ${N_COMPLEXES} self-docking, ${N_CROSS_PAIRS} cross-docking pairs"
echo "  docking runs implied     : ${N_RUNS}"
echo "  projected data footprint : ${PROJ_GB} GB, of which ${MB_TOOLS} MB is tool binaries"
echo "                             (per-complex defaults, not yet measured here)"
echo

if (( PROJ_GB > DISK_GB )); then
    SHORT=$(( PROJ_GB - DISK_GB ))
    echo "[error] the requested arms project to ${PROJ_GB} GB and the budget is ${DISK_GB} GB." >&2
    echo "        shortfall ${SHORT} GB. Raise --disk, point --data-dir at a larger" >&2
    echo "        filesystem, or drop an arm or a method. Refusing to write a" >&2
    echo "        project.conf that would die two thirds of the way through." >&2
    exit 4
fi
if (( PROJ_GB > DET_FREE_GB )); then
    SHORT=$(( PROJ_GB - DET_FREE_GB ))
    echo "[error] projection is ${PROJ_GB} GB but only ${DET_FREE_GB} GB is free at" >&2
    echo "        ${DATA_DIR}. Shortfall ${SHORT} GB." >&2
    exit 4
fi
if (( REPO_FREE_GB < 2 )); then
    echo "[warn] only ${REPO_FREE_GB} GB free on the repository filesystem. Tracked"
    echo "       output (tables, figures, report) stays under 50 MB, so this is"
    echo "       survivable, but the conda environment must not land here."
fi

# -----------------------------------------------------------------------------
# 4. Locate conda. 01_install.sh puts the Python side in an environment;
#    record where conda lives so later scripts can activate without a login
#    shell. GNINA is a downloaded binary rather than a conda package.
# -----------------------------------------------------------------------------
find_conda_sh() {
    local c root
    for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)"; do
        [[ -n "$c" && -x "$c" ]] || continue
        root="$(dirname "$(dirname "$c")")"
        [[ -r "${root}/etc/profile.d/conda.sh" ]] && { echo "${root}/etc/profile.d/conda.sh"; return; }
    done
    for c in "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/anaconda3" /opt/conda; do
        [[ -r "${c}/etc/profile.d/conda.sh" ]] && { echo "${c}/etc/profile.d/conda.sh"; return; }
    done
    echo ""
}
CONDA_SH="$(find_conda_sh)"
[[ -n "$CONDA_SH" ]] || echo "[warn] conda not found. 01_install.sh needs it; install miniforge first."

cat > "$CONF_TMP" <<CONF_EOF
# =============================================================================
#  project.conf - written by scripts/00_configure.sh on $(date -Iseconds)
#  Re-run scripts/00_configure.sh to change these. Do not edit by hand: every
#  later stage sources this file and a hand edit will not survive the next run.
# =============================================================================

# ---- hardware, as measured on this machine --------------------------------
THREADS=${THREADS}
RAM_GB=${RAM_GB}
DISK_GB=${DISK_GB}                 # budget for DATA_DIR, enforced in 06_dock.sh
JOBS=${JOBS}                       # concurrent docking processes; 1 is serial

# ---- what to run ----------------------------------------------------------
ARMS="${ARMS}"
METHODS="${METHODS}"
DATASETS="${DATASETS}"

# ---- reproducibility ------------------------------------------------------
# One seed for everything stochastic: ETKDG conformer generation, the Vina,
# smina and GNINA Monte Carlo searches, the random orientations in the null
# arm. 06_dock.sh also re-docks one complex five times with the seeds below,
# because run-to-run spread in top-1 RMSD is part of the measurement and is
# almost never reported.
SEED=20260922
SEED_REPLICATES="20260922 20260923 20260924 20260925 20260926"

# ---- docking protocol -----------------------------------------------------
# 25 A cube on the geometric centre of the crystal ligand's heavy atoms. This
# is the box Buttenschoen et al. used for Vina in the PoseBusters paper and it
# is copied here so the Vina number in results/ is comparable with the 58 per
# cent they report. It is a large box for Vina and it is arbitrary in the sense
# that no pocket is 25 A wide; it was chosen for comparability, not because it
# is right.
BOX_SIZE=25
# Search effort. Not the programs' default of 8, and this is the one place where
# this pipeline departs from the protocol it is comparing against.
#
# On a 25 Angstrom cube, exhaustiveness 8 does not converge. Complex 1HQ2_PH2,
# which has one rotatable bond and whose generated conformer is already within
# 0.21 Angstrom of the crystal conformer, returned a top-1 RMSD of 4.57 Angstrom
# on seed 20260922 and 0.33 Angstrom on seed 7, with the receptor, the ligand,
# the box and the exhaustiveness all identical. At exhaustiveness 64 the same
# complex and seed gave 0.37 Angstrom. A search whose verdict is decided by its
# seed cannot measure what an arm comparison needs it to measure.
#
# The level here was chosen from the grid in results/convergence.tsv, which is
# complexes crossed with exhaustiveness and seed, reproducible with
# scripts/06_dock.sh --convergence. The choice is the cheapest level at which
# the seeds stop disagreeing about the 2 Angstrom verdict, not the level that
# makes the headline look best: it was fixed before any success rate was read.
EXHAUSTIVENESS=${EXHAUSTIVENESS}
# Poses kept per run. Vina's default num_modes is 9; 5 is enough for the top-5
# metric and cuts the pose archive by nearly half. Anything past rank 5 is
# never read by 07_score_poses.py, so writing it would be waste.
TOP_N_POSES=5
# GNINA CNN settings. --cnn_scoring=rescore runs the CNN on the final poses the
# empirical score produced; --cnn_scoring=refinement runs it inside the search
# and is far slower on CPU. This machine has no CUDA device, so rescore is the
# only setting that finishes, and the difference between the two is large
# enough that the README states which was used. CNN_MODEL is recorded rather
# than left implicit because GNINA's default ensemble has changed between
# releases.
GNINA_CNN_SCORING=rescore
GNINA_CNN_MODEL=default
# Success thresholds in Angstrom. 2.0 is the number the field quotes and it is
# arbitrary; 1.0 is reported next to it so a reader can see whether the ranking
# of the methods survives the stricter cut.
RMSD_PASS=2.0
RMSD_STRICT=1.0

# ---- cross-docking set construction --------------------------------------
# Cutoffs for 03_build_crossdock_set.py. Both are conventional rather than
# derived: 2.5 A is the usual line for a structure good enough to dock into,
# and 95 per cent sequence identity is the usual line for calling two chains
# the same protein. Tightening either shrinks the set.
XDOCK_MAX_RESOLUTION=2.5
XDOCK_MIN_SEQ_IDENTITY=0.95
XDOCK_MAX_PARTNERS=1               # one partner per query, the highest resolution

# ---- pocket detection for the box-free arm --------------------------------
# A3 needs a box that never saw the crystal ligand. fpocket is the usual
# choice; where it is unavailable the fallback is the geometric detector in
# 05_define_boxes.py, which is named in the README. Whichever ran is recorded
# per complex in results/boxes/box_provenance.tsv.
POCKET_METHOD=auto

# ---- paths ----------------------------------------------------------------
REPO_DIR="${REPO_DIR}"
CONFIG_DIR="${REPO_DIR}/config"
SCRIPTS_DIR="${REPO_DIR}/scripts"
DATA_DIR="${DATA_DIR}"
RESULTS_DIR="${REPO_DIR}/results"
FIGURES_DIR="${REPO_DIR}/figures"
LOG_DIR="${REPO_DIR}/logs"

# ---- tools ----------------------------------------------------------------
CONDA_SH="${CONDA_SH}"
# Three environments, not one. The conda-forge smina build pins libboost 1.82
# and AutoDock Vina 1.2.7 requires libboost 1.86, so they cannot share an
# environment; asking for both gets you Vina 1.2.5 with no warning. Asking for
# Vina and PoseBusters together additionally drags numpy below 2, which drags
# RDKit to 2023.09 and PoseBusters to 0.3.1. config/env_analysis.yml has the
# full account. 06_dock.sh calls each docking program by absolute path.
CONDA_ENV_NAME=vgb_bench           # preparation, measurement, figures, report
CONDA_ENV_VINA=vgb_vina            # AutoDock Vina only
CONDA_ENV_SMINA=vgb_smina          # smina and fpocket only
# GNINA's release binary is described as static and is not: it is dynamically
# linked against CUDA 12.8 and cuDNN 9 and exits 127 before printing its version
# on a machine that has neither, GPU or no GPU. This environment holds nothing
# but those libraries. It is never activated; its lib directory goes on
# LD_LIBRARY_PATH for GNINA calls only. config/env_gnina_runtime.yml explains
# which sonames and why the CUDA version is pinned.
CONDA_ENV_GNINA_RT=vgb_gnina
# GNINA has no conda package. 01_install.sh fetches the release binary and
# records its version and sha256 in results/environment/versions.tsv.
GNINA_BIN="${DATA_DIR}/tools/gnina"
CONF_EOF

mv "$CONF_TMP" "$CONF"
echo "  wrote ${CONF}"
echo "  threads=${THREADS} ram=${RAM_GB}GB disk_budget=${DISK_GB}GB jobs=${JOBS}"
echo "  data_dir=${DATA_DIR}"
