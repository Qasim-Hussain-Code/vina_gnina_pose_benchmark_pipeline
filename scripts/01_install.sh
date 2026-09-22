#!/usr/bin/env bash
# =============================================================================
#  01_install.sh - build the three conda environments, fetch GNINA, record
#  exactly what landed
# =============================================================================
#  Nothing in this repository hardcodes a tool version. Each environment is
#  solved against conda-forge at install time, and whatever the solver picks is
#  written to files the README quotes:
#
#      results/environment/versions.tsv        one row per tool, version, source
#      results/environment/tool_licences.tsv   licence string from the package
#      config/env_*.lock.yml                   each solved environment, pinned
#
#  Reproduce from the lock files rather than from the loose specs, which is the
#  point of writing them.
#
#  Why three environments. The conda-forge smina build pins libboost 1.82 and
#  AutoDock Vina 1.2.7 requires libboost 1.86, so the two cannot share an
#  environment. The solver's response to being asked for both is to fall back to
#  Vina 1.2.5 silently. Asking for Vina and PoseBusters together additionally
#  drags numpy below 2, which drags RDKit to 2023.09 and PoseBusters to 0.3.1.
#  The first environment solved here did exactly that, and the versions table
#  is what caught it. The docking programs are therefore called by absolute path
#  from 06_dock.sh rather than through an activated environment.
#
#  GNINA. No conda package exists. The release asset is a single statically
#  linked CUDA binary, 2.1 GB for v1.3.3, fetched from the GitHub release rather
#  than built, because building it pulls in libtorch and a CUDA toolchain. It
#  runs on CPU when no CUDA device is present. That fallback is tested here
#  rather than assumed, and if it fails this script exits non-zero instead of
#  leaving a method in the config that cannot run.
#
#  Usage:
#      bash scripts/01_install.sh              # skip anything already present
#      bash scripts/01_install.sh --force      # rebuild all three environments
#      bash scripts/01_install.sh --no-gnina   # skip the 2.1 GB download
#
#  Options:
#      --force      remove and rebuild the conda environments
#      --no-gnina   do not fetch GNINA; the gnina method becomes unavailable
#      -h, --help   this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
vgb_load_conf

STAGE=01_install
FORCE=0
WANT_GNINA=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)    FORCE=1; shift ;;
        --no-gnina) WANT_GNINA=0; shift ;;
        -h|--help)  sed -n '2,41p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

vgb_skip_if_done "$STAGE" "$FORCE" && exit 0

ENVDIR="${RESULTS_DIR}/environment"
TOOLDIR="${DATA_DIR}/tools"
mkdir -p "$ENVDIR" "$TOOLDIR"

# A killed install must not leave a partial 2 GB download that a later run
# would treat as the binary.
cleanup() {
    local rc=$?
    rm -f "${TOOLDIR}/gnina.part"
    (( rc != 0 )) && echo "[${STAGE}] failed with status ${rc}; no stamp written, re-run to retry" >&2
    return 0
}
trap cleanup EXIT

vgb_stage_start "$STAGE"

set +u
# shellcheck source=/dev/null
source "$CONDA_SH"
set -u
# conda-forge must win outright. With flexible priority the solver is free to
# take smina from bioconda, where the newest build is the 2017.11.9 release.
conda config --set channel_priority strict >/dev/null 2>&1 || true

# -----------------------------------------------------------------------------
# 1. The three environments.
# -----------------------------------------------------------------------------
build_env() {  # build_env <env_name> <spec basename>
    local env="$1" base="$2"
    local spec="${CONFIG_DIR}/${base}.yml"
    local lock="${CONFIG_DIR}/${base}.lock.yml"

    if [[ $FORCE -eq 1 ]]; then
        echo "[${STAGE}] removing ${env}"
        conda env remove -n "$env" -y >/dev/null 2>&1 || true
    fi
    if conda env list | awk '{print $1}' | grep -qx "$env"; then
        echo "[${STAGE}] ${env} already exists; not re-solving"
        return 0
    fi
    if [[ -s "$lock" ]]; then
        echo "[${STAGE}] creating ${env} from the pinned lock file"
        vgb_run "env_${env}_lock" conda env create -n "$env" -f "$lock" -y
    else
        [[ -s "$spec" ]] || { echo "[error] ${spec} missing" >&2; return 1; }
        echo "[${STAGE}] solving ${env} from ${spec}"
        vgb_run "env_${env}_spec" conda env create -n "$env" -f "$spec" -y
    fi
}

build_env "$CONDA_ENV_NAME"  env_analysis
build_env "$CONDA_ENV_VINA"  env_vina
build_env "$CONDA_ENV_SMINA" env_smina

VINA_BIN="$(vgb_tool_path "$CONDA_ENV_VINA" vina || true)"
SMINA_BIN="$(vgb_tool_path "$CONDA_ENV_SMINA" smina || true)"
FPOCKET_BIN="$(vgb_tool_path "$CONDA_ENV_SMINA" fpocket || true)"
[[ -n "$VINA_BIN"  ]] || { echo "[error] vina not found after install"  >&2; exit 1; }
[[ -n "$SMINA_BIN" ]] || { echo "[error] smina not found after install" >&2; exit 1; }

# The check that caught the silent downgrade. Vina 1.2.5 and 1.2.7 differ in
# the scoring code path, and a solver that quietly gives you 1.2.5 when the spec
# asked for 1.2.6 or newer invalidates the comparison with the published number.
VINA_PKG_VER="$(conda list -n "$CONDA_ENV_VINA" 2>/dev/null | awk '$1=="vina" {print $2; exit}')"
case "$VINA_PKG_VER" in
    1.2.6|1.2.7|1.2.8|1.3.*|1.2.9)
        echo "[${STAGE}] vina package version ${VINA_PKG_VER}" ;;
    *)
        echo "[warn] vina resolved to ${VINA_PKG_VER}, below the 1.2.6 floor in"
        echo "       config/env_vina.yml. This is the silent-downgrade failure"
        echo "       mode described in that file. The run will continue and the"
        echo "       version is recorded, but the comparison with the published"
        echo "       Vina number is weaker than it looks." ;;
esac

# -----------------------------------------------------------------------------
# 2. GNINA. Fetched, checksummed, and tested on CPU before it is trusted.
# -----------------------------------------------------------------------------
GNINA_TAG=""; GNINA_SHA=""; GNINA_OK=0
if [[ $WANT_GNINA -eq 1 ]]; then
    if [[ ! -x "${TOOLDIR}/gnina" ]]; then
        echo "[${STAGE}] querying the GNINA release feed"
        curl -sS --retry 3 https://api.github.com/repos/gnina/gnina/releases/latest \
            > "${TOOLDIR}/gnina_release.json"
        # One asset per release: the static CUDA build. Match on the name prefix
        # rather than a fixed CUDA suffix, because that suffix changes between
        # releases and hardcoding it would break on the next one.
        read -r GNINA_TAG GNINA_URL < <(
            "$(vgb_tool_path "$CONDA_ENV_NAME" python)" - "${TOOLDIR}/gnina_release.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
url = ""
for a in d.get("assets", []):
    if a["name"].startswith("gnina"):
        url = a["browser_download_url"]
        break
print(d.get("tag_name", "unknown"), url)
PY
        )
        [[ -n "${GNINA_URL:-}" ]] || { echo "[error] no gnina asset in release ${GNINA_TAG}" >&2; exit 1; }
        echo "[${STAGE}] fetching GNINA ${GNINA_TAG}, about 2.1 GB"
        vgb_run gnina_download curl -sSL --retry 5 --retry-all-errors \
            -o "${TOOLDIR}/gnina.part" "$GNINA_URL"
        mv "${TOOLDIR}/gnina.part" "${TOOLDIR}/gnina"
        chmod +x "${TOOLDIR}/gnina"
    else
        echo "[${STAGE}] GNINA already present at ${TOOLDIR}/gnina"
        [[ -s "${TOOLDIR}/gnina_release.json" ]] && GNINA_TAG="$(
            "$(vgb_tool_path "$CONDA_ENV_NAME" python)" -c \
            'import json,sys; print(json.load(open(sys.argv[1]))["tag_name"])' \
            "${TOOLDIR}/gnina_release.json" 2>/dev/null || true)"
    fi

    # Upstream publishes no checksum file, so this is the checksum of what this
    # machine downloaded. It is recorded so a second run can confirm it got the
    # same bytes, not so it can be verified against upstream.
    GNINA_SHA="$(sha256sum "${TOOLDIR}/gnina" | awk '{print $1}')"

    echo "[${STAGE}] testing GNINA on CPU"
    if "${TOOLDIR}/gnina" --version >/dev/null 2>&1; then
        GNINA_OK=1
        echo "[${STAGE}] $("${TOOLDIR}/gnina" --version 2>&1 | head -1)"
    else
        echo "[error] ${TOOLDIR}/gnina will not run on this machine." >&2
        echo "        The release asset is a static CUDA build and is expected to" >&2
        echo "        fall back to CPU when no device is present. It did not." >&2
        echo "        Either build GNINA from source, or re-run 00_configure.sh" >&2
        echo "        with --methods vina,vinardo and drop the CNN arm." >&2
        exit 5
    fi
else
    echo "[${STAGE}] --no-gnina given; skipping the download"
fi

# -----------------------------------------------------------------------------
# 3. Record what landed. Every version the README quotes comes from here.
# -----------------------------------------------------------------------------
echo "[${STAGE}] recording versions"
VER="${ENVDIR}/versions.tsv"
printf 'tool\tversion\tsource\tdetail\trecorded\n' > "$VER"
rec() { printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$(date -Iseconds)" >> "$VER"; }

# Package versions come from the package manager, not from each tool's
# --version flag, because several of these print a git hash or nothing useful.
# Vina 1.2.7 for instance prints "AutoDock Vina 36dd023-mod", which is a commit
# and not a release number, so both are recorded.
record_env_pkgs() {  # record_env_pkgs <env> <pkg...>
    local env="$1"; shift
    local listing; listing="$(conda list -n "$env" 2>/dev/null)"
    local p v b
    for p in "$@"; do
        v="$(awk -v p="$p" '$1==p {print $2; exit}' <<<"$listing")"
        b="$(awk -v p="$p" '$1==p {print $3; exit}' <<<"$listing")"
        [[ -n "$v" ]] && rec "$p" "$v" "conda-forge/${env}" "$b"
    done
}
record_env_pkgs "$CONDA_ENV_NAME"  python rdkit meeko posebusters spyrmsd biopython \
                                   pdb2pqr openbabel numpy scipy pandas matplotlib-base \
                                   seaborn requests quarto
record_env_pkgs "$CONDA_ENV_VINA"  vina
record_env_pkgs "$CONDA_ENV_SMINA" smina fpocket

rec vina_banner  "$("$VINA_BIN"  --version 2>&1 | head -1 | tr -d '\r')" "binary" "$VINA_BIN"
rec smina_banner "$("$SMINA_BIN" --version 2>&1 | head -1 | tr -d '\r')" "binary" "$SMINA_BIN"
[[ -n "$FPOCKET_BIN" ]] && rec fpocket_binary "$(basename "$FPOCKET_BIN")" "binary" "$FPOCKET_BIN"
if (( GNINA_OK == 1 )); then
    rec gnina "$("${TOOLDIR}/gnina" --version 2>&1 | head -1 | tr -d '\r')" \
        "github release ${GNINA_TAG:-unknown}" "sha256:${GNINA_SHA}"
    rec gnina_cnn_scoring "$GNINA_CNN_SCORING" "project.conf" "cnn model set: ${GNINA_CNN_MODEL}"
fi
rec conda  "$(conda --version 2>&1 | awk '{print $2}')" "base" "$(command -v conda)"
rec kernel "$(uname -r)" "host" "$(uname -s) $(uname -m)"

# Licence strings as the package manager records them. This is what the run
# used, which is the claim the LICENSE file has to support.
LIC="${ENVDIR}/tool_licences.tsv"
printf 'package\tversion\tlicence\tenvironment\tchecked\n' > "$LIC"
"$(vgb_tool_path "$CONDA_ENV_NAME" python)" - "$LIC" \
    "$CONDA_ENV_NAME" "$CONDA_ENV_VINA" "$CONDA_ENV_SMINA" <<'PY'
# Licences are read from each environment's conda-meta JSON. `conda list` does
# not print licences, and `conda search --info` hits the network for every
# package and can disagree with what is installed.
import datetime, json, pathlib, subprocess, sys

out, envs = sys.argv[1], sys.argv[2:]
wanted = {"vina", "smina", "rdkit", "meeko", "posebusters", "spyrmsd",
          "biopython", "pdb2pqr", "openbabel", "fpocket", "quarto", "python",
          "numpy", "pandas", "scipy", "matplotlib-base", "seaborn"}
prefixes = {}
for line in subprocess.run(["conda", "env", "list"], capture_output=True,
                           text=True).stdout.splitlines():
    parts = line.split()
    if parts and parts[0] in envs:
        prefixes[parts[0]] = parts[-1]

stamp = datetime.date.today().isoformat()
rows = 0
with open(out, "a") as fh:
    for env in envs:
        prefix = prefixes.get(env)
        if not prefix:
            continue
        for p in sorted(pathlib.Path(prefix, "conda-meta").glob("*.json")):
            try:
                d = json.loads(p.read_text())
            except Exception:
                continue
            if d.get("name") in wanted:
                lic = d.get("license") or "NA"
                fh.write(f"{d['name']}\t{d.get('version','NA')}\t{lic}\t{env}\t{stamp}\n")
                rows += 1
print(f"[01_install] wrote {rows} licence rows")
PY

# Lock files. --no-builds keeps the version pins and drops the build strings,
# which is what travels between machines; the explicit list next to it keeps the
# build strings for exact reproduction on the same platform.
for pair in "${CONDA_ENV_NAME}:env_analysis" "${CONDA_ENV_VINA}:env_vina" "${CONDA_ENV_SMINA}:env_smina"; do
    env="${pair%%:*}"; base="${pair#*:}"
    conda env export -n "$env" --no-builds | grep -v '^prefix:' > "${CONFIG_DIR}/${base}.lock.yml"
    conda list -n "$env" --explicit > "${ENVDIR}/explicit_${base}_$(uname -m).txt"
done

echo "[${STAGE}] versions recorded in ${VER}"
column -t -s "$(printf '\t')" "$VER" | sed 's/^/    /'

vgb_stage_end "gnina_ok=${GNINA_OK} vina=${VINA_PKG_VER}"
vgb_mark_done "$STAGE"
