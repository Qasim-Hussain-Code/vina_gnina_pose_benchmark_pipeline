#!/usr/bin/env bash
# =============================================================================
#  02_fetch_benchmarks.sh - PoseBusters Benchmark set and Astex Diverse set,
#  downloaded, checksum-verified, and turned into two manifests
# =============================================================================
#  What this stage produces:
#
#      config/dataset_posebusters.tsv   308 rows, one per complex
#      config/dataset_astex.tsv          85 rows, one per complex
#      config/dataset_excluded.tsv       what was dropped here and why
#      results/data_provenance.tsv       URLs, checksums, counts, dates
#
#  The coordinates come from the PoseBusters authors' Zenodo archive rather than
#  from the PDB, so that the structures are byte-identical to the ones other
#  people report numbers on. Reassembling them from the PDB would mean
#  reproducing their solvent stripping and cofactor handling, and any difference
#  there moves the success rate without being visible in the results table.
#
#  The archive holds 428 PoseBusters complexes, which is the preprint version.
#  The published paper reports on 308, the reduction having been made during
#  peer review to drop complexes whose ligand sits against a crystal contact.
#  That 308-member list is not in the archive. config/sources.tsv says where it
#  is taken from and this script verifies, every run, that it has 308 entries
#  and is a strict subset of the 428. If either check fails the stage stops,
#  because a benchmark set that quietly changed size is not a benchmark set.
#
#  One thing the archive gives away for free and this pipeline deliberately does
#  not use: each folder ships a PDB_CCD_ligand_start_conf.sdf, a conformer the
#  authors generated with ETKDGv3 and a UFF relaxation. 04_prepare.py generates
#  its own instead, from the ligand SMILES, with the seed from project.conf
#  written into the output. The authors' file has no seed I can record, and an
#  unrecorded seed in the arm whose whole purpose is to remove crystal
#  information would undercut the comparison.
#
#  Usage:
#      bash scripts/02_fetch_benchmarks.sh
#      bash scripts/02_fetch_benchmarks.sh --force     # re-download and re-verify
#
#  Options:
#      --force      ignore the stage stamp, re-download, re-verify, rewrite
#      -h, --help   this text
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source-path=SCRIPTDIR
# shellcheck source=lib_common.sh
source "${SCRIPT_DIR}/lib_common.sh"
vgb_load_conf

STAGE=02_fetch_benchmarks
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --force)   FORCE=1; shift ;;
        -h|--help) sed -n '2,44p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "[error] unknown option: $1" >&2; exit 1 ;;
    esac
done

vgb_skip_if_done "$STAGE" "$FORCE" && exit 0

RAW="${DATA_DIR}/raw"
BENCH="${DATA_DIR}/benchmarks"
mkdir -p "$RAW" "$BENCH"

# Partial downloads must not be mistaken for finished ones on the next run.
cleanup() {
    local rc=$?
    rm -f "${RAW}"/*.part
    (( rc != 0 )) && echo "[${STAGE}] failed with status ${rc}; no stamp written" >&2
    return 0
}
trap cleanup EXIT

vgb_stage_start "$STAGE"
vgb_activate
vgb_need curl
vgb_need unzip

PY="$(command -v python)"

# -----------------------------------------------------------------------------
# 1. The Zenodo archive. The checksum is not hardcoded: the Zenodo API is asked
#    for it and the downloaded file is compared against what the API says. A
#    checksum I generated myself and pasted in would verify only that the file
#    has not changed since I looked at it, which is not the same claim.
# -----------------------------------------------------------------------------
ZEN_RECORD="$(awk -F'\t' '$1=="posebusters_structures" {print $3}' "${CONFIG_DIR}/sources.tsv")"
[[ -n "$ZEN_RECORD" ]] || { echo "[error] posebusters_structures row missing from config/sources.tsv" >&2; exit 1; }
echo "[${STAGE}] asking Zenodo about record ${ZEN_RECORD}"
curl -sS --retry 3 "https://zenodo.org/api/records/${ZEN_RECORD}" > "${RAW}/zenodo_record.json"

ZIP_NAME=posebusters_paper_data.zip
read -r ZIP_URL ZIP_SIZE ZIP_MD5 ZEN_DATE ZEN_LICENCE < <("$PY" - \
    "${RAW}/zenodo_record.json" "$ZIP_NAME" <<'PYEOF'
import json, sys
rec = json.load(open(sys.argv[1]))
want = sys.argv[2]
for f in rec.get("files", []):
    if f["key"] == want:
        # Zenodo reports "md5:<hex>"; strip the algorithm prefix.
        chk = f.get("checksum", "")
        md5 = chk.split(":", 1)[1] if ":" in chk else chk
        print(f["links"]["self"], f["size"], md5,
              rec["metadata"]["publication_date"],
              (rec["metadata"].get("license") or {}).get("id", "NA"))
        break
else:
    sys.exit("file %s not in record" % want)
PYEOF
)
echo "[${STAGE}] ${ZIP_NAME}: ${ZIP_SIZE} bytes, upstream md5 ${ZIP_MD5}"

ZIP="${RAW}/${ZIP_NAME}"
verify_md5() { [[ -s "$1" ]] && [[ "$(md5sum "$1" | awk '{print $1}')" == "$2" ]]; }

if ! verify_md5 "$ZIP" "$ZIP_MD5"; then
    echo "[${STAGE}] downloading ${ZIP_NAME}"
    # --retry-all-errors matters here: the first attempt at this file truncated
    # at 21 MB of 53.7 MB and exited zero, which is how a silent partial
    # download gets into a results table. The md5 check below is the real guard.
    rm -f "${ZIP}.part"
    vgb_run zenodo_zip curl -sSL --retry 5 --retry-all-errors --max-time 1800 \
        -o "${ZIP}.part" "$ZIP_URL"
    mv "${ZIP}.part" "$ZIP"
    if ! verify_md5 "$ZIP" "$ZIP_MD5"; then
        echo "[error] md5 mismatch on ${ZIP_NAME}." >&2
        echo "        expected ${ZIP_MD5}" >&2
        echo "        got      $(md5sum "$ZIP" | awk '{print $1}')" >&2
        echo "        got $(stat -c %s "$ZIP") bytes, expected ${ZIP_SIZE}" >&2
        exit 6
    fi
else
    echo "[${STAGE}] ${ZIP_NAME} already present and md5 matches"
fi
ZIP_SHA256="$(sha256sum "$ZIP" | awk '{print $1}')"

# -----------------------------------------------------------------------------
# 2. Extract. unzip -n rather than -o so a re-run does not rewrite 2,500 files
#    for nothing, which on a spinning disk is the difference between four
#    seconds and four minutes.
# -----------------------------------------------------------------------------
if [[ ! -d "${BENCH}/posebusters_benchmark_set" || $FORCE -eq 1 ]]; then
    echo "[${STAGE}] extracting"
    vgb_run unzip_archive unzip -q -n "$ZIP" -d "$BENCH"
else
    echo "[${STAGE}] archive already extracted"
fi
N_PB_DIRS="$(find "${BENCH}/posebusters_benchmark_set" -mindepth 1 -maxdepth 1 -type d | wc -l)"
N_AX_DIRS="$(find "${BENCH}/astex_diverse_set"          -mindepth 1 -maxdepth 1 -type d | wc -l)"
echo "[${STAGE}] extracted ${N_PB_DIRS} PoseBusters folders, ${N_AX_DIRS} Astex folders"

# -----------------------------------------------------------------------------
# 3. The 308-member list, pinned to a commit and verified against the archive.
# -----------------------------------------------------------------------------
IDS_LOCATOR="$(awk -F'\t' '$1=="posebusters_308_ids" {print $3}' "${CONFIG_DIR}/sources.tsv")"
IDS_COMMIT="$(awk -F'\t'  '$1=="posebusters_308_ids" {print $4}' "${CONFIG_DIR}/sources.tsv")"
IDS_REPO="${IDS_LOCATOR%%:*}"
IDS_PATH="${IDS_LOCATOR#*:}"
IDS_URL="https://raw.githubusercontent.com/${IDS_REPO}/${IDS_COMMIT}/${IDS_PATH}"
IDS_FILE="${RAW}/posebusters_308_ids.txt"
if [[ ! -s "$IDS_FILE" || $FORCE -eq 1 ]]; then
    echo "[${STAGE}] fetching the 308 list at ${IDS_COMMIT:0:10}"
    vgb_run ids_download curl -sSL --retry 3 -o "$IDS_FILE" "$IDS_URL"
fi
IDS_SHA256="$(sha256sum "$IDS_FILE" | awk '{print $1}')"

# -----------------------------------------------------------------------------
# 4. Manifests. Written by python because the subset check, the per-complex file
#    existence check and the TSV writing all want the same data in hand.
# -----------------------------------------------------------------------------
echo "[${STAGE}] writing dataset manifests"
# SCRIPTS_DIR is set in project.conf, which vgb_load_conf sourced above.
# shellcheck disable=SC2153
vgb_run write_manifests "$PY" "${SCRIPTS_DIR}/lib_manifest.py" \
    --bench-dir "$BENCH" \
    --ids-308 "$IDS_FILE" \
    --config-dir "$CONFIG_DIR" \
    --results-dir "$RESULTS_DIR"

# -----------------------------------------------------------------------------
# 5. Provenance. Everything a reader needs to get the same bytes.
# -----------------------------------------------------------------------------
PROV="${RESULTS_DIR}/data_provenance.tsv"
{
    printf 'item\tvalue\n'
    printf 'zenodo_record\t%s\n' "$ZEN_RECORD"
    printf 'zenodo_publication_date\t%s\n' "$ZEN_DATE"
    printf 'zenodo_licence\t%s\n' "$ZEN_LICENCE"
    printf 'archive_file\t%s\n' "$ZIP_NAME"
    printf 'archive_url\t%s\n' "$ZIP_URL"
    printf 'archive_bytes\t%s\n' "$ZIP_SIZE"
    printf 'archive_md5_upstream\t%s\n' "$ZIP_MD5"
    printf 'archive_sha256_local\t%s\n' "$ZIP_SHA256"
    printf 'archive_downloaded\t%s\n' "$(date -Iseconds)"
    printf 'posebusters_folders_in_archive\t%s\n' "$N_PB_DIRS"
    printf 'astex_folders_in_archive\t%s\n' "$N_AX_DIRS"
    printf 'ids_308_url\t%s\n' "$IDS_URL"
    printf 'ids_308_commit\t%s\n' "$IDS_COMMIT"
    printf 'ids_308_sha256\t%s\n' "$IDS_SHA256"
} > "$PROV"
echo "[${STAGE}] provenance written to ${PROV}"

vgb_stage_end "pb_folders=${N_PB_DIRS} astex_folders=${N_AX_DIRS}"
vgb_mark_done "$STAGE"
