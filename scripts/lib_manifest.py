#!/usr/bin/env python3
"""Turn the extracted benchmark archive into per-dataset manifests.

Called by scripts/02_fetch_benchmarks.sh. Not meant to be run by hand, though
it will work if you give it the same arguments.

Three jobs:

1. Verify that the 308-member list is exactly 308 entries and a strict subset of
   what the PoseBusters authors' archive contains. The archive is the 428-complex
   preprint version; the published paper reports on 308. If the list ever stops
   being a subset, the coordinates and the selection have come apart and the
   stage must fail rather than quietly benchmark a different set.

2. Check that every complex has the four files the archive README promises. A
   complex missing its protein or its crystal ligand cannot be scored, so it is
   excluded here with the reason recorded rather than failing at stage 7 with a
   stack trace.

3. Write one manifest row per complex with absolute paths, so no later stage has
   to know the archive's directory naming.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
from pathlib import Path

# The four files the archive README promises per complex folder. The suffix is
# appended to the folder name, which is "<PDB>_<CCD>".
#   protein.pdb            receptor, solvent stripped, cofactors kept, ligand of
#                          interest removed
#   ligand.sdf             one instance of the ligand of interest, crystal pose
#   ligands.sdf            every instance of that ligand in the asymmetric unit
#   ligand_start_conf.sdf  a conformer the authors generated. Recorded in the
#                          manifest for completeness and not used: 04_prepare.py
#                          generates its own with a recorded seed.
SUFFIXES = {
    "protein_pdb": "_protein.pdb",
    "ligand_sdf": "_ligand.sdf",
    "ligands_sdf": "_ligands.sdf",
    "start_conf_sdf": "_ligand_start_conf.sdf",
}

# A complex cannot be used without these two. start_conf is optional because
# nothing here reads it, and ligands.sdf is optional because the symmetry-mate
# handling in 07_score_poses.py falls back to the single instance.
REQUIRED = ("protein_pdb", "ligand_sdf")

MANIFEST_COLUMNS = [
    "dataset",
    "complex_id",
    "pdb_id",
    "ccd_id",
    "protein_pdb",
    "ligand_sdf",
    "ligands_sdf",
    "start_conf_sdf",
]


def read_ids(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def scan_set(set_dir: Path, wanted: list[str] | None):
    """Return (rows, excluded) for one benchmark set directory.

    wanted restricts the set to a list of complex ids; None takes every folder.
    """
    rows, excluded = [], []
    present = sorted(p.name for p in set_dir.iterdir() if p.is_dir())
    ids = wanted if wanted is not None else present

    present_set = set(present)
    for cid in ids:
        if cid not in present_set:
            excluded.append((cid, "folder absent from the extracted archive"))
            continue
        folder = set_dir / cid
        paths, missing = {}, []
        for key, suffix in SUFFIXES.items():
            p = folder / f"{cid}{suffix}"
            if p.is_file() and p.stat().st_size > 0:
                paths[key] = str(p)
            else:
                paths[key] = ""
                missing.append(key)
        hard = [k for k in missing if k in REQUIRED]
        if hard:
            excluded.append((cid, "missing " + ", ".join(sorted(hard))))
            continue
        # Folder names are "<PDB>_<CCD>". A few CCD ids contain no underscore
        # themselves, so a single split from the left is correct.
        pdb_id, _, ccd_id = cid.partition("_")
        rows.append(
            {
                "complex_id": cid,
                "pdb_id": pdb_id,
                "ccd_id": ccd_id,
                **paths,
            }
        )
    return rows, excluded


def write_manifest(path: Path, dataset: str, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({"dataset": dataset, **r})


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench-dir", required=True, type=Path)
    ap.add_argument("--ids-308", required=True, type=Path)
    ap.add_argument("--config-dir", required=True, type=Path)
    ap.add_argument("--results-dir", required=True, type=Path)
    args = ap.parse_args()

    pb_dir = args.bench_dir / "posebusters_benchmark_set"
    ax_dir = args.bench_dir / "astex_diverse_set"
    for d in (pb_dir, ax_dir):
        if not d.is_dir():
            print(f"[lib_manifest] {d} is missing; did the archive extract?", file=sys.stderr)
            return 1

    archive_ids = sorted(p.name for p in pb_dir.iterdir() if p.is_dir())
    ids_308 = read_ids(args.ids_308)

    # --- the two checks that make this set the published one -----------------
    if len(ids_308) != 308:
        print(f"[lib_manifest] the pinned id list has {len(ids_308)} entries, expected 308.",
              file=sys.stderr)
        print("               config/sources.tsv points at the wrong commit, or the",
              file=sys.stderr)
        print("               upstream file changed. Refusing to build a manifest.",
              file=sys.stderr)
        return 1
    stray = sorted(set(ids_308) - set(archive_ids))
    if stray:
        print(f"[lib_manifest] {len(stray)} ids in the pinned list are not in the",
              file=sys.stderr)
        print("               authors' archive, so the coordinates and the subset",
              file=sys.stderr)
        print(f"               selection have come apart: {', '.join(stray[:10])}",
              file=sys.stderr)
        return 1

    pb_rows, pb_excluded = scan_set(pb_dir, ids_308)
    ax_rows, ax_excluded = scan_set(ax_dir, None)

    write_manifest(args.config_dir / "dataset_posebusters.tsv", "posebusters", pb_rows)
    write_manifest(args.config_dir / "dataset_astex.tsv", "astex", ax_rows)

    # The 120 complexes the authors dropped in peer review are recorded here as
    # an explicit, counted exclusion rather than left implicit in the arithmetic
    # of 428 minus 308.
    dropped_in_review = sorted(set(archive_ids) - set(ids_308))
    excl_path = args.config_dir / "dataset_excluded.tsv"
    stamp = datetime.date.today().isoformat()
    with excl_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["dataset", "complex_id", "stage", "reason", "recorded"])
        for cid in dropped_in_review:
            w.writerow(["posebusters", cid, "02_fetch_benchmarks",
                        "in the 428-complex preprint archive but not in the "
                        "308-complex published set; removed by the authors "
                        "during peer review for crystal contacts", stamp])
        for cid, why in pb_excluded:
            w.writerow(["posebusters", cid, "02_fetch_benchmarks", why, stamp])
        for cid, why in ax_excluded:
            w.writerow(["astex", cid, "02_fetch_benchmarks", why, stamp])

    summary = args.results_dir / "dataset_counts.tsv"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["dataset", "in_archive", "in_manifest", "excluded_here", "recorded"])
        w.writerow(["posebusters", len(archive_ids), len(pb_rows),
                    len(dropped_in_review) + len(pb_excluded), stamp])
        w.writerow(["astex", len(list(ax_dir.iterdir())), len(ax_rows),
                    len(ax_excluded), stamp])

    print(f"[lib_manifest] posebusters: {len(archive_ids)} in archive, "
          f"{len(ids_308)} in the published list, {len(pb_rows)} usable")
    print(f"[lib_manifest] astex: {len(ax_rows)} usable")
    print(f"[lib_manifest] {len(dropped_in_review)} dropped in peer review, "
          f"{len(pb_excluded) + len(ax_excluded)} dropped here for missing files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
