#!/usr/bin/env python3
"""Define the search boxes for every arm.

A search box centred on the crystallographic ligand tells the docking program
where the pocket is. That is the standard protocol, it is what the PoseBusters
paper used for Vina, and it is a strong hint. On a 25 Angstrom cube the hint is
worth less than on a 12 Angstrom one, but it still rules out every other pocket
on the protein before the search starts, and a screening campaign against a
novel target does not have it.

So two kinds of box are written here.

  ligand_centred   a cube of BOX_SIZE Angstrom on the geometric centre of the
                   crystal ligand's heavy atoms. Arms A0, A1 and A2 use it.
  detected         a cube of the same size on the centre of the top-ranked
                   fpocket pocket, computed from the receptor alone. Arm A3
                   uses it.
  partner_ligand   for the cross-docking arm, a cube on the centre of the
                   ligand that was bound in the partner structure. This is the
                   realistic case: you have a holo structure of the target with
                   some other compound in it, and you know where that compound
                   sat.

The detected box must not be chosen using the crystal ligand, or the arm
measures nothing. fpocket's rank 1 pocket is taken every time, with no
inspection and no second choice. The distance from that pocket centre to the
true ligand centroid is recorded, because it answers a question worth asking on
its own: how often does pocket detection find the right pocket at all. It is a
diagnostic column and is never used to pick a pocket.

Box size is 25 Angstrom on a side for every arm. That is the box Buttenschoen
et al. used for Vina, copied so the numbers here are comparable with theirs, and
it is arbitrary in the sense that no pocket is 25 Angstrom wide. A smaller box
would raise every success rate and make the comparison with the published figure
meaningless. Holding it constant across arms is what makes the arms comparable
to each other, which is the point.

Usage
-----
    python scripts/05_define_boxes.py --config project.conf
    python scripts/05_define_boxes.py --config project.conf --dataset astex --jobs 8
"""

from __future__ import annotations

import argparse
import datetime
import multiprocessing as mp
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

BOX_COLUMNS = [
    "arm_box_kind", "dataset", "key", "status", "reason",
    "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
    "method", "pocket_rank", "pocket_score", "pocket_alpha_spheres",
    "distance_to_crystal_centroid", "ligand_fits_in_box", "recorded",
]


def centroid_of_sdf(path: Path):
    import numpy as np

    mol = L.load_heavy(path)
    if mol is None:
        return None, None
    pos = mol.GetConformer().GetPositions()
    return np.asarray(pos.mean(axis=0), dtype=float), np.asarray(pos, dtype=float)


def ligand_extent(pos) -> float:
    """Largest span of the ligand along any axis, in Angstrom."""
    return float((pos.max(axis=0) - pos.min(axis=0)).max())


def fpocket_pockets(receptor_pdb: Path, fpocket_bin: str, workdir: Path):
    """Run fpocket and return pockets sorted by rank, each with centre and score.

    fpocket writes next to its input, so the input is copied into a scratch
    directory first. Running it in place would scatter a _out directory beside
    every prepared receptor and those directories are large.
    """
    import numpy as np

    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / "rec.pdb"
    shutil.copyfile(receptor_pdb, local)
    p = subprocess.run([fpocket_bin, "-f", str(local)],
                       capture_output=True, text=True, timeout=3600)
    out_dir = workdir / "rec_out"
    if p.returncode != 0 or not out_dir.is_dir():
        tail = (p.stderr or p.stdout or "").strip().splitlines()[-1:] or [""]
        return [], f"fpocket failed: {tail[0]}"[:200]

    # Scores, in rank order, from the info file.
    scores: dict[int, dict] = {}
    info = out_dir / "rec_info.txt"
    if info.is_file():
        cur = None
        for line in info.read_text(errors="replace").splitlines():
            s = line.strip()
            if s.lower().startswith("pocket") and ":" in s:
                try:
                    cur = int(s.split()[1])
                except (IndexError, ValueError):
                    cur = None
                if cur is not None:
                    scores[cur] = {}
            elif cur is not None and ":" in s:
                k, _, v = s.partition(":")
                try:
                    scores[cur][k.strip().lower()] = float(v.strip())
                except ValueError:
                    pass

    pockets = []
    pdir = out_dir / "pockets"
    if pdir.is_dir():
        for f in pdir.glob("pocket*_vert.pqr"):
            try:
                rank = int(f.name.split("pocket")[1].split("_")[0])
            except (IndexError, ValueError):
                continue
            pts = []
            for line in f.read_text(errors="replace").splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        pts.append([float(line[30:38]), float(line[38:46]),
                                    float(line[46:54])])
                    except ValueError:
                        continue
            if not pts:
                continue
            arr = np.asarray(pts, dtype=float)
            sc = scores.get(rank, {})
            pockets.append({
                "rank": rank,
                "center": arr.mean(axis=0),
                "n_alpha": len(pts),
                "score": sc.get("score"),
            })
    pockets.sort(key=lambda d: d["rank"])
    return pockets, ""


_CTX: dict = {}


def _init(ctx):
    global _CTX
    _CTX = ctx


def _box_job(item):
    """One receptor: ligand-centred box and, if asked, the detected box."""
    kind, dataset, key, receptor_dir, ligand_sdf, explicit_centre = item
    size = float(_CTX["box_size"])
    stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    rows = []

    # For a self-docking complex the box centre is the crystal ligand's
    # centroid, read from the reference SDF. For a cross-docking pair it is the
    # centroid of the ligand that was bound in the partner structure, which
    # 04_prepare.py recorded while it was stripping that ligand out.
    if explicit_centre is not None:
        cen = np.asarray(explicit_centre, dtype=float)
        pos = None
    else:
        cen, pos = centroid_of_sdf(Path(ligand_sdf)) if ligand_sdf else (None, None)

    if kind in ("self", "crossdock"):
        label = "ligand_centred" if kind == "self" else "partner_ligand"
        if cen is None:
            rows.append({"arm_box_kind": label, "dataset": dataset, "key": key,
                         "status": "failed", "reason": "no box centre available",
                         "method": "centroid_of_reference_ligand", "recorded": stamp})
        else:
            fits = ligand_extent(pos) <= size if pos is not None else ""
            rows.append({
                "arm_box_kind": label, "dataset": dataset, "key": key,
                "status": "ok", "reason": "",
                "center_x": round(float(cen[0]), 3),
                "center_y": round(float(cen[1]), 3),
                "center_z": round(float(cen[2]), 3),
                "size_x": size, "size_y": size, "size_z": size,
                "method": ("centroid_of_crystal_ligand" if kind == "self"
                           else "centroid_of_partner_bound_ligand"),
                "distance_to_crystal_centroid": 0.0 if kind == "self" else "",
                "ligand_fits_in_box": int(fits) if fits != "" else "",
                "recorded": stamp,
            })

    if kind == "self" and _CTX["do_detected"]:
        rec = Path(receptor_dir) / "clean.pdb"
        row = {"arm_box_kind": "detected", "dataset": dataset, "key": key,
               "method": "fpocket_rank1", "recorded": stamp,
               "size_x": size, "size_y": size, "size_z": size}
        if not rec.is_file():
            row.update({"status": "failed", "reason": "no clean.pdb for this receptor"})
            rows.append(row)
            return rows
        with tempfile.TemporaryDirectory(prefix="fpocket_") as td:
            pockets, err = fpocket_pockets(rec, _CTX["fpocket"], Path(td))
        if not pockets:
            row.update({"status": "failed", "reason": err or "fpocket found no pockets"})
            rows.append(row)
            return rows
        # Rank 1, always. No inspection, no alternative.
        top = pockets[0]
        c = top["center"]
        row.update({
            "status": "ok", "reason": "",
            "center_x": round(float(c[0]), 3),
            "center_y": round(float(c[1]), 3),
            "center_z": round(float(c[2]), 3),
            "pocket_rank": top["rank"], "pocket_score": top["score"],
            "pocket_alpha_spheres": top["n_alpha"],
        })
        if cen is not None and pos is not None:
            d = float(np.linalg.norm(c - cen))
            row["distance_to_crystal_centroid"] = round(d, 3)
            # Does the true ligand still lie inside the detected box? If not,
            # the arm cannot possibly succeed on this complex and the failure is
            # a pocket-detection failure rather than a docking failure. Counting
            # these separately is the only way the A3 number means anything.
            half = size / 2.0
            inside = bool(np.all(np.abs(pos - c) <= half))
            row["ligand_fits_in_box"] = int(inside)
        rows.append(row)

    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--dataset", default="")
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-detected", action="store_true",
                    help="skip the fpocket arm, which is the slow part")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    config_dir = Path(conf["CONFIG_DIR"])
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"
    size = L.conf_float(conf, "BOX_SIZE", 25.0)
    jobs = args.jobs or L.conf_int(conf, "JOBS", 1)
    arms = conf.get("ARMS", "")

    out_path = config_dir / "boxes.tsv"
    if out_path.is_file() and not args.force:
        print(f"[05_define_boxes] {out_path.name} exists; skipping (--force to rebuild).")
        return 0

    do_detected = ("A3_genconf_detbox" in arms) and not args.no_detected
    fpocket = shutil.which("fpocket") or ""
    if do_detected and not fpocket:
        # fpocket lives in the smina environment, which is not the active one.
        root = Path(conf["CONDA_SH"]).parent.parent.parent
        cand = root / "envs" / conf.get("CONDA_ENV_SMINA", "vgb_smina") / "bin" / "fpocket"
        fpocket = str(cand) if cand.is_file() else ""
    if do_detected and not fpocket:
        L.eprint("[error] fpocket not found and arm A3 is in ARMS. Install it "
                 "(conda env vgb_smina) or re-run 00_configure.sh without A3.")
        return 3

    datasets = [args.dataset] if args.dataset else ["posebusters", "astex"]
    items = []
    for ds in datasets:
        p = config_dir / f"dataset_{ds}.tsv"
        if not p.is_file():
            L.eprint(f"[error] {p} missing; run 02_fetch_benchmarks.sh first")
            return 1
        rows = L.read_tsv(p)
        if args.limit:
            rows = rows[:args.limit]
        for r in rows:
            items.append(("self", ds, r["complex_id"],
                          str(prepared / "receptors" / ds / r["complex_id"]),
                          str(prepared / "reference" / ds / f"{r['complex_id']}.sdf"),
                          None))

    # The cross-docking box sits on the ligand the partner structure was solved
    # with, not on the transformed reference pose. Using the reference would
    # hand the arm the answer; the partner's own bound ligand is what an
    # experimenter actually has in front of them. 04_prepare.py recorded that
    # centroid into crossdock_align.tsv while stripping the ligand out.
    align_path = results_dir / "preparation" / "crossdock_align.tsv"
    if align_path.is_file():
        n_skipped = 0
        for a in L.read_tsv(align_path):
            if a.get("status") != "ok":
                n_skipped += 1
                continue
            try:
                centre = (float(a["partner_ligand_center_x"]),
                          float(a["partner_ligand_center_y"]),
                          float(a["partner_ligand_center_z"]))
            except (KeyError, ValueError):
                n_skipped += 1
                continue
            items.append(("crossdock", "crossdock", a["pair_id"],
                          str(prepared / "receptors" / "crossdock" / a["receptor_pdb_id"]),
                          "", centre))
            if args.limit and sum(1 for i in items if i[0] == "crossdock") >= args.limit:
                break
        if n_skipped:
            print(f"[05_define_boxes] {n_skipped} cross-docking pairs have no usable "
                  f"box centre and are skipped; see {align_path.name}")
    else:
        print("[05_define_boxes] no crossdock_align.tsv; run 04_prepare.py --what "
              "crossdock first if the A4 arm is wanted")

    ctx = {"box_size": size, "do_detected": do_detected, "fpocket": fpocket}
    print(f"[05_define_boxes] {len(items)} receptors, box {size:g} A cube, "
          f"detected arm {'on' if do_detected else 'off'}, {jobs} job(s)")

    t0 = time.time()
    rows = []
    if jobs > 1:
        with mp.Pool(jobs, initializer=_init, initargs=(ctx,)) as pool:
            for i, got in enumerate(pool.imap_unordered(_box_job, items, chunksize=1), 1):
                rows.extend(got)
                if i % 50 == 0 or i == len(items):
                    print(f"  boxes: {i}/{len(items)}", flush=True)
    else:
        _init(ctx)
        for i, it in enumerate(items, 1):
            rows.extend(_box_job(it))
            if i % 50 == 0 or i == len(items):
                print(f"  boxes: {i}/{len(items)}", flush=True)

    L.write_tsv(out_path, BOX_COLUMNS, rows)
    prov = results_dir / "boxes"
    prov.mkdir(parents=True, exist_ok=True)
    L.write_tsv(prov / "box_provenance.tsv", BOX_COLUMNS, rows)

    # --- the two numbers this stage exists to produce ----------------------
    det = [r for r in rows if r["arm_box_kind"] == "detected" and r["status"] == "ok"]
    lig = [r for r in rows if r["arm_box_kind"] == "ligand_centred" and r["status"] == "ok"]
    bad = [r for r in rows if r["status"] == "failed"]
    print(f"[05_define_boxes] {len(lig)} ligand-centred boxes, {len(det)} detected "
          f"boxes, {len(bad)} failed")
    if det:
        near = [r for r in det if isinstance(r.get("distance_to_crystal_centroid"), float)]
        if near:
            d = sorted(r["distance_to_crystal_centroid"] for r in near)
            within5 = sum(1 for x in d if x <= 5.0)
            contains = sum(1 for r in det if r.get("ligand_fits_in_box") == 1)
            print(f"[05_define_boxes] fpocket rank 1 centre to crystal ligand centroid: "
                  f"median {d[len(d) // 2]:.1f} A")
            print(f"[05_define_boxes] within 5 A: {within5}/{len(d)} "
                  f"({100.0 * within5 / len(d):.1f} per cent)")
            print(f"[05_define_boxes] crystal ligand lies inside the detected box: "
                  f"{contains}/{len(det)} ({100.0 * contains / len(det):.1f} per cent). "
                  f"The rest cannot succeed in arm A3 whatever the scoring function does.")
    print(f"[05_define_boxes] stage took {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
