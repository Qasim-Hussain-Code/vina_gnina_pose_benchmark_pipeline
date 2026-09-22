#!/usr/bin/env python3
"""Per-complex docking, pose extraction and the disk gate. Driven by 06_dock.sh.

Not meant to be run by hand. 06_dock.sh owns the stage stamp, the resource log
and the temporary directory; this owns the loop.

Three things here are worth reading before trusting the numbers.

Pose extraction. Every method is asked for PDBQT output and every pose is read
back through Meeko against the prepared ligand. The obvious alternative is to
take smina's SDF output directly, which it will happily write. That output has
lost bond orders: RDKit reads the first pose of a test complex as
[C]N([N][C]c1[c][c][c][c]c1... and neither RDKit's CalcRMS nor PoseBusters'
robust_rmsd can find a substructure match against the crystal ligand, so the
RMSD comes back NaN and the complex silently drops out of the success rate.
Meeko's PDBQT reader maps the coordinates back onto the molecule the pipeline
prepared, which keeps the bond orders it started with.

The null arm. A0 runs no search. It takes the generated conformer, moves its
centroid onto the box centre and applies a random rotation drawn from the seed
in project.conf. That is the floor: it is what you get from knowing the pocket
and nothing else. Any method whose success rate is not well clear of it has not
earned the compute.

Scores. Vina, Vinardo and the GNINA affinity head all report a number in
kcal/mol. It is the output of an empirical function fitted to a training set,
not a measured or calculated free energy of binding, and the column is named
score_kcal_per_mol_<function> so that the function is attached to the number
everywhere it appears.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import math
import multiprocessing as mp
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

RUN_COLUMNS = [
    # run_kind separates the arm runs from the seed-variance replicates. Without
    # it the five repeat-seed runs of one complex land in the same
    # arm/method/dataset cell as the single production run and that complex is
    # counted six times in the success rate. They share an arm name on purpose,
    # because they are the same protocol; only their role differs.
    "run_kind", "arm", "method", "dataset", "key", "seed", "status", "reason",
    "n_poses", "top_score", "score_function", "score_units",
    "box_center_x", "box_center_y", "box_center_z", "box_size", "exhaustiveness",
    "conformer_source", "receptor", "elapsed_s", "peak_rss_mb", "jobs",
    "pose_archive", "pose_archive_bytes", "recorded",
]

# Which conformer and which box each arm uses. Everything else about an arm is
# the same, which is the point: the arms differ in what information the search
# is given, not in how the search is run.
ARM_SPEC = {
    "A0_null":           {"conformer": "genconf", "box": "ligand_centred", "search": False},
    "A1_refconf_refbox": {"conformer": "refconf", "box": "ligand_centred", "search": True},
    "A2_genconf_refbox": {"conformer": "genconf", "box": "ligand_centred", "search": True},
    "A3_genconf_detbox": {"conformer": "genconf", "box": "detected",       "search": True},
    "A4_crossdock":      {"conformer": "genconf", "box": "partner_ligand", "search": True},
}

SCORE_FUNCTION = {"vina": "vina", "vinardo": "vinardo", "gnina": "gnina_cnn"}


# ---------------------------------------------------------------------------
def peak_rss_mb() -> float:
    """Peak resident memory of this process and its finished children, in MB."""
    try:
        import resource
        r = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        # Linux reports kilobytes; macOS reports bytes.
        return round(r / 1024.0 if sys.platform != "darwin" else r / 1048576.0, 1)
    except Exception:
        return 0.0


def build_command(method: str, cfg: dict, receptor: Path, ligand: Path,
                  out_pdbqt: Path, centre, size: float, seed: int) -> list[str]:
    cx, cy, cz = centre
    common = [
        "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
        "--size_x", f"{size:g}", "--size_y", f"{size:g}", "--size_z", f"{size:g}",
        "--exhaustiveness", str(cfg["exhaustiveness"]),
        "--num_modes", str(cfg["top_n"]),
        "--seed", str(seed),
        # One core per docking process. Vina's own threading and the process
        # pool would otherwise multiply, and a per-complex wall clock measured
        # under nested parallelism means nothing.
        "--cpu", "1",
        "--out", str(out_pdbqt),
    ]
    if method == "vina":
        return [cfg["vina_bin"]] + common
    if method == "vinardo":
        return [cfg["smina_bin"], "--scoring", "vinardo"] + common
    if method == "gnina":
        return [cfg["gnina_bin"], "--cnn_scoring", cfg["gnina_cnn_scoring"]] + common
    raise ValueError(method)


def parse_top_score(method: str, stdout: str, out_pdbqt: Path):
    """Best score from the program's own table, falling back to the PDBQT.

    Vina and smina print a mode table; GNINA prints its own columns. Reading the
    value out of the output file is the fallback because the table format has
    changed between releases and the REMARK line has not.
    """
    best = None
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                v = float(parts[1])
            except ValueError:
                continue
            best = v
            break
    if best is None and out_pdbqt.is_file():
        for line in out_pdbqt.read_text(errors="replace").splitlines():
            up = line.upper()
            if up.startswith("REMARK") and ("VINA RESULT" in up or "MINIMIZEDAFFINITY" in up
                                            or "CNNAFFINITY" in up):
                for tok in line.replace(":", " ").split():
                    try:
                        best = float(tok)
                        break
                    except ValueError:
                        continue
                if best is not None:
                    break
    return best


def poses_to_sdf_gz(out_pdbqt: Path, dest: Path, top_n: int) -> int:
    """Read docked poses back with Meeko and write top_n to a gzipped SDF.

    Returns the number of poses written. Meeko rebuilds the molecule from the
    prepared ligand's own atom typing, so bond orders survive; see the module
    docstring for what happens without that.
    """
    from rdkit import Chem, RDLogger
    from meeko import PDBQTMolecule, RDKitMolCreate

    RDLogger.DisableLog("rdApp.*")
    pm = PDBQTMolecule.from_file(str(out_pdbqt), skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    mols = [m for m in mols if m is not None]
    if not mols:
        return 0
    mol = mols[0]
    n = min(mol.GetNumConformers(), top_n)
    if n == 0:
        return 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Gzipped SDF rather than SDF. Across every arm, method and dataset the pose
    # archive is the only thing here that grows without bound, and SDF text
    # compresses by roughly four to one.
    with gzip.open(dest, "wt") as fh:
        w = Chem.SDWriter(fh)
        for i in range(n):
            single = Chem.Mol(mol)
            single.RemoveAllConformers()
            single.AddConformer(mol.GetConformer(i), assignId=True)
            single.SetProp("_Name", dest.stem)
            single.SetIntProp("pose_rank", i + 1)
            w.write(single)
        w.close()
    return n


def null_pose(ligand_pdbqt: Path, dest: Path, centre, seed: int) -> int:
    """The floor: the generated conformer, centred in the box, randomly rotated.

    No search, no scoring. This is what knowing the pocket is worth on its own.
    """
    import numpy as np
    from rdkit import Chem, RDLogger
    from rdkit.Geometry import Point3D
    from meeko import PDBQTMolecule, RDKitMolCreate

    RDLogger.DisableLog("rdApp.*")
    pm = PDBQTMolecule.from_file(str(ligand_pdbqt), skip_typing=True)
    mols = [m for m in RDKitMolCreate.from_pdbqt_mol(pm) if m is not None]
    if not mols:
        return 0
    mol = mols[0]
    conf = mol.GetConformer(0)
    pos = conf.GetPositions()

    rng = np.random.default_rng(seed)
    # Uniform random rotation via a QR decomposition of a Gaussian matrix, with
    # the sign fix that makes it uniform on SO(3) rather than O(3). Sampling
    # three Euler angles uniformly would bias towards the poles.
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    Q = Q @ np.diag(np.sign(np.diag(R)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    moved = (Q @ (pos - pos.mean(axis=0)).T).T + np.asarray(centre, dtype=float)
    for i in range(mol.GetNumAtoms()):
        x, y, z = moved[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    single = Chem.Mol(mol)
    single.RemoveAllConformers()
    single.AddConformer(conf, assignId=True)
    single.SetProp("_Name", dest.stem)
    single.SetIntProp("pose_rank", 1)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dest, "wt") as fh:
        w = Chem.SDWriter(fh)
        w.write(single)
        w.close()
    return 1


# ---------------------------------------------------------------------------
_CFG: dict = {}


def _init(cfg):
    global _CFG
    _CFG = cfg


def _dock_one(job: dict) -> dict:
    cfg = _CFG
    t0 = time.time()
    stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    row = {
        "run_kind": job.get("run_kind", "arm"),
        "arm": job["arm"], "method": cfg["method"], "dataset": job["dataset"],
        "key": job["key"], "seed": job["seed"], "status": "ok", "reason": "",
        "score_function": SCORE_FUNCTION[cfg["method"]],
        "score_units": "kcal/mol (empirical score, not a measured free energy)",
        "box_center_x": job["centre"][0], "box_center_y": job["centre"][1],
        "box_center_z": job["centre"][2], "box_size": cfg["box_size"],
        "exhaustiveness": job.get("exhaustiveness") or cfg["exhaustiveness"],
        "conformer_source": job["conformer"], "receptor": job["receptor"],
        "jobs": cfg["jobs"], "recorded": stamp,
    }
    archive = Path(job["archive"])
    if archive.is_file() and archive.stat().st_size > 0 and not cfg["force"]:
        row.update({"status": "cached", "elapsed_s": 0.0,
                    "pose_archive": str(archive),
                    "pose_archive_bytes": archive.stat().st_size})
        return row

    ligand = Path(job["ligand"])
    receptor = Path(job["receptor"])
    if not ligand.is_file():
        row.update({"status": "failed", "reason": "prepared ligand missing",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row
    if not receptor.is_file():
        row.update({"status": "failed", "reason": "prepared receptor missing",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    # --- the null arm: no search -------------------------------------------
    if not ARM_SPEC[job["arm"]]["search"]:
        try:
            n = null_pose(ligand, archive, job["centre"], job["seed"])
        except Exception as exc:
            row.update({"status": "failed", "reason": f"null pose failed: {exc!r}"[:200],
                        "elapsed_s": round(time.time() - t0, 2)})
            return row
        row.update({"n_poses": n, "top_score": "", "elapsed_s": round(time.time() - t0, 3),
                    "pose_archive": str(archive),
                    "pose_archive_bytes": archive.stat().st_size if archive.is_file() else 0,
                    "status": "ok" if n else "failed",
                    "reason": "" if n else "null pose produced no molecule"})
        return row

    # --- search -------------------------------------------------------------
    # The work directory name has to be unique across every job in the pool, not
    # just across keys. The convergence grid runs the same key and seed at four
    # exhaustiveness levels, and two of those running concurrently in one
    # directory would write the same out.pdbqt and read each other's poses.
    work = (Path(cfg["work_dir"]) /
            f"{job['key']}_s{job['seed']}_e{job.get('exhaustiveness') or cfg['exhaustiveness']}")
    work.mkdir(parents=True, exist_ok=True)
    out_pdbqt = work / "out.pdbqt"
    # The convergence grid varies exhaustiveness per job; every other run takes
    # the single value from project.conf.
    cfg_run = dict(cfg)
    if job.get("exhaustiveness"):
        cfg_run["exhaustiveness"] = int(job["exhaustiveness"])
    cmd = build_command(cfg["method"], cfg_run, receptor, ligand, out_pdbqt,
                        job["centre"], cfg["box_size"], job["seed"])
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=cfg["timeout_s"])
    except subprocess.TimeoutExpired:
        shutil.rmtree(work, ignore_errors=True)
        row.update({"status": "failed",
                    "reason": f"timed out after {cfg['timeout_s']} s",
                    "elapsed_s": round(time.time() - t0, 2)})
        return row
    if p.returncode != 0 or not out_pdbqt.is_file() or out_pdbqt.stat().st_size == 0:
        tail = (p.stderr or p.stdout or "").strip().splitlines()[-1:] or [""]
        shutil.rmtree(work, ignore_errors=True)
        row.update({"status": "failed", "reason": f"{cfg['method']}: {tail[0]}"[:200],
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    row["top_score"] = parse_top_score(cfg["method"], p.stdout, out_pdbqt)
    try:
        n = poses_to_sdf_gz(out_pdbqt, archive, cfg["top_n"])
    except Exception as exc:
        shutil.rmtree(work, ignore_errors=True)
        row.update({"status": "failed",
                    "reason": f"pose extraction failed: {exc!r}"[:200],
                    "elapsed_s": round(time.time() - t0, 2)})
        return row

    # The PDBQT goes as soon as the poses are in the archive. Across 4,500 runs
    # the intermediates are larger than everything else this pipeline keeps.
    shutil.rmtree(work, ignore_errors=True)

    row.update({"n_poses": n, "elapsed_s": round(time.time() - t0, 2),
                "pose_archive": str(archive),
                "pose_archive_bytes": archive.stat().st_size if archive.is_file() else 0,
                "peak_rss_mb": peak_rss_mb()})
    if n == 0:
        row.update({"status": "failed", "reason": "no pose survived extraction"})
    return row


# ---------------------------------------------------------------------------
def load_boxes(config_dir: Path) -> dict:
    """(kind, key) -> (cx, cy, cz) for every box 05_define_boxes.py wrote."""
    path = config_dir / "boxes.tsv"
    if not path.is_file():
        L.eprint("[error] config/boxes.tsv missing; run 05_define_boxes.py first")
        sys.exit(1)
    out = {}
    for r in L.read_tsv(path):
        if r.get("status") != "ok":
            continue
        try:
            out[(r["arm_box_kind"], r["key"])] = (
                float(r["center_x"]), float(r["center_y"]), float(r["center_z"]))
        except (KeyError, ValueError):
            continue
    return out


def build_jobs(args, conf, boxes) -> list[dict]:
    config_dir = Path(conf["CONFIG_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"
    poses = data_dir / "poses"
    spec = ARM_SPEC[args.arm]
    seed = L.conf_int(conf, "SEED", 20260922)
    jobs = []

    if args.arm == "A4_crossdock":
        align = Path(conf["RESULTS_DIR"]) / "preparation" / "crossdock_align.tsv"
        if not align.is_file():
            L.eprint("[error] crossdock_align.tsv missing; run 04_prepare.py --what crossdock")
            sys.exit(1)
        pairs = {r["pair_id"]: r for r in L.read_tsv(align) if r.get("status") == "ok"}
        cdmap = {r["pair_id"]: r for r in L.read_tsv(config_dir / "crossdock_pairs.tsv")}
        for pid, a in sorted(pairs.items()):
            centre = boxes.get(("partner_ligand", pid))
            if centre is None:
                continue
            meta = cdmap.get(pid, {})
            qid = meta.get("query_complex_id", "")
            # The ligand is the query's generated conformer. The receptor is the
            # partner structure. The reference pose, written by 04, is the query
            # ligand carried into the partner's frame.
            lig = prepared / "ligands" / "genconf" / "posebusters" / f"{qid}.pdbqt"
            jobs.append({
                "arm": args.arm, "dataset": "crossdock", "key": pid, "seed": seed,
                "ligand": str(lig),
                "receptor": str(prepared / "receptors" / "crossdock" /
                                a["receptor_pdb_id"] / "receptor.pdbqt"),
                "centre": centre, "conformer": "genconf",
                "archive": str(poses / args.arm / args.method / "crossdock" / f"{pid}.sdf.gz"),
            })
        return jobs

    rows = L.read_tsv(config_dir / f"dataset_{args.dataset}.tsv")
    for r in rows:
        key = r["complex_id"]
        centre = boxes.get((spec["box"], key))
        if centre is None:
            continue
        jobs.append({
            "arm": args.arm, "dataset": args.dataset, "key": key, "seed": seed,
            "ligand": str(prepared / "ligands" / spec["conformer"] / args.dataset / f"{key}.pdbqt"),
            "receptor": str(prepared / "receptors" / args.dataset / key / "receptor.pdbqt"),
            "centre": centre, "conformer": spec["conformer"],
            "archive": str(poses / args.arm / args.method / args.dataset / f"{key}.sdf.gz"),
        })
    return jobs


def build_convergence_jobs(conf, boxes, method) -> list[dict]:
    """A grid of complexes, exhaustiveness levels and seeds on the A2 arm.

    This is the experiment the exhaustiveness in project.conf was chosen from,
    and it is in the repository rather than in a comment because a protocol
    choice that cannot be reproduced is an opinion. Nothing downstream reads its
    output except 08_analyse.py, which summarises it into
    results/convergence.tsv.

    The complexes are taken at even intervals through the Astex manifest rather
    than as the first N, so the sample is not the alphabetical head of the set.
    Astex rather than PoseBusters because it is the smaller and easier set: if
    sampling is the limit there, it is the limit everywhere.
    """
    config_dir = Path(conf["CONFIG_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"
    poses = data_dir / "poses"
    arm = "A2_genconf_refbox"
    spec = ARM_SPEC[arm]
    levels = [int(x) for x in (conf.get("CONVERGENCE_LEVELS") or "8 16 32 64").split()]
    seeds = (conf.get("SEED_REPLICATES") or "1").split()[:3]
    n_want = L.conf_int(conf, "CONVERGENCE_N_COMPLEXES", 14)

    usable = []
    for r in L.read_tsv(config_dir / "dataset_astex.tsv"):
        key = r["complex_id"]
        lig = prepared / "ligands" / spec["conformer"] / "astex" / f"{key}.pdbqt"
        rec = prepared / "receptors" / "astex" / key / "receptor.pdbqt"
        if lig.is_file() and rec.is_file() and boxes.get((spec["box"], key)):
            usable.append((key, lig, rec, boxes[(spec["box"], key)]))
    if not usable:
        L.eprint("[error] nothing prepared for the convergence grid")
        sys.exit(1)
    step = max(1, len(usable) // n_want)
    chosen = usable[::step][:n_want]

    out = []
    for key, lig, rec, centre in chosen:
        for exh in levels:
            for sd in seeds:
                out.append({
                    "run_kind": "convergence", "arm": arm, "dataset": "astex",
                    "key": key, "seed": int(sd), "exhaustiveness": exh,
                    "ligand": str(lig), "receptor": str(rec), "centre": centre,
                    "conformer": spec["conformer"],
                    "archive": str(poses / "convergence" / method /
                                   f"{key}_e{exh}_s{sd}.sdf.gz"),
                })
    return out


def build_seed_variance_jobs(conf, boxes, method) -> list[dict]:
    """One complex, one arm, once per seed in SEED_REPLICATES.

    Run-to-run spread in top-1 RMSD is part of the measurement and is almost
    never reported. The complex is the first row of the PoseBusters manifest
    that has a prepared ligand, a prepared receptor and a box, so the choice is
    made by the data rather than by me picking a well-behaved one.
    """
    config_dir = Path(conf["CONFIG_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"
    poses = data_dir / "poses"
    arm = "A2_genconf_refbox"
    spec = ARM_SPEC[arm]
    seeds = (conf.get("SEED_REPLICATES") or str(L.conf_int(conf, "SEED", 1))).split()

    chosen = None
    for r in L.read_tsv(config_dir / "dataset_posebusters.tsv"):
        key = r["complex_id"]
        lig = prepared / "ligands" / spec["conformer"] / "posebusters" / f"{key}.pdbqt"
        rec = prepared / "receptors" / "posebusters" / key / "receptor.pdbqt"
        if lig.is_file() and rec.is_file() and boxes.get((spec["box"], key)):
            chosen = (key, lig, rec, boxes[(spec["box"], key)])
            break
    if chosen is None:
        L.eprint("[error] no complex with a prepared ligand, receptor and box; "
                 "run 04_prepare.py and 05_define_boxes.py first")
        sys.exit(1)

    key, lig, rec, centre = chosen
    out = []
    for s in seeds:
        out.append({
            "run_kind": "seed_variance",
            "arm": arm, "dataset": "posebusters", "key": key, "seed": int(s),
            "ligand": str(lig), "receptor": str(rec), "centre": centre,
            "conformer": spec["conformer"],
            "archive": str(poses / "seed_variance" / method / f"{key}_seed{s}.sdf.gz"),
        })
    return out


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--arm", default="")
    ap.add_argument("--method", required=True, choices=["vina", "vinardo", "gnina"])
    ap.add_argument("--dataset", default="")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--stage", default="06_dock")
    ap.add_argument("--seed-variance", action="store_true")
    ap.add_argument("--convergence", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--vina-bin", default="")
    ap.add_argument("--smina-bin", default="")
    ap.add_argument("--gnina-bin", default="")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results_dir = Path(conf["RESULTS_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    boxes = load_boxes(Path(conf["CONFIG_DIR"]))

    cfg = {
        "method": args.method,
        "vina_bin": args.vina_bin, "smina_bin": args.smina_bin,
        "gnina_bin": args.gnina_bin,
        "gnina_cnn_scoring": conf.get("GNINA_CNN_SCORING", "rescore"),
        "box_size": L.conf_float(conf, "BOX_SIZE", 25.0),
        "exhaustiveness": L.conf_int(conf, "EXHAUSTIVENESS", 8),
        "top_n": L.conf_int(conf, "TOP_N_POSES", 5),
        "work_dir": str(args.work_dir),
        "jobs": args.jobs,
        "force": args.force,
        # A single complex should never hold the queue. GNINA running its CNN on
        # CPU is the slow one; 45 minutes is generous for it and is a hard stop
        # for a search that has gone wrong rather than a performance target.
        "timeout_s": 2700,
    }

    if args.convergence:
        jobs = build_convergence_jobs(conf, boxes, args.method)
        out_path = results_dir / "runs" / f"convergence_{args.method}.tsv"
    elif args.seed_variance:
        jobs = build_seed_variance_jobs(conf, boxes, args.method)
        out_path = results_dir / "runs" / f"seed_variance_{args.method}.tsv"
    else:
        if not args.arm:
            L.eprint("[error] --arm required unless --seed-variance")
            return 1
        jobs = build_jobs(args, conf, boxes)
        out_path = results_dir / "runs" / f"{args.arm}_{args.method}_{args.dataset}.tsv"
    if args.limit:
        jobs = jobs[:args.limit]
    if not jobs:
        L.eprint("[error] no jobs built. Check that 04_prepare.py and "
                 "05_define_boxes.py have run for this arm and dataset.")
        return 1

    print(f"[{args.stage}] {len(jobs)} runs, {args.jobs} job(s), "
          f"box {cfg['box_size']:g} A, exhaustiveness {cfg['exhaustiveness']}, "
          f"top {cfg['top_n']} poses kept")
    if args.method == "gnina":
        print(f"[{args.stage}] gnina --cnn_scoring {cfg['gnina_cnn_scoring']}")

    # --- run, with the disk gate after the first ten -----------------------
    poses_root = data_dir / "poses"
    disk_budget_gb = L.conf_int(conf, "DISK_GB", 50)
    rows: list[dict] = []
    t0 = time.time()
    gate_done = False

    def check_gate(done: int):
        """Measure what the finished runs cost, project the rest, refuse if it
        will not fit. Called once, after the first ten have completed."""
        made = [r for r in rows if r.get("pose_archive_bytes")]
        if not made:
            return
        mean_bytes = sum(int(r["pose_archive_bytes"]) for r in made) / len(made)
        remaining = len(jobs) - done
        # This arm's remainder, then the arms and methods still to come. Those
        # are read from project.conf rather than assumed, so dropping an arm
        # from ARMS relaxes the gate as it should.
        n_arms = len([a for a in conf.get("ARMS", "").split(",") if a.strip()])
        n_methods = len([m for m in conf.get("METHODS", "").split(",") if m.strip()])
        projected = mean_bytes * (remaining + len(jobs) * (n_arms * n_methods - 1))
        used = L.du_kb(data_dir) * 1024
        free = shutil.disk_usage(data_dir).free
        proj_gb = projected / 1024 ** 3
        used_gb = used / 1024 ** 3
        free_gb = free / 1024 ** 3
        print(f"[{args.stage}] disk gate: {mean_bytes / 1024:.1f} kB per run measured "
              f"over {len(made)} runs")
        print(f"[{args.stage}] projected {proj_gb:.2f} GB for the remaining arms; "
              f"{used_gb:.2f} GB used, {free_gb:.1f} GB free, budget {disk_budget_gb} GB")
        short = (used_gb + proj_gb) - disk_budget_gb
        if short > 0 or proj_gb > free_gb - 2:
            L.eprint(f"[error] refusing to continue: the remaining arms project to "
                     f"{proj_gb:.2f} GB on top of {used_gb:.2f} GB already used.")
            L.eprint(f"        Shortfall against the {disk_budget_gb} GB budget: "
                     f"{max(short, 0.0):.2f} GB. Free disk: {free_gb:.1f} GB.")
            L.eprint("        Raise --disk in 00_configure.sh, point --data-dir at a "
                     "larger filesystem, or drop an arm or a method.")
            L.write_tsv(out_path, RUN_COLUMNS, rows)
            sys.exit(4)

    if args.jobs > 1:
        with mp.Pool(args.jobs, initializer=_init, initargs=(cfg,)) as pool:
            for i, r in enumerate(pool.imap_unordered(_dock_one, jobs, chunksize=1), 1):
                rows.append(r)
                if i == 10 and not gate_done:
                    gate_done = True
                    check_gate(i)
                if i % 25 == 0 or i == len(jobs):
                    ok = sum(1 for x in rows if x["status"] in ("ok", "cached"))
                    rate = (time.time() - t0) / max(i, 1)
                    print(f"  {i}/{len(jobs)} runs, {ok} ok, "
                          f"{rate:.1f} s/run, eta {(len(jobs) - i) * rate / 60:.0f} min",
                          flush=True)
    else:
        _init(cfg)
        for i, j in enumerate(jobs, 1):
            rows.append(_dock_one(j))
            if i == 10 and not gate_done:
                gate_done = True
                check_gate(i)
            if i % 25 == 0 or i == len(jobs):
                ok = sum(1 for x in rows if x["status"] in ("ok", "cached"))
                rate = (time.time() - t0) / max(i, 1)
                print(f"  {i}/{len(jobs)} runs, {ok} ok, "
                      f"{rate:.1f} s/run, eta {(len(jobs) - i) * rate / 60:.0f} min",
                      flush=True)

    L.write_tsv(out_path, RUN_COLUMNS, rows)
    ok = [r for r in rows if r["status"] in ("ok", "cached")]
    bad = [r for r in rows if r["status"] == "failed"]
    times = sorted(float(r["elapsed_s"]) for r in rows
                   if r["status"] == "ok" and r.get("elapsed_s"))
    print(f"[{args.stage}] {len(ok)} runs produced poses, {len(bad)} failed")
    if times:
        def q(p):
            return times[min(int(p * len(times)), len(times) - 1)]
        print(f"[{args.stage}] per-complex wall clock: median {q(0.5):.1f} s, "
              f"IQR {q(0.25):.1f} to {q(0.75):.1f} s, max {times[-1]:.1f} s "
              f"(jobs={args.jobs})")
    if bad:
        from collections import Counter
        for reason, n in Counter(r["reason"][:60] for r in bad).most_common(5):
            print(f"    {n:4d}  {reason}")
    print(f"[{args.stage}] wrote {out_path}")
    print(f"[{args.stage}] pose archive for this arm: "
          f"{L.du_kb(poses_root) / 1024:.0f} MB total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
