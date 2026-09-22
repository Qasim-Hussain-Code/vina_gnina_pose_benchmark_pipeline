#!/usr/bin/env python3
"""Score every pose: symmetry-corrected RMSD and the PoseBusters checks.

One row per pose in results/poses/, one row per run in results/scored/. Stage 8
reads those and computes the success rates; nothing is aggregated here.

What is measured per pose
-------------------------
rmsd                    PoseBusters' own robust_rmsd, in place, heavy atoms,
                        symmetry corrected. This is the success metric and it is
                        their function rather than a local reimplementation, so
                        the number is computed the same way as the number in the
                        paper being compared against.
kabsch_rmsd             the same after optimal superposition. A pose with a low
                        kabsch and a high rmsd has the right conformation in the
                        wrong place, which is a search failure and not a
                        conformer failure.
rmsd_spyrmsd            an independent graph-isomorphism implementation. Recorded
                        so that disagreement between two libraries is a column
                        rather than an assumption.
rmsd_first_match        one arbitrary valid atom correspondence instead of the
                        best one. The gap to rmsd is what symmetry correction is
                        worth, and stage 8 counts how often it moves a pose
                        across the 2 Angstrom line.
pb_* (27 columns)       the PoseBusters checks, verbatim, renamed to ASCII.
pb_valid                every check passed.

Why the validity reference is the archive's protein file and not the prepared
receptor: PoseBusters judges whether a pose is physically possible in the
structure, and the archive file is what the paper's own numbers were computed
against. It contains the cofactors and the crystallisation additives that
04_prepare.py strips before docking, so a pose can in principle be marked as
clashing with a glycerol the docking program never saw. That is a real penalty
and it is the same penalty the reference numbers carry, which is why it is kept.

Usage
-----
    python scripts/07_score_poses.py --config project.conf
    python scripts/07_score_poses.py --config project.conf --arm A2_genconf_refbox
    python scripts/07_score_poses.py --config project.conf --jobs 8 --force
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

POSE_COLUMNS_FIXED = [
    "arm", "method", "dataset", "key", "seed", "pose_rank", "status", "reason",
    "rmsd", "kabsch_rmsd", "rmsd_spyrmsd", "rmsd_first_match",
    "rmsd_implementations_disagree_by", "pb_valid", "pb_checks_failed",
    "recorded",
]

RUN_COLUMNS = [
    "arm", "method", "dataset", "key", "seed", "status", "reason", "n_poses",
    "top_score", "score_function", "score_units", "elapsed_s", "jobs",
    "top1_rmsd", "top1_kabsch_rmsd", "top1_pb_valid", "top1_pass_2a",
    "top1_pass_1a", "top1_success_2a", "top1_success_1a",
    "top5_best_rmsd", "top5_pass_2a", "top5_success_2a",
    "top1_rmsd_first_match", "symmetry_changed_verdict_2a",
    "conformer_source", "box_kind", "recorded",
]

ARM_BOX = {
    "A0_null": "ligand_centred", "A1_refconf_refbox": "ligand_centred",
    "A2_genconf_refbox": "ligand_centred", "A3_genconf_detbox": "detected",
    "A4_crossdock": "partner_ligand",
}


def ascii_col(name: str) -> str:
    """PoseBusters names one column 'rmsd_<=_2a' using non-ASCII characters.

    Every output file in this repository is ASCII, so the check names are
    transliterated rather than passed through. The mapping is mechanical and the
    original order is preserved.
    """
    out = (name.replace("≤", "_le_").replace("≥", "_ge_")
               .replace("å", "a").replace("Å", "a"))
    out = re.sub(r"[^0-9a-zA-Z]+", "_", out).strip("_").lower()
    return "pb_" + out


_CTX: dict = {}
_PB = None


def _init(ctx):
    global _CTX, _PB
    _CTX = ctx
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from posebusters import PoseBusters
    # One PoseBusters object per worker, reused across complexes. Constructing
    # it per pose costs more than the checks do.
    _PB = PoseBusters(config="redock")


def read_poses(archive: Path):
    """Read a gzipped SDF of poses into single-conformer heavy-atom molecules.

    Binary mode, not text. ForwardSDMolSupplier wraps a C++ stream and rejects a
    text-mode handle with "Need a binary mode file object"; opening the gzip with
    "rt" fails on every archive, which is how eight of eight smoke-test runs came
    back unscored with the poses sitting on disk perfectly intact.
    """
    from rdkit import Chem

    out = []
    with gzip.open(archive, "rb") as fh:
        for mol in Chem.ForwardSDMolSupplier(fh, removeHs=True):
            if mol is not None:
                out.append(mol)
    return out


def _score_run(job: dict) -> tuple[dict, list[dict]]:
    stamp = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    base = {k: job[k] for k in ("arm", "method", "dataset", "key", "seed")}
    run = dict(base)
    run.update({
        "status": "ok", "reason": "", "recorded": stamp,
        "top_score": job.get("top_score", ""),
        "score_function": job.get("score_function", ""),
        "score_units": job.get("score_units", ""),
        "elapsed_s": job.get("elapsed_s", ""), "jobs": job.get("jobs", ""),
        "conformer_source": job.get("conformer_source", ""),
        "box_kind": ARM_BOX.get(job["arm"], ""),
    })

    archive = Path(job["archive"])
    ref_path = Path(job["reference"])
    prot_path = Path(job["protein"])
    if not archive.is_file():
        run.update({"status": "failed", "reason": "pose archive missing"})
        return run, []
    if not ref_path.is_file():
        run.update({"status": "failed", "reason": "reference ligand missing"})
        return run, []

    try:
        poses = read_poses(archive)
    except Exception as exc:
        run.update({"status": "failed", "reason": f"pose archive unreadable: {exc!r}"[:200]})
        return run, []
    if not poses:
        run.update({"status": "failed", "reason": "pose archive held no molecule"})
        return run, []
    run["n_poses"] = len(poses)

    ref = L.load_heavy(ref_path)
    if ref is None:
        run.update({"status": "failed", "reason": "reference ligand would not parse"})
        return run, []

    # --- PoseBusters, one call for all poses of this run -------------------
    # Batching amortises loading the protein, which is the expensive part:
    # 1.16 s for one pose against 2.70 s for five.
    pb_rows: list[dict] = []
    pb_error = ""
    if prot_path.is_file():
        try:
            df = _PB.bust(poses, mol_true=str(ref_path), mol_cond=str(prot_path),
                          full_report=False)
            pb_rows = df.to_dict("records")
        except Exception as exc:
            # PoseBusters can raise from inside RDKit, for instance a
            # KekulizeException out of its tautomer canonicalisation fallback.
            # A run that hits that keeps its RMSD, which is computed here
            # independently, and is marked as validity-unknown rather than
            # silently counted as invalid or dropped.
            pb_error = f"{type(exc).__name__}: {exc}"[:180]
    else:
        pb_error = "protein file missing, validity not assessed"

    pose_rows = []
    for i, pose in enumerate(poses):
        pr = dict(base)
        pr.update({"pose_rank": i + 1, "status": "ok", "reason": "", "recorded": stamp})

        pbr = pb_rows[i] if i < len(pb_rows) else {}
        rmsd = pbr.get("rmsd")
        kabsch = pbr.get("kabsch_rmsd")
        if rmsd is None:
            # No PoseBusters row, or it did not return an RMSD. Compute it here
            # with the same function they use, then plain CalcRMS.
            rmsd, kabsch = L.symm_rmsd_posebusters(pose, ref)
            if rmsd is None:
                rmsd = L.symm_rmsd_rdkit(pose, ref)
        pr["rmsd"] = round(float(rmsd), 4) if _num(rmsd) else ""
        pr["kabsch_rmsd"] = round(float(kabsch), 4) if _num(kabsch) else ""

        spy = L.symm_rmsd_spyrmsd(pose, ref)
        pr["rmsd_spyrmsd"] = round(spy, 4) if spy is not None else ""
        if _num(rmsd) and spy is not None:
            pr["rmsd_implementations_disagree_by"] = round(abs(float(rmsd) - spy), 4)

        fm = L.first_match_rmsd(pose, ref)
        pr["rmsd_first_match"] = round(fm, 4) if fm is not None else ""

        checks = {}
        for k, v in pbr.items():
            if k in ("rmsd", "kabsch_rmsd", "centroid_distance", "file", "molecule"):
                continue
            if isinstance(v, bool):
                checks[ascii_col(k)] = int(v)
        pr.update(checks)
        if checks:
            failed = [k for k, v in checks.items() if v == 0]
            pr["pb_valid"] = int(not failed)
            pr["pb_checks_failed"] = ";".join(sorted(failed))
        else:
            pr["pb_valid"] = ""
            pr["pb_checks_failed"] = ""
            pr["status"] = "validity_unknown"
            pr["reason"] = pb_error
        pose_rows.append(pr)

    # --- the run-level summary --------------------------------------------
    pass_2a = L.conf_float(_CTX["conf"], "RMSD_PASS", 2.0)
    pass_1a = L.conf_float(_CTX["conf"], "RMSD_STRICT", 1.0)

    top = pose_rows[0]
    r1 = top["rmsd"] if top["rmsd"] != "" else None
    run["top1_rmsd"] = top["rmsd"]
    run["top1_kabsch_rmsd"] = top["kabsch_rmsd"]
    run["top1_pb_valid"] = top["pb_valid"]
    run["top1_rmsd_first_match"] = top["rmsd_first_match"]
    if r1 is not None:
        run["top1_pass_2a"] = int(r1 <= pass_2a)
        run["top1_pass_1a"] = int(r1 <= pass_1a)
        if top["pb_valid"] != "":
            run["top1_success_2a"] = int(r1 <= pass_2a and top["pb_valid"] == 1)
            run["top1_success_1a"] = int(r1 <= pass_1a and top["pb_valid"] == 1)
        # Did symmetry correction move this pose across the line?
        if top["rmsd_first_match"] != "":
            run["symmetry_changed_verdict_2a"] = int(
                (r1 <= pass_2a) != (float(top["rmsd_first_match"]) <= pass_2a))

    # Top-5 takes the best-RMSD pose among those kept, and asks whether that
    # pose is also valid. That is the right question for top-5: it is the best
    # pose a person could have picked out of the list.
    with_rmsd = [p for p in pose_rows if p["rmsd"] != ""]
    if with_rmsd:
        best = min(with_rmsd, key=lambda p: float(p["rmsd"]))
        run["top5_best_rmsd"] = best["rmsd"]
        run["top5_pass_2a"] = int(float(best["rmsd"]) <= pass_2a)
        if best["pb_valid"] != "":
            run["top5_success_2a"] = int(float(best["rmsd"]) <= pass_2a
                                         and best["pb_valid"] == 1)
    if not with_rmsd:
        run.update({"status": "failed",
                    "reason": "no pose yielded an RMSD; " + (pb_error or "unknown")})
    return run, pose_rows


def _num(v) -> bool:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return False
    return f == f    # excludes NaN


# ---------------------------------------------------------------------------
def build_jobs(conf, args) -> list[dict]:
    """One job per docking run, read from the tables 06_dock.sh wrote."""
    results_dir = Path(conf["RESULTS_DIR"])
    config_dir = Path(conf["CONFIG_DIR"])
    data_dir = Path(conf["DATA_DIR"])
    prepared = data_dir / "prepared"

    refs: dict[tuple[str, str], tuple[str, str]] = {}
    for ds in ("posebusters", "astex"):
        p = config_dir / f"dataset_{ds}.tsv"
        if p.is_file():
            for r in L.read_tsv(p):
                refs[(ds, r["complex_id"])] = (
                    str(prepared / "reference" / ds / f"{r['complex_id']}.sdf"),
                    r["protein_pdb"])

    # Cross-docking references were written into the partner's frame by
    # 04_prepare.py, and the protein for the validity checks is the partner
    # receptor rather than the query's.
    cd_align = results_dir / "preparation" / "crossdock_align.tsv"
    if cd_align.is_file():
        for a in L.read_tsv(cd_align):
            if a.get("status") != "ok":
                continue
            refs[("crossdock", a["pair_id"])] = (
                str(prepared / "reference" / "crossdock" / f"{a['pair_id']}.sdf"),
                str(data_dir / "crossdock" / "receptor_pdb" / f"{a['receptor_pdb_id']}.pdb"))

    jobs = []
    run_dir = results_dir / "runs"
    if not run_dir.is_dir():
        L.eprint("[error] results/runs is empty; run 06_dock.sh first")
        sys.exit(1)
    for tsv in sorted(run_dir.glob("*.tsv")):
        for r in L.read_tsv(tsv):
            if r.get("status") not in ("ok", "cached"):
                continue
            if args.arm and r["arm"] != args.arm:
                continue
            if args.method and r["method"] != args.method:
                continue
            if args.dataset and r["dataset"] != args.dataset:
                continue
            ref = refs.get((r["dataset"], r["key"]))
            if ref is None:
                continue
            jobs.append({
                "arm": r["arm"], "method": r["method"], "dataset": r["dataset"],
                "key": r["key"], "seed": r.get("seed", ""),
                "archive": r["pose_archive"], "reference": ref[0], "protein": ref[1],
                "top_score": r.get("top_score", ""),
                "score_function": r.get("score_function", ""),
                "score_units": r.get("score_units", ""),
                "elapsed_s": r.get("elapsed_s", ""), "jobs": r.get("jobs", ""),
                "conformer_source": r.get("conformer_source", ""),
                "source_table": tsv.name,
            })
    return jobs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--arm", default="")
    ap.add_argument("--method", default="")
    ap.add_argument("--dataset", default="")
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results_dir = Path(conf["RESULTS_DIR"])
    jobs_n = args.jobs or L.conf_int(conf, "JOBS", 1)

    out_pose = results_dir / "poses" / "pose_scores.tsv"
    out_run = results_dir / "scored" / "run_scores.tsv"
    if out_run.is_file() and not args.force:
        print(f"[07_score_poses] {out_run} exists; skipping (--force to redo).")
        return 0

    jobs = build_jobs(conf, args)
    if args.limit:
        jobs = jobs[:args.limit]
    if not jobs:
        L.eprint("[error] no completed docking runs matched. Run 06_dock.sh first.")
        return 1

    print(f"[07_score_poses] {len(jobs)} runs to score, {jobs_n} job(s)")
    ctx = {"conf": conf}
    t0 = time.time()
    run_rows, pose_rows = [], []

    if jobs_n > 1:
        with mp.Pool(jobs_n, initializer=_init, initargs=(ctx,)) as pool:
            for i, (rr, pp) in enumerate(pool.imap_unordered(_score_run, jobs, chunksize=1), 1):
                run_rows.append(rr)
                pose_rows.extend(pp)
                if i % 100 == 0 or i == len(jobs):
                    rate = (time.time() - t0) / i
                    print(f"  {i}/{len(jobs)} runs, {rate:.2f} s/run, "
                          f"eta {(len(jobs) - i) * rate / 60:.0f} min", flush=True)
    else:
        _init(ctx)
        for i, j in enumerate(jobs, 1):
            rr, pp = _score_run(j)
            run_rows.append(rr)
            pose_rows.extend(pp)
            if i % 100 == 0 or i == len(jobs):
                rate = (time.time() - t0) / i
                print(f"  {i}/{len(jobs)} runs, {rate:.2f} s/run, "
                      f"eta {(len(jobs) - i) * rate / 60:.0f} min", flush=True)

    # The PoseBusters check columns are discovered from the data rather than
    # hardcoded, because the set of checks has changed between releases and a
    # hardcoded list would silently drop a new one.
    extra = sorted({k for p in pose_rows for k in p
                    if k.startswith("pb_") and k not in ("pb_valid", "pb_checks_failed")})
    L.write_tsv(out_pose, POSE_COLUMNS_FIXED + extra, pose_rows)
    L.write_tsv(out_run, RUN_COLUMNS, run_rows)

    ok = [r for r in run_rows if r["status"] == "ok"]
    unknown = [p for p in pose_rows if p["status"] == "validity_unknown"]
    print(f"[07_score_poses] scored {len(ok)} runs, "
          f"{len(run_rows) - len(ok)} could not be scored")
    if unknown:
        print(f"[07_score_poses] {len(unknown)} poses have an RMSD but no validity "
              f"verdict, because PoseBusters raised on them. They are counted as "
              f"validity-unknown and excluded from the pb_valid rates, not as invalid.")
    print(f"[07_score_poses] wrote {out_run} and {out_pose} in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
