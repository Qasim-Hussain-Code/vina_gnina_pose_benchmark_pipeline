#!/usr/bin/env python3
"""Aggregate the scored runs into the tables the README quotes.

Reads results/scored/run_scores.tsv and results/poses/pose_scores.tsv and writes:

    results/success_rates.tsv         the headline table, one row per
                                      arm x method x dataset
    results/arm_comparison.tsv        what removing each piece of crystal
                                      information costs, in percentage points
    results/rmsd_distributions.tsv    median and IQR, because a fraction below a
                                      cutoff hides the shape
    results/timing_distributions.tsv  per-complex wall clock as a distribution,
                                      since a few ligands dominate the mean
    results/pb_failure_waterfall.tsv  which checks fail most often, per method
    results/symmetry_effect.tsv       how often symmetry correction moved a pose
                                      across the 2 Angstrom line
    results/seed_variance.tsv         spread of top-1 RMSD over repeat seeds
    results/excluded.tsv              every excluded complex, every stage, with
                                      its reason
    results/headline.tsv              the handful of numbers the summary
                                      paragraph of the README uses

Denominators
------------
Three different ones are possible and they give different numbers, so all three
are in the table.

    n_attempted   runs that were built for this arm, method and dataset
    n_scored      runs that produced at least one pose with an RMSD
    n_assessed    runs whose top pose also got a validity verdict

Rates for RMSD use n_scored. Rates involving validity use n_assessed. A run that
failed to dock at all is a failure of the method and counts against it in the
n_attempted column, which is why that column is reported next to the others
rather than dropped.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

ARM_LABEL = {
    "A0_null": "no search, conformer dropped in the box",
    "A1_refconf_refbox": "crystal conformer, ligand-centred box",
    "A2_genconf_refbox": "generated conformer, ligand-centred box",
    "A3_genconf_detbox": "generated conformer, detected box",
    "A4_crossdock": "cross-docked into another structure of the same protein",
}


def wilson(k: int, n: int, z: float = 1.959963985) -> tuple[float, float]:
    """Wilson score interval for a proportion.

    Wilson rather than the normal approximation because several cells here have
    fewer than thirty runs and a few have a proportion near zero, where the
    normal interval runs below zero and stops meaning anything.
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def quantiles(vals: list[float]):
    if not vals:
        return {}
    s = sorted(vals)

    def q(p):
        if len(s) == 1:
            return s[0]
        i = p * (len(s) - 1)
        lo, hi = int(math.floor(i)), int(math.ceil(i))
        return s[lo] + (s[hi] - s[lo]) * (i - lo)

    return {"n": len(s), "min": s[0], "q25": q(0.25), "median": q(0.5),
            "q75": q(0.75), "q90": q(0.90), "max": s[-1],
            "mean": sum(s) / len(s)}


def f(x, nd=4):
    return "" if x is None or (isinstance(x, float) and x != x) else round(float(x), nd)


def pct(k, n):
    return "" if not n else round(100.0 * k / n, 2)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results_dir = Path(conf["RESULTS_DIR"])
    config_dir = Path(conf["CONFIG_DIR"])
    pass_2a = L.conf_float(conf, "RMSD_PASS", 2.0)
    pass_1a = L.conf_float(conf, "RMSD_STRICT", 1.0)
    stamp = datetime.date.today().isoformat()

    run_path = results_dir / "scored" / "run_scores.tsv"
    if not run_path.is_file():
        L.eprint("[error] results/scored/run_scores.tsv missing; run 07_score_poses.py")
        return 1
    runs = L.read_tsv(run_path)
    poses = L.read_tsv(results_dir / "poses" / "pose_scores.tsv") \
        if (results_dir / "poses" / "pose_scores.tsv").is_file() else []

    def fv(r, k):
        v = r.get(k, "")
        if v == "" or v is None:
            return None
        try:
            x = float(v)
        except ValueError:
            return None
        return None if x != x else x

    # -----------------------------------------------------------------------
    # 1. success rates
    # -----------------------------------------------------------------------
    # Only the production arm runs. The seed-variance replicates share an arm
    # name with A2 and would otherwise be counted as five extra complexes.
    arm_runs = [r for r in runs if (r.get("run_kind") or "arm") == "arm"]
    cells = defaultdict(list)
    for r in arm_runs:
        cells[(r["arm"], r["method"], r["dataset"])].append(r)

    sr_cols = ["arm", "arm_description", "method", "dataset",
               "n_attempted", "n_scored", "n_assessed",
               "rmsd_le_2a_n", "rmsd_le_2a_pct", "rmsd_le_2a_ci_low", "rmsd_le_2a_ci_high",
               "pb_valid_n", "pb_valid_pct",
               "success_2a_n", "success_2a_pct", "success_2a_ci_low", "success_2a_ci_high",
               "rmsd_le_1a_n", "rmsd_le_1a_pct",
               "success_1a_n", "success_1a_pct",
               "top5_rmsd_le_2a_n", "top5_rmsd_le_2a_pct",
               "top5_success_2a_n", "top5_success_2a_pct",
               "top1_minus_top5_gap_pp", "rmsd_threshold_a", "strict_threshold_a",
               "recorded"]
    sr_rows = []
    for (arm, method, dataset), rs in sorted(cells.items()):
        scored = [r for r in rs if fv(r, "top1_rmsd") is not None]
        assessed = [r for r in scored if r.get("top1_pb_valid") not in ("", None)]
        k2 = sum(1 for r in scored if fv(r, "top1_rmsd") <= pass_2a)
        k1 = sum(1 for r in scored if fv(r, "top1_rmsd") <= pass_1a)
        kv = sum(1 for r in assessed if str(r.get("top1_pb_valid")) == "1")
        ks2 = sum(1 for r in assessed if str(r.get("top1_success_2a")) == "1")
        ks1 = sum(1 for r in assessed if str(r.get("top1_success_1a")) == "1")
        t5 = [r for r in rs if fv(r, "top5_best_rmsd") is not None]
        k5 = sum(1 for r in t5 if fv(r, "top5_best_rmsd") <= pass_2a)
        t5a = [r for r in t5 if r.get("top5_success_2a") not in ("", None)]
        k5s = sum(1 for r in t5a if str(r.get("top5_success_2a")) == "1")
        lo2, hi2 = wilson(k2, len(scored))
        los, his = wilson(ks2, len(assessed))
        gap = ""
        if len(scored) and len(t5):
            gap = round(100.0 * k5 / len(t5) - 100.0 * k2 / len(scored), 2)
        sr_rows.append({
            "arm": arm, "arm_description": ARM_LABEL.get(arm, ""), "method": method,
            "dataset": dataset, "n_attempted": len(rs), "n_scored": len(scored),
            "n_assessed": len(assessed),
            "rmsd_le_2a_n": k2, "rmsd_le_2a_pct": pct(k2, len(scored)),
            "rmsd_le_2a_ci_low": f(100 * lo2, 2), "rmsd_le_2a_ci_high": f(100 * hi2, 2),
            "pb_valid_n": kv, "pb_valid_pct": pct(kv, len(assessed)),
            "success_2a_n": ks2, "success_2a_pct": pct(ks2, len(assessed)),
            "success_2a_ci_low": f(100 * los, 2), "success_2a_ci_high": f(100 * his, 2),
            "rmsd_le_1a_n": k1, "rmsd_le_1a_pct": pct(k1, len(scored)),
            "success_1a_n": ks1, "success_1a_pct": pct(ks1, len(assessed)),
            "top5_rmsd_le_2a_n": k5, "top5_rmsd_le_2a_pct": pct(k5, len(t5)),
            "top5_success_2a_n": k5s, "top5_success_2a_pct": pct(k5s, len(t5a)),
            "top1_minus_top5_gap_pp": gap,
            "rmsd_threshold_a": pass_2a, "strict_threshold_a": pass_1a,
            "recorded": stamp,
        })
    L.write_tsv(results_dir / "success_rates.tsv", sr_cols, sr_rows)

    # -----------------------------------------------------------------------
    # 2. arm comparison: what each piece of crystal information is worth
    # -----------------------------------------------------------------------
    by = {(r["arm"], r["method"], r["dataset"]): r for r in sr_rows}
    comparisons = [
        ("A1_refconf_refbox", "A2_genconf_refbox",
         "cost of generating the starting conformer instead of using the crystal one"),
        ("A2_genconf_refbox", "A3_genconf_detbox",
         "cost of taking the box from pocket detection instead of the crystal ligand"),
        ("A1_refconf_refbox", "A3_genconf_detbox",
         "cost of removing both, which is the gap between the quoted protocol and a blind one"),
        ("A2_genconf_refbox", "A4_crossdock",
         "cost of docking into a structure solved with a different ligand"),
        ("A0_null", "A2_genconf_refbox",
         "what the search buys over dropping the conformer in the box"),
    ]
    ac_cols = ["from_arm", "to_arm", "meaning", "method", "dataset",
               "from_success_2a_pct", "to_success_2a_pct", "delta_pp",
               "from_rmsd_le_2a_pct", "to_rmsd_le_2a_pct", "delta_rmsd_pp",
               "from_n", "to_n", "recorded"]
    ac_rows = []
    for a_from, a_to, meaning in comparisons:
        for (arm, method, dataset) in sorted(by):
            if arm != a_from:
                continue
            tgt = by.get((a_to, method, dataset))
            if tgt is None:
                # A4 is only defined on the crossdock dataset, so the pairing
                # with a posebusters or astex row does not exist. Compare it
                # against the posebusters A2 row instead, which is the set the
                # cross-docking queries came from.
                if a_to == "A4_crossdock":
                    tgt = by.get(("A4_crossdock", method, "crossdock"))
                if tgt is None:
                    continue
            src = by[(arm, method, dataset)]

            def d(x, y):
                try:
                    return round(float(y) - float(x), 2)
                except (TypeError, ValueError):
                    return ""
            ac_rows.append({
                "from_arm": a_from, "to_arm": a_to, "meaning": meaning,
                "method": method, "dataset": dataset,
                "from_success_2a_pct": src["success_2a_pct"],
                "to_success_2a_pct": tgt["success_2a_pct"],
                "delta_pp": d(src["success_2a_pct"], tgt["success_2a_pct"]),
                "from_rmsd_le_2a_pct": src["rmsd_le_2a_pct"],
                "to_rmsd_le_2a_pct": tgt["rmsd_le_2a_pct"],
                "delta_rmsd_pp": d(src["rmsd_le_2a_pct"], tgt["rmsd_le_2a_pct"]),
                "from_n": src["n_assessed"], "to_n": tgt["n_assessed"],
                "recorded": stamp,
            })
    L.write_tsv(results_dir / "arm_comparison.tsv", ac_cols, ac_rows)

    # -----------------------------------------------------------------------
    # 3. RMSD distributions
    # -----------------------------------------------------------------------
    rd_cols = ["arm", "method", "dataset", "quantity", "n", "min", "q25",
               "median", "q75", "q90", "max", "mean", "recorded"]
    rd_rows = []
    for (arm, method, dataset), rs in sorted(cells.items()):
        for quantity, key in (("top1_rmsd", "top1_rmsd"),
                              ("top5_best_rmsd", "top5_best_rmsd"),
                              ("top1_kabsch_rmsd", "top1_kabsch_rmsd")):
            vals = [fv(r, key) for r in rs]
            vals = [v for v in vals if v is not None]
            q = quantiles(vals)
            if not q:
                continue
            rd_rows.append({"arm": arm, "method": method, "dataset": dataset,
                            "quantity": quantity, "recorded": stamp,
                            **{k: f(v, 3) for k, v in q.items()}})
    L.write_tsv(results_dir / "rmsd_distributions.tsv", rd_cols, rd_rows)

    # -----------------------------------------------------------------------
    # 4. timing
    # -----------------------------------------------------------------------
    td_cols = ["arm", "method", "dataset", "jobs", "n", "min", "q25", "median",
               "q75", "q90", "max", "mean", "total_core_seconds", "recorded"]
    td_rows = []
    for (arm, method, dataset), rs in sorted(cells.items()):
        by_jobs = defaultdict(list)
        for r in rs:
            v = fv(r, "elapsed_s")
            if v is not None and v > 0:
                by_jobs[str(r.get("jobs", ""))].append(v)
        for jb, vals in sorted(by_jobs.items()):
            q = quantiles(vals)
            td_rows.append({"arm": arm, "method": method, "dataset": dataset,
                            "jobs": jb, "recorded": stamp,
                            "total_core_seconds": round(sum(vals), 1),
                            **{k: f(v, 2) for k, v in q.items()}})
    L.write_tsv(results_dir / "timing_distributions.tsv", td_cols, td_rows)

    # -----------------------------------------------------------------------
    # 5. PoseBusters failure waterfall
    # -----------------------------------------------------------------------
    wf_cols = ["arm", "method", "dataset", "check", "n_top1_assessed",
               "n_failed", "pct_failed", "recorded"]
    wf_rows = []
    top1 = [p for p in poses if str(p.get("pose_rank")) == "1"
            and (p.get("run_kind") or "arm") == "arm"]
    # Discovered from the data rather than hardcoded: the set of PoseBusters
    # checks has changed between releases and a fixed list would silently drop a
    # new one.
    check_names = sorted({k for p in poses for k in p
                          if k.startswith("pb_")
                          and k not in ("pb_valid", "pb_checks_failed")})
    groups = defaultdict(list)
    for p in top1:
        groups[(p["arm"], p["method"], p["dataset"])].append(p)
    for keyt, ps in sorted(groups.items()):
        assessed = [p for p in ps if p.get("pb_valid") not in ("", None)]
        if not assessed:
            continue
        for chk in check_names:
            nf = sum(1 for p in assessed if str(p.get(chk)) == "0")
            if nf == 0:
                continue
            wf_rows.append({"arm": keyt[0], "method": keyt[1], "dataset": keyt[2],
                            "check": chk, "n_top1_assessed": len(assessed),
                            "n_failed": nf, "pct_failed": pct(nf, len(assessed)),
                            "recorded": stamp})
    wf_rows.sort(key=lambda r: (r["arm"], r["method"], r["dataset"], -r["n_failed"]))
    L.write_tsv(results_dir / "pb_failure_waterfall.tsv", wf_cols, wf_rows)

    # -----------------------------------------------------------------------
    # 6. symmetry correction
    # -----------------------------------------------------------------------
    sy_cols = ["arm", "method", "dataset", "n_with_both", "n_verdict_changed",
               "pct_verdict_changed", "median_abs_difference",
               "max_abs_difference", "n_spyrmsd_disagrees_over_0p1a",
               "median_spyrmsd_disagreement", "recorded"]
    sy_rows = []
    for keyt, ps in sorted(groups.items()):
        both = [p for p in ps
                if p.get("rmsd") not in ("", None) and p.get("rmsd_first_match") not in ("", None)]
        if not both:
            continue
        diffs, changed = [], 0
        for p in both:
            a, b = float(p["rmsd"]), float(p["rmsd_first_match"])
            diffs.append(abs(a - b))
            if (a <= pass_2a) != (b <= pass_2a):
                changed += 1
        spy = [abs(float(p["rmsd_implementations_disagree_by"])) for p in ps
               if p.get("rmsd_implementations_disagree_by") not in ("", None)]
        sy_rows.append({
            "arm": keyt[0], "method": keyt[1], "dataset": keyt[2],
            "n_with_both": len(both), "n_verdict_changed": changed,
            "pct_verdict_changed": pct(changed, len(both)),
            "median_abs_difference": f(quantiles(diffs).get("median"), 4),
            "max_abs_difference": f(quantiles(diffs).get("max"), 4),
            "n_spyrmsd_disagrees_over_0p1a": sum(1 for x in spy if x > 0.1),
            "median_spyrmsd_disagreement": f(quantiles(spy).get("median"), 5),
            "recorded": stamp,
        })
    L.write_tsv(results_dir / "symmetry_effect.tsv", sy_cols, sy_rows)

    # -----------------------------------------------------------------------
    # 6b. the convergence grid that EXHAUSTIVENESS was chosen from
    # -----------------------------------------------------------------------
    cv_cols = ["method", "exhaustiveness", "n_runs", "n_complexes",
               "rmsd_le_2a_pct", "median_rmsd", "median_elapsed_s",
               "complexes_where_seeds_disagree", "pct_seeds_disagree",
               "mean_seed_range_a", "max_seed_range_a", "recorded"]
    cv_rows = []
    conv = [r for r in runs if (r.get("run_kind") or "") == "convergence"]
    if conv:
        by_m_e = defaultdict(list)
        for r in conv:
            by_m_e[(r["method"], r.get("exhaustiveness", ""))].append(r)
        for (method, exh), rs in sorted(by_m_e.items(),
                                        key=lambda kv: (kv[0][0], int(kv[0][1] or 0))):
            vals = [(r["key"], fv(r, "top1_rmsd")) for r in rs]
            vals = [(k, v) for k, v in vals if v is not None]
            if not vals:
                continue
            per_key = defaultdict(list)
            for k, v in vals:
                per_key[k].append(v)
            disagree, ranges = 0, []
            for k, vs in per_key.items():
                if len(vs) < 2:
                    continue
                if len({v <= pass_2a for v in vs}) > 1:
                    disagree += 1
                ranges.append(max(vs) - min(vs))
            times = [fv(r, "elapsed_s") for r in rs]
            times = sorted(t for t in times if t)
            k2 = sum(1 for _, v in vals if v <= pass_2a)
            cv_rows.append({
                "method": method, "exhaustiveness": exh, "n_runs": len(vals),
                "n_complexes": len(per_key),
                "rmsd_le_2a_pct": pct(k2, len(vals)),
                "median_rmsd": f(quantiles([v for _, v in vals]).get("median"), 3),
                "median_elapsed_s": f(quantiles(times).get("median"), 1) if times else "",
                "complexes_where_seeds_disagree": disagree,
                "pct_seeds_disagree": pct(disagree, max(len(ranges), 1)),
                "mean_seed_range_a": f(sum(ranges) / len(ranges), 3) if ranges else "",
                "max_seed_range_a": f(max(ranges), 3) if ranges else "",
                "recorded": stamp,
            })
    L.write_tsv(results_dir / "convergence.tsv", cv_cols, cv_rows)

    # -----------------------------------------------------------------------
    # 7. seed variance
    # -----------------------------------------------------------------------
    sv_cols = ["method", "complex_id", "arm", "n_seeds", "seeds", "rmsd_values",
               "min_rmsd", "max_rmsd", "range_rmsd", "median_rmsd", "stdev_rmsd",
               "verdicts_at_2a", "verdict_unanimous", "recorded"]
    sv_rows = []
    sv = defaultdict(list)
    for r in runs:
        if (r.get("run_kind") or "arm") == "seed_variance":
            sv[(r["method"], r["key"], r["arm"])].append(r)
    for keyt, rs in sorted(sv.items()):
        seeds = sorted({str(r["seed"]) for r in rs})
        if len(seeds) < 2:
            continue
        vals, verdicts = [], []
        for s in seeds:
            got = [fv(r, "top1_rmsd") for r in rs if str(r["seed"]) == s]
            got = [g for g in got if g is not None]
            if got:
                vals.append(got[0])
                verdicts.append(1 if got[0] <= pass_2a else 0)
        if len(vals) < 2:
            continue
        mean = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (len(vals) - 1))
        sv_rows.append({
            "method": keyt[0], "complex_id": keyt[1], "arm": keyt[2],
            "n_seeds": len(vals), "seeds": ",".join(seeds),
            "rmsd_values": ",".join(f"{v:.3f}" for v in vals),
            "min_rmsd": f(min(vals), 3), "max_rmsd": f(max(vals), 3),
            "range_rmsd": f(max(vals) - min(vals), 3),
            "median_rmsd": f(quantiles(vals)["median"], 3),
            "stdev_rmsd": f(sd, 4),
            "verdicts_at_2a": ",".join(str(v) for v in verdicts),
            "verdict_unanimous": int(len(set(verdicts)) == 1),
            "recorded": stamp,
        })
    L.write_tsv(results_dir / "seed_variance.tsv", sv_cols, sv_rows)

    # -----------------------------------------------------------------------
    # 8. consolidated exclusions. Every complex dropped anywhere, with the
    #    stage that dropped it. Silent exclusion is the most common way a
    #    benchmark number gets inflated, so this table is built from every
    #    source rather than maintained by hand.
    # -----------------------------------------------------------------------
    ex_cols = ["dataset", "key", "stage", "reason", "recorded"]
    ex_rows = []
    for src, ds_default in ((config_dir / "dataset_excluded.tsv", None),
                            (config_dir / "crossdock_excluded.tsv", "crossdock")):
        if src.is_file():
            for r in L.read_tsv(src):
                ex_rows.append({"dataset": r.get("dataset", ds_default or ""),
                                "key": r.get("complex_id", ""),
                                "stage": r.get("stage", src.stem),
                                "reason": r.get("reason", ""),
                                "recorded": r.get("recorded", stamp)})
    prep = results_dir / "preparation"
    for name, keycol, dscol in (("receptor_prep.tsv", "key", "scope"),
                                ("ligand_prep.tsv", "complex_id", "dataset"),
                                ("crossdock_align.tsv", "pair_id", None)):
        p = prep / name
        if not p.is_file():
            continue
        for r in L.read_tsv(p):
            if r.get("status") == "failed":
                ex_rows.append({"dataset": r.get(dscol, "crossdock") if dscol else "crossdock",
                                "key": r.get(keycol, ""), "stage": f"04_prepare/{p.stem}",
                                "reason": r.get("reason", ""), "recorded": stamp})
    # Docking failures come from the raw run tables, not from the scored ones.
    # 07_score_poses only scores runs that produced a pose, so a complex the
    # docking program refused never reaches run_scores.tsv and would vanish from
    # this table entirely. That is exactly the silent exclusion this file exists
    # to prevent, so the per-arm tables are read directly here.
    run_dir = results_dir / "runs"
    if run_dir.is_dir():
        for tsv in sorted(run_dir.glob("*.tsv")):
            for r in L.read_tsv(tsv):
                if r.get("status") != "failed":
                    continue
                if (r.get("run_kind") or "arm") != "arm":
                    continue
                ex_rows.append({"dataset": r.get("dataset", ""), "key": r.get("key", ""),
                                "stage": f"06_dock/{r.get('arm','')}/{r.get('method','')}",
                                "reason": r.get("reason", ""), "recorded": stamp})
    L.write_tsv(results_dir / "excluded.tsv", ex_cols, ex_rows)

    # -----------------------------------------------------------------------
    # 9. headline numbers, so the README's summary paragraph has one file to
    #    cite rather than five.
    # -----------------------------------------------------------------------
    hl = []

    def add(k, v, note=""):
        hl.append({"item": k, "value": v, "note": note})

    for method in sorted({r["method"] for r in sr_rows}):
        for arm in ("A1_refconf_refbox", "A2_genconf_refbox", "A3_genconf_detbox",
                    "A0_null", "A4_crossdock"):
            for dataset in ("posebusters", "astex", "crossdock"):
                r = by.get((arm, method, dataset))
                if not r:
                    continue
                add(f"{method}.{dataset}.{arm}.rmsd_le_2a_pct", r["rmsd_le_2a_pct"],
                    f"{r['rmsd_le_2a_n']} of {r['n_scored']} scored")
                add(f"{method}.{dataset}.{arm}.success_2a_pct", r["success_2a_pct"],
                    f"{r['success_2a_n']} of {r['n_assessed']} assessed, "
                    f"RMSD at or below {pass_2a} A and every PoseBusters check passed")
    counts = results_dir / "dataset_counts.tsv"
    if counts.is_file():
        for r in L.read_tsv(counts):
            add(f"dataset.{r['dataset']}.in_manifest", r["in_manifest"],
                f"{r['excluded_here']} excluded at stage 2")
    add("exclusions.total_rows", len(ex_rows),
        "every complex dropped at any stage; see results/excluded.tsv")
    add("thresholds.pass_angstrom", pass_2a, "the number the field quotes")
    add("thresholds.strict_angstrom", pass_1a, "reported alongside")
    L.write_tsv(results_dir / "headline.tsv", ["item", "value", "note"], hl)

    # -----------------------------------------------------------------------
    print(f"[08_analyse] {len(sr_rows)} arm/method/dataset cells")
    print(f"[08_analyse] {len(ex_rows)} exclusion rows consolidated")
    hdr = f"{'arm':20s} {'method':8s} {'dataset':12s} {'n':>5s} {'<=2A%':>7s} {'valid%':>7s} {'succ%':>7s} {'<=1A%':>7s} {'top5%':>7s}"
    print(hdr)
    print("-" * len(hdr))
    for r in sr_rows:
        print(f"{r['arm']:20s} {r['method']:8s} {r['dataset']:12s} "
              f"{r['n_scored']:>5} {str(r['rmsd_le_2a_pct']):>7} "
              f"{str(r['pb_valid_pct']):>7} {str(r['success_2a_pct']):>7} "
              f"{str(r['rmsd_le_1a_pct']):>7} {str(r['top5_rmsd_le_2a_pct']):>7}")
    if cv_rows:
        print()
        print(f"{'method':9s} {'exh':>5s} {'n':>5s} {'<=2A%':>7s} {'med_s':>7s} "
              f"{'seed_disagree':>14s} {'mean_range':>11s}")
        for r in cv_rows:
            print(f"{r['method']:9s} {str(r['exhaustiveness']):>5} {r['n_runs']:>5} "
                  f"{str(r['rmsd_le_2a_pct']):>7} {str(r['median_elapsed_s']):>7} "
                  f"{r['complexes_where_seeds_disagree']:>6}/{r['n_complexes']:<7} "
                  f"{str(r['mean_seed_range_a']):>11}")
    if sv_rows:
        print()
        for r in sv_rows:
            print(f"[08_analyse] seed variance, {r['method']} on {r['complex_id']}: "
                  f"top-1 RMSD {r['rmsd_values']} A, range {r['range_rmsd']} A, "
                  f"verdict unanimous: {'yes' if r['verdict_unanimous'] else 'no'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
