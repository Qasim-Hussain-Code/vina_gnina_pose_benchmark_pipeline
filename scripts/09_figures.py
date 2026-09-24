#!/usr/bin/env python3
"""The figures the README shows. Reads only files in results/, writes figures/.

Every figure here is produced from a table stage 8 wrote, so any number visible
in a figure is traceable to a TSV. No figure recomputes anything.

Design notes, since these are the plots a reader will judge the work by.

Colour. Three methods and five arms, so two categorical palettes taken in fixed
slot order and never cycled: blue, orange, aqua for the methods, with yellow and
magenta added for the arms. Both were checked for colour-vision separation before
use (worst adjacent CVD delta E 9.1, worst normal-vision 19.6, both above the
floor). Three of the five sit below 3:1 contrast against the surface, which
obliges a visible label on every bar rather than relying on the fill, so every
bar carries its number. That suits the brief anyway.

Form. Magnitude comparisons are horizontal bars, sorted, because the category
names are long and a reader compares lengths from a common baseline more reliably
than angles or areas. Distributions are empirical cumulative distributions rather
than box plots: the question about RMSD is where the mass sits either side of
2 Angstrom, and a box plot hides a bimodal distribution behind a median line.
No figure has two y-axes.

Theme. Light surface, fixed. These are PNG files in a README, so they cannot
follow a reader's dark mode; a transparent background would invert the axis text
into illegibility on a dark page, so the surface is painted explicitly.

Usage
-----
    python scripts/09_figures.py --config project.conf
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import PercentFormatter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib_vgb as L  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d9d8d3"
# Fixed slot order. A fourth method would take slot 4, not a generated hue.
METHOD_COLOUR = {"vina": "#2a78d6", "vinardo": "#eb6834", "gnina": "#1baf7a"}
ARM_COLOUR = {
    "A1_refconf_refbox": "#2a78d6",
    "A2_genconf_refbox": "#eb6834",
    "A3_genconf_detbox": "#1baf7a",
    "A4_crossdock": "#eda100",
    "A0_null": "#e87ba4",
}
ARM_SHORT = {
    "A0_null": "A0 no search (floor)",
    "A1_refconf_refbox": "A1 crystal conformer, ligand box",
    "A2_genconf_refbox": "A2 generated conformer, ligand box",
    "A3_genconf_detbox": "A3 generated conformer, detected box",
    "A4_crossdock": "A4 cross-docked",
}
ARM_ORDER = ["A1_refconf_refbox", "A2_genconf_refbox", "A3_genconf_detbox",
             "A4_crossdock", "A0_null"]
METHOD_ORDER = ["vina", "vinardo", "gnina"]


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": GRID,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "legend.frameon": False,
        "figure.dpi": 150,
    })


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def save(fig, path: Path, note: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"[09_figures] {path.name}  {note}")


# ---------------------------------------------------------------------------
def fig_arm_ladder(sr, out: Path, dataset="posebusters"):
    """The centrepiece: success rate per arm, grouped by method, with the floor."""
    rows = [r for r in sr if r["dataset"] in (dataset, "crossdock")]
    arms = [a for a in ARM_ORDER if any(r["arm"] == a for r in rows)]
    methods = [m for m in METHOD_ORDER if any(r["method"] == m for r in rows)]
    if not arms or not methods:
        return False

    fig, ax = plt.subplots(figsize=(8.4, 0.62 * len(arms) * len(methods) + 1.6))
    h = 0.8 / len(methods)
    ypos, ylabels = [], []
    for i, arm in enumerate(arms):
        for j, m in enumerate(methods):
            r = next((x for x in rows if x["arm"] == arm and x["method"] == m), None)
            if r is None:
                continue
            v = num(r["success_2a_pct"])
            if v is None:
                continue
            y = i + (j - (len(methods) - 1) / 2) * h
            ax.barh(y, v, height=h * 0.86, color=METHOD_COLOUR[m],
                    label=m if i == 0 else None, zorder=3)
            # Direct label on every bar: three of the five palette slots sit
            # below 3:1 against the surface, so the number carries the value and
            # the fill only carries identity.
            ax.text(v + 0.7, y, f"{v:.1f}", va="center", ha="left",
                    fontsize=8.2, color=INK_2, zorder=4)
        ypos.append(i)
        ylabels.append(ARM_SHORT.get(arm, arm))

    # The null floor as a reference line, taken from A0 rather than assumed.
    floor = [num(r["success_2a_pct"]) for r in rows if r["arm"] == "A0_null"]
    floor = [x for x in floor if x is not None]
    if floor:
        fl = max(floor)
        ax.axvline(fl, color=INK_2, lw=1.0, ls=(0, (4, 3)), zorder=2)
        # Above the top bar, horizontally. Rotated text hanging below the axis
        # made tight_layout reserve a band of empty figure as tall as the plot.
        ax.annotate(f"null floor {fl:.1f} per cent",
                    xy=(fl, -0.75), xytext=(fl + 2.0, -0.75),
                    fontsize=8, color=INK_2, va="center", ha="left",
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=INK_2, lw=0.8))

    ax.set_yticks(ypos, ylabels)
    ax.invert_yaxis()
    ax.set_xlabel("top-ranked pose within 2 A of the crystal pose and passing every "
                  "PoseBusters check (per cent)")
    ax.set_title("What each piece of crystal information is worth", loc="left", pad=10)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", ncols=len(methods), fontsize=9)
    save(fig, out, "arm ladder, success = RMSD and validity together")
    return True


def fig_rmsd_vs_validity(sr, out: Path):
    """The gap between the number the field quotes and the number it does not."""
    rows = [r for r in sr if r["arm"] in ("A1_refconf_refbox", "A2_genconf_refbox",
                                          "A3_genconf_detbox", "A4_crossdock")]
    keys = [(r["arm"], r["method"], r["dataset"]) for r in rows]
    keys = [k for k in keys if k[2] in ("posebusters", "crossdock")]
    keys.sort(key=lambda k: (ARM_ORDER.index(k[0]) if k[0] in ARM_ORDER else 9,
                             METHOD_ORDER.index(k[1]) if k[1] in METHOD_ORDER else 9))
    if not keys:
        return False

    fig, ax = plt.subplots(figsize=(8.4, 0.42 * len(keys) + 1.9))
    labels = []
    for i, k in enumerate(keys):
        r = next(x for x in rows if (x["arm"], x["method"], x["dataset"]) == k)
        a, b = num(r["rmsd_le_2a_pct"]), num(r["success_2a_pct"])
        if a is None or b is None:
            labels.append("")
            continue
        # One bar to the RMSD-only figure, a darker segment to the part that also
        # passes validity. Same hue, two steps: this is one quantity being
        # reduced, not two categories.
        ax.barh(i, a, height=0.62, color=METHOD_COLOUR[r["method"]], alpha=0.38,
                zorder=3, label="RMSD at or below 2 A" if i == 0 else None)
        ax.barh(i, b, height=0.62, color=METHOD_COLOUR[r["method"]], zorder=4,
                label="and physically valid" if i == 0 else None)
        ax.text(a + 0.7, i, f"{a:.1f}", va="center", fontsize=8, color=INK_2, zorder=5)
        if a - b > 0.05:
            ax.text(b - 0.8, i, f"{b:.1f}", va="center", ha="right", fontsize=8,
                    color=SURFACE, zorder=6)
        labels.append(f"{k[0].split('_')[0]} {k[1]}")

    ax.set_yticks(range(len(keys)), labels)
    ax.invert_yaxis()
    ax.set_xlabel("per cent of scored complexes")
    ax.set_title("Requiring the pose to be physically possible, as well as close",
                 loc="left", pad=10)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=9)
    save(fig, out, "RMSD-only against RMSD plus validity")
    return True


def fig_benchmark_gap(sr, out: Path):
    """Astex against PoseBusters: how much the choice of set moves the headline."""
    pairs = []
    for arm in ("A1_refconf_refbox", "A2_genconf_refbox", "A3_genconf_detbox"):
        for m in METHOD_ORDER:
            a = next((r for r in sr if r["arm"] == arm and r["method"] == m
                      and r["dataset"] == "astex"), None)
            p = next((r for r in sr if r["arm"] == arm and r["method"] == m
                      and r["dataset"] == "posebusters"), None)
            if a and p and num(a["success_2a_pct"]) is not None \
                    and num(p["success_2a_pct"]) is not None:
                pairs.append((arm, m, num(a["success_2a_pct"]), num(p["success_2a_pct"])))
    if not pairs:
        return False

    fig, ax = plt.subplots(figsize=(8.0, 0.46 * len(pairs) + 1.9))
    for i, (arm, m, av, pv) in enumerate(pairs):
        ax.plot([pv, av], [i, i], color=GRID, lw=2.4, zorder=2, solid_capstyle="round")
        ax.plot(pv, i, "o", ms=8, color=METHOD_COLOUR[m], zorder=4,
                label="PoseBusters Benchmark set" if i == 0 else None,
                markeredgecolor=SURFACE, markeredgewidth=1.4)
        ax.plot(av, i, "D", ms=7, color=METHOD_COLOUR[m], zorder=4,
                markerfacecolor=SURFACE, markeredgewidth=1.8,
                label="Astex Diverse set" if i == 0 else None)
        ax.text(max(av, pv) + 1.0, i, f"{av - pv:+.1f} pp", va="center",
                fontsize=8, color=INK_2)
    ax.set_yticks(range(len(pairs)),
                  [f"{a.split('_')[0]} {m}" for a, m, _, _ in pairs])
    ax.invert_yaxis()
    ax.set_xlabel("success (2 A and physically valid), per cent")
    ax.set_title("The same protocol on an older, easier set", loc="left", pad=10)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=9)
    save(fig, out, "Astex against PoseBusters, per arm and method")
    return True


def fig_rmsd_ecdf(runs, out: Path, dataset="posebusters"):
    """Where the RMSD mass actually sits, per arm."""
    series = defaultdict(list)
    for r in runs:
        if r["dataset"] not in (dataset, "crossdock"):
            continue
        v = num(r.get("top1_rmsd"))
        if v is not None:
            series[(r["arm"], r["method"])].append(v)
    if not series:
        return False

    methods = [m for m in METHOD_ORDER if any(k[1] == m for k in series)]
    fig, axes = plt.subplots(1, len(methods), figsize=(3.5 * len(methods), 3.5),
                             sharey=True)
    if len(methods) == 1:
        axes = [axes]
    for ax, m in zip(axes, methods):
        for arm in ARM_ORDER:
            vals = sorted(series.get((arm, m), []))
            if not vals:
                continue
            y = [(i + 1) / len(vals) * 100 for i in range(len(vals))]
            ax.step(vals, y, where="post", lw=2.0, color=ARM_COLOUR[arm],
                    label=arm.split("_")[0], zorder=3)
        for x, lab in ((1.0, "1 A"), (2.0, "2 A")):
            ax.axvline(x, color=INK_2, lw=0.9, ls=(0, (3, 3)), zorder=2)
            ax.text(x, 101, lab, fontsize=7.5, color=INK_2, ha="center")
        ax.set_xscale("log")
        ax.set_xlim(0.2, 40)
        ax.set_ylim(0, 100)
        ax.set_title(m, loc="left")
        ax.set_xlabel("top-1 symmetry-corrected RMSD (A)")
        ax.set_axisbelow(True)
    axes[0].set_ylabel("cumulative per cent of complexes")
    axes[0].yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    axes[-1].legend(loc="lower right", fontsize=8.2, title="arm", title_fontsize=8.2)
    fig.suptitle("A success rate is one point on this curve", x=0.01, ha="left",
                 fontsize=10.5)
    save(fig, out, "empirical cumulative distribution of top-1 RMSD")
    return True


def fig_pb_waterfall(wf, out: Path, arm="A2_genconf_refbox"):
    """Which physical checks fail, and how often, per method."""
    rows = [r for r in wf if r["arm"] == arm and r["dataset"] == "posebusters"]
    if not rows:
        rows = [r for r in wf if r["dataset"] == "posebusters"]
    if not rows:
        return False
    tot = defaultdict(float)
    for r in rows:
        v = num(r["pct_failed"])
        if v:
            tot[r["check"]] += v
    checks = [c for c, _ in sorted(tot.items(), key=lambda kv: -kv[1])][:12]
    if not checks:
        return False
    methods = [m for m in METHOD_ORDER if any(r["method"] == m for r in rows)]

    fig, ax = plt.subplots(figsize=(8.2, 0.44 * len(checks) + 1.9))
    h = 0.8 / max(len(methods), 1)
    for i, chk in enumerate(checks):
        for j, m in enumerate(methods):
            r = next((x for x in rows if x["check"] == chk and x["method"] == m), None)
            v = num(r["pct_failed"]) if r else 0.0
            v = v or 0.0
            y = i + (j - (len(methods) - 1) / 2) * h
            ax.barh(y, v, height=h * 0.86, color=METHOD_COLOUR[m],
                    label=m if i == 0 else None, zorder=3)
            if v > 0:
                ax.text(v + 0.15, y, f"{v:.1f}", va="center", fontsize=7.6,
                        color=INK_2, zorder=4)
    ax.set_yticks(range(len(checks)),
                  [c.replace("pb_", "").replace("_", " ") for c in checks])
    ax.invert_yaxis()
    ax.set_xlabel("per cent of top-ranked poses failing the check")
    ax.set_title(f"Which physical checks the top pose fails ({arm.split('_')[0]})",
                 loc="left", pad=10)
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=9)
    save(fig, out, "PoseBusters failure waterfall")
    return True


def fig_timing(runs, out: Path):
    """Per-complex wall clock as a distribution, because the mean is a lie here."""
    series = defaultdict(list)
    for r in runs:
        if r["arm"] == "A0_null":
            continue
        v = num(r.get("elapsed_s"))
        if v and v > 0:
            series[r["method"]].append(v)
    series = {k: sorted(v) for k, v in series.items() if v}
    if not series:
        return False

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for m in METHOD_ORDER:
        vals = series.get(m)
        if not vals:
            continue
        y = [(i + 1) / len(vals) * 100 for i in range(len(vals))]
        ax.step(vals, y, where="post", lw=2.2, color=METHOD_COLOUR[m], zorder=3,
                label=f"{m} (median {vals[len(vals) // 2]:.0f} s, max {vals[-1]:.0f} s)")
    ax.set_xscale("log")
    ax.set_xlabel("wall clock per complex (s, one core)")
    ax.set_ylabel("cumulative per cent of runs")
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_title("A few ligands dominate the mean", loc="left", pad=10)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8.6)
    save(fig, out, "per-complex wall clock distribution")
    return True


def fig_top1_top5(sr, out: Path):
    """Where the search found the pose and the ranking threw it away."""
    rows = [r for r in sr if r["arm"] != "A0_null"
            and num(r["rmsd_le_2a_pct"]) is not None
            and num(r["top5_rmsd_le_2a_pct"]) is not None]
    if not rows:
        return False
    fig, ax = plt.subplots(figsize=(5.6, 5.4))
    lim = 100
    ax.plot([0, lim], [0, lim], color=GRID, lw=1.2, zorder=1)
    for r in rows:
        x, y = num(r["rmsd_le_2a_pct"]), num(r["top5_rmsd_le_2a_pct"])
        ax.plot(x, y, "o", ms=9, color=METHOD_COLOUR.get(r["method"], INK_2),
                markeredgecolor=SURFACE, markeredgewidth=1.4, zorder=3)
        ax.annotate(f"{r['arm'].split('_')[0]}", (x, y), textcoords="offset points",
                    xytext=(7, -3), fontsize=7.6, color=INK_2)
    for m in METHOD_ORDER:
        if any(r["method"] == m for r in rows):
            ax.plot([], [], "o", color=METHOD_COLOUR[m], label=m, ms=8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("top-1 within 2 A (per cent)")
    ax.set_ylabel("best of top 5 within 2 A (per cent)")
    ax.set_title("Distance above the line is a ranking failure, not a search failure",
                 loc="left", pad=10, fontsize=9.6)
    ax.xaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=9)
    save(fig, out, "top-1 against top-5")
    return True


def fig_seed_variance(sv, out: Path):
    """Run-to-run spread, which almost nobody reports."""
    if not sv:
        return False
    fig, ax = plt.subplots(figsize=(7.0, 0.7 * len(sv) + 1.9))
    for i, r in enumerate(sv):
        vals = [num(x) for x in str(r["rmsd_values"]).split(",")]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        c = METHOD_COLOUR.get(r["method"], INK_2)
        ax.plot([min(vals), max(vals)], [i, i], color=c, lw=2.6, zorder=3,
                solid_capstyle="round", alpha=0.45)
        for v in vals:
            ax.plot(v, i, "o", ms=8, color=c, markeredgecolor=SURFACE,
                    markeredgewidth=1.4, zorder=4)
        ax.text(max(vals) + 0.05, i,
                f"range {max(vals) - min(vals):.2f} A over {len(vals)} seeds",
                va="center", fontsize=8, color=INK_2)
    ax.axvline(2.0, color=INK_2, lw=1.0, ls=(0, (3, 3)), zorder=2)
    ax.text(2.0, len(sv) - 0.4, " 2 A", fontsize=8, color=INK_2, va="top")
    ax.set_yticks(range(len(sv)), [f"{r['method']} on {r['complex_id']}" for r in sv])
    ax.invert_yaxis()
    ax.set_xlabel("top-1 symmetry-corrected RMSD (A), one complex, repeat seeds")
    ax.set_title("The same complex, the same protocol, a different seed",
                 loc="left", pad=10)
    ax.set_axisbelow(True)
    ax.grid(axis="y", visible=False)
    save(fig, out, "seed variance")
    return True


def fig_box_detection(boxes, out: Path):
    """Whether pocket detection finds the site at all, which bounds arm A3."""
    det = [r for r in boxes if r.get("arm_box_kind") == "detected"
           and r.get("status") == "ok"
           and num(r.get("distance_to_crystal_centroid")) is not None]
    if not det:
        return False
    vals = sorted(num(r["distance_to_crystal_centroid"]) for r in det)
    inside = sum(1 for r in det if str(r.get("ligand_fits_in_box")) == "1")
    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    y = [(i + 1) / len(vals) * 100 for i in range(len(vals))]
    ax.step(vals, y, where="post", lw=2.2, color=ARM_COLOUR["A3_genconf_detbox"],
            zorder=3)
    # Log x, because the distribution has most of its mass under 30 Angstrom and
    # a tail out past 100 on the receptors where fpocket ranked a surface groove
    # on the far side of the protein first. A linear axis spends four fifths of
    # its width on that tail.
    ax.set_xscale("log")
    half = 12.5
    ax.axvline(half, color=INK_2, lw=1.0, ls=(0, (3, 3)), zorder=2)
    ax.text(half * 1.06, 52, "half the box width\nbeyond here the ligand\n"
            "cannot be inside the box", fontsize=7.8, color=INK_2,
            ha="left", va="center")
    ax.set_xlabel("distance from fpocket rank 1 centre to the crystal ligand centroid (A)")
    ax.set_ylabel("cumulative per cent of receptors")
    ax.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    ax.set_title(f"The crystal ligand lies inside the detected box in "
                 f"{100.0 * inside / len(det):.0f} per cent of receptors "
                 f"({inside} of {len(det)})", loc="left", pad=10, fontsize=9.6)
    ax.set_axisbelow(True)
    save(fig, out, "pocket detection against the true site")
    return True


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    args = ap.parse_args()

    conf = L.load_conf(args.config)
    results = Path(conf["RESULTS_DIR"])
    figures = Path(conf["FIGURES_DIR"])
    config_dir = Path(conf["CONFIG_DIR"])
    figures.mkdir(parents=True, exist_ok=True)
    style()

    def tsv(p):
        return L.read_tsv(p) if Path(p).is_file() else []

    sr = tsv(results / "success_rates.tsv")
    runs = tsv(results / "scored" / "run_scores.tsv")
    wf = tsv(results / "pb_failure_waterfall.tsv")
    sv = tsv(results / "seed_variance.tsv")
    boxes = tsv(config_dir / "boxes.tsv")
    if not sr:
        L.eprint("[error] results/success_rates.tsv missing; run 08_analyse.py first")
        return 1

    made = 0
    made += fig_arm_ladder(sr, figures / "fig1_arm_ladder.png")
    made += fig_rmsd_vs_validity(sr, figures / "fig2_rmsd_vs_validity.png")
    made += fig_benchmark_gap(sr, figures / "fig3_benchmark_gap.png")
    made += fig_rmsd_ecdf(runs, figures / "fig4_rmsd_ecdf.png")
    made += fig_pb_waterfall(wf, figures / "fig5_pb_failure_waterfall.png")
    made += fig_top1_top5(sr, figures / "fig6_top1_vs_top5.png")
    made += fig_timing(runs, figures / "fig7_timing.png")
    made += fig_seed_variance(sv, figures / "fig8_seed_variance.png")
    made += fig_box_detection(boxes, figures / "fig9_pocket_detection.png")
    print(f"[09_figures] wrote {made} figures to {figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
