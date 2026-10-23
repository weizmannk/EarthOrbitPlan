import sys

sys.path.insert(0, sys.argv[1])
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

from paper_results import load

D, out = sys.argv[2], sys.argv[3]
us, uv = load(f"{D}/lcs.ecsv"), load(f"{D}/uvex.ecsv")
assert (us["key"] == uv["key"]).all()
w = us["weight"]
det = lambda m, a: (w * m["pdet"])[m["selected"] & a].sum()

mpl.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "font.serif": ["STIXGeneral", "DejaVu Serif"],
        "font.size": 8,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.edgecolor": "#55554f",
        "xtick.color": "#3d3d3a",
        "ytick.color": "#3d3d3a",
    }
)
COL = {"ULTRASAT": "#2a78d6", "UVEX": "#eb6834"}
panels = [
    (
        "area",
        [0, 100, 500, 1000, 2000, 5000, np.inf],
        ["<100", "100–500", "500–1k", "1k–2k", "2k–5k", ">5k"],
        r"90% credible area (deg$^2$)",
    ),
    (
        "distance",
        [0, 200, 300, 400, 500, 700, np.inf],
        ["<200", "200–300", "300–400", "400–500", "500–700", ">700"],
        "Luminosity distance (Mpc)",
    ),
]
fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.3))
rows = []
for ax, (key, edges, labels, xlabel), tag in zip(axes, panels, "ab"):
    x = np.arange(len(labels))
    gw = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        a = (us[key] >= lo) & (us[key] < hi)
        gw.append(w[a].sum())
        for j, (name, m) in enumerate((("ULTRASAT", us), ("UVEX", uv))):
            v = det(m, a)
            ax.bar(
                x[i] + (j - 0.5) * 0.38,
                v,
                width=0.36,
                color=COL[name],
                edgecolor="white",
                linewidth=0.8,
                label=name if i == 0 else None,
                zorder=3,
            )
        rows.append(
            (
                key,
                labels[i],
                round(gw[-1], 1),
                round(det(us, a), 2),
                round(det(uv, a), 2),
            )
        )
    ax.set_xticks(x, labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Detected KNe over O5")
    ax.grid(axis="y", color="#e6e6e3", linewidth=0.6, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2.5)
    top = ax.get_ylim()[1]
    for xi, g in zip(x, gw):
        ax.text(
            xi,
            top * 1.0,
            f"{g:.0f}",
            ha="center",
            va="bottom",
            fontsize=6.5,
            color="#6e6e69",
        )
    ax.text(
        -0.75, top * 1.0, "GW:", ha="right", va="bottom", fontsize=6.5, color="#6e6e69"
    )
    ax.set_ylim(0, top * 1.12)
    ax.text(0.99, 0.80, f"({tag})", transform=ax.transAxes, ha="right", fontsize=8)
axes[0].legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.78), fontsize=7)
fig.tight_layout(h_pad=1.2)
fig.savefig(f"{out}/detections_by_area_distance.pdf")
fig.savefig(f"{out}/detections_by_area_distance.png", dpi=200)
for r in rows:
    print(r)
