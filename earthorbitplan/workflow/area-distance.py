from pathlib import Path

import numpy as np
import synphot
from astropy import units as u
from astropy.coordinates import ICRS, EarthLocation
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable
from astropy.time import Time
from astropy.visualization import quantity_support
from astropy_healpix import HEALPix
from m4opt import missions
from m4opt.synphot import observing
from matplotlib import patheffects
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats


def customize_style(columns=1):
    if columns == 1:
        target_width = 3.5  # ApJ column size in inches
    else:
        target_width = 7.25  # ApJ two-column text width in inches
    width, height = plt.rcParams["figure.figsize"]
    plt.style.use("seaborn-v0_8-paper")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["figure.figsize"] = (target_width, height * target_width / width)


quantity_support()

customize_style()

main_table = QTable.read("../data/events.ecsv")


plan_args = QTable.read(
    Path("data") / main_table[0]["run"] / f"{main_table[0]['coinc_event_id']}.ecsv"
).meta["args"]
plan_args.pop("skymap")
skymap_area_cl = 90
hpx = HEALPix(nside=plan_args["nside"], frame=ICRS(), order="nested")
mission = getattr(missions, plan_args["mission"])
cutoff = plan_args["cutoff"]

with observing(
    observer_location=EarthLocation(0 * u.m, 0 * u.m, 0 * u.m),
    target_coord=hpx.healpix_to_skycoord(np.arange(hpx.npix)),
    obstime=Time("2025-01-01"),
):
    limmag = mission.detector.get_limmag(
        plan_args["snr"],
        min(plan_args["deadline"] - plan_args["delay"], plan_args["exptime_max"]),
        synphot.SourceSpectrum(synphot.ConstFlux1D, amplitude=0 * u.ABmag),
        plan_args["bandpass"],
    ).max()

skymap_area_cl = 90
min_area = (mission.fov.width * mission.fov.height).to(u.deg**2)

chisq_ppf = stats.chi2(df=2).ppf
area_factor = chisq_ppf(skymap_area_cl / 100) / chisq_ppf(plan_args["cutoff"])
max_area = (
    area_factor
    * min_area
    * (plan_args["deadline"] - plan_args["delay"])
    / (plan_args["visits"] * plan_args["exptime_min"])
).to(u.deg**2)
max_distance = (
    10
    ** (
        0.2
        * (
            limmag.to_value(u.mag)
            - (
                plan_args["absmag_mean"]
                - stats.norm.ppf(1 - plan_args["cutoff"]) * plan_args["absmag_stdev"]
            )
            - 25
        )
    )
    * u.Mpc
)
crossover_distance = (
    max_distance * (min_area / max_area).to_value(u.dimensionless_unscaled) ** 0.25
)

for run in ["O5", "O6"]:
    table = main_table[main_table["run"] == run]
    selected_table = table[table["objective_value"] >= cutoff]

    z = z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc).to_value(
        u.dimensionless_unscaled
    )
    m2 = table["mass2"] / (1 + z)
    table = table[(m2 <= 3)]

    width, height = plt.rcParams["figure.figsize"]
    default_fig_width_height_ratio = width / height
    fig_width_height_ratio = 0.7
    width = 7
    fig = plt.figure(figsize=(width, width * fig_width_height_ratio))

    left, bottom, width, height = (
        0.2,
        0.065,
        0.475,
        0.475 / default_fig_width_height_ratio,
    )
    depth = 0.2
    sep = 0.025
    ax_joint = fig.add_subplot(
        (left, bottom / fig_width_height_ratio, width, height / fig_width_height_ratio),
        aspect=0.25,
        xlim=(50 * u.Mpc, 4000 * u.Mpc),
        ylim=(
            (
                10 ** (np.log10(50 / 4000) * 4 / default_fig_width_height_ratio)
                * u.spat
            ).to(u.deg**2),
            (1 * u.spat).to(u.deg**2),
        ),
        xscale="log",
        yscale="log",
    )
    ax_joint.set_xlabel(r"Luminosity distance, $d_\mathrm{L}$ (Mpc)")
    ax_joint.set_ylabel(
        rf"{skymap_area_cl}% credible area, $A_{{{skymap_area_cl}\%}}$ (deg$^2$)"
    )
    ax_x = fig.add_subplot(
        (
            left,
            (bottom + height + sep) / fig_width_height_ratio,
            width,
            depth / fig_width_height_ratio,
        ),
        sharex=ax_joint,
    )
    ax_y = fig.add_subplot(
        (
            left + width + sep,
            bottom / fig_width_height_ratio,
            depth,
            height / fig_width_height_ratio,
        ),
        sharey=ax_joint,
    )

    ax_joint.annotate(
        "Relative\nfrequency",
        (
            left + width + sep,
            (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        ),
        (
            left + width + sep + 0.5 * depth,
            (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        ),
        xycoords="figure fraction",
        textcoords="figure fraction",
        ha="center",
        va="center",
        arrowprops=dict(
            facecolor="k",
            edgecolor="none",
            linewidth=0,
            arrowstyle="simple",
        ),
        fontsize=plt.rcParams["axes.labelsize"],
    )
    ax_joint.annotate(
        "Relative\nfrequency",
        (
            left + width + sep + 0.5 * depth,
            (bottom + height + sep) / fig_width_height_ratio,
        ),
        (
            left + width + sep + 0.5 * depth,
            (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        ),
        xycoords="figure fraction",
        textcoords="figure fraction",
        ha="center",
        va="center",
        arrowprops=dict(
            facecolor="k",
            edgecolor="none",
            linewidth=0,
            arrowstyle="simple",
        ),
        fontsize=plt.rcParams["axes.labelsize"],
        color="none",
    )

    ax_joint.fill_between(
        [
            ax_joint.get_xlim()[0] * u.Mpc,
            crossover_distance,
            max_distance,
            max_distance,
            ax_joint.get_xlim()[1] * u.Mpc,
        ],
        [
            max_area,
            max_area,
            min_area,
            ax_joint.get_ylim()[0] * u.deg**2,
            ax_joint.get_ylim()[0] * u.deg**2,
        ],
        [ax_joint.get_ylim()[1] * u.deg**2] * 5,
        color="gainsboro",
    )

    color = "tab:blue"
    ax_joint.plot(
        u.Quantity(
            [
                ax_joint.get_xlim()[0] * u.Mpc,
                crossover_distance,
                max_distance,
                max_distance,
            ]
        ),
        u.Quantity([max_area, max_area, min_area, ax_joint.get_ylim()[0] * u.deg**2]),
        color=color,
    )
    kwargs = dict(
        color=color,
        ha="center",
        va="bottom",
        rotation_mode="anchor",
        linespacing=0.1,
        path_effects=[patheffects.withStroke(linewidth=2, foreground="white")],
        fontsize=plt.rcParams["ytick.labelsize"],
    )
    ax_joint.text(
        np.sqrt(ax_joint.get_xlim()[0] * u.Mpc * crossover_distance),
        max_area,
        "Max area\n",
        **kwargs,
    )
    ax_joint.text(
        max_distance,
        np.sqrt(ax_joint.get_ylim()[0] * u.deg**2 * min_area),
        "Max distance\n",
        rotation=-90,
        **kwargs,
    )
    ax_joint.text(
        np.sqrt(crossover_distance * max_distance),
        np.sqrt(min_area * max_area),
        "Area $\propto$ distance$^{-4}$\n",
        rotation=-45,
        **kwargs,
    )

    ax_joint.scatter(
        "distance",
        f"area({skymap_area_cl})",
        s=2,
        facecolor="silver",
        edgecolor="none",
        data=table,
    )

    cmap = plt.get_cmap("cool")
    cmap = LinearSegmentedColormap.from_list(
        "truncated_cool", cmap(np.linspace(1 / 3, 1))
    )
    marker_scale = 30
    scatter = ax_joint.scatter(
        "distance",
        f"area({skymap_area_cl})",
        s=selected_table["objective_value"] * marker_scale,
        c=selected_table["detection_probability_known_position"],
        cmap=cmap,
        vmin=0,
        vmax=1,
        data=selected_table,
    )

    ax_legend = fig.add_subplot(
        (
            0.075,
            (bottom + height + sep + 0.5 * depth - 0.08) / fig_width_height_ratio,
            0.08,
            0.08 / fig_width_height_ratio,
        ),
        aspect=1,
    )
    ax_legend.margins(0.2)
    dx = 1 / 4
    x = np.arange(0, 1 + dx, dx)
    ax_legend.set_xticks([0, 0.5, 1])
    ax_legend.set_yticks([0, 0.5, 1])
    x, y = (a.ravel() for a in np.meshgrid(x, x))
    ax_legend.scatter(
        x,
        y,
        s=2,
        facecolor="silver",
        edgecolor="none",
        clip_on=False,
    )
    ax_legend.scatter(
        x, y, s=marker_scale * x, c=y, vmin=0, vmax=1, cmap=cmap, clip_on=False
    )
    twin = ax_legend.twiny()
    twin.set_frame_on(False)
    twin.set_ylim(ax_legend.get_ylim())
    twin.set_xticks([])
    twin.set_xlabel("Objective val.")
    ax_legend.set_ylabel("Detection prob.")
    bbox_cbar, bbox_joint = [
        ax.get_window_extent().transformed(fig.transFigure.inverted())
        for ax in [ax_legend, ax_joint]
    ]
    ax_legend.annotate(
        "",
        (bbox_joint.xmin, bbox_joint.ymax),
        (bbox_cbar.xmax, bbox_cbar.ymin),
        "figure fraction",
        "figure fraction",
        clip_on=False,
        arrowprops=dict(
            facecolor="k",
            edgecolor="none",
            linewidth=0,
            arrowstyle="simple",
            shrinkA=6,
            shrinkB=6,
        ),
    )
    ax_legend.spines[["top", "right"]].set_visible(False)

    ticks = [
        np.quantile(
            selected_table["distance"],
            [0.9],
            method="inverted_cdf",
        ).item(),
        np.quantile(
            selected_table["distance"],
            [0.9],
            weights=selected_table["detection_probability_known_position"],
            method="inverted_cdf",
        ).item(),
    ]
    bins = np.linspace(*np.log(ax_joint.get_xlim()), 16)
    color = "silver"
    values, _ = np.histogram(np.log(table["distance"]), bins=bins)
    ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=0)
    ax_x.set_ylim(0, values.max() * 1.1)
    color = cmap(0)
    values, _ = np.histogram(
        np.log(selected_table["distance"]),
        bins=bins,
    )
    ax_x.axvline(
        ticks[0],
        linewidth=plt.rcParams["xtick.major.width"],
        color=plt.rcParams["xtick.color"],
        zorder=1,
    )
    ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=2)
    color = cmap(np.inf)
    values, _ = np.histogram(
        np.log(selected_table["distance"]),
        weights=selected_table["detection_probability_known_position"],
        bins=bins,
    )
    ax_x.axvline(
        ticks[1],
        linewidth=plt.rcParams["xtick.major.width"],
        color=plt.rcParams["xtick.color"],
        zorder=3,
    )
    ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=4)
    twin = ax_x.twiny()
    twin.set_frame_on(False)
    twin.set_xlim(*ax_joint.get_xlim())
    twin.set_xscale(ax_joint.get_xscale())
    twin.set_xticks(
        [*ticks, np.prod(np.asarray(ax_joint.get_xlim()) ** [0.6, 0.4])],
        [*(f"{np.round(tick):g}\nMpc" for tick in ticks), "90th\npercentile"],
    )
    twin.xaxis.minorticks_off()
    tick = twin.xaxis.get_major_ticks()[2]
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)
    tick = twin.xaxis.get_major_ticks()[0]
    tick.label2.set_ha("left")
    tick = twin.xaxis.get_major_ticks()[1]
    tick.label2.set_ha("right")
    tick = twin.xaxis.get_major_ticks()[2]
    tick.label2.set_ha("right")

    ticks = [
        np.quantile(
            selected_table[f"area({skymap_area_cl})"],
            [0.9],
            method="inverted_cdf",
        ).item(),
        np.quantile(
            selected_table[f"area({skymap_area_cl})"],
            [0.9],
            weights=selected_table["detection_probability_known_position"],
            method="inverted_cdf",
        ).item(),
    ]
    bins = np.linspace(*np.log(ax_joint.get_ylim()), 16)
    color = "silver"
    values, _ = np.histogram(np.log(table[f"area({skymap_area_cl})"]), bins=bins)
    ax_y.stairs(
        values, np.exp(bins), color=color, fill=True, orientation="horizontal", zorder=0
    )
    ax_y.set_xlim(0, values.max() * 1.1)
    color = cmap(0)
    values, _ = np.histogram(
        np.log(selected_table[f"area({skymap_area_cl})"]),
        bins=bins,
    )
    ax_y.axhline(
        ticks[0],
        linewidth=plt.rcParams["ytick.major.width"],
        color=plt.rcParams["ytick.color"],
        zorder=1,
    )
    artist = ax_y.stairs(
        values, np.exp(bins), color=color, fill=True, orientation="horizontal", zorder=2
    )
    color = cmap(np.inf)
    values, _ = np.histogram(
        np.log(selected_table[f"area({skymap_area_cl})"]),
        weights=selected_table["detection_probability_known_position"],
        bins=bins,
    )
    ax_y.axhline(
        ticks[1],
        linewidth=plt.rcParams["ytick.major.width"],
        color=plt.rcParams["ytick.color"],
        zorder=3,
    )
    ax_y.stairs(
        values, np.exp(bins), color=color, fill=True, orientation="horizontal", zorder=4
    )
    twin = ax_y.twinx()
    twin.set_frame_on(False)
    twin.set_ylim(*ax_joint.get_ylim())
    twin.set_yscale(ax_joint.get_yscale())
    twin.set_yticks(
        [*ticks, np.prod(np.asarray(ax_joint.get_ylim()) ** [0.1, 0.9])],
        [*(f"{np.round(tick):g} deg$^2$" for tick in ticks), "90th\npercentile"],
    )
    twin.yaxis.minorticks_off()
    tick = twin.yaxis.get_major_ticks()[0]
    tick.label2.set_va("bottom")
    tick = twin.yaxis.get_major_ticks()[1]
    tick.label2.set_va("top")
    tick = twin.yaxis.get_major_ticks()[2]
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)

    kwargs = dict(
        # color="black",
        transform=ax_x.transAxes,
        fontsize=plt.rcParams["legend.fontsize"],
        zorder=5,
        ha="left",
        va="top",
    )
    ax_x.text(0.05, 0.9, "All events", color="dimgray", **kwargs)
    ax_x.text(0.05, 0.78, "Triggered", color="tab:blue", **kwargs)
    ax_x.text(0.05, 0.66, "Detected", color="magenta", **kwargs)

    plt.setp(ax_x.get_xticklabels() + ax_y.get_yticklabels(), visible=False)
    ax_x.set_yticks([])
    ax_y.set_xticks([])

    fig.text(
        0.9,
        0.9,
        run,
        fontsize="x-large",
        zorder=1000,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "boxstyle": "round"},
    )
    fig.savefig(f"figures/area-distance-{run}.pdf")
    plt.close(fig)
