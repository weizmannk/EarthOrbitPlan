import warnings

import numpy as np
import synphot
from astropy import units as u
from astropy.coordinates import ICRS, EarthLocation
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

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
warnings.filterwarnings("ignore", ".*dubious year.*")
warnings.filterwarnings(
    "ignore", "Tried to get polar motions for times after IERS data is valid.*"
)


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


def load_events(events_file):
    """
    Load and prepare area-distance data from an ECSV file.

    Parameters
    ----------
    events_file : str or Path
        Path to the ECSV file containing the event data.

    Returns
    -------
    main_table : QTable
        The full event table.
    constants : dict
        Dictionary of constant simulation parameters.
    """
    main_table = QTable.read(events_file)

    # Define columns that must remain constant
    constant_cols = [
        "cutoff",
        "mission",
        "nside",
        "snr",
        "deadline",
        "delay",
        "exptime_min",
        "exptime_max",
        "bandpass",
        "absmag_mean",
        "absmag_stdev",
        "visits",
        "skydir",
    ]
    constants = {}
    for col in constant_cols:
        if col in main_table.colnames:
            if not np.all(main_table[col] == main_table[0][col]):
                raise ValueError(
                    f"Column '{col}' has multiple values; expected a constant."
                )
            constants[col] = main_table[0][col]

    return main_table, constants


def preprocess_events(main_table, constants, ns_max=3.0):
    r"""
    Preprocess the event table and select events for each run.

    This function processes the simulated events table, grouping by run and
    selecting events with an objective value above the cutoff threshold.
    It is intended for use with simulations that only include binary neutron
    star (BNS) and neutron star–black hole (NSBH) events, i.e., systems with
    secondary mass below a specified upper limit (default: 3.0, i.e., neutron star upper bound, \( M_\odot \)).
    This mass cut is consistent with the original simulation setup, and is applied here for completeness.

    This function selects only BNS and NSBH systems (i.e., systems with secondary mass <= ns_max),
    and provides both the full filtered sample and the sub-sample passing the objective_value cutoff.

    Parameters
    ----------
    main_table : `~astropy.table.QTable`
        The full event table as returned by :func:`load_events`.
    constants : dict
        Dictionary of constant simulation parameters as returned by :func:`load_events`.
    ns_max : float, optional
        Maximum secondary mass (in solar masses, default: 3.0 :math:`M_\odot`).

    tables : dict
        {run_name: table of BNS/NSBH events}
    selected_tables : dict
        {run_name: table of BNS/NSBH events with objective_value >= cutoff}
        and objective_value >= cutoff.
    """

    runs = np.unique(main_table["run"])
    main_table_dict = {}
    selected_tables = {}
    for run in runs:
        table = main_table[main_table["run"] == run]
        # z = z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc).to_value(
        #     u.dimensionless_unscaled
        # )
        # m2 = table["mass2"] / (1 + z)

        # FIXME: This is redundant because the table should only contain BNS/NSBH events.
        # However, we keep this check for completeness and to ensure compatibility with other datasets.
        # Note: The input table should already be filtered to include only BNS/NSBH events during the simulation setup
        # when reading the zip file and selecting the relevant events.
        # Ideally selected_tables[run] = table[table["objective_value"] >= constants["cutoff"]] and  main_table_dict[run] = table

        # Filter: secondary mass <= ns_max
        # keep = m2 <= ns_max
        # main_table_dict[run] = table[(m2 <= 3)]
        main_table_dict[run] = table
        selected_tables[run] = table[table["objective_value"] >= constants["cutoff"]]

        # Filter objective_value >= cutoff (and mass)
        # selected_tables[run] = table[
        #     keep & (table["objective_value"] >= constants["cutoff"])
        # ]
    return (
        runs,
        main_table_dict,
        selected_tables,
    )


def area_distance_limits(constants, skymap_area_cl=90):
    """
    Compute min/max area and limiting distances for a survey.

    Parameters
    ----------
    constants : dict
        Survey configuration and simulation constantsas returned by :func:`load_events`.
    skymap_area_cl : float, optional
        Percent credible level for the area (default: 90).

    Returns
    -------
    min_area : `~astropy.units.Quantity`
        Minimum area covered (telescope field of view, FoV, deg^2).
    max_area : `~astropy.units.Quantity`
        Maximum credible sky area  that can be covered given the observing constraints (deg^2).
    max_distance : `~astropy.units.Quantity`
        Maximum luminosity distance (Mpc) at which a source can be detected (given limiting magnitude).
    crossover_distance : `~astropy.units.Quantity`
        Crossover distance for area–distance scaling (Mpc).
    """
    # Get Healpix grid for sky coordinates
    hpx = HEALPix(nside=constants["nside"], frame=ICRS(), order="nested")
    mission = getattr(missions, constants["mission"])

    # Compute the limiting magnitude for the configuration
    with observing(
        observer_location=EarthLocation(0 * u.m, 0 * u.m, 0 * u.m),
        target_coord=hpx.healpix_to_skycoord(np.arange(hpx.npix)),
        obstime=Time("2025-01-01"),
    ):
        limmag = mission.detector.get_limmag(
            constants["snr"],
            min(constants["deadline"] - constants["delay"], constants["exptime_max"]),
            synphot.SourceSpectrum(synphot.ConstFlux1D, amplitude=0 * u.ABmag),
            constants["bandpass"],
        ).max()

    # Minimum area: single telescope field of view (FoV)
    min_area = (mission.fov.width * mission.fov.height).to(u.deg**2)

    # Area factor using chi-square inverse cumulative distribution
    chisq_ppf = stats.chi2(df=2).ppf
    area_factor = chisq_ppf(skymap_area_cl / 100) / chisq_ppf(constants["cutoff"])

    # Maximum credible area that can be covered within observing constraints
    max_area = (
        area_factor
        * min_area
        * (constants["deadline"] - constants["delay"])
        / (constants["visits"] * constants["exptime_min"])
    ).to(u.deg**2)

    # Maximum distance (depth) for detection at limiting magnitude
    max_distance = (
        10
        ** (
            0.2
            * (
                limmag.to_value(u.mag)
                - (
                    constants["absmag_mean"]
                    - stats.norm.ppf(1 - constants["cutoff"])
                    * constants["absmag_stdev"]
                )
                - 25
            )
        )
        * u.Mpc
    )

    # Crossover distance: transition from FoV-limited to depth-limited regime
    # At distances below, area is FoV-limited; above, area is depth-limited.
    crossover_distance = (
        max_distance * (min_area / max_area).to_value(u.dimensionless_unscaled) ** 0.25
    )

    print(crossover_distance, max_distance, min_area, max_area, limmag)
    return min_area, max_area, max_distance, crossover_distance


def plot_area_distance(
    runs,
    tables,
    selected_tables,
    min_area,
    max_area,
    max_distance,
    crossover_distance,
    skymap_area_cl=90,
):
    """
    Plot 2D scatter of credible area vs luminosity distance for each run,
    with marginal histograms.

    Parameters
    ----------
    runs : list
        List of run names/labels to process.
    tables : dict
        Dictionary of all event tables, keyed by run as process in :func:`preprocess_event`.
    selected_tables : dict
        Dictionary of selected event tables (e.g. triggered/detected), keyed by run as created by :func:`preprocess_events`.
    min_area : float
        Minimum credible area (deg^2) to plot threshold, from :func:`area_distance_limits`.
    max_area : float
        Maximum credible area (deg^2) to plot threshold, from :func:`area_distance_limits`.
    max_distance : float
        Maximum distance (Mpc) to plot threshold, from :func:`area_distance_limits`.
    crossover_distance : float
        Distance at which area threshold transitions, from :func:`area_distance_limits`.
    skymap_area_cl : int, optional
        Credible level for sky localization area (default is 90).
    """

    quantity_support()

    customize_style()

    for run in runs:
        # === Data preparation ===
        table = tables[run]
        selected_table = selected_tables[run]

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
            (
                left,
                bottom / fig_width_height_ratio,
                width,
                height / fig_width_height_ratio,
            ),
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
            u.Quantity(
                [max_area, max_area, min_area, ax_joint.get_ylim()[0] * u.deg**2]
            ),
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
        ax_joint.scatter(
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
            values,
            np.exp(bins),
            color=color,
            fill=True,
            orientation="horizontal",
            zorder=0,
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
        ax_y.stairs(
            values,
            np.exp(bins),
            color=color,
            fill=True,
            orientation="horizontal",
            zorder=2,
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
            values,
            np.exp(bins),
            color=color,
            fill=True,
            orientation="horizontal",
            zorder=4,
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

        # === Figure and axes layout ===

        # width, height = plt.rcParams["figure.figsize"]
        # default_fig_width_height_ratio = width / height
        # fig_width_height_ratio = 0.7
        # width = 7
        # fig = plt.figure(figsize=(width, width * fig_width_height_ratio))

        # # Axes geometry
        # left, bottom, width, height = 0.2, 0.065, 0.475, 0.475 /  default_fig_width_height_ratio
        # depth, sep = 0.2, 0.025

        # # Main scatter axis (joint)
        # ax_joint = fig.add_subplot(
        #     (
        #         left,
        #         bottom /  default_fig_width_height_ratio,
        #         width,
        #         height /  default_fig_width_height_ratio,
        #     ),
        #     aspect=0.25,
        #     xlim=(50 * u.Mpc, 4000 * u.Mpc),
        #     ylim=(
        #         (
        #             10 ** (np.log10(50 / 4000) * 4 /  default_fig_width_height_ratio)
        #             * u.spat
        #         ).to(u.deg**2),
        #         (1 * u.spat).to(u.deg**2),
        #     ),
        #     xscale="log",
        #     yscale="log",
        # )
        # ax_joint.set_xlabel(r"Luminosity distance, $d_\mathrm{L}$ (Mpc)")
        # ax_joint.set_ylabel(
        #     rf"{skymap_area_cl}% credible area, $A_{{{skymap_area_cl}\%}}$ (deg$^2$)"
        # )

        # # Top histogram (distance)
        # ax_x = fig.add_subplot(
        #     (
        #         left,
        #         (bottom + height + sep) / fig_width_height_ratio,
        #         width,
        #         depth / fig_width_height_ratio,
        #     ),
        #     sharex=ax_joint,
        # )

        # # Right histogram (area)
        # ax_y = fig.add_subplot(
        #     (
        #         left + width + sep,
        #         bottom / fig_width_height_ratio,
        #         depth,
        #         height / fig_width_height_ratio,
        #     ),
        #     sharey=ax_joint,
        # )

        # # === Annotate marginal axes (arrows, labels) ===

        # ax_joint.annotate(
        #     "Relative\nfrequency",
        #     (
        #         left + width + sep,
        #         (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        #     ),
        #     (
        #         left + width + sep + 0.5 * depth,
        #         (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        #     ),
        #     xycoords="figure fraction",
        #     textcoords="figure fraction",
        #     ha="center",
        #     va="center",
        #     arrowprops=dict(
        #         facecolor="k",
        #         edgecolor="none",
        #         linewidth=0,
        #         arrowstyle="simple",
        #     ),
        #     fontsize=plt.rcParams["axes.labelsize"],
        # )
        # ax_joint.annotate(
        #     "Relative\nfrequency",
        #     (
        #         left + width + sep + 0.5 * depth,
        #         (bottom + height + sep) / fig_width_height_ratio,
        #     ),
        #     (
        #         left + width + sep + 0.5 * depth,
        #         (bottom + height + sep + 0.5 * depth) / fig_width_height_ratio,
        #     ),
        #     xycoords="figure fraction",
        #     textcoords="figure fraction",
        #     ha="center",
        #     va="center",
        #     arrowprops=dict(
        #         facecolor="k",
        #         edgecolor="none",
        #         linewidth=0,
        #         arrowstyle="simple",
        #     ),
        #     fontsize=plt.rcParams["axes.labelsize"],
        #     color="none",
        # )

        # # === Plot shaded regions for area/distance thresholds ===
        # ax_joint.fill_between(
        #     [
        #         ax_joint.get_xlim()[0] * u.Mpc,
        #         crossover_distance,
        #         max_distance,
        #         max_distance,
        #         ax_joint.get_xlim()[1] * u.Mpc,
        #     ],
        #     [
        #         max_area,
        #         max_area,
        #         min_area,
        #         ax_joint.get_ylim()[0] * u.deg**2,
        #         ax_joint.get_ylim()[0] * u.deg**2,
        #     ],
        #     [ax_joint.get_ylim()[1] * u.deg**2] * 5,
        #     color="gainsboro",
        # )

        # # === Add "max area" and "max distance" boundaries ===
        # color = "tab:blue"
        # ax_joint.plot(
        #     u.Quantity(
        #         [
        #             ax_joint.get_xlim()[0] * u.Mpc,
        #             crossover_distance,
        #             max_distance,
        #             max_distance,
        #         ]
        #     ),
        #     u.Quantity(
        #         [max_area, max_area, min_area, ax_joint.get_ylim()[0] * u.deg**2]
        #     ),
        #     color=color,
        # )

        # kwargs = dict(
        #     color=color,
        #     ha="center",
        #     va="bottom",
        #     rotation_mode="anchor",
        #     linespacing=0.1,
        #     path_effects=[patheffects.withStroke(linewidth=2, foreground="white")],
        #     fontsize=plt.rcParams["ytick.labelsize"],
        # )
        # ax_joint.text(
        #     np.sqrt(ax_joint.get_xlim()[0] * u.Mpc * crossover_distance),
        #     max_area,
        #     "Max area\n",
        #     **kwargs,
        # )
        # ax_joint.text(
        #     max_distance,
        #     np.sqrt(ax_joint.get_ylim()[0] * u.deg**2 * min_area),
        #     "Max distance\n",
        #     rotation=-90,
        #     **kwargs,
        # )
        # ax_joint.text(
        #     np.sqrt(crossover_distance * max_distance),
        #     np.sqrt(min_area * max_area),
        #     "Area $\propto$ distance$^{-4}$\n",
        #     rotation=-45,
        #     **kwargs,
        # )

        # # === Scatter plots: all events (background, grey), selected events (color/size) ===
        # ax_joint.scatter(
        #     "distance",
        #     f"area({skymap_area_cl})",
        #     s=2,
        #     facecolor="silver",
        #     edgecolor="none",
        #     data=table,
        # )

        # # Custom colormap for selected points
        # base_cmap = plt.get_cmap("cool")
        # cmap = LinearSegmentedColormap.from_list(
        #     "truncated_cool", base_cmap(np.linspace(1 / 3, 1))
        # )
        # marker_scale = 30
        # ax_joint.scatter(
        #     "distance",
        #     f"area({skymap_area_cl})",
        #     s=selected_table["objective_value"] * marker_scale,
        #     c=selected_table["detection_probability_known_position"],
        #     cmap=cmap,
        #     vmin=0,
        #     vmax=1,
        #     data=selected_table,
        # )

        # # === Mini legend for objective value and detection probability ===
        # ax_legend = fig.add_subplot(
        #     (
        #         0.075,
        #         (bottom + height + sep + 0.5 * depth - 0.08) /  default_fig_width_height_ratio,
        #         0.08,
        #         0.08 /  default_fig_width_height_ratio,
        #     ),
        #     aspect=1,
        # )
        # ax_legend.margins(0.2)
        # dx = 1 / 4
        # x = np.arange(0, 1 + dx, dx)
        # ax_legend.set_xticks([0, 0.5, 1])
        # ax_legend.set_yticks([0, 0.5, 1])
        # x, y = (a.ravel() for a in np.meshgrid(x, x))
        # ax_legend.scatter(
        #     x,
        #     y,
        #     s=2,
        #     facecolor="silver",
        #     edgecolor="none",
        #     clip_on=False,
        # )
        # ax_legend.scatter(
        #     x, y, s=marker_scale * x, c=y, vmin=0, vmax=1, cmap=cmap, clip_on=False
        # )
        # twin = ax_legend.twiny()
        # twin.set_frame_on(False)
        # twin.set_ylim(ax_legend.get_ylim())
        # twin.set_xticks([])
        # twin.set_xlabel("Objective val.")
        # ax_legend.set_ylabel("Detection prob.")
        # bbox_cbar, bbox_joint = [
        #     ax.get_window_extent().transformed(fig.transFigure.inverted())
        #     for ax in [ax_legend, ax_joint]
        # ]
        # ax_legend.annotate(
        #     "",
        #     (bbox_joint.xmin, bbox_joint.ymax),
        #     (bbox_cbar.xmax, bbox_cbar.ymin),
        #     "figure fraction",
        #     "figure fraction",
        #     clip_on=False,
        #     arrowprops=dict(
        #         facecolor="k",
        #         edgecolor="none",
        #         linewidth=0,
        #         arrowstyle="simple",
        #         shrinkA=6,
        #         shrinkB=6,
        #     ),
        # )
        # ax_legend.spines[["top", "right"]].set_visible(False)

        # ticks = [
        #     np.quantile(
        #         selected_table["distance"],
        #         [0.9],
        #         method="inverted_cdf",
        #     ).item(),
        #     np.quantile(
        #         selected_table["distance"],
        #         [0.9],
        #         weights=selected_table["detection_probability_known_position"],
        #         method="inverted_cdf",
        #     ).item(),
        # ]

        # # === Marginal histograms (distance/area) ===

        # ticks = [
        #     np.quantile(
        #         selected_table["distance"],
        #         [0.9],
        #         method="inverted_cdf",
        #     ).item(),
        #     np.quantile(
        #         selected_table["distance"],
        #         [0.9],
        #         weights=selected_table["detection_probability_known_position"],
        #         method="inverted_cdf",
        #     ).item(),
        # ]
        # bins = np.linspace(*np.log(ax_joint.get_xlim()), 16)
        # color = "silver"
        # values, _ = np.histogram(np.log(table["distance"]), bins=bins)
        # ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=0)
        # ax_x.set_ylim(0, values.max() * 1.1)
        # color = cmap(0)
        # values, _ = np.histogram(
        #     np.log(selected_table["distance"]),
        #     bins=bins,
        # )
        # ax_x.axvline(
        #     ticks[0],
        #     linewidth=plt.rcParams["xtick.major.width"],
        #     color=plt.rcParams["xtick.color"],
        #     zorder=1,
        # )
        # ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=2)
        # color = cmap(np.inf)
        # values, _ = np.histogram(
        #     np.log(selected_table["distance"]),
        #     weights=selected_table["detection_probability_known_position"],
        #     bins=bins,
        # )
        # ax_x.axvline(
        #     ticks[1],
        #     linewidth=plt.rcParams["xtick.major.width"],
        #     color=plt.rcParams["xtick.color"],
        #     zorder=3,
        # )
        # ax_x.stairs(values, np.exp(bins), color=color, fill=True, zorder=4)
        # twin = ax_x.twiny()
        # twin.set_frame_on(False)
        # twin.set_xlim(*ax_joint.get_xlim())
        # twin.set_xscale(ax_joint.get_xscale())
        # twin.set_xticks(
        #     [*ticks, np.prod(np.asarray(ax_joint.get_xlim()) ** [0.6, 0.4])],
        #     [*(f"{np.round(tick):g}\nMpc" for tick in ticks), "90th\npercentile"],
        # )
        # twin.xaxis.minorticks_off()
        # tick = twin.xaxis.get_major_ticks()[2]
        # tick.tick1line.set_visible(False)
        # tick.tick2line.set_visible(False)
        # tick = twin.xaxis.get_major_ticks()[0]
        # tick.label2.set_ha("left")
        # tick = twin.xaxis.get_major_ticks()[1]
        # tick.label2.set_ha("right")
        # tick = twin.xaxis.get_major_ticks()[2]
        # tick.label2.set_ha("right")

        # ticks = [
        #     np.quantile(
        #         selected_table[f"area({skymap_area_cl})"],
        #         [0.9],
        #         method="inverted_cdf",
        #     ).item(),
        #     np.quantile(
        #         selected_table[f"area({skymap_area_cl})"],
        #         [0.9],
        #         weights=selected_table["detection_probability_known_position"],
        #         method="inverted_cdf",
        #     ).item(),
        # ]
        # bins = np.linspace(*np.log(ax_joint.get_ylim()), 16)
        # color = "silver"
        # values, _ = np.histogram(np.log(table[f"area({skymap_area_cl})"]), bins=bins)
        # ax_y.stairs(
        #     values,
        #     np.exp(bins),
        #     color=color,
        #     fill=True,
        #     orientation="horizontal",
        #     zorder=0,
        # )
        # ax_y.set_xlim(0, values.max() * 1.1)
        # color = cmap(0)
        # values, _ = np.histogram(
        #     np.log(selected_table[f"area({skymap_area_cl})"]),
        #     bins=bins,
        # )
        # ax_y.axhline(
        #     ticks[0],
        #     linewidth=plt.rcParams["ytick.major.width"],
        #     color=plt.rcParams["ytick.color"],
        #     zorder=1,
        # )
        # ax_y.stairs(
        #     values,
        #     np.exp(bins),
        #     color=color,
        #     fill=True,
        #     orientation="horizontal",
        #     zorder=2,
        # )
        # color = cmap(np.inf)
        # values, _ = np.histogram(
        #     np.log(selected_table[f"area({skymap_area_cl})"]),
        #     weights=selected_table["detection_probability_known_position"],
        #     bins=bins,
        # )
        # ax_y.axhline(
        #     ticks[1],
        #     linewidth=plt.rcParams["ytick.major.width"],
        #     color=plt.rcParams["ytick.color"],
        #     zorder=3,
        # )
        # ax_y.stairs(
        #     values,
        #     np.exp(bins),
        #     color=color,
        #     fill=True,
        #     orientation="horizontal",
        #     zorder=4,
        # )
        # twin = ax_y.twinx()
        # twin.set_frame_on(False)
        # twin.set_ylim(*ax_joint.get_ylim())
        # twin.set_yscale(ax_joint.get_yscale())
        # twin.set_yticks(
        #     [*ticks, np.prod(np.asarray(ax_joint.get_ylim()) ** [0.1, 0.9])],
        #     [*(f"{np.round(tick):g} deg$^2$" for tick in ticks), "90th\npercentile"],
        # )
        # twin.yaxis.minorticks_off()
        # tick = twin.yaxis.get_major_ticks()[0]
        # tick.label2.set_va("bottom")
        # tick = twin.yaxis.get_major_ticks()[1]
        # tick.label2.set_va("top")
        # tick = twin.yaxis.get_major_ticks()[2]
        # tick.tick1line.set_visible(False)
        # tick.tick2line.set_visible(False)

        # kwargs = dict(
        #     # color="black",
        #     transform=ax_x.transAxes,
        #     fontsize=plt.rcParams["legend.fontsize"],
        #     zorder=5,
        #     ha="left",
        #     va="top",
        # )
        # ax_x.text(0.05, 0.9, "All events", color="dimgray", **kwargs)
        # ax_x.text(0.05, 0.78, "Triggered", color="tab:blue", **kwargs)
        # ax_x.text(0.05, 0.66, "Detected", color="magenta", **kwargs)

        # plt.setp(ax_x.get_xticklabels() + ax_y.get_yticklabels(), visible=False)
        # ax_x.set_yticks([])
        # ax_y.set_xticks([])

        # fig.text(
        #     0.9,
        #     0.9,
        #     run,
        #     fontsize="x-large",
        #     zorder=1000,
        #     ha="right",
        #     va="top",
        #     bbox={"facecolor": "white", "boxstyle": "round"},
        # )

        # === Save to PDF ===
        fig.savefig(f"area-distance-{run}.pdf", bbox_inches="tight")
        plt.show()
        plt.close(fig)


# if __name__ == "__main__":


#     ns_max=3.0
#     skymap_area_cl = 90

#     # Load events and constants
#     root = get_project_root()
#     events_file = root / "data" / "events.ecsv"

#     main_table, constants = load_events(events_file)

#     # process events and filter the objective values then masses and return the runs, main_table_dict and selected_tables
#     runs, main_table_dict, selected_tables = preprocess_events(main_table, constants, ns_max=ns_max)
#     min_area, max_area, max_distance, crossover_distance = area_distance_limits(constants=constants, skymap_area_cl=skymap_area_cl)

#     # Plot area vs distance for each run
#     plot_area_distance(runs, main_table_dict, selected_tables, min_area, max_area, max_distance, crossover_distance, skymap_area_cl=90,)
