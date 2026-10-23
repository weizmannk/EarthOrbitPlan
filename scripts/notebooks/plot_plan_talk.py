#!/usr/bin/env python
"""Talk-style figure of one M4OPT observing plan: zoomed sky map + timeline.

Reproduces the layout of m4opt-paper/figures-for-talk/3628.pdf (UVEX) for any
schedule written by ``m4opt schedule`` (here: ULTRASAT, event 721, O5c).

Top panel   : required exposure time (min) over the zoomed 90% region,
              observed fields (one shade per visit), field being observed now
              outlined in black, 90% credible contour.
Bottom panel: cumulative probability (blue) and area (orange) covered;
              light = covered at least once, dark = covered at least twice.

Example
-------
    python plot_plan_talk.py ultrasat-non-overlap/schedules/O5c/721.ecsv \
        --skymap ultrasat-non-overlap/skymaps/O5c/721.fits \
        --center "0d -14d" --radius 47deg --aspect 0.66 --still-time 4.75h \
        --still ultrasat_721_plan.pdf --movie ultrasat_721_plan.mp4

Exposure-time map, footprints and coverage are computed with the same M4OPT
functions as ``m4opt animate``.
"""

import argparse
from itertools import accumulate

import ligo.skymap.plot  # noqa: F401  (registers the "astro zoom" projection)
import numpy as np
import synphot
from astropy import units as u
from astropy.coordinates import ICRS, Distance, SkyCoord
from astropy.table import QTable
from astropy.time import Time
from astropy_healpix import HEALPix
from ligo.skymap.bayestar import rasterize
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import find_greedy_credible_levels
from m4opt import missions
from m4opt.fov import footprint, footprint_healpix
from m4opt.synphot import observing
from m4opt.synphot.extinction import DustExtinction
from matplotlib import pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.colors import to_rgb
from matplotlib.transforms import blended_transform_factory
from regions import Regions

PROB_COLOR = "#1f77ff"  # probability / 90% contour
AREA_COLOR = "#ff4d00"  # area / footprints


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("schedule", help="schedule .ecsv from m4opt schedule")
    p.add_argument(
        "--skymap", help="sky map .fits (default: path stored in the schedule)"
    )
    p.add_argument("--center", required=True, help='zoom centre, e.g. "4.5d -15d"')
    p.add_argument("--radius", default="36deg", help="zoom radius")
    p.add_argument(
        "--still-time", default="4.75h", help="time of the still frame since the event"
    )
    p.add_argument("--still", help="output still frame (.pdf or .png)")
    p.add_argument("--movie", help="output animation (.mp4 or .gif)")
    p.add_argument("--time-step", default="2.5min", help="animation time step")
    p.add_argument(
        "--aspect",
        type=float,
        default=1.0,
        help="height/width of the map panel (<1 crops in Dec)",
    )
    p.add_argument(
        "--vmax", type=float, default=None, help="colour-bar maximum in minutes"
    )
    p.add_argument("--dpi", type=float, default=300)
    return p.parse_args(argv)


def wrap_ra(ra_deg, center_deg):
    """Keep RA continuous around the zoom centre (handles fields crossing 0h)."""
    return (ra_deg - center_deg + 180) % 360 - 180 + center_deg


def main(argv=None):
    args = parse_args(argv)
    center = SkyCoord(args.center)
    radius = u.Quantity(args.radius)
    still_time = u.Quantity(args.still_time)
    time_step = u.Quantity(args.time_step)

    # ---- schedule --------------------------------------------------------
    table = QTable.read(args.schedule, format="ascii.ecsv")
    meta = table.meta["args"]
    table = table[table["action"] == "observe"]
    table["end_time"] = table["start_time"] + table["duration"]
    mission = getattr(missions, meta["mission"])
    nside, visits = meta["nside"], meta["visits"]
    deadline, delay = meta["deadline"], meta["delay"]
    bandpass = meta["bandpass"]
    bandpass = bandpass if isinstance(bandpass, str) else bandpass[0]
    event_time = Time(meta["event_time"])

    # ---- sky map -----------------------------------------------------------
    hpx = HEALPix(nside, frame=ICRS(), order="nested")
    skymap_moc = read_sky_map(args.skymap or meta["skymap"], moc=True)
    probs = rasterize(skymap_moc["UNIQ", "PROBDENSITY"], hpx.level)["PROB"]
    cls = find_greedy_credible_levels(probs)

    # ---- required exposure time (as in m4opt animate) ------------------------
    distmod = Distance(skymap_moc.meta["distmean"] * u.Mpc).distmod
    t0 = event_time + delay
    with observing(
        observer_location=mission.observer_location(t0),
        target_coord=hpx.healpix_to_skycoord(np.arange(hpx.npix)),
        obstime=t0,
    ):
        exptime = mission.detector.get_exptime(
            meta["snr"],
            synphot.SourceSpectrum(
                synphot.ConstFlux1D(meta["absmag_mean"] * u.ABmag + distmod)
            )
            * DustExtinction(),
            bandpass,
        ).to_value(u.min)

    # ---- cumulative coverage after each observation ------------------------
    ivisit = np.arange(visits)
    area = np.empty((len(table), visits))
    prob = np.empty((len(table), visits))
    for k, pix in enumerate(
        accumulate(
            footprint_healpix(hpx, mission.fov, table["target_coord"], table["roll"]),
            lambda *a: np.concatenate(a),
        )
    ):
        visited = np.bincount(pix, minlength=hpx.npix) > ivisit[:, np.newaxis]
        area[k] = np.count_nonzero(visited, axis=1) * hpx.pixel_area.to_value(u.deg**2)
        prob[k] = np.sum(probs * visited, axis=1)

    # ---- figure ---------------------------------------------------------------
    plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(4.2, 4.7), dpi=args.dpi)
    ax = fig.add_axes(
        (0.12, 0.415, 0.84, 0.52),
        projection="astro degrees zoom",
        center=center,
        radius=radius,
    )
    ax.grid(color="0.85", lw=0.5)
    for key in ("ra", "dec"):
        ax.coords[key].set_major_formatter("d")
    ax.coords["ra"].set_ticks_position("t")
    ax.coords["ra"].set_ticklabel_position("t")
    ax.coords["dec"].set_ticks_position("l")
    ax.coords["dec"].set_ticklabel_position("l")
    ax.coords["ra"].set_axislabel("")
    ax.coords["dec"].set_axislabel("")

    vmax = args.vmax or 10 * np.ceil(table["duration"].max().to_value(u.min) / 10)
    if (
        args.aspect < 1
    ):  # crop in Dec to get a wide panel (fields span more RA than Dec)
        y0, y1 = ax.get_ylim()
        half = 0.5 * (y1 - y0) * args.aspect
        ym = 0.5 * (y0 + y1)
        ax.set_ylim(ym - half, ym + half)
        ax.set_anchor("S")

    im = ax.imshow_hpx(exptime, nested=True, cmap="viridis", vmin=0, vmax=vmax)
    ax.contour_hpx(cls, levels=[0.9], colors=[PROB_COLOR], nested=True, linewidths=1.2)

    world = ax.get_transform("world")
    regions = footprint(mission.fov, table["target_coord"], table["roll"])
    patches, outlines = [], []
    for region in regions:
        verts = [
            np.column_stack(
                (wrap_ra(r.vertices.ra.deg, center.ra.deg), r.vertices.dec.deg)
            )
            for r in (region.regions if isinstance(region, Regions) else [region])
        ]
        patches.append(
            [
                ax.add_patch(
                    plt.Polygon(
                        v,
                        transform=world,
                        facecolor=AREA_COLOR,
                        alpha=1 / visits,
                        lw=0,
                        visible=False,
                    )
                )
                for v in verts
            ]
        )
        outlines.append(
            [
                ax.add_patch(
                    plt.Polygon(
                        v,
                        transform=world,
                        facecolor="none",
                        edgecolor="black",
                        lw=1.2,
                        visible=False,
                    )
                )
                for v in verts
            ]
        )

    cax = fig.add_axes((0.12, 0.315, 0.84, 0.02))
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Required exposure time (min)")
    cb.ax.xaxis.set_label_position("top")

    axt = fig.add_axes((0.12, 0.075, 0.72, 0.14))
    axa = axt.twinx()
    hours = lambda t: (t - event_time).to_value(u.hour)  # noqa: E731
    for i in range(visits):
        alpha = (i + 1) / visits
        pc, ac = (
            np.asarray(to_rgb(c)) * alpha + (1 - alpha)
            for c in (PROB_COLOR, AREA_COLOR)
        )
        axt.hlines(
            prob[:, i], hours(table["start_time"]), hours(table["end_time"]), pc, lw=1.5
        )
        axa.hlines(
            area[:, i], hours(table["start_time"]), hours(table["end_time"]), ac, lw=1.5
        )
    axt.axvspan(0, delay.to_value(u.hour), color="0.88", lw=0)
    axt.set_xlim(0, deadline.to_value(u.hour))
    axt.set_ylim(0, 1)
    axa.set_ylim(0, None)
    axt.set_xlabel("time in hours since event")
    axt.set_ylabel("probability", color=PROB_COLOR)
    axa.set_ylabel("area (deg$^2$)", color=AREA_COLOR)
    axt.tick_params(axis="y", colors=PROB_COLOR)
    axa.tick_params(axis="y", colors=AREA_COLOR)
    now_line = axt.axvline(0, color="black", lw=1)
    now_label = axt.text(
        0,
        1.02,
        "now",
        ha="center",
        va="bottom",
        transform=blended_transform_factory(axt.transData, axt.transAxes),
    )

    def draw(time):
        h = hours(time)
        now_line.set_xdata([h])
        now_label.set_x(h)
        started = table["start_time"] <= time
        running = started & (table["end_time"] > time)
        for ps, os_, s, r in zip(patches, outlines, started, running):
            for p_ in ps:
                p_.set_visible(bool(s))
            for o in os_:
                o.set_visible(bool(r))
        return []

    if args.movie:
        times = (
            event_time
            + np.arange(
                0, (deadline + time_step).to_value(u.s), time_step.to_value(u.s)
            )
            * u.s
        )
        writer = (
            FFMpegWriter(fps=12, bitrate=4000)
            if args.movie.endswith(".mp4")
            else PillowWriter(fps=12)
        )
        FuncAnimation(fig, draw, frames=list(times)).save(
            args.movie, writer=writer, dpi=min(args.dpi, 200)
        )

    if args.still:
        draw(event_time + still_time)
        fig.savefig(args.still, dpi=args.dpi)


if __name__ == "__main__":
    main()
