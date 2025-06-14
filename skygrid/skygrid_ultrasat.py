#!/usr/bin/env python3
"""
Generate a Mollweide projection map showing UV extinction (A_U) for V45==1,
with overplotted ULTRASAT footprints.
Usage:
     python skygrid_ultrasat.py  ULTRASAT_LCS_nonoverlapping_V45.pdf
"""

from pathlib import Path
from typing import Annotated

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import typer
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from ligo.skymap.plot.poly import cut_prime_meridian
from m4opt._cli.core import app
from m4opt.fov import footprint
from m4opt.missions import ultrasat
from m4opt.utils.console import progress


def read_skygrid():
    """Read the CSV file and return SkyCoord and UV_extinction."""
    table = Table.read("LCS_nonoverlapping_grid.csv", format="ascii.csv")
    mask = table["V45"] == 1
    skygrid = SkyCoord(table["RA"][mask], table["Dec"][mask], unit=u.deg)
    UV_extinction = table["AU"][mask]
    print(f"Max UV extinction: {max(UV_extinction):.2f} mag")
    return skygrid, UV_extinction


@app.command()
@progress()
def make_plot(
    still: Annotated[
        typer.FileBinaryWrite,
        typer.Argument(help="Output file (e.g., skygrid_ultrasat.pdf)"),
    ],
    dpi: Annotated[
        float | None,
        typer.Option(help="Figure DPI (dots per inch)", show_default=False),
    ] = 300,
):
    """Generate static Mollweide UV extinction map with ULTRASAT footprints."""
    mission = ultrasat
    skygrid, UV_ext = read_skygrid()

    fig_width_inch = 7
    fig_height_inch = fig_width_inch / 2
    fig = plt.figure(
        figsize=(fig_width_inch, fig_height_inch), dpi=dpi, constrained_layout=True
    )
    ax = fig.add_subplot(111, projection="astro hours mollweide")
    ax.set_facecolor("#999999")
    ax.grid(color="white", linestyle=":", linewidth=0.3, alpha=0.3)

    lon = skygrid.ra.wrap_at(180 * u.deg).radian
    lat = skygrid.dec.radian
    vmin, vmax = UV_ext.min(), UV_ext.max()

    sc = ax.scatter(
        lon,
        lat,
        c=UV_ext,
        vmin=vmin,
        vmax=vmax,
        cmap="plasma",
        s=0,
        alpha=0.7,
        edgecolors="none",
        transform=ax.get_transform("world"),
    )

    cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.05, aspect=40)
    cbar.ax.xaxis.set_label_position("bottom")
    cbar.set_label(r"UV mean extinction (AB mag)", y=1)

    for i, region in enumerate(footprint(mission.fov, skygrid, 0 * u.deg)):
        frac = (UV_ext[i] - vmin) / (vmax - vmin)
        color = cm.plasma(frac)
        verts = np.column_stack(
            [
                region.vertices.ra.wrap_at(180 * u.deg).radian,
                region.vertices.dec.radian,
            ]
        )
        for seg in cut_prime_meridian(verts):
            poly = plt.Polygon(
                np.rad2deg(seg),
                transform=ax.get_transform("world"),
                facecolor="none",
                edgecolor=color,
                linewidth=0.3,
                alpha=0.7,
            )
            ax.add_patch(poly)

    # os.makedirs("../figures", exist_ok=True)

    fig.savefig(still, bbox_inches="tight", format=Path(still.name).suffix.lstrip("."))
    plt.show()


if __name__ == "__main__":
    typer.run(make_plot)
