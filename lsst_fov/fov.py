import ligo.skymap.plot  # noqa: F401
import numpy as np
import regions
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy_healpix import HEALPix
from m4opt.fov import footprint_healpix
from m4opt.missions import rubin
from matplotlib import pyplot as plt

# Map of available missions
mission_map = {
    # "ultrasat": ultrasat,
    # "uvex": uvex,
    "rubin": rubin,
    # "ztf": ztf,
}

mission_name = "rubin"
mission = mission_map[mission_name]


def customize_style(columns=1):
    """Customize Matplotlib style for publication."""
    if columns == 1:
        target_width = 3.5  # ApJ column size in inches
    else:
        target_width = 7.25  # ApJ two-column width in inches
    width, height = plt.rcParams["figure.figsize"]
    plt.style.use("seaborn-v0_8-paper")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["figure.figsize"] = (target_width, height * target_width / width)


def plot_boundaries(ipix, ax, hpx, transform, **kwargs):
    """Plot HEALPix boundaries for given pixel indices."""
    lons, lats = hpx.boundaries_lonlat(ipix, 1)
    coords = np.moveaxis([lons.to_value(u.deg), lats.to_value(u.deg)], 0, -1)
    for coord in coords:
        ax.add_patch(plt.Polygon(coord, clip_on=True, transform=transform, **kwargs))


if __name__ == "__main__":
    customize_style()
    plt.rcParams["figure.figsize"][1] = plt.rcParams["figure.figsize"][0]

    center = SkyCoord(0 * u.deg, 0 * u.deg)
    hpx = HEALPix(nside=128, frame=ICRS())

    fig = plt.figure(tight_layout=True)
    ax = fig.add_subplot(projection="astro zoom", center=center, radius=2.5 * u.deg)
    ax.coords.frame.set_color("none")
    transform = ax.get_transform("world")

    # Plot the plot limit boundaries
    plot_limit_ipix = footprint_healpix(
        hpx,
        regions.PolygonSkyRegion(SkyCoord(*ax.wcs.calc_footprint().T, unit=u.deg)),
        center,
    )
    plot_boundaries(
        plot_limit_ipix,
        ax,
        hpx,
        transform,
        facecolor="none",
        edgecolor=plt.rcParams["grid.color"],
        linewidth=plt.rcParams["grid.linewidth"],
    )

    # Plot the mission FOV boundaries
    fov_boundaries = footprint_healpix(hpx, mission.fov, center)
    plot_boundaries(
        fov_boundaries,
        ax,
        hpx,
        transform,
        facecolor="#9400D3",  # deep violet
        edgecolor="black",
        linewidth=plt.rcParams["lines.linewidth"],
    )

    # Optionally overlay the FOV shape as a Matplotlib artist
    fov_mission = mission.fov
    if hasattr(fov_mission, "__len__") and len(fov_mission) > 1:
        # For LSST
        for sub_fov in fov_mission:
            ax.add_artist(
                sub_fov.to_pixel(ax.wcs).as_artist(
                    edgecolor="black", linewidth=plt.rcParams["axes.linewidth"]
                )
            )
    else:
        ax.add_artist(
            fov_mission.to_pixel(ax.wcs).as_artist(
                edgecolor="black", linewidth=plt.rcParams["axes.linewidth"]
            )
        )

    # Remove axis labels and ticks
    for coord in ax.coords:
        coord.set_axislabel(None)
        coord.set_ticks_visible(False)
        coord.set_ticklabel_visible(False)

    # Save the figure
    fig.savefig(f"fov_{mission_name}.pdf")
    plt.show()
