import ligo.skymap.plot  # noqa: F401
import numpy as np
import regions
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy_healpix import HEALPix
from matplotlib import pyplot as plt

from m4opt.fov import footprint_healpix
from m4opt.missions import ultrasat


def customize_style(columns=1):
    target_width = 3.5 if columns == 1 else 7.25
    width, height = plt.rcParams["figure.figsize"]
    plt.style.use("seaborn-v0_8-paper")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["figure.figsize"] = (target_width, height * target_width / width)


customize_style()
plt.rcParams["figure.figsize"] = [5.0, 5.0]

# -- couleurs ----------------------------------------------
FOV_FACE = "#E8B84C"
FOV_EDGE = "#A87A1A"
GRID_COLOR = "#CCCCCC"
BG_COLOR = "#FDF9F0"

center = SkyCoord(0 * u.deg, 0 * u.deg)
hpx = HEALPix(nside=128, frame=ICRS())

fig = plt.figure(tight_layout=True, facecolor="white")
ax = fig.add_subplot(
    projection="astro zoom",
    center=center,
    radius=9 * u.deg,
)
ax.set_facecolor(BG_COLOR)
transform = ax.get_transform("world")


def plot_boundaries(ipix, **kwargs):
    lons, lats = hpx.boundaries_lonlat(ipix, 1)
    coords = np.moveaxis([lons.to_value(u.deg), lats.to_value(u.deg)], 0, -1)
    for coord in coords:
        ax.add_patch(plt.Polygon(coord, transform=transform, **kwargs))


# -- grille HEALPix ----------------------------------------
plot_limit_ipix = footprint_healpix(
    hpx,
    regions.PolygonSkyRegion(SkyCoord(*ax.wcs.calc_footprint().T, unit=u.deg)),
    center,
)
plot_boundaries(
    plot_limit_ipix,
    facecolor="none",
    edgecolor=GRID_COLOR,
    linewidth=0.4,
    zorder=1,
)

# -- pixels FOV -------------------------------------------
fov_ipix = footprint_healpix(hpx, ultrasat.fov, center)
plot_boundaries(
    fov_ipix,
    facecolor=FOV_FACE,
    edgecolor=FOV_FACE,
    alpha=0.45,
    linewidth=0.0,
    zorder=2,
)

# -- contour FOV ------------------------------------------
ax.add_artist(
    ultrasat.fov.to_pixel(ax.wcs).as_artist(
        facecolor="none",
        edgecolor=FOV_EDGE,
        linewidth=2,
        zorder=3,
    )
)

# -- boresight --------------------------------------------
ax.plot(
    center.ra.deg,
    center.dec.deg,
    "+",
    color=FOV_EDGE,
    markersize=8,
    markeredgewidth=1.4,
    transform=transform,
    zorder=4,
)

# -- annotation -------------------------------------------
ax.text(
    0.50,
    0.02,
    r"FOV $= 204\,\mathrm{deg}^2$",
    transform=ax.transAxes,
    ha="center",
    va="bottom",
    fontsize=10,
    color=FOV_EDGE,
    bbox=dict(
        boxstyle="round,pad=0.25", fc="white", ec=FOV_EDGE, alpha=0.85, linewidth=0.8
    ),
)

# -- axes RA/Dec en degrés ---------------------------------
ra = ax.coords["ra"]
dec = ax.coords["dec"]

# format degrés au lieu de heures
ra.set_major_formatter("d.d")
dec.set_major_formatter("d.d")

ra.set_ticks(spacing=5 * u.deg)
dec.set_ticks(spacing=5 * u.deg)

ra.set_axislabel("R.A. (deg)", fontsize=9, minpad=0.5)
dec.set_axislabel("Dec. (deg)", fontsize=9, minpad=0.5)

ra.set_ticklabel(size=8, color="#555555")
dec.set_ticklabel(size=8, color="#555555")

ra.set_ticks_visible(True)
dec.set_ticks_visible(True)
ra.set_ticklabel_visible(True)
dec.set_ticklabel_visible(True)

ax.coords.frame.set_linewidth(0.8)
ax.coords.frame.set_color("#AAAAAA")

fig.savefig("fov_ultrasat_204.pdf", bbox_inches="tight", dpi=300)
fig.savefig("fov_ultrasat_204.png", bbox_inches="tight", dpi=300)
plt.close()
