import ligo.skymap.plot  # noqa: F401
import numpy as np
from astropy import units as u
from m4opt.fov import footprint
from m4opt.missions import lsst
from matplotlib import pyplot as plt
from tqdm.auto import tqdm


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


fig = plt.figure(dpi=300)
ax = fig.add_subplot(projection="astro globe")
for coord in ax.coords:
    coord.set_ticks_visible(False)
    coord.set_ticklabel_visible(False)
fig.tight_layout()

frames = np.arange(0, 90, 5)
artists_to_delete = []
with tqdm(total=len(frames) + 3, desc="drawing", unit="frame") as progress:

    def func(rotation):
        result = []
        while artists_to_delete:
            artist = artists_to_delete.pop()
            result.append(artist)
            artist.remove()

        # Optionally overlay the FOV shape as a Matplotlib artist
        fov_mission = lsst.fov
        # if hasattr(fov_mission, "__len__") and len(fov_mission) > 1:
        for sub_region in fov_mission:
            for region in footprint(sub_region, lsst.skygrid, rotation * u.deg):
                artist = region.to_pixel(ax.wcs).as_artist(linewidth=0.25)
                artists_to_delete.append(artist)
                result.append(artist)
                ax.add_patch(artist)
            progress.update()
        return result

    func(45)
    fig.savefig("skygrid_lsst.pdf")
