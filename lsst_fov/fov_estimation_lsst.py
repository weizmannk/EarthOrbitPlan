#!/usr/bin/env python
"""

Constructs a plus-shaped arrangement of LSST Field of Vied in celestial coordinates.
Each raft is defined as a 3×3 block of CCDs, placed in a 5×5 grid
with the four corner positions removed. This yields 21 total rafts, creating
the characteristic plus layout.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy_healpix import HEALPix
from m4opt.fov import footprint_healpix
from regions import PolygonSkyRegion, RectangleSkyRegion, Regions


def convert_rectangle_to_polygon_skyregion(
    rect: RectangleSkyRegion,
) -> PolygonSkyRegion:
    """
    Convert a RectangleSkyRegion to a PolygonSkyRegion, preserving RA/Dec coordinates.

    Parameters
    ----------
    rect : RectangleSkyRegion
        A rectangle on the celestial sphere, defined by a center (RA, Dec)
        and width/height in degrees.

    Returns
    -------
    PolygonSkyRegion
        Polygonal region whose corners match the rectangle's corners in RA/Dec.
    """
    ra_center = rect.center.ra.deg
    dec_center = rect.center.dec.deg

    half_width = rect.width.to_value(u.deg) / 2.0
    half_height = rect.height.to_value(u.deg) / 2.0

    corners_ra = [
        ra_center - half_width,
        ra_center + half_width,
        ra_center + half_width,
        ra_center - half_width,
    ]
    corners_dec = [
        dec_center - half_height,
        dec_center - half_height,
        dec_center + half_height,
        dec_center + half_height,
    ]

    return PolygonSkyRegion(
        vertices=SkyCoord(ra=corners_ra * u.deg, dec=corners_dec * u.deg)
    )


def make_plus_shaped_fov() -> Regions:
    """
    Create a plus-shaped LSST focal plane layout in RA/Dec coordinates.

    Each raft is modeled as a 3×3 array of CCDs. The rafts are arranged in
    a 5×5 grid, where the four corner rafts are omitted to form a plus shape.
    This results in 21 rafts.

    Returns
    -------
    Regions
        A collection of PolygonSkyRegion objects, each representing one raft.
    """
    #  The side length of a single LSST CCD in degree
    ccd_side_deg = np.sqrt(3.5 / 189)

    # Each raft is a 3×3 block of CCDs
    raft_width_ccd = 3
    raft_height_ccd = 3

    raft_regions = []
    nrows, ncols = 5, 5

    for row in range(nrows):
        for col in range(ncols):
            # Skip the four corners (top=0, bottom=4 => keep columns=1..3)
            if (row in [0, 4]) and (col not in [1, 2, 3]):
                continue

            # Determine the center of this raft in RA/Dec
            ra_center_deg = (col - (ncols // 2)) * raft_width_ccd * ccd_side_deg
            dec_center_deg = ((nrows // 2) - row) * raft_height_ccd * ccd_side_deg

            center_coord = SkyCoord(
                ra_center_deg * u.deg, dec_center_deg * u.deg, frame="icrs"
            )

            # Convert the raft (as a rectangle) to a polygon
            width_deg = raft_width_ccd * ccd_side_deg
            height_deg = raft_height_ccd * ccd_side_deg
            rect = RectangleSkyRegion(
                center=center_coord, width=width_deg * u.deg, height=height_deg * u.deg
            )
            poly = convert_rectangle_to_polygon_skyregion(rect)
            raft_regions.append(poly)

    return Regions(raft_regions)


plus_sign = make_plus_shaped_fov()


# Initialize a HEALPix object (using ICRS frame).
hpx = HEALPix(nside=512, frame=ICRS())

# Used footprint_healpix with these regions
target_coords = SkyCoord(*(np.meshgrid([-15, 0, 15], [-15, 0, 15]) * u.deg))
pixels = np.unique(
    np.concatenate(footprint_healpix(hpx, plus_sign, target_coords).ravel())
)
