"""
Create the Field of View (FOV) for LSST science detectors from `lsst_science_detectors.csv`.
Converts detector offsets to RA/Dec, defines sky regions, and maps them onto a HEALPix grid.
Offsets are in millimeters (mm), converted using a plate scale of 0.2 arcsec/pixel.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import Table
from astropy_healpix import HEALPix
from m4opt.fov import footprint_healpix
from regions import PolygonSkyRegion, RectangleSkyRegion, Regions


def read_csv(file_path):
    """Read a CSV file containing detector information."""
    return Table.read(file_path, format="csv")


def get_bbox_size(physical_type):
    """Return the detector's bounding box size (width, height in pixels)."""
    return {"E2V": [4095, 4003], "ITL": [4071, 3999]}.get(physical_type, None)


def convert_rectangle_to_polygon_skyregion(rect):
    """Convert a RectangleSkyRegion to a PolygonSkyRegion (RA/Dec)."""
    ra, dec = rect.center.ra.deg, rect.center.dec.deg
    half_w, half_h = rect.width.to_value(u.deg) / 2, rect.height.to_value(u.deg) / 2
    corners_ra = [ra - half_w, ra + half_w, ra + half_w, ra - half_w]
    corners_dec = [dec - half_h, dec - half_h, dec + half_h, dec + half_h]
    return PolygonSkyRegion(
        vertices=SkyCoord(ra=corners_ra * u.deg, dec=corners_dec * u.deg)
    )


def make_fov(detectors):
    """Generate FOV regions in RA/Dec from detector positions."""
    plate_scale = 0.2 * u.arcsec  # arcsec/pixel
    pixel_size = 0.01  # mm/pixel
    mm_to_arcsec = plate_scale / pixel_size
    fov_regions = []

    for det in detectors:
        bbox_size = get_bbox_size(det["physical_type"])
        if not bbox_size:
            continue

        x_arcsec, y_arcsec = (
            det["x_offset"] * mm_to_arcsec,
            det["y_offset"] * mm_to_arcsec,
        )
        width, height = bbox_size * plate_scale
        center_coord = SkyCoord(ra=x_arcsec, dec=y_arcsec, frame="icrs")
        rect = RectangleSkyRegion(center=center_coord, width=width, height=height)
        fov_regions.append(convert_rectangle_to_polygon_skyregion(rect))

    return Regions(fov_regions)


# Load data and compute FOV
detectors = read_csv("lsst_science_detectors.csv")
lsst_fov = make_fov(detectors)

# Map FOV onto HEALPix grid
hpx = HEALPix(nside=512, frame=ICRS())
target_coords = SkyCoord(*(np.meshgrid([-15, 0, 15], [-15, 0, 15]) * u.deg))
pixels = np.unique(
    np.concatenate(footprint_healpix(hpx, lsst_fov, target_coords).ravel())
)
