import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
from astropy_healpix import HEALPix
from regions import RectangleSkyRegion, PolygonSkyRegion, Regions

from m4opt.fov import footprint_healpix

# Initialize a HEALPix object (using ICRS frame).
hpx = HEALPix(nside=64, frame=ICRS())

def convert_rectangle_to_polygon_skyregion(rect: RectangleSkyRegion) -> PolygonSkyRegion:
    """
    Convert a RectangleSkyRegion to a PolygonSkyRegion, preserving RA/Dec coordinates.
    
    This function computes the corner points of the rectangle (in degrees)
    and uses them to create an equivalent polygonal region.
    
    Parameters
    ----------
    rect : RectangleSkyRegion
        The rectangular region to convert, including its center and dimensions.
    
    Returns
    -------
    PolygonSkyRegion
        PolygonSkyRegion that matches the original rectangle’s footprint.
    """
    ra_center = rect.center.ra.deg
    dec_center = rect.center.dec.deg

    half_width = rect.width.to_value(u.deg) / 2.0
    half_height = rect.height.to_value(u.deg) / 2.0

    corners_ra = [
        ra_center - half_width,
        ra_center + half_width,
        ra_center + half_width,
        ra_center - half_width
    ]
    corners_dec = [
        dec_center - half_height,
        dec_center - half_height,
        dec_center + half_height,
        dec_center + half_height
    ]

    return PolygonSkyRegion(
        vertices=SkyCoord(ra=corners_ra*u.deg, dec=corners_dec*u.deg)
    )

# Example parameters
ccd_side_deg = np.sqrt(9.6 / 189)

# Define three rectangular regions forming a "plus" shape
center_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, 0*u.deg),
    width=15 * ccd_side_deg * u.deg,
    height=9 * ccd_side_deg * u.deg
)

top_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, 6 * ccd_side_deg * u.deg),
    width=9 * ccd_side_deg * u.deg,
    height=3 * ccd_side_deg * u.deg
)

bottom_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, -6 * ccd_side_deg * u.deg),
    width=9 * ccd_side_deg * u.deg,
    height=3 * ccd_side_deg * u.deg
)


# Convert rectangles into polygons and combine them into a single "plus sign" region
plus_sign = Regions([
    convert_rectangle_to_polygon_skyregion(r) 
    for r in [top_rect, center_rect, bottom_rect]
])


# Used footprint_healpix with these regions
target_coords = SkyCoord(*(np.meshgrid([-15, 0, 15], [-15, 0, 15]) * u.deg))
pixels = np.unique(
    np.concatenate(
        footprint_healpix(hpx, plus_sign, target_coords).ravel()
    )
)

print("Unique HEALPix pixels intersecting the plus_sign footprint:", pixels)
