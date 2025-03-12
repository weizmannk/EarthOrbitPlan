import numpy as np
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord
from astropy_healpix import HEALPix
from regions import RectangleSkyRegion, PolygonSkyRegion
from shapely.geometry import Polygon 

from m4opt.fov import footprint_healpix

hpx = HEALPix(nside=64, frame=ICRS())

def rectangle_to_shapely(rect: RectangleSkyRegion):
    """Convert a RectangleSkyRegion to a shapely Polygon (RA/Dec degrees)."""
    ra_center = rect.center.ra.deg
    dec_center = rect.center.dec.deg
    half_width = rect.width.value / 2
    half_height = rect.height.value / 2

    corners = [
        (ra_center - half_width, dec_center - half_height),
        (ra_center + half_width, dec_center - half_height),
        (ra_center + half_width, dec_center + half_height),
        (ra_center - half_width, dec_center + half_height)
    ]

    return Polygon(corners)


# Parameters for rectangles
#this is the lenght side of each CCD
# FOV = 9.6 with a total number of 189 CDDs inside the Plus sign
ccd_side_deg = np.sqrt(9.6 / 189)

# Define three rectangles forming a plus shape
center_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, 0*u.deg),
    width=15*ccd_side_deg*u.deg,
    height=9*ccd_side_deg*u.deg
)

top_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, 6*ccd_side_deg*u.deg),
    width=9*ccd_side_deg*u.deg,
    height=3*ccd_side_deg*u.deg
)

bottom_rect = RectangleSkyRegion(
    center=SkyCoord(0*u.deg, -6*ccd_side_deg*u.deg),
    width=9*ccd_side_deg*u.deg,
    height=3*ccd_side_deg*u.deg
)

# Convert rectangles to shapely polygons
shapely_polys = [rectangle_to_shapely(rect) for rect in [top_rect, center_rect, bottom_rect]]

# union of polygons
union_poly = shapely_polys[0]
for poly in shapely_polys[1:]:
    union_poly = union_poly.union(poly)
# Extract outer boundary
exterior_coords = list(union_poly.exterior.coords)

# Create PolygonSkyRegion
ra_vals, dec_vals = zip(*exterior_coords)
plus_sign = PolygonSkyRegion(SkyCoord(ra=ra_vals*u.deg, dec=dec_vals*u.deg))


# footprint with healpix
target_coords = SkyCoord(*(np.meshgrid([-15, 0, 15], [-15, 0, 15]) * u.deg))
pixels = np.unique(
        np.concatenate(footprint_healpix(hpx, plus_sign, target_coords).ravel())
    )