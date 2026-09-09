"""
Build a follow-up skymap for UVEX from an ULTRASAT schedule.

Strategy
--------
ULTRASAT observes the full GW skymap for 6 hours with its wide FOV.
UVEX then follows up on the *same region* covered by ULTRASAT, going
deeper to detect sources that ULTRASAT may have missed due to its
shallower magnitude limit.

The output skymap retains only the pixels covered by ULTRASAT,
renormalized to sum to 1, so that the UVEX MILP scheduler focuses
exclusively on those pixels.
"""

import numpy as np
from astropy.coordinates import ICRS
from astropy.table import QTable, unique
from astropy_healpix import HEALPix
from ligo.skymap.bayestar import derasterize, rasterize
from ligo.skymap.io import read_sky_map, write_sky_map
from m4opt import missions
from m4opt.fov import footprint_healpix

from earthorbitplan.utils.path import get_project_root

root = get_project_root()
ultrasat_schedule = root / "data" / "ultrasat" / "1000.ecsv"
ultrasat_skymap = root / "data" / "ultrasat" / "1000.fits"
output_skymap = root / "data" / "ultrasat" / "1000_uvex_followup.fits"

# ------------------------------------------------------------------
# 1. Load the ULTRASAT schedule and extract metadata
# ------------------------------------------------------------------
plan = QTable.read(ultrasat_schedule)
args = plan.meta["args"]

nside = args["nside"]  # 128
mission = getattr(missions, args["mission"])  # ultrasat

hpx = HEALPix(nside=nside, frame=ICRS(), order="nested")

# ------------------------------------------------------------------
# 2. Load the original GW skymap
# ------------------------------------------------------------------
skymap_moc = read_sky_map(ultrasat_skymap, moc=True)
skymap_flat = rasterize(skymap_moc, hpx.level)
probs = skymap_flat["PROB"].copy()

# ------------------------------------------------------------------
# 3. Find the unique fields observed by ULTRASAT
#    (visits=2 means each field appears twice — keep one per position)
# ------------------------------------------------------------------
observations = plan[plan["action"] == "observe"].filled()

coords = observations["target_coord"].to_table()
coords["i"] = np.arange(len(coords))
i = np.sort(unique(coords, keys=["ra", "dec"])["i"])
fields = observations[i]

# ------------------------------------------------------------------
# 4. Compute the HEALPix pixels covered by ULTRASAT
# ------------------------------------------------------------------
covered_pixels = np.unique(
    np.concatenate(
        [
            footprint
            for footprint in footprint_healpix(
                hpx, mission.fov, fields["target_coord"], fields["roll"]
            )
        ]
    )
)

print(f"Unique fields    : {len(fields)}")
print(f"Covered pixels   : {len(covered_pixels)}")
print(f"Covered prob     : {probs[covered_pixels].sum():.4f}")

# ------------------------------------------------------------------
# 5. Build the follow-up skymap
#    Keep ONLY the pixels covered by ULTRASAT and renormalize.
#    UVEX will then schedule observations on this reduced skymap,
#    going deeper than ULTRASAT on the same sky region.
# ------------------------------------------------------------------
probs_followup = np.zeros(hpx.npix)
probs_followup[covered_pixels] = probs[covered_pixels]

total = probs_followup.sum()
print(f"Follow-up prob   : {total:.4f}")

if total > 0:
    probs_followup /= total

# ------------------------------------------------------------------
# 6. Reconstruct the MOC and write to disk
# ------------------------------------------------------------------
skymap_flat_followup = skymap_flat.copy()
skymap_flat_followup["PROB"] = probs_followup

skymap_followup_moc = derasterize(skymap_flat_followup)
skymap_followup_moc.meta = skymap_moc.meta  # preserve distmean, diststd, gps_time

write_sky_map(output_skymap, skymap_followup_moc, moc=True)
print(f"Done => {output_skymap}")
