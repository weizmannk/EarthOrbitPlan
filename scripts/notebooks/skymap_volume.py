## Using a command line  :
#   ligo-skymap-stats -o 14.dat  ./14.fits --cosmology --contour 20 50 90 -j


from typing import Iterable, List, Optional

import numpy as np
from astropy.table import Table
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import crossmatch


def statistic_skymaps_no_db(
    fitsfiles: Iterable[str],
    contours_percent: Optional[List[float]] = None,
    areas: Optional[List[float]] = None,
    modes: bool = False,
    cosmology: bool = False,
) -> Table:
    """
    Compute the skymap statistics (no injection database).
    no injection database , means no true location (this reserved for the simulation), no SNR ...


    Parameters
    ----------
    fitsfiles : iterable of str
        Paths to input sky maps (FITS or FITS.gz). Globs should be expanded by the caller.
    contours_percent : list of float, optional
        Credible levels P (in percent) for which to report area(P), dist(P), vol(P).
        Example: [20, 50, 90]. Default: [] (no contour-based columns).
    areas : list of float, optional
        Fixed sky areas A (in deg^2) for which to report prob(A).
        Example: [10, 50, 100]. Default: [] (no area-based probability columns).
    modes : bool, default False
        If True, report modes(P): number of disjoint figures in the P% credible region.
    cosmology : bool, default False
        If True, report volumes as comoving volumes.

    Returns
    -------
    pandas.DataFrame
        A table with one row per input sky map and columns similar to ligo-skymap-stats
        (excluding DB-dependent columns).
    """

    true_coord = None
    contours_percent = contours_percent or []
    areas = areas or []
    pvalues = 0.01 * np.asarray(contours_percent, dtype=float)  # convert % => fraction

    # Column headers
    colnames = [
        "coinc_event_id",
        "runtime",
        "distmean",
        "diststd",
        "log_bci",
        "log_bsn",
    ]
    colnames += [f"area({p:g})" for p in contours_percent]
    colnames += [f"prob({a:g})" for a in areas]
    colnames += [f"dist({p:g})" for p in contours_percent]
    colnames += [f"vol({p:g})" for p in contours_percent]
    if modes:
        colnames += [f"modes({p:g})" for p in contours_percent]

    rows = []

    for f in fitsfiles:
        sky_map = read_sky_map(str(f), moc=True)

        coinc_event_id = sky_map.meta.get("objid", np.nan)
        runtime = sky_map.meta.get("runtime", np.nan)
        distmean = sky_map.meta.get("distmean", np.nan)
        diststd = sky_map.meta.get("diststd", np.nan)
        log_bci = sky_map.meta.get("log_bci", np.nan)
        log_bsn = sky_map.meta.get("log_bsn", np.nan)

        # No database => no true sky/distance => pass true_coord=None
        (
            searched_area,
            searched_prob,
            offset,
            searched_modes,
            contour_areas,
            area_probs,
            contour_modes,
            searched_prob_dist,
            contour_dists,
            searched_vol,
            searched_prob_vol,
            contour_vols,
            probdensity,
            probdensity_vol,
        ) = crossmatch(
            sky_map,
            true_coord,
            contours=pvalues,
            areas=areas,
            modes=modes,
            cosmology=cosmology,
        )

        row = [
            coinc_event_id,
            runtime,
            distmean,
            diststd,
            log_bci,
            log_bsn,
        ]
        row += contour_areas
        row += area_probs
        row += contour_dists
        row += contour_vols
        if modes:
            row += contour_modes

        rows.append(row)

    table = Table(rows=rows, names=colnames)

    return table


if __name__ == "__main__":
    example_files = ["./14.fits"]
    if example_files:
        table = statistic_skymaps_no_db(
            example_files,
            contours_percent=[20, 50, 90],
            areas=[],
            modes=False,
            cosmology=False,
        )

        table.write("allsky_no_db.dat", format="ascii.fast_tab", overwrite=True)

        print(table)
