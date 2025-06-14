import numpy as np
import synphot
from astropy import units as u
from astropy.coordinates import ICRS, Distance
from astropy.table import unique
from astropy_healpix import HEALPix, nside_to_level
from ligo.skymap import distance
from ligo.skymap.bayestar import rasterize
from m4opt import missions
from m4opt.fov import footprint_healpix
from m4opt.synphot import observing
from m4opt.synphot.extinction import DustExtinction
from scipy import stats


def get_detection_probability_known_position(plan, event_row, plan_args):
    if len(plan) == 0:
        return 0

    hpx = HEALPix(nside=plan_args["nside"], order="nested", frame=ICRS())
    mission = getattr(missions, plan_args["mission"])

    observations = plan[plan["action"] == "observe"].filled()
    coords = observations["target_coord"].to_table()
    coords["i"] = np.arange(len(coords))
    i = np.sort(unique(coords, keys=["ra", "dec"])["i"])
    fields = observations[i]

    target_ipix = hpx.lonlat_to_healpix(
        event_row["longitude"] * u.rad, event_row["latitude"] * u.rad
    )
    target_in_field = [
        target_ipix in footprint
        for footprint in footprint_healpix(
            hpx, mission.fov, fields["target_coord"], fields["roll"]
        )
    ]
    fields = fields[target_in_field]
    if len(fields) == 0:
        return 0

    with observing(
        observer_location=fields["observer_location"],
        target_coord=fields["target_coord"],
        obstime=(fields["start_time"] + 0.5 * fields["duration"]),
    ):
        spectrum = synphot.SourceSpectrum(
            synphot.ConstFlux1D, amplitude=0 * u.ABmag
        ) * synphot.SpectralElement(DustExtinction())
        limmag = mission.detector.get_limmag(
            plan_args["snr"], fields["duration"], spectrum, plan_args["bandpass"]
        ).max()
    lim_absmag = limmag - Distance(event_row["distance"] * u.Mpc).distmod
    return stats.norm(
        loc=plan_args["absmag_mean"], scale=plan_args["absmag_stdev"]
    ).cdf(lim_absmag.to_value(u.mag))


def get_detection_probability_unknown_position(plan, skymap_moc, plan_args):
    if len(plan) == 0:
        return 0

    hpx = HEALPix(nside=plan_args["nside"], order="nested", frame=ICRS())
    mission = getattr(missions, plan_args["mission"])

    skymap = rasterize(skymap_moc, order=nside_to_level(plan_args["nside"]))

    observations = plan[plan["action"] == "observe"].filled()
    coords = observations["target_coord"].to_table()
    coords["i"] = np.arange(len(coords))
    i = np.sort(unique(coords, keys=["ra", "dec"])["i"])
    fields = observations[i]

    durations = np.zeros(hpx.npix)
    for ipix, duration in zip(
        footprint_healpix(hpx, mission.fov, fields["target_coord"], fields["roll"]),
        fields["duration"].to_value(u.s),
    ):
        durations[ipix] = np.maximum(durations[ipix], duration)
    skymap["duration"] = durations * u.s
    skymap["ipix"] = np.arange(hpx.npix)

    skymap = skymap[durations > 0]

    with observing(
        observer_location=plan["observer_location"][0],
        target_coord=hpx.healpix_to_skycoord(skymap["ipix"]),
        obstime=plan["start_time"][0],
    ):
        spectrum = synphot.SourceSpectrum(
            synphot.ConstFlux1D, amplitude=0 * u.ABmag
        ) * synphot.SpectralElement(DustExtinction())
        skymap["limmag"] = mission.detector.get_limmag(
            plan_args["snr"], skymap["duration"], spectrum, plan_args["bandpass"]
        )
    skymap["limmag"][np.isnan(skymap["limmag"])] = -np.inf * u.mag

    distmean, diststd, distnorm = distance.parameters_to_moments(
        skymap["DISTMU"], skymap["DISTSIGMA"]
    )
    sigma2_log = np.log1p(np.square(diststd / distmean))
    logdistsigma = np.sqrt(sigma2_log)
    logdistmu = np.log(distmean) - 0.5 * sigma2_log

    absmagmu = plan_args["absmag_mean"]
    absmagsigma = plan_args["absmag_stdev"]
    a = 5 / np.log(10)
    appmagmu = absmagmu + a * logdistmu + 25
    appmagsigma = np.sqrt(np.square(absmagsigma) + np.square(a * logdistsigma))
    skymap["appmagmu"] = appmagmu
    skymap["appmagsigma"] = appmagsigma
    return (
        skymap["PROB"]
        * stats.norm(loc=skymap["appmagmu"], scale=skymap["appmagsigma"]).cdf(
            skymap["limmag"]
        )
    ).sum()
