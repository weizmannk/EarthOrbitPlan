# This probability is from https://github.com/m4opt/m4opt-paper

import warnings

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
from m4opt.synphot.background import update_missions
from scipy import stats
from tqdm import tqdm

warnings.filterwarnings("ignore", ".*Wswiglal-redir-stdio.*")
warnings.filterwarnings("ignore", ".*dubious year.*")
warnings.filterwarnings("ignore", ".*polar motions.*")


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

    # with observing(
    #     observer_location=fields["observer_location"],
    #     target_coord=fields["target_coord"],
    #     obstime=(fields["start_time"] + 0.5 * fields["duration"]),
    # ):
        
    #     # Update mission parameters with contextual background information
    #     # Compute limmag per field: the Cerenkov background depends on the satellite's
    #     # position in the radiation belts (AE8 model), which varies along the orbit.
    #     update_missions(
    #         mission, 
    #         fields["observer_location"][0], 
    #         (fields["start_time"] + 0.5 * fields["duration"])[0],
    #     )
        
        
    #     spectrum = synphot.SourceSpectrum(
    #         synphot.ConstFlux1D, amplitude=0 * u.ABmag
    #     ) * synphot.SpectralElement(DustExtinction())
    #     limmag = mission.detector.get_limmag(
    #         plan_args["snr"], fields["duration"], spectrum, plan_args["bandpass"]
    #     ).max()

    def _limmag_obs_time(field, mission, plan_args) -> float:
        """
        Evaluate limiting magnitude at a single orbital position.
        
        Takes into account:
        - Observer location in orbit (affects Cerenkov background via AE8 model)
        - Observation time (affects SNR accumulation)
        - Dust extinction along line of sight
        """
        with observing(
            observer_location=field["observer_location"],
            target_coord=field["target_coord"],
            obstime=(field["start_time"] + 0.5 * field["duration"]),
        ):
            # Update mission parameters based on orbital position
            # (Cerenkov background varies with radiation belt exposure)
            update_missions(
                mission, 
                field["observer_location"], 
                (field["start_time"] + 0.5 * field["duration"]),
            )
        
            # Create spectrum with dust extinction
            spectrum = synphot.SourceSpectrum(
                synphot.ConstFlux1D, amplitude=0 * u.ABmag
            ) * synphot.SpectralElement(DustExtinction())
            
            # Compute limiting magnitude for this field
            limmag = mission.detector.get_limmag(
                plan_args["snr"], 
                field["duration"],
                spectrum, 
                plan_args["bandpass"]
            )
        return limmag

    # Compute limiting magnitude for all fields covering the true position
    # Take the law  detectable limiting magnitude across all observations
    limmag = max(
        _limmag_obs_time(field, mission, plan_args) 
        for field in tqdm(
            fields,
            desc=f"Computing limmag ({mission.name})",
            total=len(fields),
            unit="field",
        )
    )

    # Convert to absolute magnitude limit at source distance
    lim_absmag =  limmag - Distance(event_row["distance"] * u.Mpc).distmod
    # Return probability that intrinsic magnitude is brighter than limit
    # (based on Gaussian distribution of kilonova absolute magnitudes)
    return stats.norm(
        loc=plan_args["absmag_mean"], scale=plan_args["absmag_stdev"]
    ).cdf(lim_absmag.to_value(u.mag))
    

def get_detection_probability_unknown_position(plan, skymap_moc, plan_args):
    """
    Estimate the detection probability of a transient with an uncertain sky position,
    described by a gravitational wave probability skymap.

    This function computes the overall detection probability by integrating over all
    possible sky positions and distances from the skymap. For each sky location, it
    considers the total exposure from the observation plan, computes the limiting
    magnitude, and estimates the chance that a transient with an assumed absolute
    magnitude distribution would be detectable at that distance. The final probability
    is the sum of the detection probabilities weighted by the skymap's localization
    probabilities.

    Parameters
    ----------
    plan : astropy.table.Table
        Table of scheduled observations, including pointing coordinates, times, and durations.
    skymap_moc :
        Gravitational  wave localization skymap, encoding the sky probability
        and distance estimates per pixel.
    plan_args : dict
        Dictionary of observation and mission parameters (e.g., nside, mission name, SNR threshold,
        bandpass, mean and standard deviation of absolute magnitude).

    Returns
    -------
    probability : float
        Integrated probability that the transient would be detected by the observation plan,
        marginalized over the skymap localization.
    """

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
        # Update mission parameters with contextual background information
        # Compute limmag per field: the Cerenkov background depends on the satellite's
        # position in the radiation belts (AE8 model), which varies along the orbit.
        update_missions(mission, plan["observer_location"][0], plan["start_time"][0])
        
        
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
