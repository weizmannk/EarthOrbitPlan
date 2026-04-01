#!/usr/bin/env python
"""
mission_limmag.py — Compute per-field limmag for any M4OPT mission
-------------------------------------------------------------------

Usage
-----
    python mission_limmag.py --config earthorbitplan/config/params_ultrasat.ini

Configuration file example
--------------------------
    [params]
    input_table = ultrasat_bns.ecsv
    plan_dir    = data/O5
    output_dir  = add_limmag
"""

import argparse
import configparser
import logging
import sys
from pathlib import Path

import numpy as np
import synphot
from astropy import units as u
from astropy.coordinates import ICRS, EarthLocation
from astropy.table import QTable
from astropy.utils.masked import Masked
from astropy_healpix import HEALPix
from ligo.skymap.bayestar import rasterize
from ligo.skymap.distance import parameters_to_marginal_moments, parameters_to_moments
from ligo.skymap.io import read_sky_map
from m4opt import missions
from m4opt.fov import footprint_healpix
from m4opt.synphot import observing
from m4opt.synphot._math import countrate
from m4opt.synphot.background import update_missions
from m4opt.synphot.extinction import DustExtinction
from scipy import stats


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
        force=True,
    )


# ---------------------------------------------------------------------------
# Argument / config parsing
# ---------------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compute per-field limmag for any M4OPT mission."
    )
    parser.add_argument("--config", type=str, help="Path to .ini configuration file")
    args, remaining = parser.parse_known_args()

    if args.config:
        cfg_parser = configparser.ConfigParser()
        cfg_parser.read(args.config)
        cfg = cfg_parser["params"]
        return argparse.Namespace(
            input_table=cfg.get("input_table", fallback="bns.ecsv"),
            plan_dir=cfg.get("plan_dir", fallback="data/O5"),
            output_dir=cfg.get("output_dir", fallback="add_limmag"),
        )

    parser.add_argument("--input-table", type=str, default="bns.ecsv")
    parser.add_argument("--plan-dir",    type=str, default="data/O5")
    parser.add_argument("--output-dir",  type=str, default="add_limmag")
    return parser.parse_args(remaining)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
def _masked_col(values, n_obs):
    """Interleave observe/slew rows; mask the slew rows."""
    arr      = np.asarray(values, dtype=float)
    repeated = np.repeat(arr, 2)
    mask     = np.tile([False, True], n_obs)
    return Masked(repeated, mask=mask)[:-1]


# ---------------------------------------------------------------------------
# Single process function
# ---------------------------------------------------------------------------
def process(plan_args, mission, hpx, skymap, observations, footprints):
    """
    Compute limmag and detection metrics for every observation field.

    update_missions() is always called — it handles Cerenkov internally
    (injects AE8-based spectrum for ULTRASAT, does nothing for others).
    """
    snr              = plan_args["snr"]
    bandpass         = plan_args["bandpass"]
    bandpass_spectrum = mission.detector.bandpasses[bandpass]
    absmagmu         = plan_args["absmag_mean"]
    absmagsigma      = plan_args["absmag_stdev"]
    a                = 5 / np.log(10)

    flat_spectrum        = synphot.SourceSpectrum(synphot.ConstFlux1D, amplitude=0 * u.ABmag)
    dusty_flat_spectrum  = flat_spectrum * DustExtinction()
    observation_midtimes = observations["start_time"] + 0.5 * observations["duration"]

    distmean_all, diststd_all, _ = parameters_to_moments(
        skymap["DISTMU"], skymap["DISTSIGMA"]
    )

    rows = []

    for i_obs, (obs, obs_time, footprint) in enumerate(
        zip(observations, observation_midtimes, footprints)
    ):
        obs_location  = obs["observer_location"]
        n_pix         = len(footprint)
        duration_arr  = np.full(n_pix, obs["duration"].to_value(u.s)) * u.s

        obs_loc_tiled = EarthLocation.from_geocentric(
            *(
                np.tile(
                    np.array([coord.value for coord in obs_location.geocentric]),
                    (n_pix, 1),
                ).T
                * obs_location.geocentric[0].unit
            )
        )

        try:
            with observing(
                observer_location=obs_loc_tiled,
                target_coord=hpx.healpix_to_skycoord(footprint),
                obstime=obs_time + np.zeros(n_pix) * u.s,
            ):
                # Always call update_missions - handles Cerenkov automatically
                # (ULTRASAT: injects AE8 spectrum / other missions: no-op)
                update_missions(mission, obs_location, obs_time)

                limmag_no_dust = mission.detector.get_limmag(
                    snr, duration_arr, flat_spectrum, bandpass
                ).to_value(u.mag)

                limmag_dust = mission.detector.get_limmag(
                    snr, duration_arr, dusty_flat_spectrum, bandpass
                )

                sky_bg = (
                    countrate(mission.detector.background, bandpass_spectrum)
                    / countrate(flat_spectrum, bandpass_spectrum)
                ).to_value(u.mag(u.dimensionless_unscaled))

                dust_ext = (
                    countrate(DustExtinction() * flat_spectrum, bandpass_spectrum)
                    / countrate(flat_spectrum, bandpass_spectrum)
                ).to_value(u.mag(u.dimensionless_unscaled))

        except ValueError as e:
            logging.warning(f"Skipping obs {i_obs}: {e}")
            return None

        # Detection probability per pixel
        dm           = distmean_all[footprint]
        ds           = diststd_all[footprint]
        sigma2_log   = np.log1p(np.square(ds / dm))
        logdistsigma = np.sqrt(sigma2_log)
        logdistmu    = np.log(dm) - 0.5 * sigma2_log
        appmagmu_pix    = absmagmu + a * logdistmu + 25
        appmagsigma_pix = np.sqrt(
            np.square(absmagsigma) + np.square(a * logdistsigma)
        )
        det_prob_pix = skymap["PROB"][footprint] * stats.norm(
            loc=appmagmu_pix, scale=appmagsigma_pix
        ).cdf(limmag_dust)

        rows.append({
            "limmag_no_dust":         float(np.median(limmag_no_dust)),
            "limmag_dust":    float(np.median(limmag_dust)),
            "sky_background":         float(sky_bg),
            "dust":           float(dust_ext),
            "prob":           float(skymap["PROB"][footprint].sum()),
            "dist":           float(
                parameters_to_marginal_moments(
                    skymap[footprint]["PROB"] / skymap[footprint]["PROB"].sum(),
                    skymap[footprint]["DISTMU"],
                    skymap[footprint]["DISTSIGMA"],
                )[0]
            ),
            "detection_prob": float(np.sum(det_prob_pix)),
        })

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    setup_logging()
    args = parse_arguments()

    plan_dir   = Path(args.plan_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"Reading simulation IDs from: {args.input_table}")
    simulation_ids = QTable.read(args.input_table)["simulation_id"]
    logging.info(f"Found {len(simulation_ids)} simulations.")

    for sim_id in simulation_ids:
        logging.info(f"Processing simulation {sim_id}")

        plan_file = plan_dir / f"{sim_id}.ecsv"
        if not plan_file.exists():
            logging.warning(f"Missing plan file: {plan_file}, skipping.")
            continue

        plan      = QTable.read(plan_file)
        plan_args = plan.meta["args"]

        hpx     = HEALPix(nside=plan_args["nside"], order="nested", frame=ICRS())
        mission = getattr(missions, plan_args["mission"])

        skymap_moc = read_sky_map(plan_args["skymap"], moc=True)
        skymap     = rasterize(skymap_moc, order=hpx.level)

        observations = plan[plan["action"] == "observe"].filled()
        if len(observations) == 0:
            logging.warning(f"No observations in plan {sim_id}, skipping.")
            continue

        footprints = footprint_healpix(
            hpx, mission.fov, observations["target_coord"], observations["roll"]
        )

        rows = process(plan_args, mission, hpx, skymap, observations, footprints)

        if not rows:
            logging.warning(f"No valid rows for simulation {sim_id}, skipping.")
            continue

        n_obs = len(rows)
        for key in ["limmag_no_dust", "limmag_dust" "sky_background", "dust", "prob", "dist", "detection_prob"]:
            plan[key] = _masked_col([r[key] for r in rows], n_obs)

        plan["start_time"].precision = 0
        out_path = output_dir / f"{sim_id}_limmag.ecsv"
        plan.write(out_path, format="ascii.ecsv", overwrite=True)
        logging.info(f"  Saved → {out_path}")

    logging.info("Done.")


if __name__ == "__main__":
    main()
