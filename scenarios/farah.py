"""
Farah GWTC-3 Distribution Processor
===================================

This script converts the Farah GWTC-3 population synthesis distribution into a format
suitable for `bayestar-inject`, computes astrophysical merger rates, and downloads
the required file if missing.

Source:
    LIGO-T2100512/public/ (https://dcc.ligo.org/LIGO-T2100512)

Usage
-----
Run directly from the command line to generate a processed HDF5 file and compute merger rates.

Example:
    $ python farah_processor.py --outdir ./scenarios
"""

import argparse
import logging
import os
import sys

import numpy as np
import requests
import scipy.stats as stats
from astropy.table import Table
from scipy import integrate, optimize, special
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("farah_processing.log"),
        logging.StreamHandler(sys.stdout),
    ],
    force=True,
)


def parse_arguments():
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments with attributes like `outdir`.
    """
    parser = argparse.ArgumentParser(
        description="Run the ULTRASAT workflow and submit it using HTCondor."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        default="./scenarios",
        help="Path to the params file.",
    )
    return parser.parse_args()


def betabinom_k_n(k, n):
    """
    Compute the beta-binomial distribution given k and n.

    Parameters
    ----------
    k : int
        Number of observed successes.
    n : int
        Total number of trials.

    Returns
    -------
    scipy.stats._distn_infrastructure.rv_frozen
        A frozen beta-binomial distribution instance.
    """
    return stats.betabinom(n, k + 1, n - k + 1)


@np.vectorize
def poisson_lognormal_rate_cdf(k, mu, sigma):
    """
    Compute the marginal CDF of a Poisson distribution with a log-normal prior on the rate.

    Parameters
    ----------
    k : float
        Count value.
    mu : float
        Mean of the log of the Poisson rate.
    sigma : float
        Standard deviation of the log of the Poisson rate.

    Returns
    -------
    float
        The CDF evaluated at k.
    """
    lognorm_pdf = stats.lognorm(s=sigma, scale=np.exp(mu)).pdf

    def func(lam):
        prior = lognorm_pdf(lam)
        # poisson_pdf = np.exp(special.xlogy(k, lam) - special.gammaln(k + 1) - lam)
        poisson_cdf = special.gammaincc(k + 1, lam)
        return poisson_cdf * prior

    cdf, _ = integrate.quad(func, 0, np.inf, epsabs=0)
    return cdf


@np.vectorize
def poisson_lognormal_rate_quantiles(p, mu, sigma):
    """
    Find quantiles of a Poisson distribution with a log-normal prior on the rate.

    Parameters
    ----------
    p : float
        Quantile level (e.g., 0.05, 0.5, 0.95).
    mu : float
        Mean of the log of the Poisson rate.
    sigma : float
        Standard deviation of the log of the Poisson rate.

    Returns
    -------
    float
        The event count at the specified quantile.

    Notes
    -----
    Uses `scipy.optimize.root_scalar` to solve for quantiles.
    """

    def func(k):
        return poisson_lognormal_rate_cdf(k, mu, sigma) - p

    if func(0) >= 0:
        return 0

    result = optimize.root_scalar(func, bracket=[0, 1e6])
    return result.root


def merger_rate(farah_file, ns_max_mass=3, quantiles=[0.05, 0.5, 0.95]):
    """
    Compute astrophysical merger rates from a processed Farah GWTC-3 file.

    Parameters
    ----------
    farah_file : str
        Path to HDF5 file containing mass1 and mass2.
    ns_max_mass : float, optional
        Threshold mass to define neutron stars (default is 3).
    quantiles : list of float, optional
        Quantile levels to compute (default is [0.05, 0.5, 0.95]).

    Returns
    -------
    astropy.table.Table
        Table with populations (BNS, NSBH, BBH) and their log-normal parameters.
    tuple of float
        Lower, median, and upper rate quantiles.
    """
    rates_table = Table(
        [
            {"population": "BNS", "lower": 100.0, "mid": 240.0, "upper": 510.0},
            {"population": "NSBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
            {"population": "BBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
        ]
    )

    table = Table.read(farah_file)
    source_mass1 = table["mass1"]
    source_mass2 = table["mass2"]
    rates_table["mass_fraction"] = np.array(
        [
            np.sum((source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)),
            np.sum((source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)),
            np.sum((source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)),
        ]
    ) / len(table)

    for key in ["lower", "mid", "upper"]:
        rates_table[key] *= rates_table["mass_fraction"]

    standard_90pct_interval = np.diff(stats.norm.interval(0.9))[0]
    rates_table["mu"] = np.log(rates_table["mid"])
    rates_table["sigma"] = (
        np.log(rates_table["upper"]) - np.log(rates_table["lower"])
    ) / standard_90pct_interval

    lo, mid, hi = poisson_lognormal_rate_quantiles(
        quantiles, rates_table["mu"], rates_table["sigma"]
    )

    return rates_table, (lo, mid, hi)


def download_file(file_url, file_name):
    """
    Download a file from the specified URL with a progress bar.

    Parameters
    ----------
    file_url : str
        URL to download the file from.
    file_name : str
        Local file path to save the downloaded file.

    Returns
    -------
    str or None
        Path to the downloaded file or None if the download failed.
    """
    if os.path.exists(file_name):
        logging.info(f"File already exists: {file_name}")
        return file_name

    try:
        response = requests.get(file_url, stream=True)
        response.raise_for_status()
        file_size = int(response.headers.get("content-length", 0))

        with (
            open(file_name, "wb") as file,
            tqdm(
                total=file_size, unit="B", unit_scale=True, desc="Downloading"
            ) as progress,
        ):
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
                progress.update(len(chunk))

        logging.info(f"File downloaded successfully to {file_name}.")
        return file_name
    except requests.RequestException as e:
        logging.error(f"Error downloading the file: {e}")
        return None


if __name__ == "__main__":
    args = parse_arguments()
    output_dir = os.path.abspath(args.outdir)
    os.makedirs(output_dir, exist_ok=True)
    farah_file = os.path.join(output_dir, "farah.h5")

    if not os.path.exists(farah_file):
        file_url = "https://dcc.ligo.org/LIGO-T2100512/public/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"
        file_name = os.path.join(output_dir, file_url.split("/")[-1])
        input_file = download_file(file_url, file_name)

        if input_file is None:
            logging.error("Failed to download the required file.")
            sys.exit(1)

        data = Table.read(input_file)
        Table(
            {
                "mass1": data["mass_1"],
                "mass2": data["mass_2"],
                "spin1z": data["a_1"] * data["cos_tilt_1"],
                "spin2z": data["a_2"] * data["cos_tilt_2"],
            }
        ).write(farah_file, overwrite=True)

        os.remove(input_file)
        logging.info(f"Removed temporary file: {input_file}")
        logging.info(
            f"number of BNS: {len(data[data['mass_1'] < 3])}, number of NSBH: "
            f"{len(data[(data['mass_1'] >= 3) & (data['mass_2'] < 3)])}, number of BBH: "
            f"{len(data[data['mass_2'] >= 3])}"
        )
        logging.info(f"Processed file saved at: {farah_file}")

    ns_max_mass = 3
    quantiles = [0.05, 0.5, 0.95]
    logging.info("Computing Astrophysical Merger Rate...")
    rate_table, poison_stats = merger_rate(farah_file, ns_max_mass, quantiles)
    logging.info(f"Merger rates computed:\n{rate_table}")

    print(poison_stats)
