# """
# Farah GWTC-3 Distribution Processor
# ===================================

# This script converts the Farah GWTC-3 population synthesis distribution into a format
# suitable for `bayestar-inject`, computes astrophysical merger rates, and downloads
# the required file if missing.

# Source:
#     LIGO-T2100512/public/ (https://dcc.ligo.org/LIGO-T2100512)

# Usage
# -----
# Run directly from the command line to generate a processed HDF5 file and compute merger rates.

# Example:
#     $ python farah_processor.py --outdir ./scenarios
# """

# import argparse
# import logging
# import os
# import sys

# import numpy as np
# import requests
# import scipy.stats as stats
# from astropy.table import Table
# from rate_stats import format_with_errorbars, poisson_lognormal_rate_quantiles
# from tqdm import tqdm

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(message)s",
#     handlers=[
#         logging.FileHandler("farah_processing.log"),
#         logging.StreamHandler(sys.stdout),
#     ],
#     force=True,
# )


# def parse_arguments():
#     parser = argparse.ArgumentParser(description="Process Farah GWTC-3 distribution.")
#     parser.add_argument(
#         "-o",
#         "--outdir",
#         type=str,
#         default="./scenarios",
#         help="Output directory to store processed results.",
#     )
#     return parser.parse_args()


# def merger_rate(farah_file, ns_max_mass=3, quantiles=[0.05, 0.5, 0.95]):
#     """
#     Compute astrophysical merger rates from a processed Farah GWTC-3 file.

#     Parameters
#     ----------
#     farah_file : str
#         Path to HDF5 file containing mass1 and mass2.
#     ns_max_mass : float
#         Maximum mass to consider an object a neutron star.
#     quantiles : list
#         Quantiles to compute.

#     Returns
#     -------
#     astropy.table.Table
#         Table with log-normal parameters for each population.
#     tuple of str
#         Formatted rate estimates with uncertainties.
#     """
#     sim_rate = 2.712359951521142e3  # Gpc^-3 yr^-1
#     run_duration = 1.0  # years

#     rates_table = Table(
#         [
#             {"population": "BNS", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "NSBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "BBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#         ]
#     )

#     table = Table.read(farah_file)
#     m1 = table["mass1"]
#     m2 = table["mass2"]
#     rates_table["mass_fraction"] = np.array(
#         [
#             np.sum((m1 < ns_max_mass) & (m2 < ns_max_mass)),
#             np.sum((m1 >= ns_max_mass) & (m2 < ns_max_mass)),
#             np.sum((m1 >= ns_max_mass) & (m2 >= ns_max_mass)),
#         ]
#     ) / len(table)

#     for key in ["lower", "mid", "upper"]:
#         rates_table[key] *= rates_table["mass_fraction"]

#     delta = np.diff(stats.norm.interval(0.9))[0]
#     rates_table["mu"] = (
#         np.log(rates_table["mid"]) + np.log(run_duration) - np.log(sim_rate)
#     )
#     rates_table["sigma"] = (
#         np.log(rates_table["upper"]) - np.log(rates_table["lower"])
#     ) / delta

#     lo_n, mid_n, hi_n = poisson_lognormal_rate_quantiles(
#         quantiles, rates_table["mu"], rates_table["sigma"]
#     )

#     lo = lo_n * sim_rate / run_duration
#     mid = mid_n * sim_rate / run_duration
#     hi = hi_n * sim_rate / run_duration

#     return rates_table, format_with_errorbars(mid, lo, hi)


# def download_file(file_url, file_name):
#     if os.path.exists(file_name):
#         logging.info(f"File already exists: {file_name}")
#         return file_name

#     try:
#         response = requests.get(file_url, stream=True)
#         response.raise_for_status()
#         total = int(response.headers.get("content-length", 0))
#         with (
#             open(file_name, "wb") as f,
#             tqdm(
#                 total=total, unit="B", unit_scale=True, desc="Downloading"
#             ) as progress,
#         ):
#             for chunk in response.iter_content(chunk_size=8192):
#                 f.write(chunk)
#                 progress.update(len(chunk))
#         return file_name
#     except requests.RequestException as e:
#         logging.error(f"Error downloading the file: {e}")
#         return None


# if __name__ == "__main__":
#     args = parse_arguments()
#     output_dir = os.path.abspath(args.outdir)
#     os.makedirs(output_dir, exist_ok=True)
#     farah_file = os.path.join(output_dir, "farah.h5")

#     if not os.path.exists(farah_file):
#         url = "https://dcc.ligo.org/LIGO-T2100512/public/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"
#         tmp_file = download_file(url, os.path.join(output_dir, url.split("/")[-1]))
#         if tmp_file is None:
#             sys.exit("Download failed.")
#         data = Table.read(tmp_file)
#         Table(
#             {
#                 "mass1": data["mass_1"],
#                 "mass2": data["mass_2"],
#                 "spin1z": data["a_1"] * data["cos_tilt_1"],
#                 "spin2z": data["a_2"] * data["cos_tilt_2"],
#             }
#         ).write(farah_file, overwrite=True)
#         os.remove(tmp_file)

#     rate_table, rate_stats = merger_rate(farah_file)
#     logging.info(f"Merger rate stats: {rate_stats}")


# """
# Farah GWTC-3 Distribution Processor
# ===================================

# This script converts the Farah GWTC-3 population synthesis distribution into a format
# suitable for `bayestar-inject`, computes astrophysical merger rates, and downloads
# the required file if missing.

# Source:
#     LIGO-T2100512/public/ (https://dcc.ligo.org/LIGO-T2100512)

# Usage
# -----
# Run directly from the command line to generate a processed HDF5 file and compute merger rates.

# Example:
#     $ python farah_processor.py --outdir ./scenarios
# """

# import logging
# import os
# import sys

# from astropy.table import Table

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(message)s",
#     handlers=[
#         logging.FileHandler("farah_processing.log"),
#         logging.StreamHandler(sys.stdout),
#     ],
#     force=True,
# )


# def parse_arguments():
#     """
#     Parse command-line arguments.

#     Returns
#     -------
#     argparse.Namespace
#         Parsed command-line arguments with attributes like `outdir`.
#     """
#     parser = argparse.ArgumentParser(
#         description="Run the ULTRASAT workflow and submit it using HTCondor."
#     )
#     parser.add_argument(
#         "-o",
#         "--outdir",
#         type=str,
#         default="./scenarios",
#         help="Path to the params file.",
#     )
#     return parser.parse_args()


# def merger_rate(
#     farah_file,
#     ns_max_mass=3,
#     quantiles=[0.05, 0.5, 0.95],
#     duration_yr=1.0,
#     rate_sim=None,
# ):
#     """
#     Compute astrophysical merger rates from a processed Farah GWTC-3 file.

#     Parameters
#     ----------
#     farah_file : str
#         Path to HDF5 file containing mass1 and mass2.
#     ns_max_mass : float, optional
#         Threshold mass to define neutron stars (default is 3).
#     quantiles : list of float, optional
#         Quantile levels to compute (default is [0.05, 0.5, 0.95]).
#     duration_yr : float
#         Duration of the observing run in years. Default is 1.
#     rate_sim : float or None
#         Simulated effective rate in Gpc^-3 yr^-1. If provided, will convert counts to rates.

#     Returns
#     -------
#     astropy.table.Table
#         Table with populations (BNS, NSBH, BBH) and their log-normal parameters.
#     tuple of float or tuple of str
#         Quantiles of detected counts, or converted to rates if rate_sim is given.
#     """
#     rates_table = Table(
#         [
#             {"population": "BNS", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "NSBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "BBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#         ]
#     )

#     table = Table.read(farah_file)
#     source_mass1 = table["mass1"]
#     source_mass2 = table["mass2"]
#     rates_table["mass_fraction"] = np.array(
#         [
#             np.sum((source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)),
#             np.sum((source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)),
#             np.sum((source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)),
#         ]
#     ) / len(table)

#     for key in ["lower", "mid", "upper"]:
#         rates_table[key] *= rates_table["mass_fraction"]

#     standard_90pct_interval = np.diff(stats.norm.interval(0.9))[0]
#     rates_table["mu"] = np.log(rates_table["mid"])
#     rates_table["sigma"] = (
#         np.log(rates_table["upper"]) - np.log(rates_table["lower"])
#     ) / standard_90pct_interval

#     lo_n, mid_n, hi_n = poisson_lognormal_rate_quantiles(
#         quantiles, rates_table["mu"], rates_table["sigma"]
#     )

#     if rate_sim is not None:
#         lo_r = lo_n * rate_sim / duration_yr
#         mid_r = mid_n * rate_sim / duration_yr
#         hi_r = hi_n * rate_sim / duration_yr
#         return rates_table, (lo_r, mid_r, hi_r)

#     return rates_table, (lo_n, mid_n, hi_n)

#     # Compute astrophysical merger rates from a processed Farah GWTC-3 file.

#     # O3 R&P paper Table II row 1 last column:
#     # 5%, 50%, and 95% quantiles of the total merger rate
#     # in Gpc^-3 yr^-1.
#     # See https://doi.org/10.1103/PhysRevX.13.011048

#     # === Simulated BNS Merger Rates ===
#     # Use this command to retrieve comments:
#     # 1- sqlite3 events.sqlite
#     # 2- select comment from process;
#     # Simulated BNS merger rate in yr^-1 Gpc^-3 for O5 and O6-HLVK configuration (SNR = 10)
#     # From kiendrebeogo et al. 2023 the simulated rate is given by the by (yr^-1 Mpc^-3)
#     # so this need to be convert in yr^-1 Gpc^-3, before add use it here.
#     # simulation rate  sim_rate =2.712359951521142e3 (u.Gpc**-3 * u.yr**-1)

#     rates_table = Table(
#         [
#             {"population": "BNS", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "NSBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#             {"population": "BBH", "lower": 100.0, "mid": 240.0, "upper": 510.0},
#         ]
#     )

#     # lo = 100
#     # mid = 240
#     # hi = 510

#     # (standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))
#     # log_target_rate_mu = np.log(mid)
#     # log_target_rate_sigma = np.log(hi / lo) / standard_90pct_interval
#     # log_target_rate_mu, log_target_rate_sigma

#     # log_simulation_effective_rate_by_run = {
#     #     key: np.log(value.to_value(u.Gpc**-3 * u.yr**-1))
#     #     for key, value in main_table.meta["effective_rate"].items()
#     # }

#     table = Table.read(farah_file)
#     source_mass1 = table["mass1"]
#     source_mass2 = table["mass2"]
#     rates_table["mass_fraction"] = np.array(
#         [
#             np.sum((source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)),
#             np.sum((source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)),
#             np.sum((source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass)),
#         ]
#     ) / len(table)

#     for key in ["lower", "mid", "upper"]:
#         rates_table[key] *= rates_table["mass_fraction"]

#     (standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))
#     rates_table["mu"] = np.log(
#         rates_table["mid"]
#     )  # + np.log(run_duration) - np.log(sim_rate)
#     rates_table["sigma"] = (
#         np.log(rates_table["upper"]) - np.log(rates_table["lower"])
#     ) / standard_90pct_interval

#     lo, mid, hi = poisson_lognormal_rate_quantiles(
#         quantiles, rates_table["mu"], rates_table["sigma"]
#     )
#     mid, lo, hi = format_with_errorbars(mid, lo, hi)

#     return rates_table, (lo, mid, hi)


# def download_file(file_url, file_name):
#     """
#     Download a file from the specified URL with a progress bar.

#     Parameters
#     ----------
#     file_url : str
#         URL to download the file from.
#     file_name : str
#         Local file path to save the downloaded file.

#     Returns
#     -------
#     str or None
#         Path to the downloaded file or None if the download failed.
#     """
#     if os.path.exists(file_name):
#         logging.info(f"File already exists: {file_name}")
#         return file_name

#     try:
#         response = requests.get(file_url, stream=True)
#         response.raise_for_status()
#         file_size = int(response.headers.get("content-length", 0))

#         with (
#             open(file_name, "wb") as file,
#             tqdm(
#                 total=file_size, unit="B", unit_scale=True, desc="Downloading"
#             ) as progress,
#         ):
#             for chunk in response.iter_content(chunk_size=8192):
#                 file.write(chunk)
#                 progress.update(len(chunk))

#         logging.info(f"File downloaded successfully to {file_name}.")
#         return file_name
#     except requests.RequestException as e:
#         logging.error(f"Error downloading the file: {e}")
#         return None


# if __name__ == "__main__":
#     args = parse_arguments()
#     output_dir = os.path.abspath(args.outdir)
#     os.makedirs(output_dir, exist_ok=True)
#     farah_file = os.path.join(output_dir, "farah.h5")

#     if not os.path.exists(farah_file):
#         file_url = "https://dcc.ligo.org/LIGO-T2100512/public/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"
#         file_name = os.path.join(output_dir, file_url.split("/")[-1])
#         input_file = download_file(file_url, file_name)

#         if input_file is None:
#             logging.error("Failed to download the required file.")
#             sys.exit(1)

#         data = Table.read(input_file)
#         Table(
#             {
#                 "mass1": data["mass_1"],
#                 "mass2": data["mass_2"],
#                 "spin1z": data["a_1"] * data["cos_tilt_1"],
#                 "spin2z": data["a_2"] * data["cos_tilt_2"],
#             }
#         ).write(farah_file, overwrite=True)

#         os.remove(input_file)
#         logging.info(f"Removed temporary file: {input_file}")
#         logging.info(
#             f"number of BNS: {len(data[data['mass_1'] < 3])}, number of NSBH: "
#             f"{len(data[(data['mass_1'] >= 3) & (data['mass_2'] < 3)])}, number of BBH: "
#             f"{len(data[data['mass_2'] >= 3])}"
#         )
#         logging.info(f"Processed file saved at: {farah_file}")

#     ns_max_mass = 3
#     quantiles = [0.05, 0.5, 0.95]
#     logging.info("Computing Astrophysical Merger Rate...")
#     rate_table, poison_stats = merger_rate(farah_file, ns_max_mass, quantiles)
#     logging.info(f"Merger rates computed:\n{rate_table}")

#     print(poison_stats)
