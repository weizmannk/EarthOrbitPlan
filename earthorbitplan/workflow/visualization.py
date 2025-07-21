import warnings

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable
from matplotlib import pyplot as plt

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
warnings.filterwarnings("ignore", ".*dubious year.*")
warnings.filterwarnings(
    "ignore", "Tried to get polar motions for times after IERS data is valid.*"
)


def customize_style(columns=1):
    if columns == 1:
        target_width = 3.5  # ApJ column size in inches
    else:
        target_width = 7.25  # ApJ two-column text width in inches
    width, height = plt.rcParams["figure.figsize"]
    plt.style.use("seaborn-v0_8-paper")
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["figure.figsize"] = (target_width, height * target_width / width)


def load_events(events_file):
    """
    Load and prepare area-distance data from an ECSV file.

    Parameters
    ----------
    events_file : str or Path
        Path to the ECSV file containing the event data.

    Returns
    -------
    table : QTable
        The full event table.
    constants : dict
        Dictionary of constant simulation parameters.
    """
    table = QTable.read(events_file)

    # Define columns that must remain constant
    constant_cols = [
        "cutoff",
        "mission",
        "nside",
        "snr",
        "deadline",
        "delay",
        "exptime_min",
        "exptime_max",
        "bandpass",
        "absmag_mean",
        "absmag_stdev",
        "visits",
        "skydir",
    ]
    constants = {}
    for col in constant_cols:
        if col in table.colnames:
            if not np.all(table[col] == table[0][col]):
                raise ValueError(
                    f"Column '{col}' has multiple values; expected a constant."
                )
            constants[col] = table[0][col]

    return table, constants


def preprocess_events(main_table, constants, ns_max=3.0):
    r"""
    Preprocess the event table and select events for each run.

    This function processes the simulated events table, grouping by run and
    selecting events with an objective value above the cutoff threshold.
    It is intended for use with simulations that only include binary neutron
    star (BNS) and neutron star–black hole (NSBH) events, i.e., systems with
    secondary mass below a specified upper limit (default: 3.0, i.e., neutron star upper bound, \( M_\odot \)).
    This mass cut is consistent with the original simulation setup, and is applied here for completeness.

    This function selects only BNS and NSBH systems (i.e., systems with secondary mass <= ns_max),
    and provides both the full filtered sample and the sub-sample passing the objective_value cutoff.

    Parameters
    ----------
    main_table : `~astropy.table.QTable`
        The full event table as returned by :func:`load_events`.
    constants : dict
        Dictionary of constant simulation parameters as returned by :func:`load_events`.
    ns_max : float, optional
        Maximum secondary mass (in solar masses, default: 3.0 :math:`M_\odot`).

    tables : dict
        {run_name: table of BNS/NSBH events}
    selected_tables : dict
        {run_name: table of BNS/NSBH events with objective_value >= cutoff}
        and objective_value >= cutoff.
    """
    runs = np.unique(main_table["run"])
    tables = {}
    selected_tables = {}
    for run in runs:
        table = main_table[main_table["run"] == run]
        z = z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc).to_value(
            u.dimensionless_unscaled
        )
        m2 = table["mass2"] / (1 + z)
        # Filter: secondary mass <= ns_max
        keep = m2 <= ns_max
        tables[run] = table

        # Filter objective_value >= cutoff (and mass)
        selected_tables[run] = table[
            keep & (table["objective_value"] >= constants["cutoff"])
        ]
    return runs, tables, selected_tables
