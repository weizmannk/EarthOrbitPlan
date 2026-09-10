# """
# Zenodo GW Injection Data Unpacker
# =================================

# This module automates the unpacking, filtering, and conversion of injection datasets
# (GWTC-5.0, FullPop or PixelPop population model) from Zenodo ZIP archives. It
# processes event tables and associated localization files for specific observing
# runs (e.g., O5a, IR1HL), and outputs filtered ECSV tables and organized FITS files.

# Usage
# -----
# Run from the command line:

#     python unpacker.py --zip runs.zip --subdir runs --pop fullpop --runs IR1HL IR1HLV O5a O5b O5c  --data-dir ./data --mass-threshold 3 --skymap-dir skymaps

# Or use a config file (must contain a [params] section):

#     python earthorbitplan/workflow/unpacker.py --config  earthorbitplan/scenarios/params_scenarios.ini

# Or import and call `process_zip()` in your Python code.

# Source data (concept DOIs, always resolve to the latest version):
#     FullPop  : https://doi.org/10.5281/zenodo.22550047
#     PixelPop : https://doi.org/10.5281/zenodo.22555948
# """

import argparse
import configparser
import pathlib
import sqlite3
import zipfile
from functools import reduce
from shutil import copyfileobj
from tempfile import NamedTemporaryFile

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable, join, vstack
from tqdm.auto import tqdm


def parse_arguments():
    """
    Parse command-line arguments or ini config file.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Unpack and process GW injection ZIP data.",
        allow_abbrev=False,
    )
    parser.add_argument("--config", type=str, help="Path to .ini config file")
    args, remaining_args = parser.parse_known_args()

    if args.config:
        config = configparser.ConfigParser(inline_comment_prefixes=("#",))
        config.read(args.config)

        def get(key, fallback=None, sections=("scenarios", "paths", "params")):
            """Read a key from the first section that defines it.

            The mission configs split settings across [scenarios] and [paths];
            params_scenarios.ini keeps everything in a flat [params] section.
            """
            for section in sections:
                v = config.get(section, key, fallback=None)
                if v is not None:
                    return v
            return fallback

        return argparse.Namespace(
            zip=get("zip"),
            subdir=get("subdir", "runs"),
            pop=get("pop", "fullpop"),
            runs=get("runs", "IR1HL IR1HLV O5a O5b O5c").split(),
            detectors=get("detectors", ""),
            data_dir=get("data_dir", "data"),
            skymap_dir=get("skymap_dir", "skymaps"),
            mass_threshold=float(get("mass_threshold", 3.0)),
        )

    parser.add_argument(
        "--zip", type=str, required=True, help="Path to the ZIP archive."
    )
    parser.add_argument(
        "--subdir",
        type=str,
        default="runs",
        help="Subdirectory inside the ZIP archive.",
    )
    parser.add_argument(
        "--pop",
        type=str,
        choices=["fullpop", "pixelpop"],
        default="fullpop",
        help="Population-model subdirectory inside each run (default: fullpop).",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        default=["IR1HL", "IR1HLV", "O5a", "O5b", "O5c"],
        help="Observation runs to process (folder names inside <subdir>/).",
    )
    parser.add_argument(
        "--detectors",
        type=str,
        default="",
        help="Detector combination tag (e.g., HLVK).",
    )
    parser.add_argument(
        "--data-dir", type=str, default="data", help="Directory to store output."
    )
    parser.add_argument(
        "--skymap-dir", type=str, default="skymaps", help="Directory to store output."
    )
    parser.add_argument(
        "--mass-threshold",
        type=float,
        default=3.0,
        help="Maximum secondary mass for filtering.",
    )
    return parser.parse_args(remaining_args)


def process_zip(
    zip_path,
    runs,
    outdir,
    skymap_dir,
    max_ns_mass=3.0,
    subdir="runs",
    detectors="None",
    pop="fullpop",
):
    """
    Extract and filter GW injection tables from a Zenodo-style ZIP archive.

    Parameters
    ----------
    zip_path : str or Path
        Path to the input ZIP archive.
    runs : list of str
        Observing run labels to extract (e.g., ["IR1HL", "O5a"]).
    outdir : str or Path
        Destination directory for output files.
    skymap_dir : str or Path
        Directory for output skymap FITS files.
    max_ns_mass : float, optional
        Maximum mass for secondary object (used to filter BNS/NSBH).
    subdir : str, optional
        Name of the root folder inside the ZIP archive.
    detectors : str, optional
        Detector label used in file path construction (e.g., HLVK).
    pop : str, optional
        Population-model folder inside each run ("fullpop" or "pixelpop").

    Outputs
    -------
    - ECSV summary table at <outdir>/observing-scenarios.ecsv
    - FITS localization files under <outdir>/<skymap_dir>/<run>/ directories
    """

    out_root = pathlib.Path(outdir)
    skymap_root = out_root / skymap_dir
    for run in runs:
        (skymap_root / run).mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as archive:
        in_root = zipfile.Path(archive) / subdir

        tables = []
        for run in tqdm(runs, desc="Reading summary tables"):
            in_run = in_root / f"{run}{detectors}" / pop
            table = reduce(
                join,
                (
                    QTable.read((in_run / filename).read_text(), format="ascii")
                    for filename in ["coincs.dat", "allsky.dat", "injections.dat"]
                ),
            )

            table["run"] = run
            table.meta.clear()

            with (
                (in_run / "events.sqlite").open("rb") as in_file,
                NamedTemporaryFile() as out_file,
            ):
                copyfileobj(in_file, out_file)
                out_file.flush()
                with sqlite3.connect(f"file:{out_file.name}?mode=ro", uri=True) as db:
                    ((comment,),) = db.execute(
                        "SELECT comment FROM process WHERE program = 'bayestar-inject'"
                    )
            table.meta["effective_rate"] = {run: u.Quantity(comment)}

            tables.append(table)

        table = vstack(tables)
        del tables

        # filter the BNS and NSBH event with NS max of 3 Sun mass
        table["redshift"] = z_at_value(
            cosmo.luminosity_distance, table["distance"] * u.Mpc
        ).to_value(u.dimensionless_unscaled)

        table = table[table["mass2"] <= max_ns_mass * (1 + table["redshift"])]
        table["source_class"] = np.where(
            table["mass1"] <= max_ns_mass * (1 + table["redshift"]), "BNS", "NSBH"
        )

        table.write(f"{out_root}/observing-scenarios.ecsv", overwrite=True)

        # Copy FITS skymaps
        for row in tqdm(table, desc="Copying FITS files"):
            filename = f"{row['coinc_event_id']}.fits"
            in_path = in_root / f"{row['run']}{detectors}" / pop / "allsky" / filename
            out_path = skymap_root / row["run"] / filename
            with in_path.open("rb") as in_file, out_path.open("wb") as out_file:
                copyfileobj(in_file, out_file)


if __name__ == "__main__":
    args = parse_arguments()
    process_zip(
        zip_path=args.zip,
        runs=args.runs,
        outdir=args.data_dir,
        skymap_dir=args.skymap_dir,
        max_ns_mass=args.mass_threshold,
        subdir=args.subdir,
        detectors=args.detectors,
        pop=args.pop,
    )
