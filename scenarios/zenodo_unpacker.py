"""
Zenodo GW Injection Data Unpacker
=================================

This module automates the unpacking, filtering, and conversion of gravitational-wave injection datasets
(e.g., Farah / GWTC-3) from Zenodo ZIP archives. It processes event tables and associated
localization files for specific observing runs (e.g., O5, O6), and outputs
filtered ECSV tables and organized FITS localization files.

This tool is adapted from:
    https://github.com/m4opt/m4opt-paper/blob/main/scripts/unpack-observing-scenarios.py

Usage
-----

Run directly from the command line by specifying parameters manually:

    $ python zenodo_unpacker.py --zip runs_SNR-10.zip --subdir runs_SNR-10 --runs O5 O6 --detectors HLVK --outdir ./data  --outfile observing-scenarios.ecsv --mass-threshold 3

Or run using a configuration file:

    $ python zenodo_unpacker.py --config params.ini

Or integrate into a larger Python pipeline by calling the `process_zip()` function directly.

Source
------

Injection datasets are available at:
    https://zenodo.org/records/14585837
"""

import argparse
import configparser
import pathlib
import sqlite3
import zipfile
from functools import reduce
from shutil import copyfileobj
from tempfile import NamedTemporaryFile

from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable, join, vstack
from tqdm.auto import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Unpack and process GW injection ZIP data.",
        allow_abbrev=False,
    )
    parser.add_argument("--config", type=str, help="Path to .ini config file")
    args, remaining_args = parser.parse_known_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)
        cfg = config["params"]

        return argparse.Namespace(
            zip=cfg.get("zip"),
            subdir=cfg.get("subdir", fallback="runs_SNR-10"),
            runs=cfg.get("runs", fallback="O5 O6").split(),
            detectors=cfg.get("detectors", fallback="HLVK"),
            outdir=cfg.get("outdir", fallback="data"),
            outfile=cfg.get("outfile", fallback="observing-scenarios.ecsv"),
            mass_threshold=cfg.getfloat("mass_threshold", fallback=3.0),
        )

    # If no config file, define arguments normally
    parser.add_argument(
        "--zip", type=str, required=True, help="Path to the ZIP archive."
    )
    parser.add_argument(
        "--subdir",
        type=str,
        default="runs_SNR-10",
        help="Subdirectory inside the ZIP archive.",
    )
    parser.add_argument(
        "--runs", nargs="+", default=["O5", "O6"], help="Observation runs to process."
    )
    parser.add_argument(
        "--detectors",
        type=str,
        default="HLVK",
        help="Detector combination tag (e.g., HLVK).",
    )
    parser.add_argument(
        "--outdir", type=str, default="data", help="Directory to store output."
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default="observing-scenarios.ecsv",
        help="Name of the output ECSV summary table.",
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
    outfile,
    max_mass2=3.0,
    subdir="runs_SNR-10",
    detectors="HLVK",
):
    """
    Extract and filter GW injection tables from a Zenodo-style ZIP archive.

    Parameters
    ----------
    zip_path : str or Path
        Path to the input ZIP archive.
    runs : list of str
        Observing run labels to extract (e.g., ["O5", "O6"]).
    outdir : str or Path
        Destination directory for output files.
    max_mass2 : float
        Maximum mass for secondary object (used to filter BNS/NSBH).
    subdir : str
        Name of the root folder inside the ZIP archive.
    detectors : str
        Detector label used in file path construction (e.g., HLVK).

    Outputs
    -------
    - ECSV summary table at <outdir>/observing-scenarios.ecsv
    - FITS localization files under <outdir>/<run>/ directories
    """
    out_root = pathlib.Path(outdir)
    for run in runs:
        (out_root / run).mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(zip_path) as archive:
        in_root = zipfile.Path(archive) / subdir

        tables = []
        for run in tqdm(runs, desc="Reading summary tables"):
            in_run = in_root / f"{run}{detectors}" / "farah"
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

        z = z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc).to_value()
        zp1 = 1 + z
        source_mass2 = table["mass2"] / zp1

        table = table[source_mass2 <= max_mass2]

        table.write(f"{out_root}/{outfile}", overwrite=True)

        for row in tqdm(table, desc="Copying FITS files"):
            filename = f"{row['coinc_event_id']}.fits"
            in_path = (
                in_root / f"{row['run']}{detectors}" / "farah" / "allsky" / filename
            )
            out_path = out_root / row["run"] / filename
            with in_path.open("rb") as in_file, out_path.open("wb") as out_file:
                copyfileobj(in_file, out_file)


if __name__ == "__main__":
    args = parse_arguments()
    process_zip(
        zip_path=args.zip,
        runs=args.runs,
        outdir=args.outdir,
        outfile=args.outfile,
        max_mass2=args.mass_threshold,
        subdir=args.subdir,
        detectors=args.detectors,
    )
