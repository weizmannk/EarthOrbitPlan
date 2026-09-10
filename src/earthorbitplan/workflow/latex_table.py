"""Render the events summary table as LaTeX rows (``events.tex``)."""

import argparse

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable


def build_table(events_file):
    table = QTable.read(events_file)
    zp1 = 1 + z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc)
    table["mass1"] /= zp1
    table["mass2"] /= zp1
    return table[table["mass2"] <= 3]


def write_tex(table, out_file, n=8):
    with open(out_file, "w") as f:
        for row in table[:n]:
            print(
                row["run"],
                row["coinc_event_id"],
                np.format_float_positional(row["mass1"], 3, fractional=True),
                np.format_float_positional(row["mass2"], 3, fractional=True),
                np.format_float_positional(
                    np.rad2deg(row["longitude"]), 4, fractional=True
                ),
                np.format_float_positional(
                    np.rad2deg(row["latitude"]), 4, fractional=True, sign=True
                ),
                np.format_float_positional(
                    row["distance"], 0, trim="-", fractional=True
                ),
                np.format_float_positional(
                    row["area(90)"], 0, trim="-", fractional=True
                ),
                (
                    r"\phantom{$<$}"
                    + np.format_float_positional(
                        row["objective_value"],
                        2,
                        min_digits=2,
                        fractional=True,
                        trim="k",
                    )
                    if row["objective_value"] >= 0.1
                    else r"$<$0.10"
                ),
                (
                    np.format_float_positional(
                        row["detection_probability_known_position"],
                        2,
                        min_digits=2,
                        fractional=True,
                        trim="k",
                    )
                    if row["objective_value"] >= 0.1
                    else "---"
                ),
                sep=" & ",
                end=" \\\\\n",
                file=f,
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", default="events.ecsv", help="Input events ECSV")
    parser.add_argument("--out", default="events.tex", help="Output LaTeX fragment")
    parser.add_argument("-n", type=int, default=8, help="Number of rows")
    args = parser.parse_args()
    write_tex(build_table(args.events), args.out, args.n)


if __name__ == "__main__":
    main()
