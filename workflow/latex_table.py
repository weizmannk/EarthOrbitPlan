import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable

table = QTable.read("../events.ecsv")
zp1 = 1 + z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc)
table["mass1"] /= zp1
table["mass2"] /= zp1
table = table[table["mass2"] <= 3]


def maybe_zero(s: str):
    if s == "0.":
        s = "0.00"
    return s


n = 8
with open("tables/events.tex", "w") as f:
    for i, row in enumerate(table[:n]):
        print(
            "O5",
            row["coinc_event_id"],
            np.format_float_positional(row["mass1"], 3, fractional=True),
            np.format_float_positional(row["mass2"], 3, fractional=True),
            np.format_float_positional(
                np.rad2deg(row["longitude"]), 4, fractional=True
            ),
            np.format_float_positional(
                np.rad2deg(row["latitude"]), 4, fractional=True, sign=True
            ),
            np.format_float_positional(row["distance"], 0, trim="-", fractional=True),
            np.format_float_positional(row["area(90)"], 0, trim="-", fractional=True),
            "\phantom{$<$}"
            + np.format_float_positional(
                row["objective_value"], 2, min_digits=2, fractional=True, trim="k"
            )
            if row["objective_value"] >= 0.1
            else r"$<$0.10",
            np.format_float_positional(
                row["detection_probability_known_position"],
                2,
                min_digits=2,
                fractional=True,
                trim="k",
            )
            if row["objective_value"] >= 0.1
            else "---",
            sep=" & ",
            end=" \\\\\n",
            file=f,
        )
