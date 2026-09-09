import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import QTable

from earthorbitplan.utils.path import get_project_root


def format_row_latex(row):
    """Return a LaTeX table row for a given event."""

    # Format each field as required
    mass1 = np.format_float_positional(row["mass1"], 3, fractional=True)
    mass2 = np.format_float_positional(row["mass2"], 3, fractional=True)
    lon = np.format_float_positional(np.rad2deg(row["longitude"]), 4, fractional=True)
    lat = np.format_float_positional(
        np.rad2deg(row["latitude"]), 4, fractional=True, sign=True
    )
    distance = np.format_float_positional(row["distance"], 0, trim="-", fractional=True)
    area90 = np.format_float_positional(row["area(90)"], 0, trim="-", fractional=True)

    if row["objective_value"] >= 0.1:
        objective = "\\phantom{$<$}" + np.format_float_positional(
            row["objective_value"], 2, min_digits=2, fractional=True, trim="k"
        )
        detection = np.format_float_positional(
            row["detection_probability_known_position"],
            2,
            min_digits=2,
            fractional=True,
            trim="k",
        )
    else:
        objective = r"$<$0.10"
        detection = "---"
    # Construct the LaTeX table row
    return (
        f"{row['run']} & {row['coinc_event_id']} & {mass1} & {mass2} & "
        f"{lon} & {lat} & {distance} & {area90} & {objective} & {detection} \\\\"
    )


# --- Main processing ---
root = get_project_root()
events_file = root / "data" / "events.ecsv"
table = QTable.read(events_file)

zp1 = 1 + z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc)
table["mass1"] /= zp1
table["mass2"] /= zp1
table = table[table["mass2"] <= 3]

n_show = 4  # Number of first/last events to show
rows = [format_row_latex(row) for row in table]

# Write full table to LaTeX file
with open(root / "data" / "events.tex", "w") as f:
    for row in rows:
        print(row, file=f)

# Print first and last 4 rows
print("% The first 4 events:")
for r in rows[:n_show]:
    print(r)
print("\n% The last 4 events:")
for r in rows[-n_show:]:
    print(r)
