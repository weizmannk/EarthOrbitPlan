from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import QTable
from scipy import stats

from ..probability.rate import poisson_lognormal_rate_quantiles

base_path = Path("runs_SNR-10")
runs = ["O5", "O6"]

main_table = QTable.read("../events.ecsv")

runs = main_table.unique("run")
# Throw away events with objective value less than cutoff
plan_args = QTable.read(
    f"data/{main_table[0]['run']}/{main_table[0]['coinc_event_id']}.ecsv"
).meta["args"]
cutoff = plan_args["cutoff"]
main_table = main_table[main_table["objective_value"] >= cutoff]

event_tables_by_run = {run: main_table[main_table["run"] == run] for run in runs}


# O3 R&P paper Table II row 1 last column:
# 5%, 50%, and 95% quantiles of the total merger rate
# in Gpc^-3 yr^-1.
# See https://doi.org/10.1103/PhysRevX.13.011048
lo = 100
mid = 240
hi = 510

(standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))
log_target_rate_mu = np.log(mid)
log_target_rate_sigma = np.log(hi / lo) / standard_90pct_interval
log_target_rate_mu, log_target_rate_sigma

log_simulation_effective_rate_by_run = {
    key: np.log(value.to_value(u.Gpc**-3 * u.yr**-1))
    for key, value in main_table.meta["effective_rate"].items()
}

prob_quantiles = np.asarray([0.5, 0.05, 0.95])
run_duration = 1.5  # years
mu = np.asarray(
    [
        log_target_rate_mu
        + np.log(run_duration)
        - log_simulation_effective_rate_by_run[run]
        + np.log(
            [
                np.sum(_)
                for _ in [
                    np.ones_like(event_tables_by_run[run]["objective_value"]),
                    event_tables_by_run[run]["detection_probability_known_position"],
                ]
            ]
        )
        for run in runs
    ]
)

rate_quantiles = poisson_lognormal_rate_quantiles(
    prob_quantiles[np.newaxis, np.newaxis, :],
    mu.T[:, :, np.newaxis],
    log_target_rate_sigma,
)

with open("tables/selected-detected.tex", "w") as f:
    for i, (label, row) in enumerate(
        zip(["Number of events selected", "Number of events detected"], rate_quantiles)
    ):
        print(
            label,
            *(
                "${}_{{-{}}}^{{+{}}}$".format(
                    *np.rint([mid, mid - lo, hi - mid]).astype(int)
                )
                for mid, lo, hi in row
            ),
            sep=" & ",
            end=" \\\\\n" if i < len(runs) - 1 else "\n",
            file=f,
        )
