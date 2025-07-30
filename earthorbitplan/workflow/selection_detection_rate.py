import numpy as np
from astropy import units as u
from astropy.table import QTable
from scipy import stats


def summarize_selected_detected_events(
    events_file,
    quantiles=(0.5, 0.05, 0.95),
    merger_rate_lo=100,
    merger_rate_mid=240,
    merger_rate_hi=510,
    run_duration=1.5,
    poisson_lognormal_rate_quantiles=None,
    verbose=True,
):
    """
    Summarize selected and detected events by run, returning a LaTeX table of rate quantiles.

    Computes the 5%, 50%, and 95% quantiles of the merger rate for each run using a Poisson-lognormal model,
    following the O3 Rate & Population (R&P) methodology
    (see `Phys. Rev. X 13, 011048 <https://doi.org/10.1103/PhysRevX.13.011048>`_).

    The log-expected rate for each run is:
    .. math::

        \mu = \log(\mathrm{target median}) + \log(\mathrm{run duration}) - \log(\mathrm{simulation effective rate}) + \log(N)

    where :math:`N` is the number of selected or detected events.

    Parameters
    ----------
    events_file : str or Path
        Path to the ECSV table of candidate events.
    quantiles : sequence of float, optional
        Probability quantiles to compute (default: (0.5, 0.05, 0.95)).
    merger_rate_lo : float, optional
        Lower bound of target merger rate (Gpc\\ :sup:`-3`\\ yr\\ :sup:`-1`).
    merger_rate_mid : float, optional
        Median of target merger rate (Gpc\\ :sup:`-3`\\ yr\\ :sup:`-1`).
    merger_rate_hi : float, optional
        Upper bound of target merger rate (Gpc\\ :sup:`-3`\\ yr\\ :sup:`-1`).
    run_duration : float, optional
        Duration of the observing run in years (default: 1.5).
    poisson_lognormal_rate_quantiles : callable
        Function to calculate Poisson-lognormal rate quantiles.
    verbose : bool, optional
        If True, prints the LaTeX table.

    Returns
    -------
    str
        LaTeX-formatted table as a string, with one row for the number of events selected and one row for the expected number of events detected, 
        along with their 90% credible intervals for each run. This output can be directly included in documentation or reports.


    Notes
    -----
    The merger rate quantiles are standardized using the O3 Rate & Population (R&P) reference values
    (Table II, row 1, last column :footcite:`2023PhRvX..13a1048A`), but can be changed for other runs or models.

    Example
    -------
    >>> from earthorbitplan.workflow.selection_detection_rate import summarize_selected_detected_events
    >>> from earthorbitplan.probability.rate import poisson_lognormal_rate_quantiles
    >>> from earthorbitplan.utils.path import get_project_root
    >>> root = get_project_root()
    >>> events_file = root / "data" / "events.ecsv"
    >>> latex_table = summarize_selected_detected_events(
    ...     events_file, poisson_lognormal_rate_quantiles=poisson_lognormal_rate_quantiles
    ... )
    >>> print(latex_table)
    Run & O5 & O6 \\
    Number of events selected & $43_{-26}^{+56}$ & $55_{-33}^{+72}$ \\
    Number of events detected & $19_{-12}^{+26}$ & $25_{-16}^{+34}$

    References
    ----------
    .. footbibliography::
    """

    # Load main event table
    main_table = QTable.read(events_file)
    runs = np.unique(main_table["run"])

    # Get cutoff value use for the simulation
    cutoff_values = main_table["cutoff"]
    if not np.all(cutoff_values == cutoff_values[0]):
        raise ValueError("Multiple cutoff values found in the table.")
    cutoff = cutoff_values[0]

    # Apply cutoff filter
    main_table = main_table[main_table["objective_value"] >= cutoff]

    # Check that at least one event remains
    if len(main_table) == 0:
        raise RuntimeError("No events passed the cutoff filter.")

    event_tables_by_run = {run: main_table[main_table["run"] == run] for run in runs}

    # 90% confidence interval width for standard normal (used to scale log-normal sigma)
    (standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))

    # Log-mean (mu) and log-standard deviation (sigma) for the target merger rate distribution
    log_target_rate_mu = np.log(merger_rate_mid)
    log_target_rate_sigma = (
        np.log(merger_rate_hi / merger_rate_lo) / standard_90pct_interval
    )

    # Log-effective of observing scenarios simulation rates per run, from metadata
    log_sim_effective_rate_by_run = {
        key: np.log(value.to_value(u.Gpc**-3 * u.yr**-1))
        for key, value in main_table.meta["effective_rate"].items()
    }

    # Prepare mu for each run
    prob_quantiles = np.array(quantiles)
    mu = []
    for run in runs:
        obj_vals = event_tables_by_run[run]["objective_value"]
        det_probs = event_tables_by_run[run]["detection_probability_known_position"]
        n_selected = len(obj_vals)
        n_detected = np.sum(det_probs)
        log_effective_rate = log_sim_effective_rate_by_run[run]
        mu_run = (
            log_target_rate_mu
            + np.log(run_duration)
            - log_effective_rate
            + np.log([n_selected, n_detected])
        )
        mu.append(mu_run)
    mu = np.array(mu).T  # shape (2, n_runs)

    # Compute quantiles
    rate_quantiles = poisson_lognormal_rate_quantiles(
        prob_quantiles[np.newaxis, :],  # shape (1, 3)
        mu[:, :, np.newaxis],  # shape (2, n_runs, 1)
        log_target_rate_sigma,
    )  # shape (2, n_runs, 3)

    latex_rows = []
    header = "Run & " + " & ".join(list(runs)) + r" \\"
    latex_rows.append(header)

    labels = ["Number of events selected", "Number of events detected"]

    for i, (label, row) in enumerate(zip(labels, rate_quantiles)):
        formatted = [
            "${}_{{-{}}}^{{+{}}}$".format(
                *np.rint([mid, mid - lo, hi - mid]).astype(int)
            )
            for mid, lo, hi in row
        ]

        if i < len(labels) - 1:
            line = " & ".join([label] + formatted) + r" \\"
        else:
            line = " & ".join([label] + formatted)
        latex_rows.append(line)

    latex_table = "\n".join(latex_rows)

    # Display table (print LaTeX-ready lines for any number of runs)
    if verbose:
        print(latex_table)

    return latex_table
