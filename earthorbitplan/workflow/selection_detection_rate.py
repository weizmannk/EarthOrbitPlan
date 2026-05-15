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
    Summarize selected and detected events by class an run, returning a LaTeX table of rate quantiles.

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
        LaTeX-formatted table as a string, with one row for the number of 
        events selected and one row for the expected number of events detected,
        broken down by source class (BNS, NSBH, All) and observing run,
        along with their 90% credible intervals.

    Notes
    -----
    Source classes: BNS (both components < 3 Sun-Mass), NSBH (one component > 3 Sun-Mass), 
    All (no mass cut). The merger rate quantiles are standardized using the O3 
    Rate & Population (R&P) reference values ...

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
    Run & O5a (BNS) & O5a (NSBH) & O5a (All) & O5b (BNS) & O5b (NSBH) & O5b (All) \\
    Number of events selected & $43_{-26}^{+56}$ & $12_{-7}^{+15}$ & $55_{-33}^{+72}$ & ... \\
    Number of events detected & $19_{-12}^{+26}$ & $8_{-5}^{+11}$ & $27_{-16}^{+34}$ & ...
    
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

    # Derive source class from component masses before filtering
    main_table["source_class"] = np.where(
        (main_table["mass1"] <= 3) & (main_table["mass2"] <= 3),
        "BNS",
        np.where(
            (main_table["mass1"] > 3) | (main_table["mass2"] > 3),
            "NSBH",
            "Other",
        ),
    )

    is_class_by_category = {
        "BNS": lambda table: table["source_class"] == "BNS",
        "NSBH": lambda table: table["source_class"] == "NSBH",
        "All": lambda table: np.ones(len(table), dtype=bool),
    }

    event_tables_by_run = {run: main_table[main_table["run"] == run] for run in runs}

    event_tables_by_run_and_class = {
        run: {
            cls: table[is_class(table)]
            for cls, is_class in is_class_by_category.items()
        }
        for run, table in event_tables_by_run.items()
    }

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

    # Prepare mu for each run and class
    mu = []
    for run in runs:
        mu_run = []
        for table_class in event_tables_by_run_and_class[run].values():
            # for cls, table_class in event_tables_by_run_and_class[run].items():
            n_selected = len(table_class)
            n_detected = np.sum(table_class["detection_probability_known_position"])
            mu_run.append(
                log_target_rate_mu
                + np.log(run_duration)
                - log_sim_effective_rate_by_run[run]
                + np.log([n_selected, n_detected])
            )
        mu.append(mu_run)

    mu = np.moveaxis(np.array(mu), 2, 0)

    # Compute quantiles => shape (2, n_runs, n_classes, 3)
    # by_run a shape (n_runs, n_classes, 3)
    prob_quantiles = np.array(quantiles)
    rate_quantiles = poisson_lognormal_rate_quantiles(
        prob_quantiles[np.newaxis, np.newaxis, :],
        mu[:, :, :, np.newaxis],
        log_target_rate_sigma,
    )

    # Build LaTeX table
    classes = list(is_class_by_category.keys())
    n_classes = len(classes)
    col_spec = "l" + "c" * (len(runs) * n_classes)

    # Row 1: run names spanning n_classes columns each
    run_headers = " & ".join(
        rf"\multicolumn{{{n_classes}}}{{c}}{{{run}}}" for run in runs
    )

    # Row 2: class names under each run
    class_headers = " & ".join(cls for _ in runs for cls in classes)

    # Data rows
    labels = ["Number of events selected", "Number of events detected"]
    data_rows = []
    for i, (label, by_run) in enumerate(zip(labels, rate_quantiles)):
        formatted = [
            "${}_{{-{}}}^{{+{}}}$".format(
                *np.rint([mid, mid - lo, hi - mid]).astype(int)
            )
            for run_quantiles in by_run
            for mid, lo, hi in run_quantiles
        ]
        data_rows.append(" & ".join([label] + formatted) + r" \\")

    latex_table = "\n".join(
        [
            rf"\begin{{tabular}}{{{col_spec}}}",
            r"\hline",
            rf"& {run_headers} \\",
            rf"& {class_headers} \\",
            r"\hline",
            *data_rows,
            r"\hline",
            r"\end{tabular}",
        ]
    )

    if verbose:
        print(latex_table)

    return latex_table
