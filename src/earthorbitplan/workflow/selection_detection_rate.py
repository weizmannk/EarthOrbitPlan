import numpy as np
from astropy import units as u
from astropy.table import QTable
from scipy import stats


def summarize_selected_detected_events(
    events_file,
    quantiles=(0.5, 0.05, 0.95),
    merger_rate_lo=50,
    merger_rate_mid=130,
    merger_rate_hi=290,
    run_duration=1.0,
    poisson_lognormal_rate_quantiles=None,
    output_file=None,
    verbose=True,
):
    """
    Summarize selected and detected events by class and run, returning a LaTeX table of rate quantiles.

    Computes the 5%, 50%, and 95% quantiles of the merger rate for each run using a Poisson-lognormal model,
    following the O4 Rate & Population (R&P) methodology
    (see `GWTC-4 R&P <https://arxiv.org/pdf/2508.18083>`_).

    The log-expected rate for each run $r$ and source class $c$ is:
    .. math::

        \mu_{r,c} = \log(\mathrm{target\ median}) + \log(\mathrm{run\ duration})
                  - \log(\mathrm{effective\ rate}_{\mathrm{sim},r}) + \log(N_{r,c})

    where :math:`N_{r,c}` is the number of selected or detected events for class $c$ in run $r$.
    Note the effective simulation rate depends only on the run, not the class.

    Parameters
    ----------
    events_file : str or Path
        Path to the ECSV table of candidate events.
    quantiles : sequence of float, optional
        Probability quantiles to compute (default: (0.5, 0.05, 0.95)).
    merger_rate_lo : float, optional
        Lower bound of target merger rate (Gpc^-3 yr^-1).
    merger_rate_mid : float, optional
        Median of target merger rate (Gpc^-3 yr^-1).
    merger_rate_hi : float, optional
        Upper bound of target merger rate (Gpc^-3 yr^-1).
    run_duration : float, optional
        Duration of the observing run in years (default: 1.5).
    poisson_lognormal_rate_quantiles : callable
        Function to calculate Poisson-lognormal rate quantiles.
    verbose : bool, optional
        If True, displays the tabular in Jupyter via IPython.display.Latex.

    Returns
    -------
    str
        Full LaTeX table (\\begin{table}...\\end{table}) as a string.
        Rows = observing runs x {Selected, Detected}. Columns = BNS, NSBH, All.
        Style: \\Xhline thick rules, \\multirow for run names, \\textbf headers.
        Requires: multirow, makecell packages in LaTeX preamble.

    Notes
    -----
    Source classes: BNS (both components < 3 Msun), NSBH (one component > 3 Msun),
    All (no mass cut). Column order is BNS => NSBH => All, from most specific to most inclusive.

    References
    ----------
    .. footbibliography::
    """

    # Load main event table
    main_table = QTable.read(events_file)

    # Get unique run names, excluding O5a-HL
    runs = np.unique(main_table["run"])
    runs = runs[runs != "O5a-HL"]

    # Extract mission and skygrid from metadata
    mission = np.unique(main_table["mission"])[0]
    skygrid_raw = np.unique(main_table["skygrid"])[0]
    skygrid = "" if skygrid_raw == "None" else skygrid_raw

    # Map skygrid for ULTRASAT to be add in the the caption
    skygrid_label = {
        "allsky": r"\emph{All-Sky Survey} (AllSS)",
        "non-overlap": r"\emph{Low-Cadence Survey} (LCS)",
    }.get(skygrid, "")

    # Get cutoff value used for the simulation
    cutoff_values = main_table["cutoff"]
    if not np.all(cutoff_values == cutoff_values[0]):
        raise ValueError("Multiple cutoff values found in the table.")
    cutoff = cutoff_values[0]

    # Apply cutoff filter
    main_table = main_table[main_table["objective_value"] >= cutoff]

    if len(main_table) == 0:
        raise RuntimeError("No events passed the cutoff filter.")

    # Class order: BNS => NSBH => All (specific to inclusive)
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

    # 90% CI width for standard normal
    (standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))

    # Lognormal prior parameters
    log_target_rate_mu = np.log(merger_rate_mid)
    log_target_rate_sigma = (
        np.log(merger_rate_hi / merger_rate_lo) / standard_90pct_interval
    )

    # Log effective simulation rate per run (from metadata)
    log_sim_effective_rate_by_run = {
        key: np.log(value.to_value(u.Gpc**-3 * u.yr**-1))
        for key, value in main_table.meta["effective_rate"].items()
    }

    # Compute mu_{r,c} for each run and class
    mu = []
    for run in runs:
        mu_run = []
        for table_class in event_tables_by_run_and_class[run].values():
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

    # Compute Poisson-lognormal quantiles
    # rate_quantiles shape: (2, n_runs, n_classes, 3)
    #   axis 0: [selected, detected]
    #   axis 1: runs
    #   axis 2: classes (BNS, NSBH, All)
    #   axis 3: quantiles (median, lo, hi)
    prob_quantiles = np.array(quantiles)
    rate_quantiles = poisson_lognormal_rate_quantiles(
        prob_quantiles[np.newaxis, np.newaxis, :],
        mu[:, :, :, np.newaxis],
        log_target_rate_sigma,
    )

    # -------------------------------------------------------------------------
    # Build LaTeX table — article style:
    #   - \Xhline{3\arrayrulewidth} for thick top/bottom/mid rules
    #   - \multirow{2}{*} to merge run name across Selected/Detected rows
    #   - \textbf for column headers and run names
    #   - \hline between runs
    # Requires: \usepackage{multirow, makecell} in preamble
    # -------------------------------------------------------------------------
    classes = list(is_class_by_category.keys())  # ["BNS", "NSBH", "All"]
    row_labels = ["Selected", "Detected"]
    xhline = r"\Xhline{3\arrayrulewidth}"
    hline = r"\hline"

    def fmt(mid, lo, hi):
        m, lo_val, h_val = np.rint([mid, mid - lo, hi - mid]).astype(int)
        return f"${m}_{{-{lo_val}}}^{{+{h_val}}}$"

    col_spec = "llccc"
    header = r"\textbf{Run} & & \textbf{BNS} & \textbf{NSBH} & \textbf{All} \\"

    data_rows = []
    for i_run, run in enumerate(runs):
        for i_label, label in enumerate(row_labels):
            # \multirow on first sub-row only; blank on second
            run_cell = (
                rf"\multirow{{2}}{{*}}{{\textbf{{{run}}}}}" if i_label == 0 else ""
            )
            cells = [
                fmt(*rate_quantiles[i_label, i_run, i_cls, :])
                for i_cls in range(len(classes))
            ]
            data_rows.append(f"{run_cell} & {label} & " + " & ".join(cells) + r" \\")
        # thin \hline between runs, nothing after last
        if i_run < len(runs) - 1:
            data_rows.append(hline)

    # tabular only — KaTeX-compatible for Jupyter display
    tabular = "\n".join(
        [
            rf"\begin{{tabular}}{{{col_spec}}}",
            xhline,
            header,
            xhline,
            *data_rows,
            xhline,
            r"\end{tabular}",
        ]
    )
    # Caption with mission and skygrid
    skygrid_str = rf" using the {skygrid_label} strategy," if skygrid_label else ","

    latex_table = "\n".join(
        [
            r"\begin{table}",
            r"\renewcommand\arraystretch{1.3}",
            r"\setlength{\tabcolsep}{0.3cm}",
            r"\centering",
            rf"\caption{{Expected number of selected and detected events for {mission.upper()}"
            rf"{skygrid_str} per observing run and source class."
            r" BNS: both components $\leq 3\,M_\odot$;"
            r" NSBH: one component $> 3\,M_\odot$;"
            r" All: BNS $+$ NSBH combined."
            r" Values are medians with 90\% credible intervals.}",
            rf"\label{{tab:{mission}-{skygrid}-selected-detected-{run_duration}yr}}",
            tabular,
            r"\end{table}",
        ]
    )

    if verbose:

        def make_rst_table(headers, rows):
            columns = [headers] + rows
            n_cols = len(headers)
            col_widths = [
                max(len(str(row[i])) for row in columns) for i in range(n_cols)
            ]

            def sep(char="+", fill="-"):
                return char + char.join(fill * (w + 2) for w in col_widths) + char

            def fmt_row(row):
                return (
                    "| "
                    + " | ".join(str(cell).ljust(w) for cell, w in zip(row, col_widths))
                    + " |"
                )

            lines = [sep(), fmt_row(headers), sep("=", "=")]
            for row in rows:
                lines.append(fmt_row(row))
                lines.append(sep())
            return "\n".join(lines)

        headers = ["Run", "", "BNS", "NSBH", "All"]
        rst_rows = []
        for i_run, run in enumerate(runs):
            for i_label, label in enumerate(row_labels):
                run_cell = run if i_label == 0 else ""
                cells = [
                    "{}_{{-{}}}^{{+{}}}".format(
                        *np.rint([mid, mid - lo, hi - mid]).astype(int)
                    )
                    for mid, lo, hi in [
                        rate_quantiles[i_label, i_run, i_cls, :]
                        for i_cls in range(len(classes))
                    ]
                ]
                rst_rows.append([run_cell, label] + cells)

        print(make_rst_table(headers, rst_rows))

    if output_file is not None:
        from pathlib import Path

        Path(output_file).write_text(latex_table)

    return latex_table
