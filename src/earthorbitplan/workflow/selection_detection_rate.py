import numpy as np
from astropy import units as u
from astropy.table import QTable
from scipy import stats

from earthorbitplan.probability.rate import poisson_lognormal_rate_cdf
from earthorbitplan.utils.table import get_skygrid


def summarize_selected_detected_events(
    events_file,
    merger_rate_lo,
    merger_rate_mid,
    merger_rate_hi,
    quantiles=(0.5, 0.05, 0.95),
    run_duration=1.0,
    poisson_lognormal_rate_quantiles=None,
    aggregate_runs=None,
    aggregate_label="O5 total",
    include_zero_probability=False,
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
    merger_rate_lo, merger_rate_mid, merger_rate_hi : float
        5%, 50% and 95% quantiles of the target merger rate density, in
        Gpc^-3 yr^-1. Required, with no default: the rate depends on the
        population model the run was drawn from, so the caller states it
        rather than inheriting a silent default.
    run_duration : float, optional
        Duration of the observing run in years (default: 1.5).
    poisson_lognormal_rate_quantiles : callable
        Function to calculate Poisson-lognormal rate quantiles.
    aggregate_runs : sequence of str, optional
        Run names to combine into one extra block at the foot of the table,
        e.g. ``("O5a", "O5b", "O5c")``. The merger-rate prior is the same
        astrophysical rate for every run, so its uncertainty is fully
        correlated across them: the aggregate is built by summing the
        expected counts, ``lambda_total = sum_r exp(mu_r)``, and re-deriving
        the quantiles from ``log(lambda_total)`` at the *same* prior width.
        Convolving the per-run predictive distributions as if independent, or
        summing their quantiles, would both be wrong -- the median is not
        additive. ``None`` (default) adds no such block.
    aggregate_label : str, optional
        Row label for that block (default ``"O5 total"``).
    include_zero_probability : bool, optional
        If True, add a ``P(N=0)`` column per class: the probability of
        observing no event at all, from the Poisson-lognormal CDF at zero.
        It follows neither from lambda nor from the median, and it is what
        tells the reader the risk of an empty year. Default False.
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

    # Get unique run names
    runs = np.unique(main_table["run"])

    # Extract mission and skygrid from metadata
    mission = np.unique(main_table["mission"])[0]
    skygrid = get_skygrid(main_table)

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

    # Expected count, straight from mu and before the CDF inversion. The
    # table reports it because the median is uninformative here: it rounds to
    # zero for every Detected entry, including those at lambda = 1.07 and
    # 1.16, whose continuous medians are 0.41 and 0.50. Only lambda tells the
    # sub-periods apart.
    lam = np.exp(mu)

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

    # Probability of seeing nothing at all, P(N = 0), per run and class.
    zero_prob = (
        poisson_lognormal_rate_cdf(0, mu, log_target_rate_sigma)
        if include_zero_probability
        else None
    )

    # Aggregate block. The merger-rate prior is one astrophysical rate shared
    # by every run, so its uncertainty is fully correlated across them: the
    # expected counts add, and the quantiles are re-derived once from their
    # sum at the same prior width. Summing quantiles, or convolving the
    # per-run predictive distributions as if independent, would both be wrong.
    aggregate = None
    if aggregate_runs:
        run_index = {run: i for i, run in enumerate(runs)}
        missing = [run for run in aggregate_runs if run not in run_index]
        if missing:
            raise ValueError(
                f"aggregate_runs names absent from {events_file}: {missing}. "
                f"Available runs: {sorted(run_index)}."
            )
        columns = [run_index[run] for run in aggregate_runs]

        # Sum exp(mu) in full precision -- never the rounded printed values.
        lam_aggregate = lam[:, columns, :].sum(axis=1)
        mu_aggregate = np.log(lam_aggregate)
        aggregate_quantiles = poisson_lognormal_rate_quantiles(
            prob_quantiles[np.newaxis, np.newaxis, :],
            mu_aggregate[:, :, np.newaxis],
            log_target_rate_sigma,
        )
        aggregate = {
            "lam": lam_aggregate,
            "quantiles": aggregate_quantiles,
            "zero_prob": (
                poisson_lognormal_rate_cdf(0, mu_aggregate, log_target_rate_sigma)
                if include_zero_probability
                else None
            ),
        }

    # -------------------------------------------------------------------------
    # Build LaTeX table — article style:
    #   - \Xhline{3\arrayrulewidth} for thick top/bottom/mid rules
    #   - \multirow{2}{*} to merge run name across Selected/Detected rows
    #   - \textbf for column headers and run names
    #   - \hline between runs
    # Requires: \usepackage{multirow, makecell} in preamble
    # -------------------------------------------------------------------------
    # Same rule as the figure file names: skip the empty skygrid so a
    # mission without grid variants does not produce "tab:uvex--...".
    label_stem = "-".join(part for part in (str(mission), skygrid) if part)

    classes = list(is_class_by_category.keys())  # ["BNS", "NSBH", "All"]
    row_labels = ["Selected", "Detected"]
    xhline = r"\Xhline{3\arrayrulewidth}"
    hline = r"\hline"

    def interval_parts(mid, lo, hi):
        """Rounded (median, minus, plus), or None if every percentile is 0."""
        m, lo_val, h_val = np.rint([mid, mid - lo, hi - mid]).astype(int)
        if m == 0 and lo_val == 0 and h_val == 0:
            # "0^{+0}_{-0}" says nothing that the lambda column does not say
            # better, so the renderers below collapse it to "<1".
            return None
        return m, lo_val, h_val

    def fmt(mid, lo, hi):
        """LaTeX cell, AAS convention: superscript before subscript."""
        parts = interval_parts(mid, lo, hi)
        if parts is None:
            return r"$<1$"
        m, lo_val, h_val = parts
        return f"${m}^{{+{h_val}}}_{{-{lo_val}}}$"

    def fmt_plain(mid, lo, hi):
        """Same cell, for the plain-text table shown in a notebook."""
        parts = interval_parts(mid, lo, hi)
        if parts is None:
            return "<1"
        m, lo_val, h_val = parts
        return f"{m}^{{+{h_val}}}_{{-{lo_val}}}"

    # One sub-column per quantity, per class.
    sub_headers = [r"$\lambda$", r"90\% CI"]
    if include_zero_probability:
        sub_headers.append(r"$P(N{=}0)$")
    n_sub = len(sub_headers)

    col_spec = "ll" + "c" * n_sub * len(classes)
    header = "\n".join(
        [
            r"\textbf{Run} & & "
            + " & ".join(
                rf"\multicolumn{{{n_sub}}}{{c}}{{\textbf{{{cls}}}}}" for cls in classes
            )
            + r" \\",
            " & & " + " & ".join([" & ".join(sub_headers)] * len(classes)) + r" \\",
        ]
    )

    def cells_for(lam_row, quantile_row, zero_row, i_cls):
        """The sub-columns of one class, in header order."""
        out = [f"{lam_row[i_cls]:.2f}", fmt(*quantile_row[i_cls, :])]
        if include_zero_probability:
            out.append(rf"{100 * zero_row[i_cls]:.1f}\%")
        return out

    data_rows = []
    for i_run, run in enumerate(runs):
        for i_label, label in enumerate(row_labels):
            # \multirow on first sub-row only; blank on second
            run_cell = (
                rf"\multirow{{2}}{{*}}{{\textbf{{{run}}}}}" if i_label == 0 else ""
            )
            cells = []
            for i_cls in range(len(classes)):
                cells += cells_for(
                    lam[i_label, i_run],
                    rate_quantiles[i_label, i_run],
                    zero_prob[i_label, i_run] if include_zero_probability else None,
                    i_cls,
                )
            data_rows.append(f"{run_cell} & {label} & " + " & ".join(cells) + r" \\")
        # thin \hline between runs, nothing after last
        if i_run < len(runs) - 1:
            data_rows.append(hline)

    if aggregate is not None:
        data_rows.append(xhline)
        for i_label, label in enumerate(row_labels):
            run_cell = (
                rf"\multirow{{2}}{{*}}{{\textbf{{{aggregate_label}}}}}"
                if i_label == 0
                else ""
            )
            cells = []
            for i_cls in range(len(classes)):
                cells += cells_for(
                    aggregate["lam"][i_label],
                    aggregate["quantiles"][i_label],
                    (
                        aggregate["zero_prob"][i_label]
                        if include_zero_probability
                        else None
                    ),
                    i_cls,
                )
            data_rows.append(f"{run_cell} & {label} & " + " & ".join(cells) + r" \\")

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

    aggregate_caption = (
        (
            rf" The {aggregate_label} block combines "
            + ", ".join(aggregate_runs)
            + r". The merger-rate prior is a single astrophysical rate shared by"
            r" those runs, so its uncertainty is fully correlated across them:"
            r" the aggregate adds the expected counts, $\lambda_\mathrm{tot} ="
            r" \sum_r \lambda_r$, and its interval is re-derived from"
            r" $\ln \lambda_\mathrm{tot}$ at the same prior width. Medians are not"
            r" additive and must not be summed."
        )
        if aggregate is not None
        else ""
    )
    zero_caption = (
        r" $P(N{=}0)$ is the probability of observing no event at all, from the"
        r" Poisson-lognormal CDF at zero; it follows neither from $\lambda$ nor"
        r" from the median."
        if include_zero_probability
        else ""
    )

    latex_table = "\n".join(
        [
            r"\begin{table}",
            r"\renewcommand\arraystretch{1.3}",
            r"\setlength{\tabcolsep}{0.3cm}",
            r"\centering",
            rf"\caption{{Expected number of selected and detected events for {mission.upper()}"
            rf"{skygrid_str} per observing run and source class.{aggregate_caption}{zero_caption}"
            r" BNS: both components $\leq 3\,M_\odot$;"
            r" NSBH: one component $> 3\,M_\odot$;"
            r" All: BNS $+$ NSBH combined."
            r" $\lambda = \mathcal{R}_{50}\,T\,N / \mathcal{R}_{\mathrm{sim}}$ is the expected count at the fiducial (median) merger rate. The interval"
            r" is the median and 90\% credible range of the Poisson-lognormal predictive distribution, which marginalises over the rate prior; its"
            r" mean is $\lambda e^{\sigma^2/2}$, not $\lambda$. Entries marked $<1$ have every percentile below 0.5, where only $\lambda$ discriminates.}",
            rf"\label{{tab:{label_stem}-selected-detected-{run_duration}yr}}",
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

        headers = ["Run", ""]
        for cls in classes:
            headers += [f"{cls} lambda", f"{cls} 90% CI"]
            if include_zero_probability:
                headers.append(f"{cls} P(N=0)")

        def plain_cells_for(lam_row, quantile_row, zero_row, i_cls):
            out = [f"{lam_row[i_cls]:.2f}", fmt_plain(*quantile_row[i_cls, :])]
            if include_zero_probability:
                out.append(f"{100 * zero_row[i_cls]:.1f}%")
            return out

        rst_rows = []
        for i_run, run in enumerate(runs):
            for i_label, label in enumerate(row_labels):
                run_cell = run if i_label == 0 else ""
                cells = []
                for i_cls in range(len(classes)):
                    cells += plain_cells_for(
                        lam[i_label, i_run],
                        rate_quantiles[i_label, i_run],
                        zero_prob[i_label, i_run] if include_zero_probability else None,
                        i_cls,
                    )
                rst_rows.append([run_cell, label] + cells)

        if aggregate is not None:
            for i_label, label in enumerate(row_labels):
                run_cell = aggregate_label if i_label == 0 else ""
                cells = []
                for i_cls in range(len(classes)):
                    cells += plain_cells_for(
                        aggregate["lam"][i_label],
                        aggregate["quantiles"][i_label],
                        (
                            aggregate["zero_prob"][i_label]
                            if include_zero_probability
                            else None
                        ),
                        i_cls,
                    )
                rst_rows.append([run_cell, label] + cells)

        print(make_rst_table(headers, rst_rows))

    if output_file is not None:
        from pathlib import Path

        # Trailing newline: without it the end-of-file-fixer pre-commit hook
        # rewrites the file on every regeneration, showing a spurious diff.
        Path(output_file).write_text(latex_table + "\n")

    return latex_table
