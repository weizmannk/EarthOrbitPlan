.. _postprocess:


Post-processing
===============

After generating observation plans, post-processing computes detection probabilities and extracts key optimization metrics for each event.

.. dropdown:: Extract detection Probability and Optimization metrics

    .. tab-set::

        .. tab-item:: Configuration File

            Use a configuration file (e.g., :doc:`params_ultrasat.ini <../../config/params_ultrasat.ini>`) to specify all parameters:

            .. code-block:: console

                python earthorbitplan/workflow/postprocess.py --config ./earthorbitplan/config/params_ultrasat.ini

        .. tab-item:: CLI

            Run the post-processing directly with a data directory:

            .. code-block:: console

                python postprocess.py --data-dir ./data --event-table ./data/observing-scenarios.ecsv  --output-file ./data/events.ecsv --sched-dir ./data/schedules

    .. important::

        The results are aggregated and saved in an ECSV table (by default, ``events.ecsv``), ready for statistical analysis or further reporting.



    .. admonition:: Main Processing (`earthorbitplan.workflow.postprocess`)
        :class: tip


        .. .. automodapi:: earthorbitplan.workflow.postprocess
        ..     :no-inheritance-diagram:


.. dropdown:: Kilonova detection rate and Statistics

     Calculating poisson-lognormal rate quantiles and formatting

    .. tab-set::

        .. tab-item:: python


            Once activated, you can edit and run any code block below in your browser.
            Powered by Thebe and Binder.

            .. jupyter-execute::
                :raises:

                import warnings

                import numpy as np
                from astropy import units as u
                from astropy.table import QTable
                from scipy import stats

                from earthorbitplan.probability.rate import poisson_lognormal_rate_quantiles
                from earthorbitplan.utils.path import get_project_root

                # ------------------------------------------------------------------------------
                # Suppress known warnings for cleaner output
                warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
                warnings.filterwarnings("ignore", ".*dubious year.*")
                warnings.filterwarnings(
                    "ignore", "Tried to get polar motions for times after IERS data is valid.*"
                )

                # ------------------------------------------------------------------------------
                # Load simulated event data
                root = get_project_root()
                events_file = root / "data" / "events.ecsv"
                main_table = QTable.read(events_file)

                # Get unique run names
                runs = np.unique(main_table["run"])

                # Filter events by objective_value cutoff
                cutoff = main_table["cutoff"][0]
                main_table = main_table[main_table["objective_value"] >= cutoff]

                # Group events by run
                event_tables_by_run = {run: main_table[main_table["run"] == run] for run in runs}

                # ------------------------------------------------------------------------------
                # Set merger rate priors from O3 R&P Table II (last column)
                lo, mid, hi = 100, 240, 510  # In Gpc^-3 yr^-1

                # Log-normal width for 90% interval
                (standard_90pct_interval,) = np.diff(stats.norm.interval(0.9))
                log_target_rate_mu = np.log(mid)
                log_target_rate_sigma = np.log(hi / lo) / standard_90pct_interval

                # Get effective rate for each run
                log_simulation_effective_rate_by_run = {
                    key: np.log(value.to_value(u.Gpc**-3 * u.yr**-1))
                    for key, value in main_table.meta["effective_rate"].items()
                }

                # ------------------------------------------------------------------------------
                # Compute median and quantiles for each run
                prob_quantiles = np.asarray([0.5, 0.05, 0.95])  # Median, 5%, 95%
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

                # Compute Poisson-Lognormal rate quantiles for all runs
                rate_quantiles = poisson_lognormal_rate_quantiles(
                    prob_quantiles[np.newaxis, np.newaxis, :],
                    mu.T[:, :, np.newaxis],
                    log_target_rate_sigma,
                )

                # ------------------------------------------------------------------------------
                # Utility: Format a table as reStructuredText grid table
                def make_rst_table(headers, rows):
                    columns = [headers] + rows
                    n_cols = len(headers)
                    col_widths = [max(len(str(row[i])) for row in columns) for i in range(n_cols)]

                    def sep(char="+", fill="-"):
                        return char + char.join(fill * (w + 2) for w in col_widths) + char

                    def fmt_row(row):
                        return (
                            "| "
                            + " | ".join(str(cell).ljust(w) for cell, w in zip(row, col_widths))
                            + " |"
                        )

                    lines = [
                        sep(),
                        fmt_row(headers),
                        sep("=", "="),
                    ]
                    for row in rows:
                        lines.append(fmt_row(row))
                        lines.append(sep())
                    return "\n".join(lines)

                # Example: Prepare headers and format quantile results
                headers = ["Run"] + list(runs)
                labels = ["Number of events selected", "Number of events detected"]
                rst_rows = []

                for label, row in zip(labels, rate_quantiles):
                    formatted = [
                        "${}_{{-{}}}^{{+{}}}$".format(*np.rint([mid, mid - lo, hi - mid]).astype(int))
                        for mid, lo, hi in row
                    ]
                    rst_rows.append([label] + formatted)

                rst_table = make_rst_table(headers, rst_rows)

                # Print the table for RST documentation
                print(rst_table)


        .. tab-item:: Notebook

            .. seealso::

                You can explore and run the calculations in the
                - :download:`Download this notebook <../../earthorbitplan/tutorials/kilonovae_detection_rate.ipynb>`
                - `Open in Binder <https://mybinder.org/v2/gh/weizmannk/EarthOrbitPlan/HEAD?urlpath=https%3A%2F%2Fgithub.com%2Fweizmannk%2FEarthOrbitPlan%2Ftree%2Fmain%2Fearthorbitplan%2Ftutorials%2Fkilonovae_detection_rate.ipynb>`_

                Alternatively, view it rendered in the documentation if included in the toctree.
