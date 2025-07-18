#!/usr/bin/env python
"""
postproces.py: Compute Detection Probabilities and Extract Metrics
==============================================================================

This script processes observation plans generated by M4OPT, computing detection
probabilities and extracting relevant optimization metrics.

It supports both command-line arguments and `.ini` configuration files.

Usage
-----
From the command line:

    python postprocess.py  --data-dir data

Or using a config file:

    python earthorbitplan/workflow/postprocess.py --config earthorbitplan/config/params_ultrasat.ini


Description
-----------
This process computes the detection probability using detection_probability.py,
collects all relevant data from the simulation, and stores the results in an
`events.ecsv` file located in the data directory (default: ./data). You can
then use this file later for further statistical processing.
"""

import argparse
import configparser
import logging
import sys
import warnings
from pathlib import Path

from astropy.table import QTable
from ligo.skymap.util.progress import progress_map

from ..probability.detection_probability import get_detection_probability_known_position

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
warnings.filterwarnings("ignore", ".*dubious year.*")
warnings.filterwarnings(
    "ignore", "Tried to get polar motions for times after IERS data is valid.*"
)


def parse_arguments():
    """
    Parse command-line arguments or load them from a .ini configuration file.
    """
    parser = argparse.ArgumentParser(description="Post-process M4OPT ECSV plans.")
    parser.add_argument("--config", type=str, help="Path to .ini configuration file")

    args, remaining_args = parser.parse_known_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)
        cfg = config["params"]
        return argparse.Namespace(
            data_dir=cfg.get("data_dir", fallback="data"),
            event_table=cfg.get("event_table", fallback="observing-scenarios.ecsv"),
            output_file=cfg.get("output_file", fallback="events.ecsv"),
            sched_dir=cfg.get("sched_dir", fallback="schedules"),
        )

    # CLI parsing if no config
    parser.add_argument(
        "--data-dir", type=str, default="data", help="Directory containing ECSV files"
    )
    parser.add_argument(
        "--event-table",
        type=str,
        default="observing-scenarios.ecsv",
        help="Input summary table",
    )
    parser.add_argument(
        "--output-file", type=str, default="events.ecsv", help="Output filename"
    )
    parser.add_argument(
        "--sched-dir", type=str, default="schedules", help="Schedule directory"
    )

    return parser.parse_args(remaining_args)


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
        force=True,
    )


def process(row, sched_path):
    """
    Process a single event row and its associated schedule.

    Parameters
    ----------
    row : astropy.table.Row
        A row from the input event table, containing event and run information.
    sched_path : pathlib.Path
        Base path to the directory containing schedule files.

    Returns
    -------
    tuple
        A tuple containing:
            - detection probability (float or None)
            - objective value (float or None)
            - best bound (float or None)
            - solution status (str or None)
            - solution time (float or None)
            - number of observed fields (int or None)
            - cutoff value (float or None)
        Returns all None if the schedule file is missing.
    """

    run = row["run"]
    event_id = row["coinc_event_id"]
    plan_file = sched_path / run / f"{event_id}.ecsv"

    if not plan_file.exists():
        logging.warning(f"Missing schedule file: {plan_file}")
        return (None, None, None, None, None, None)

    plan = QTable.read(plan_file)
    plan_args = {**plan.meta["args"]}
    plan_args.pop("skymap", None)

    return (
        get_detection_probability_known_position(plan, row, plan_args),
        plan.meta.get("objective_value"),
        plan.meta.get("best_bound"),
        plan.meta.get("solution_status"),
        plan.meta.get("solution_time"),
        len(plan[plan["action"] == "observe"]) // plan_args["visits"],
        plan_args.get("cutoff"),
    )


def main():
    args = parse_arguments()
    setup_logging()

    base_path = Path(args.data_dir)
    input_path = base_path / args.event_table
    sched_path = base_path / args.sched_dir
    output_path = base_path / args.output_file

    logging.info(f"Reading input table from {input_path}")
    table = QTable.read(input_path)

    logging.info("Processing observation plans...")
    (
        table["detection_probability_known_position"],
        table["objective_value"],
        table["best_bound"],
        table["solution_status"],
        table["solution_time"],
        table["num_fields"],
        table["cutoff"],
    ) = zip(*progress_map(lambda row: process(row, sched_path), table))

    logging.info(f"Writing results to {output_path}")
    table.write(output_path, overwrite=True)


if __name__ == "__main__":
    main()
