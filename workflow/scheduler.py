#!/usr/bin/env python
"""
m4opt-scheduler: Batch Scheduling for Gravitational Wave Follow-up
==================================================================

This script executes batch scheduling using the M4OPT framework for a list of gravitational wave sky maps.
It supports various execution backends (HTCondor, local parallel with joblib, or Dask) and can be configured
via command-line arguments or a configuration `.ini` file.

Usage
-----
You can run the script either with Command-Line Interface arguments:

    python scheduler.py --mission ULTRASAT --bandpass NUV ...

Or with a configuration file:

    python workscheduler.py --config params.ini

In the root directory of the project:
    python -m workflow.scheduler --config params.ini

"""

import argparse
import configparser
import logging
import os
import subprocess
import sys

from astropy.table import QTable
from tqdm.auto import tqdm

try:
    # Try relative imports (if running as part of a package)
    from ..backend.condor import submit_condor_job
    from ..backend.dask import run_dask
    from ..backend.parallel import run_parallel
except (ImportError, ValueError) as e:
    logging.warning(
        f"Relative import failed ({e}). "
        "Falling back to absolute import. "
        "Tip: To use relative imports properly, run with 'python -m workflow.scheduler' from the project root."
    )
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    from backend.condor import submit_condor_job

    # from backend.dask import run_dask
    from backend.parallel import run_parallel


def setup_logging(log_dir):
    """
    Configure logging for the script.
    """
    log_file = os.path.join(log_dir, "workflow.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
        force=True,
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Submit M4OPT scheduling jobs using HTCondor, locally, or via Dask.",
        allow_abbrev=False,
    )
    parser.add_argument("--config", type=str, help="Path to .ini config file")
    args, remaining_args = parser.parse_known_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)
        cfg = config["params"]
        return argparse.Namespace(
            mission=cfg.get("mission"),
            bandpass=cfg.get("bandpass"),
            absmag_mean=cfg.getfloat("absmag_mean", fallback=-16),
            absmag_stdev=cfg.getfloat("absmag_stdev", fallback=1.3),
            exptime_min=cfg.getint("exptime_min", fallback=300),
            exptime_max=cfg.getint("exptime_max", fallback=14400),
            snr=cfg.getint("snr", fallback=10),
            delay=cfg.get("delay", fallback="15min"),
            deadline=cfg.get("deadline", fallback="24hour"),
            timelimit=cfg.get("timelimit", fallback="2hour"),
            nside=cfg.getint("nside", fallback=128),
            job_cpu=cfg.getint("job_cpu", fallback=8),
            skymap_dir=cfg.get("skymap_dir", fallback="data"),
            sched_dir=cfg.get("sched_dir", fallback="data"),
            log_dir=cfg.get("log_dir", fallback="logs"),
            prog_dir=cfg.get("prog_dir", fallback="progress"),
            event_table=cfg.get(
                "event_table", fallback="data/observing-scenarios.ecsv"
            ),
            backend=cfg.get("backend", fallback="condor"),
            n_cores=cfg.getint("n_cores", fallback=4),
        )

    parser.add_argument("--mission", type=str, required=True, help="Mission name")
    parser.add_argument("--bandpass", type=str, required=True, help="Bandpass filter")
    parser.add_argument(
        "--absmag-mean", type=float, default=-16, help="Fiducial AB magnitude mean"
    )
    parser.add_argument(
        "--absmag-stdev", type=float, default=1.3, help="Magnitude standard deviation"
    )
    parser.add_argument(
        "--exptime-min", type=int, default=300, help="Minimum exposure time (s)"
    )
    parser.add_argument(
        "--exptime-max", type=int, default=14400, help="Maximum exposure time (s)"
    )
    parser.add_argument("--snr", type=int, default=10, help="SNR threshold")
    parser.add_argument("--delay", type=str, default="0h", help="Visibility delay")
    parser.add_argument(
        "--deadline", type=str, default="24hour", help="Observation deadline"
    )
    parser.add_argument(
        "--timelimit", type=str, default="2hour", help="Observation time limit"
    )
    parser.add_argument(
        "--nside", type=int, default=128, help="HEALPix nside resolution"
    )
    parser.add_argument("--job-cpu", type=int, default=8, help="CPUs per job")
    parser.add_argument(
        "--skymap-dir", type=str, default="data", help="Input sky map directory"
    )
    parser.add_argument(
        "--sched-dir", type=str, default="data", help="Output schedule directory"
    )
    parser.add_argument("--log-dir", type=str, default="logs", help="Logging directory")
    parser.add_argument(
        "--prog-dir", type=str, default="progress", help="Progress files directory"
    )
    parser.add_argument(
        "--event-table",
        type=str,
        default="data/observing-scenarios.ecsv",
        help="Event parameters file",
    )
    parser.add_argument(
        "--backend",
        type=str,
        default="condor",
        choices=["condor", "parallel", "dask"],
        help="Execution backend",
    )
    parser.add_argument(
        "--n-cores",
        type=int,
        default=4,
        help="Number of parallel cores (parallel backend)",
    )

    return parser.parse_args(remaining_args)


def create_wrapper(run_name, event_id, args, m4opt_executable):
    os.makedirs(os.path.join(args.sched_dir, run_name), exist_ok=True)
    os.makedirs(os.path.join(args.prog_dir, run_name), exist_ok=True)

    skymap_file = os.path.join(args.skymap_dir, f"{run_name}/{event_id}.fits")
    sched_file = os.path.join(args.sched_dir, f"{run_name}/{event_id}.ecsv")
    prog_file = os.path.join(args.prog_dir, f"{run_name}/PROGRESS_{event_id}.ecsv")

    wrapper_script = os.path.join(args.log_dir, f"wrapper_{run_name}_{event_id}.sh")

    wrapper_content = f"""#!/bin/bash
{m4opt_executable} schedule \
{skymap_file} \
{sched_file} \
--mission={args.mission} \
--bandpass={args.bandpass} \
--absmag-mean={args.absmag_mean} \
--absmag-stdev={args.absmag_stdev} \
--exptime-min='{args.exptime_min} s' \
--exptime-max='{args.exptime_max} s' \
--snr={args.snr} \
--delay='{args.delay}' \
--deadline='{args.deadline}' \
--timelimit='{args.timelimit}' \
--nside={args.nside} \
--write-progress {prog_file} \
--jobs {args.job_cpu} \
--cutoff=0.1
"""
    with open(wrapper_script, "w") as f:
        f.write(wrapper_content)
    os.chmod(wrapper_script, 0o755)
    return wrapper_script


if __name__ == "__main__":
    args = parse_arguments()
    os.makedirs(args.log_dir, exist_ok=True)
    setup_logging(args.log_dir)
    os.makedirs(args.prog_dir, exist_ok=True)
    os.makedirs(args.sched_dir, exist_ok=True)

    try:
        m4opt_executable = subprocess.check_output(["which", "m4opt"]).decode().strip()
        if not os.path.exists(m4opt_executable):
            raise FileNotFoundError("m4opt executable not found.")
    except Exception:
        logging.error("m4opt executable not found in PATH.")
        sys.exit(1)

    table = QTable.read(args.event_table)

    event_ids = table["coinc_event_id"].tolist()
    run_names = table["run"].tolist()

    if args.backend == "condor":
        for run_name, event_id in tqdm(zip(run_names, event_ids), total=len(event_ids)):
            wrapper_script = create_wrapper(run_name, event_id, args, m4opt_executable)
            submit_condor_job(run_name, event_id, args, wrapper_script)

    elif args.backend == "parallel":
        wrapper_scripts = [
            create_wrapper(run_name, event_id, args, m4opt_executable)
            for run_name, event_id in zip(run_names, event_ids)
        ]
        run_parallel(args, wrapper_scripts)

    elif args.backend == "dask":
        print("Running with Dask HTCondorCluster backend")
        run_dask(run_names, event_ids, args)

    else:
        print("Unknown backend: use 'condor', 'parallel', or 'dask'")
        sys.exit(1)
