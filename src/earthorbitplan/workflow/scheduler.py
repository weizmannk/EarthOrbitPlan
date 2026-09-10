#!/usr/bin/env python

# m4opt-scheduler: Batch Scheduling for Gravitational Wave Follow-up
# ==================================================================

# This script executes batch scheduling using the M4OPT framework for a list of gravitational wave sky maps.
# It supports various execution backends (HTCondor, local parallel with joblib, or Dask) and can be configured
# via command-line arguments or a configuration `.ini` file.

# Usage
# -----
# You can run the script either with Command-Line Interface arguments:

#     python ./earthorbitplan/workflow/scheduler.py --mission ULTRASAT --bandpass NUV ...

# Or with a configuration file:

#     python ./earthorbitplan/workflow/scheduler.py --config ./earthorbitplan/config/params_ultrasat.ini

# In the root directory of the project:
#     python -m earthorbitplan.workflow.scheduler --config ./earthorbitplan/config/params_ultrasat.ini


import argparse
import configparser
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.table import QTable
from tqdm.auto import tqdm

from earthorbitplan.backend.condor import submit_condor_job
from earthorbitplan.backend.dask import run_dask
from earthorbitplan.backend.parallel import run_parallel


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


def parse_inf_int(value):
    """Parse int or 'inf'"""
    if value.lower() == "inf":
        return np.inf
    return int(value)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Submit M4OPT scheduling jobs using HTCondor, locally, or via Dask.",
        allow_abbrev=False,
    )
    parser.add_argument("--config", type=str, help="Path to .ini config file")
    parser.add_argument(
        "--backend",
        type=str,
        default="condor",
        choices=["condor", "parallel", "dask"],
        help="Execution backend",
    )

    args, remaining_args = parser.parse_known_args()

    if args.config:
        config = configparser.ConfigParser(inline_comment_prefixes=("#",))
        config.read(args.config)

        def get(section, key, fallback=None):
            return config.get(section, key, fallback=fallback)

        def getf(section, key, fallback):
            v = config.get(section, key, fallback=None)
            return fallback if v is None else float(v)

        def geti(section, key, fallback):
            v = config.get(section, key, fallback=None)
            return fallback if v is None else int(v)

        def getb(section, key, fallback):
            v = config.get(section, key, fallback=None)
            return (
                fallback
                if v is None
                else v.strip().lower() in ("1", "true", "yes", "on")
            )

        def get_opt(section, key):
            """Return the raw string, or None when the key is missing/blank or
            the literal ``none`` (so the m4opt flag is omitted and its default
            used)."""
            v = config.get(section, key, fallback=None)
            if v is None or v.strip().lower() in ("", "none"):
                return None
            return v

        return argparse.Namespace(
            # [mission]
            mission=get("mission", "mission"),
            skygrid=get_opt("mission", "skygrid"),
            bandpass=get("mission", "bandpass"),
            nside=geti("mission", "nside", 128),
            # [kilonova]
            absmag_mean=getf("kilonova", "absmag_mean", -16),
            absmag_stdev=getf("kilonova", "absmag_stdev", 1.3),
            # [scheduling]
            snr=geti("scheduling", "snr", 10),
            delay=get("scheduling", "delay", "15min"),
            deadline=get("scheduling", "deadline", "24hour"),
            exptime_min=geti("scheduling", "exptime_min", 300),
            exptime_max=getf("scheduling", "exptime_max", np.inf),
            # Optional: blank / missing -> flag omitted, m4opt default used.
            visits=get_opt("scheduling", "visits"),
            cadence=get_opt("scheduling", "cadence"),
            appmag_dist=(
                None
                if config.get("scheduling", "appmag_dist", fallback=None) in (None, "")
                else getb("scheduling", "appmag_dist", True)
            ),
            max_fields=get_opt("scheduling", "max_fields"),
            # [solver]
            timelimit=get("solver", "timelimit", "20min"),
            memory=get("solver", "memory", ""),
            jobs=geti("solver", "jobs", 0),
            # [paths]
            data_dir=get("paths", "data_dir", "data"),
            skymap_dir=get("paths", "skymap_dir", "skymaps"),
            sched_dir=get("paths", "sched_dir", "schedules"),
            log_dir=get("paths", "log_dir", "logs"),
            prog_dir=get("paths", "prog_dir", "progress"),
            event_table=get("paths", "event_table", "observing-scenarios.ecsv"),
            # [execution]
            backend=get("execution", "backend", "condor"),
            n_cores=geti("execution", "n_cores", 4),
        )

    parser.add_argument("--mission", type=str, required=True, help="Mission name")
    parser.add_argument(
        "--skygrid",
        type=str,
        default=None,
        help="Name of sky grid to use, if the mission supports multiple sky grids",
    )
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
        "--exptime-max",
        type=parse_inf_int,
        default=np.inf,
        help="Maximum exposure time (s)",
    )
    parser.add_argument("--snr", type=int, default=10, help="SNR threshold")
    parser.add_argument("--delay", type=str, default="0h", help="Visibility delay")
    parser.add_argument(
        "--deadline", type=str, default="24hour", help="Observation deadline"
    )
    # These four pass straight through to `m4opt schedule`. Left as None by
    # default so the flag is omitted and m4opt's own default applies.
    parser.add_argument(
        "--visits", type=int, default=None, help="Number of observations per field"
    )
    parser.add_argument(
        "--cadence", type=str, default=None, help="Minimum separation between visits"
    )
    parser.add_argument(
        "--appmag-dist",
        dest="appmag_dist",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Variable exposure with an absolute-magnitude distribution",
    )
    parser.add_argument(
        "--max-fields",
        type=int,
        default=None,
        help="Keep only the N most probable fields",
    )
    parser.add_argument(
        "--timelimit", type=str, default="20min", help="Time limit for MILP solver"
    )
    parser.add_argument(
        "--memory",
        type=str,
        default="",
        help="Maximum solver memory usage before terminating",
    )
    parser.add_argument(
        "--nside", type=int, default=128, help="HEALPix nside resolution"
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=0,
        help="Threads for solving one MILP (0 = all cores)",
    )

    parser.add_argument("--data-dir", type=str, default="data", help="Data directory")
    parser.add_argument(
        "--skymap-dir", type=str, default="skymaps", help="Input sky map directory"
    )
    parser.add_argument(
        "--sched-dir", type=str, default="schedules", help="Output schedule directory"
    )
    parser.add_argument("--log-dir", type=str, default="logs", help="Logging directory")
    parser.add_argument(
        "--prog-dir", type=str, default="progress", help="Progress files directory"
    )
    parser.add_argument(
        "--event-table",
        type=str,
        default="observing-scenarios.ecsv",
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


def create_wrapper(run_name, event_id, base_path, args, m4opt_executable):
    # Build directories from base_path and args
    sched_dir = base_path / args.sched_dir
    prog_dir = base_path / args.prog_dir
    log_dir = base_path / args.log_dir
    skymap_dir = base_path / args.skymap_dir

    # File paths
    skymap_file = skymap_dir / run_name / f"{event_id}.fits"
    sched_file = sched_dir / run_name / f"{event_id}.ecsv"
    prog_file = prog_dir / run_name / f"PROGRESS_{event_id}.ecsv"
    wrapper_script = log_dir / f"wrapper_{run_name}_{event_id}.sh"

    cmd = [
        m4opt_executable,
        "schedule",
        str(skymap_file),
        str(sched_file),
        f"--mission={args.mission}",
        f"--bandpass={args.bandpass}",
        f"--absmag-mean={args.absmag_mean}",
        f"--absmag-stdev={args.absmag_stdev}",
        f"--exptime-min='{args.exptime_min} s'",
        f"--exptime-max='{args.exptime_max} s'",
        f"--snr={args.snr}",
        f"--delay='{args.delay}'",
        f"--deadline='{args.deadline}'",
        f"--timelimit='{args.timelimit}'",
        f"--nside={args.nside}",
        f"--write-progress {prog_file}",
        f"--jobs {args.jobs}",
        "--cutoff=0.1",
    ]
    # Optional pass-through flags: only added when explicitly set. Otherwise
    # m4opt applies its own default (and stays correct if that default moves).
    if args.skygrid:
        cmd.insert(5, f"--skygrid={args.skygrid}")
    if args.memory:
        cmd.append(f"--memory='{args.memory}'")
    if args.visits is not None:
        cmd.append(f"--visits={args.visits}")
    if args.cadence is not None:
        cmd.append(f"--cadence='{args.cadence}'")
    if args.appmag_dist is not None:
        cmd.append("--appmag-dist" if args.appmag_dist else "--no-appmag-dist")
    if args.max_fields is not None:
        cmd.append(f"--max-fields={args.max_fields}")

    # --event-time is intentionally omitted: m4opt reads the trigger time from
    # the DATE-OBS / gps_time header of the input sky map.

    wrapper_content = "#!/bin/bash\nexport OMP_NUM_THREADS=1\n\n" + " ".join(cmd) + "\n"

    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        with open(wrapper_script, "w") as f:
            f.write(wrapper_content)
        os.chmod(wrapper_script, 0o755)
        return wrapper_script
    except Exception as e:
        logging.error(f"Failed to create wrapper script for event {event_id}: {e}")
        return None


if __name__ == "__main__":
    args = parse_arguments()

    # Set up main data directory
    base_path = Path(args.data_dir)
    base_path.mkdir(exist_ok=True)

    # Create necessary subdirectories
    log_dir = base_path / args.log_dir
    log_dir.mkdir(parents=True, exist_ok=True)

    prog_dir = base_path / args.prog_dir
    prog_dir.mkdir(parents=True, exist_ok=True)

    sched_dir = base_path / args.sched_dir
    sched_dir.mkdir(parents=True, exist_ok=True)

    skymap_dir = base_path / args.skymap_dir
    skymap_dir.mkdir(parents=True, exist_ok=True)

    # Input event table path
    event_table_path = base_path / args.event_table

    # Initialize logging
    setup_logging(str(log_dir))

    try:
        m4opt_executable = subprocess.check_output(["which", "m4opt"]).decode().strip()
        if not os.path.exists(m4opt_executable):
            raise FileNotFoundError("m4opt executable not found.")
    except Exception:
        logging.error("m4opt executable not found in PATH.")
        sys.exit(1)

    try:
        table = QTable.read(event_table_path)
    except Exception as e:
        logging.error(f"Failed to read event table {event_table_path}: {e}")
        sys.exit(1)

    event_ids = table["coinc_event_id"].tolist()
    run_names = table["run"].tolist()

    # Run HTCondor process
    if args.backend == "condor":
        for run_name, event_id in tqdm(zip(run_names, event_ids), total=len(event_ids)):
            # Create run subdirectories if needed
            (sched_dir / run_name).mkdir(parents=True, exist_ok=True)
            (prog_dir / run_name).mkdir(parents=True, exist_ok=True)

            skymap_file = skymap_dir / run_name / f"{event_id}.fits"
            if not skymap_file.exists():
                logging.warning(
                    f"Missing skymap: {skymap_file}. Skipping event {event_id}."
                )
                continue

            wrapper_script = create_wrapper(
                run_name, event_id, base_path, args, m4opt_executable
            )
            if wrapper_script is None:
                logging.error(
                    f"Failed to create wrapper script for event {event_id}, skipping submission."
                )
                continue

            submit_condor_job(run_name, event_id, log_dir, wrapper_script)

    # Parallel backend
    elif args.backend == "parallel":
        wrapper_scripts = [
            create_wrapper(run_name, event_id, base_path, args, m4opt_executable)
            for run_name, event_id in zip(run_names, event_ids)
        ]
        run_parallel(args, wrapper_scripts)

    elif args.backend == "dask":
        print("Running with Dask HTCondorCluster backend")
        run_dask(run_names, event_ids, args)

    else:
        print("Unknown backend: use 'condor', 'parallel', or 'dask'")
        sys.exit(1)
