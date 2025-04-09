#!/usr/bin/env python
"""
Documentation for m4opt-scheduler.py: Batch Scheduling Script
===========================================================

This script executes M4OPT-based scheduling over a batch of gravitational wave sky maps
using different backend methods (Condor, local parallel, or Dask). It supports command-line
configuration and .ini-based setup.
"""

import argparse
import configparser
import logging
import os
import subprocess
import sys

from astropy.table import QTable
from joblib import Parallel, delayed
from tqdm.auto import tqdm


def parse_arguments():
    """
    Parse command-line arguments or load them from a .ini configuration file.

    If the --config option is provided, the function loads parameters from the specified
    .ini file under the [params] section. Otherwise, it uses standard CLI arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments.
    """
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
            job_cpu=cfg.getint("job", fallback=8),
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

    parser.add_argument("--mission", help="Mission name", type=str, required=True)
    parser.add_argument(
        "--bandpass", help="Mission bandpass, filter", type=str, required=True
    )
    parser.add_argument(
        "--absmag-mean",
        help="Fiducial kilonova AB absolute bolometric magnitude of the source",
        type=float,
        default=-16,
    )
    parser.add_argument(
        "--absmag-stdev",
        help="Standard deviation of AB absolute magnitude of the source",
        type=float,
        default=1.3,
    )
    parser.add_argument(
        "--exptime-min",
        help="Minimum per-tile exposure time in seconds",
        type=int,
        default=300,
    )
    parser.add_argument(
        "--exptime-max",
        help="Maximum per-tile exposure time in seconds",
        type=int,
        default=14400,
    )
    parser.add_argument(
        "--snr", help="Signal-to-noise ratio threshold", type=int, default=10
    )
    parser.add_argument(
        "--delay",
        help="Time until the kilonova event becomes visible",
        type=str,
        default="0h",
    )
    parser.add_argument(
        "--deadline",
        help="Deadline for observation completion",
        type=str,
        default="24hour",
    )
    parser.add_argument("--timelimit", help="", type=str, default="2hour")
    parser.add_argument("--nside", help="HEALPix resolution", type=int, default=128)
    parser.add_argument("--job-cpu", type=int, default=8)
    parser.add_argument(
        "--skymap-dir", help="GW Sky map filename", type=str, default="data"
    )
    parser.add_argument(
        "--sched-dir",
        help="Output filename for generated schedule",
        type=str,
        default="data",
    )
    parser.add_argument(
        "--log-dir", help="Track simulations process", type=str, default="logs"
    )
    parser.add_argument(
        "--prog-dir",
        help="Save a time series of the CPLEX objective value and best bound to this file",
        type=str,
        default="progress",
    )
    parser.add_argument(
        "--event-table",
        help="Observing scenarios parameters",
        type=str,
        default="data/observing-scenarios.ecsv",
    )
    parser.add_argument(
        "--backend",
        help="Selection of the type of process",
        type=str,
        default="condor",
        choices=["condor", "parallel", "dask"],
    )
    parser.add_argument(
        "--n-cores", help="Number of jobs for parallel processing", type=int, default=4
    )

    return parser.parse_args(remaining_args)


def create_wrapper(event_id, args, m4opt_path):
    """
    Create a bash wrapper script for submitting an M4OPT schedule job.

    Parameters
    ----------
    event_id : str or int
        Unique identifier for the sky map event.
    args : argparse.Namespace
        Parsed command-line or config arguments.
    m4opt_path : str
        Path to the m4opt executable.

    Returns
    -------
    str
        Path to the created wrapper script.
    """
    skymap_file = os.path.join(args.skymap_dir, f"{event_id}.fits")
    sched_file = os.path.join(args.sched_dir, f"{event_id}.ecsv")
    prog_file = os.path.join(args.prog_dir, f"PROGRESS_{event_id}.ecsv")
    wrapper_script = os.path.join(args.log_dir, f"wrapper_{event_id}.sh")

    wrapper_content = f"""#!/bin/bash
{m4opt_path} schedule \
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


def run_script_locally(wrapper_script):
    """
    Execute a wrapper bash script locally.

    Parameters
    ----------
    wrapper_script : str
        Path to the wrapper shell script to execute.
    """
    try:
        result = subprocess.run(
            ["bash", wrapper_script], check=True, capture_output=True
        )
        print(result.stdout.decode().strip())
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing {wrapper_script}: {e.stderr.decode().strip()}")


def run_parallel(event_ids, args, m4opt_path):
    """
    Run M4OPT schedule jobs in parallel using joblib.

    Parameters
    ----------
    event_ids : list
        List of event IDs to process.
    args : argparse.Namespace
        Parsed arguments.
    m4opt_path : str
        Path to the m4opt executable.
    """
    wrapper_scripts = [
        create_wrapper(event_id, args, m4opt_path) for event_id in event_ids
    ]
    Parallel(n_jobs=args.n_cores)(
        delayed(run_script_locally)(script) for script in wrapper_scripts
    )


def run_dask(event_ids, args):
    """
    Run M4OPT scheduling tasks using Dask with an HTCondorCluster backend.

    Parameters
    ----------
    event_ids : list
        List of event IDs to process.
    args : argparse.Namespace
        Parsed arguments.
    """
    import shlex

    from dask_jobqueue import HTCondorCluster
    from distributed import Client, as_completed
    from m4opt._cli import app

    def task(event_id):
        skymap_file = os.path.join(args.skymap_dir, f"{event_id}.fits")
        sched_file = os.path.join(args.sched_dir, f"{event_id}.ecsv")
        cmdline = (
            f"schedule {skymap_file} {sched_file} --mission={args.mission} --bandpass={args.bandpass} "
            f"--absmag-mean={args.absmag_mean} --absmag-stdev={args.absmag_stdev} --exptime-min={args.exptime_min}s "
            f"--exptime-max={args.exptime_max}s --snr={args.snr} --delay='{args.delay}' --deadline='{args.deadline}' "
            f"--timelimit='{args.timelimit}' --nside={args.nside} --cutoff=0.1 --jobs={args.job_cpu}"
        )
        try:
            app(shlex.split(cmdline))
        except SystemExit as e:
            if e.code != 0:
                raise RuntimeError(
                    f"Task for event {event_id} failed with exit code {e.code}"
                )

    cluster = HTCondorCluster(
        cores=1,
        memory="4GB",
        disk="2GB",
        job_extra={"accounting_group": "ligo.dev.o4.cbc.pe.bayestar"},
        job_script_prologue=["export OMP_NUM_THREADS=1"],
    )
    client = Client(cluster)
    cluster.adapt(minimum=1, maximum=20)

    futures = client.map(task, event_ids)
    for future in tqdm(as_completed(futures), total=len(futures)):
        future.result()


def submit_condor_job(event_id, args, m4opt_path):
    """
    Submit an M4OPT job to HTCondor using a wrapper script.

    Parameters
    ----------
    event_id : str or int
        The ID of the sky map event.
    args : argparse.Namespace
        Parsed arguments.
    m4opt_path : str
        Path to the m4opt executable.
    """
    wrapper_script = create_wrapper(event_id, args, m4opt_path)
    condor_submit_script = f"""
+MaxHours = 24
universe = vanilla
accounting_group = ligo.dev.o4.cbc.pe.bayestar
envfile = {wrapper_script}
getenv = true
executable = {wrapper_script}
output = {args.log_dir}/$(Cluster)_$(Process).out
error = {args.log_dir}/$(Cluster)_$(Process).err
log = {args.log_dir}/$(Cluster)_$(Process).log
request_memory = 50000 MB
request_disk = 8000 MB
request_cpus = 1
on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
on_exit_hold_reason = (ExitBySignal == True \
    ? strcat("The job exited with signal ", ExitSignal) \
    : strcat("The job exited with code ", ExitCode))
environment = "OMP_NUM_THREADS=1"
queue 1
"""
    proc = subprocess.Popen(
        ["condor_submit"],
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = proc.communicate(input=condor_submit_script)
    if proc.returncode != 0:
        logging.error(f"Condor submit error for {event_id}: {stderr.strip()}")
    else:
        print(f"Submitted condor job for {event_id}: {stdout.strip()}")


if __name__ == "__main__":
    args = parse_arguments()
    os.makedirs(args.log_dir, exist_ok=True)
    os.makedirs(args.prog_dir, exist_ok=True)

    try:
        m4opt_path = subprocess.check_output(["which", "m4opt"]).decode().strip()
        if not os.path.exists(m4opt_path):
            raise FileNotFoundError("m4opt executable not found.")
    except Exception:
        logging.error("m4opt executable not found in PATH.")
        sys.exit(1)

    table = QTable.read(args.event_table)
    event_ids = list(table["coinc_event_id"])

    if args.backend == "condor":
        for event_id in tqdm(event_ids):
            submit_condor_job(event_id, args, m4opt_path)
    elif args.backend == "parallel":
        print(f"Running in parallel mode with {args.n_cores} cores")
        run_parallel(event_ids, args, m4opt_path)
    elif args.backend == "dask":
        print("Running with Dask HTCondorCluster backend")
        run_dask(event_ids, args)
    else:
        print("Unknown backend: use 'condor', 'parallel', or 'dask'")
        sys.exit(1)
