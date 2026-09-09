import os

from dask_jobqueue import HTCondorCluster
from distributed import Client, as_completed
from tqdm import tqdm

from m4opt._cli import app


def _schedule_argv(run_name, event_id, args):
    """Build the ``m4opt schedule`` argument list (kept in sync with
    earthorbitplan.workflow.scheduler.create_wrapper)."""
    skymap_file = os.path.join(args.skymap_dir, f"{run_name}/{event_id}.fits")
    sched_file = os.path.join(args.sched_dir, f"{run_name}/{event_id}.ecsv")
    argv = [
        "schedule",
        skymap_file,
        sched_file,
        f"--mission={args.mission}",
        f"--bandpass={args.bandpass}",
        f"--absmag-mean={args.absmag_mean}",
        f"--absmag-stdev={args.absmag_stdev}",
        f"--exptime-min={args.exptime_min} s",
        f"--exptime-max={args.exptime_max} s",
        f"--snr={args.snr}",
        f"--delay={args.delay}",
        f"--deadline={args.deadline}",
        f"--timelimit={args.timelimit}",
        f"--nside={args.nside}",
        "--cutoff=0.1",
        f"--jobs={args.jobs}",
    ]
    # Optional pass-through flags: only added when explicitly set, otherwise
    # m4opt applies its own default.
    if args.skygrid:
        argv.insert(4, f"--skygrid={args.skygrid}")
    if args.memory:
        argv.append(f"--memory={args.memory}")
    if args.visits is not None:
        argv.append(f"--visits={args.visits}")
    if args.cadence is not None:
        argv.append(f"--cadence={args.cadence}")
    if args.appmag_dist is not None:
        argv.append("--appmag-dist" if args.appmag_dist else "--no-appmag-dist")
    if args.max_fields is not None:
        argv.append(f"--max-fields={args.max_fields}")
    # --event-time is intentionally omitted: the trigger time comes from the
    # DATE-OBS / gps_time header of the sky map.
    return argv


def run_dask(run_names, event_ids, args):
    """
    Run M4OPT scheduling tasks using Dask with an HTCondorCluster backend.

    Parameters
    ----------
    run_names : list
        List of run names corresponding to each event ID.
    event_ids : list
        List of event IDs to process.
    args : argparse.Namespace
        Parsed arguments.
    """

    def task(run_name, event_id):
        try:
            app(_schedule_argv(run_name, event_id, args))
        except SystemExit as e:
            if e.code != 0:
                raise RuntimeError(
                    f"Task for event {event_id} failed with exit code {e.code}"
                )

    cluster = HTCondorCluster(
        cores=1,
        memory="50GB",
        disk="8GB",
        job_extra={"accounting_group": "ligo.dev.o4.cbc.pe.bayestar"},
        job_script_prologue=["export OMP_NUM_THREADS=1"],
    )
    client = Client(cluster)
    cluster.adapt(minimum=1, maximum=20)

    futures = client.map(task, run_names, event_ids)
    for future in tqdm(as_completed(futures), total=len(futures)):
        future.result()
