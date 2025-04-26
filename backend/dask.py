import os
import shlex

from dask_jobqueue import HTCondorCluster
from distributed import Client, as_completed
from m4opt._cli import app
from tqdm import tqdm


def run_dask(event_ids, args):
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
        skymap_file = os.path.join(args.skymap_dir, f"{run_name}/{event_id}.fits")
        sched_file = os.path.join(args.sched_dir, f"{run_name}/{event_id}.ecsv")
        cmdline = (
            f"schedule {skymap_file} {sched_file} --mission={args.mission} --bandpass={args.bandpass} "
            f"--absmag-mean={args.absmag_mean} --absmag-stdev={args.absmag_stdev} --exptime-min={args.exptime_min}s "
            f"--exptime-max={args.exptime_max}s --snr={args.snr} --delay='{args.delay}' --deadline='{args.deadline}' "
            f"--timelimit='{args.timelimit}' --memory='{args.memory}' --nside={args.nside} --cutoff=0.1 --jobs={args.jobs}"
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
        memory="50GB",
        disk="8GB",
        job_extra={"accounting_group": "ligo.dev.o4.cbc.pe.bayestar"},
        job_script_prologue=["export OMP_NUM_THREADS=1"],
    )
    client = Client(cluster)
    cluster.adapt(minimum=1, maximum=20)

    futures = client.map(task, event_ids)
    for future in tqdm(as_completed(futures), total=len(futures)):
        future.result()
