import logging
import re
import subprocess
import textwrap

from memory_profiler import profile

_UNITS = {"": 1, "KB": 1, "K": 1, "MB": 1024, "M": 1024, "GB": 1024**2, "G": 1024**2}


def _to_kb(value):
    """Convert a Condor-style size ("80000 MB", "8 GB", "4096") to KiB.

    HTCondor treats a bare ``request_memory`` as MiB and a bare
    ``request_disk`` as KiB, so expressions that mix a configured value with
    ``MemoryUsage`` / ``DiskUsage`` must agree on one unit. We normalise
    everything to KiB and re-scale at the point of use.
    """
    match = re.fullmatch(r"\s*([0-9.]+)\s*([KMG]?B?)\s*", str(value), re.IGNORECASE)
    if match is None:
        raise ValueError(f"cannot parse size: {value!r}")
    number, unit = match.groups()
    return int(float(number) * _UNITS[unit.upper()])


@profile
def submit_condor_job(
    run_name,
    event_id,
    log_dir,
    wrapper_script,
    request_memory="40000 MB",
    request_disk="8000 MB",
    max_retries=3,
    growth_factor=2.0,
    max_request_memory="200000 MB",
):
    """Submit one m4opt scheduling job to HTCondor.

    ``request_memory`` is the cgroup limit for the whole process and is
    unrelated to the ``--memory`` flag passed to m4opt, which only caps
    CPLEX's node file. The MILP model itself lives outside that cap.

    Requests are adaptive rather than fixed. ``request_memory`` and
    ``request_disk`` are the starting floors; if the job is held for going
    over either limit, ``periodic_release`` puts it back in the queue and the
    ``ifthenelse`` expressions re-request ``growth_factor`` times the usage
    actually measured on the previous run. This happens at most
    ``max_retries`` times. Starting low keeps the job matchable on many more
    slots; only the events that genuinely need a large machine escalate to
    one.

    Escalation is capped at ``max_request_memory``, which must stay below the
    RAM of the largest slot in the pool (``condor_status -af Memory | sort -n
    | tail -1``). Without that cap a job could escalate past every machine and
    sit idle forever instead of failing visibly.
    """
    mem_kb = _to_kb(request_memory)
    mem_cap_kb = _to_kb(max_request_memory)
    disk_kb = _to_kb(request_disk)
    # request_memory is expressed in MiB, request_disk in KiB.
    mem_mb = max(1, mem_kb // 1024)
    mem_cap_mb = max(mem_mb, mem_cap_kb // 1024)
    condor_submit_script = textwrap.dedent(f"""\
        +MaxHours = 24
        universe = vanilla
        accounting_group = ligo.dev.o4.cbc.pe.bayestar
        executable = {wrapper_script}
        output = {log_dir}/$(Cluster)_$(Process).out
        error = {log_dir}/$(Cluster)_$(Process).err
        log = {log_dir}/$(Cluster)_$(Process).log
        request_cpus = 1
        request_memory = ifthenelse(MemoryUsage =!= undefined, min({{{mem_cap_mb}, max({{{mem_mb}, MemoryUsage * {growth_factor}}})}}), {mem_mb})
        request_disk = ifthenelse(DiskUsage =!= undefined, max({{{disk_kb}, DiskUsage * {growth_factor}}}), {disk_kb})
        periodic_release = (HoldReasonCode == 34 || HoldReasonCode == 104) && (NumJobStarts <= {max_retries})
        on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
        on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
        on_exit_hold_reason = (ExitBySignal == True \
            ? strcat("The job exited with signal ", ExitSignal) \
            : strcat("The job exited with code ", ExitCode))
        environment = "OMP_NUM_THREADS=1"
        queue 1
        """)

    try:
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
            logging.info(f"Submitted condor job for {event_id}")

    except Exception as e:
        logging.error(f"Error submitting Condor job for {event_id}: {e}")
