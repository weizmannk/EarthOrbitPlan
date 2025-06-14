import logging
import subprocess


def submit_condor_job(run_name, event_id, log_dir, wrapper_script):
    condor_submit_script = f"""
        +MaxHours = 24
        universe = vanilla
        accounting_group = ligo.dev.o4.cbc.pe.bayestar
        getenv = true
        executable = {wrapper_script}
        output = {log_dir}/$(Cluster)_$(Process).out
        error = {log_dir}/$(Cluster)_$(Process).err
        log = {log_dir}/$(Cluster)_$(Process).log
        request_memory = 80000 MB
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
