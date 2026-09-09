import logging
import os
import subprocess


def submit_slurm(wrapper_scripts, log_dir, job_name="M4OPT_batch"):
    """
    Submit a SLURM array job where each array task runs one wrapper script.

    Parameters
    ----------
    wrapper_scripts : list of str
        Paths to the wrapper scripts to be executed.
    log_dir : str
        Directory where logs and intermediate files are saved.
    job_name : str, optional
        Name of the SLURM job (default is "M4OPT_batch").
    """
    if not wrapper_scripts:
        logging.error("No wrapper scripts provided for SLURM submission.")
        return

    # Make sure log directory exists
    os.makedirs(log_dir, exist_ok=True)

    # Save list of wrapper scripts
    wrapper_list_path = os.path.join(log_dir, "wrapper_list.txt")
    with open(wrapper_list_path, "w") as f:
        for wrapper in wrapper_scripts:
            f.write(f"{wrapper}\n")

    array_length = len(wrapper_scripts)
    slurm_array_script_path = os.path.join(log_dir, "slurm_batch_array.sh")

    # Write the SLURM array batch script
    slurm_array_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={log_dir}/slurm_%A_%a.out
#SBATCH --error={log_dir}/slurm_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --account=bcrv-delta-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --array=0-{array_length - 1}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=leggi014@umn.edu

export OMP_NUM_THREADS=1

WRAPPER=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" {wrapper_list_path})
echo "Running wrapper: $WRAPPER"
bash "$WRAPPER"
"""

    try:
        with open(slurm_array_script_path, "w") as f:
            f.write(slurm_array_script)
        os.chmod(slurm_array_script_path, 0o755)

        # Submit the SLURM array job
        subprocess.run(["sbatch", slurm_array_script_path], check=True)
        print(f"Submitted SLURM array job with {array_length} tasks.")
    except Exception as e:
        logging.error(f"Failed to submit SLURM array job: {e}")
