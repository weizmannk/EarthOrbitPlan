import logging
import subprocess

from joblib import Parallel, delayed


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
            ["bash", wrapper_script],
            check=True,
            capture_output=True,
            stderr=subprocess.STDOUT,
        )
        print(result.stdout.decode().strip())
        logging.info(f"Successfully executed {wrapper_script}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing {wrapper_script}: {e.stdout.decode().strip()}")


def run_parallel(args, wrapper_scripts):
    """
    Run M4OPT schedule jobs in parallel using joblib.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed arguments.
    wrapper_scripts : list
        List of wrapper script paths to execute.
    """
    print(f"Running in parallel mode with {args.n_cores} cores")
    Parallel(n_jobs=args.n_cores)(
        delayed(run_script_locally)(script) for script in wrapper_scripts
    )
