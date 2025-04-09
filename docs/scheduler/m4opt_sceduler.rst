.. _m4opt_scheduler:

Documentation for m4opt-scheduler.py: Batch Scheduling Script
===========================================================

This script executes M4OPT-based scheduling over a batch of gravitational wave sky maps
using different backend methods (Condor, local parallel, or Dask). It supports command-line
configuration and .ini-based setup.

You can run the script with:

.. code-block:: bash

    python m4opt-scheduler.py --config params.ini

Functions
---------

.. autofunction:: parse_arguments
   Parses arguments from the command line or .ini config file into an argparse.Namespace object.

.. autofunction:: create_wrapper
   Creates a bash wrapper script to run the M4OPT scheduling command for a given event.

.. autofunction:: run_script_locally
   Executes a bash wrapper script locally and logs output or errors.

.. autofunction:: run_parallel
   Uses Joblib to execute scheduling wrapper scripts in parallel on multiple CPU cores.

.. autofunction:: run_dask
   Schedules M4OPT tasks on a Dask HTCondor cluster dynamically using dask_jobqueue.

.. autofunction:: submit_condor_job
   Submits the bash wrapper script as a job to an HTCondor queue using condor_submit.


.. Module Reference
.. ----------------

.. .. automodule:: m4opt_scheduler
..    :members:
..    :undoc-members:
..    :show-inheritance:
