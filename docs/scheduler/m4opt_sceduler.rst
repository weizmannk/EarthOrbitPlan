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


.. automodapi:: workflows.m4opt-scheduler
   :show-inheritance:
   :members:
   :private-members:
   :undoc-members:
