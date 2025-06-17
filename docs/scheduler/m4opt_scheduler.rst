.. _m4opt_scheduler:

Documentation for m4opt-scheduler.py: Batch Scheduling Script
===========================================================

This script executes M4OPT-based scheduling over a batch of gravitational wave sky maps
using different backend methods (Condor, local parallel, or Dask). It supports command-line
configuration and .ini-based setup.

You can run the script with:

.. code-block:: bash

    python m4opt_scheduler.py --config params.ini

Functions
---------

.. literalinclude:: ../../workflow/m4opt_scheduler.py
   :language: python
   :caption: Full code of `m4opt_scheduler.py`


`m4opt_scheduler`: Function to execute M4OPT scheduling.
----------------
.. automodapi:: workflow.m4opt_scheduler
   .. :show-inheritance:
   .. :members:
   .. :private-members:
   .. :undoc-members:
