.. _m4opt_scheduler:

Documentation for m4opt-scheduler.py: Batch Scheduling Script
===========================================================

This script executes M4OPT-based scheduling over a batch of gravitational wave sky maps
using different backend methods (Condor, local parallel, or Dask). It supports command-line
configuration and .ini-based setup.

You can run the script with:

.. code-block:: bash

    python scheduler.py --config params.ini


Qucik Start
-----------
Here we show an example how to run your first test in M4OPT with ULTRASAT mission


Functions
---------

.. literalinclude:: ../../workflow/scheduler.py
   :language: python
   :caption: Full code of `scheduler.py`


`m4opt_scheduler`: Function to execute M4OPT scheduling.
----------------
.. automodapi:: workflow.scheduler
   .. :show-inheritance:
   .. :members:
   .. :private-members:
   .. :undoc-members:
