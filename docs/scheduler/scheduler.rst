.. _scheduler:

Documentation for scheduler.py: Batch Scheduling Script
=============================================================

This script executes :math:`\mathrm{M^4OPT}`-based scheduling over a batch of gravitational wave sky maps
using different backend methods (Condor, local parallel, or Dask). It supports command-line
configuration and .ini-based setup.

You can run the script with:

.. code-block:: bash

    python earthorbitplan/workflow/scheduler.py --config ./earthorbitplan/config/params_ultrasat.ini

.. Functions
.. ---------

.. .. literalinclude:: ../../workflow/scheduler.py
..    :language: python
..    :caption: Full code of `scheduler.py`


`m4opt_scheduler`: Function to execute M4OPT scheduling
-------------------------------------------------------
.. automodapi:: earthorbitplan.workflow.scheduler
   .. :show-inheritance:
   .. :members:
   .. :private-members:
   .. :undoc-members:
