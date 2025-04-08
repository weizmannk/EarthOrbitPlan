.. _params_ini:

params.ini — Configuration File Reference
=========================================

The `params.ini` file defines all user-configurable parameters for the `m4opt-scheduler.py` script.
This script allows you to run M4OPT scheduling on a batch of gravitational wave sky maps using different execution backends.

Usage
-----

To use this configuration file with the scheduler:

.. code-block:: bash

   python m4opt-scheduler.py --config params.ini

Supported Backends
------------------

The following execution backends are supported via the `backend` field:

- ``condor``: Submits each job to an HTCondor batch system.
- ``parallel``: Executes jobs locally in parallel using Joblib.
- ``dask``: Distributes jobs across an HTCondor cluster using Dask.

Example Configuration
---------------------

Below is a complete example configuration for the ULTRASAT mission during the O5 observing run.

.. literalinclude:: ../../params.ini
   :language: ini
   :caption: Example configuration for M4OPT scheduling using ULTRASAT
