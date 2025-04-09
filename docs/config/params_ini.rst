.. _params_ini:

Configuration File Reference
=============================

The `params.ini` file defines all user-configurable parameters for the `m4opts_cheduler.py` script.
This script allows you to run M4OPT scheduling on a batch of gravitational wave sky maps using different execution backends.

Usage
-----

To use this configuration file with the scheduler:

.. code-block:: bash

   python m4opt_scheduler.py --config params.ini

Supported Backends
------------------

The following execution backends are supported via the `backend` field:

- ``condor``: Submits each job to an HTCondor batch system.
- ``parallel``: Executes jobs locally in parallel using Joblib.
- ``dask``: Distributes jobs across an HTCondor cluster using Dask.

Example Configuration
---------------------

Below is a complete parameters configuration for M4OPT scheduling using ULTRASAT

.. literalinclude:: ../../params.ini
   :language: ini
   :caption: Ultrasat parameters `params.ini`
