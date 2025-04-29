.. _params_ini:

Configuration File Reference
=============================

The `params.ini` file defines all user-configurable parameters for the `scheduler.py` script
and for the `zenodo_unpacker.py` script.
These scripts allow you to run M4OPT scheduling on a batch of gravitational wave skymaps
and to process injection datasets from Zenodo archives, respectively.

Usage
-----

To use this configuration file with the scheduler:

.. code-block:: bash

   python scheduler.py --config params.ini

To use this configuration file with the unpacker:

.. code-block:: bash

   python zenodo_unpacker.py --config params.ini

Supported Backends
------------------

The following execution backends are supported via the `backend` field:

- ``condor``: Submits each scheduling job to an HTCondor batch system.
- ``parallel``: Executes jobs locally using parallel processing with Joblib.
- ``dask``: Distributes jobs dynamically across an HTCondor cluster via Dask.

Example Configuration
----------------------

Below is a complete configuration example for M4OPT scheduling using ULTRASAT parameters.

.. literalinclude:: ../../params.ini
   :language: ini
   :caption: Example `params.ini` for ULTRASAT scheduling
