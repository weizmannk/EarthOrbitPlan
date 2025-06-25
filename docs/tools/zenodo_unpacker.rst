.. _zenodo_unpacker:

Zenodo Injection Unpacker
==========================

This module automates the unpacking, filtering, and conversion of injection datasets
(e.g., Farah / GWTC-3) from Zenodo ZIP archives. It processes event tables and associated
localization files for specific observing runs (e.g., O5, O6), and outputs
filtered ECSV tables and organized FITS files.

Usage
-----

You can run this unpacker from the command line as follows:

.. code-block:: bash

   python zenodo_unpacker.py --zip runs_SNR-10.zip --subdir runs_SNR-10 --runs O5 O6 --detectors HLVK --outdir ./data --mass-threshold 3

Source
------

Zenodo Dataset: https://zenodo.org/records/14585837

.. Full Code
.. ---------

.. .. literalinclude:: ../../scenarios/zenodo_unpacker.py
..    :language: python
..    :caption: Full code of `zenodo_unpacker.py`

Module Reference
----------------

.. automodapi:: scenarios.zenodo_unpacker
   .. :show-inheritance:
   .. :members:
   .. :private-members:
   .. :undoc-members:
   .. :special-members: __init__, __call__
   .. :exclude-members: __weakref__, __dict__, __module__, __class__
