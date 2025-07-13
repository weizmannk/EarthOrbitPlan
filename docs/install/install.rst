.. _install:

This guide explains how to install the required components to run EarthOrbitPlan. The pipeline is built around
the ` M⁴OPT toolkit <https://m4opt.readthedocs.io/en/latest/>`_
and also makes use of `observing-scenarios`_  to generate realistic
gravitational wave follow-up situations. It is recommended to use a virtual environment (e.g. via `conda` or `venv`)
to avoid conflicts between dependencies.

.. note::

   EarthOrbitPlan is designed for educational use and builds directly on :math:`\mathrm{M^4OPT}`.
   It provides a guided way to run simulations based on skymaps from GW detectors and
   schedule follow-up facilities with  such as UVEX, ULTRASAT, ZTF, and Rubin.

.. dropdown:: Requirements

   You will need the following:

   - Python >= 3.11
   - `M⁴OPT <https://m4opt.readthedocs.io/en/latest/>`_
   - `observing-scenarios`_


.. dropdown:: Environment setup

   .. tab-item::

      .. tab-item:: Conda

         .. code-block:: bash

         conda create -n earthorbitplan_env python=3.11
         conda activate earthorbitplan_env

      .. tab-item:: Python venv

         .. code-block:: bash

               python3.11 -m venv earthorbitplan-env
               source earthorbitplan-env/bin/activate
               pip install --upgrade pip


.. dropdown:: Install EarthOrbitPlan

   .. tab-item::

      .. tab-item:: PyPI

         The recommended way to install EarthOrbitPlan is using pip:

         .. code-block:: bash

            pip install earthorbitplan

      .. tab-item:: From source (editable)

          Alternatively, you can clone the repository and install it in editable (development) mode:

         .. code-block:: bash

            git clone https://github.com/weizmannk/EarthOrbitPlan.git
            cd EarthOrbitPlan/
            pip install -e .
            cd ..

.. dropdown:: Install CPLEX

    CPLEX is required for optimization.
    See the official instructions:
    https://m4opt.readthedocs.io/en/latest/install/cplex.html

.. dropdown:: Check your installation

   .. code-block:: bash

      m4opt --help
      m4opt schedule --help

.. note::

   For most EarthOrbitPlan use cases, installing :math:`\mathrm{M^4OPT}` and its dependencies is sufficient.
   If you wish to simulate additional GW-only scenarios, consider installing the `observing-scenarios` package as below.


.. dropdown:: Optional: Observing Scenarios

   .. code-block:: bash

      curl -sSL https://install.python-poetry.org | python3 -
      git clone https://github.com/lpsinger/observing-scenarios-simulations.git
      cd observing-scenarios-simulations
      poetry install
      poetry shell

.. note::

   The `observing-scenarios`_ package is optional, but useful for testing standalone GW follow-up strategies without electromagnetic scheduling.
   It provides realistic skymaps and scenarios commonly used in follow-up campaigns.


.. _observing-scenarios: https://github.com/lpsinger/observing-scenarios-simulations

.. _m4opt: https://m4opt.readthedocs.io/en/latest/
