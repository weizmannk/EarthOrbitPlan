.. _install:

============
Installation
============

This guide explains how to install the required components to run EarthOrbitPlan. The pipeline is built around
the `M⁴OPT toolkit <https://m4opt.readthedocs.io/en/latest/>`_
and also makes use of `observing-scenarios`_  to generate realistic
gravitational wave follow-up situations. It is recommended to use a virtual environment (e.g. via `conda` or `venv`)
to avoid conflicts between dependencies.

.. note::

   EarthOrbitPlan is designed for educational use and builds directly on :math:`\mathrm{M^4OPT}`.
   It provides a guided way to run simulations based on skymaps from GW detectors and
   schedule follow-up facilities with  such as UVEX, ULTRASAT, ZTF, and Rubin.

.. dropdown:: Requirements

   You will need the following:

   .. button-link:: https://m4opt.readthedocs.io/en/latest/
      :color: info
      :shadow:

      :math:`\mathrm{M^4OPT}`

   .. button-link::  https://github.com/lpsinger/observing-scenarios-simulations/
      :color: info
      :shadow:

      Observing Scenarios

   :bdg-warning:`Python >= 3.11`

.. dropdown:: Environment setup

   .. tab-set::

      .. tab-item:: Conda

         .. code-block:: bash

            conda create -n earthorbitplan_env python=3.11
            conda activate earthorbitplan_env

      .. tab-item:: Python

         .. code-block:: bash

            python3.11 -m venv earthorbitplan-env
            source earthorbitplan-env/bin/activate
            pip install --upgrade pip


.. dropdown:: Install EarthOrbitPlan

   .. tab-set::

      .. tab-item:: PyPI

         .. admonition:: Recommended
            :class: success

            The recommended way to install EarthOrbitPlan is using pip:

            .. code-block:: console

               $  pip install earthorbitplan

      .. tab-item:: Dev

         .. admonition:: Development install
            :class: info

            To install EarthOrbitPlan in development mode, run:

            .. code-block:: bash

               git clone https://github.com/weizmannk/EarthOrbitPlan.git
               cd EarthOrbitPlan/
               pip install -e .
               cd ..

      .. tab-item:: Latest

         .. admonition:: Unreleased version
            :class: danger

            If you want to use the latest version of EarthOrbitPlan, you can install the current main branch directly from GitHub:

            .. code-block:: console

               $ pip install git+https://github.com/weizmannk/EarthOrbitPlan.git@main

.. dropdown:: Install CPLEX

    CPLEX is required for optimization.
    See the official instructions:

   .. button-link:: https://m4opt.readthedocs.io/en/latest/install/cplex.html
      :color: info
      :shadow:

      CPLEX Install Guide

.. dropdown:: Check your installation

   .. code-block:: bash

      m4opt --help
      m4opt schedule --help

.. important::

   For most EarthOrbitPlan use cases, installing :math:`\mathrm{M^4OPT}` and its dependencies is sufficient.
   If you wish to simulate additional GW-only scenarios, consider installing the `observing-scenarios`_ package as below.


.. dropdown:: Optional: Observing Scenarios

   .. code-block:: bash

      curl -sSL https://install.python-poetry.org | python3 -
      git clone https://github.com/lpsinger/observing-scenarios-simulations.git
      cd observing-scenarios-simulations
      poetry install
      poetry shell

.. seealso::

   The `observing-scenarios`_ package is optional, but useful for testing standalone GW follow-up strategies without electromagnetic scheduling.
   It provides realistic skymaps and scenarios commonly used in follow-up campaigns.


.. _observing-scenarios: https://github.com/lpsinger/observing-scenarios-simulations

.. _m4opt: https://m4opt.readthedocs.io/en/latest/
