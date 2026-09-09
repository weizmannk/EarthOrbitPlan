EarthOrbitPlan: An Educational Framework for Multimessenger Observing Scenarios
===============================================================================

.. image:: https://readthedocs.org/projects/EarthOrbitPlan/badge/?version=latest
   :target: https://EarthOrbitPlan.readthedocs.io/en/latest/
   :alt: Documentation Status
.. image:: https://github.com/weizmannk/EarthOrbitPlan/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/weizmannk/EarthOrbitPlan/actions/workflows/ci.yml
   :alt: CI Status
.. image:: https://img.shields.io/badge/License-BSD_3--Clause-blue.svg
   :target: https://opensource.org/licenses/BSD-3-Clause
   :alt: License: BSD-3-Clause
.. image:: https://img.shields.io/badge/python-3.12%2B-blue.svg
   :target: https://www.python.org/downloads/
   :alt: Python 3.12+
.. Replace XXXXXXX with the Zenodo record id once the release is archived.
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg
   :target: https://doi.org/10.5281/zenodo.XXXXXXX
   :alt: DOI

**EarthOrbitPlan** is an open-source educational framework designed to introduce students, educators, and early-career researchers to the challenges and opportunities of gravitational wave (GW) multi-messenger follow-up observations.

By simulating realistic observing scenarios, from the initial detection of a GW event to the coordinated search for its electromagnetic counterpart. EarthOrbitPlan offers a hands-on approach to learning about the coordination and optimization required in modern astronomy.

Install
-------

EarthOrbitPlan requires **Python 3.12+** (inherited from `M⁴OPT <https://github.com/m4opt/m4opt>`_).

.. code-block:: bash

   python3.12 -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install earthorbitplan            # from PyPI

   # or, for a development checkout:
   git clone https://github.com/weizmannk/EarthOrbitPlan.git
   cd EarthOrbitPlan
   pip install -e ".[dev,docs]"

The M⁴OPT scheduler needs the CPLEX solver; see the
`M⁴OPT CPLEX guide <https://m4opt.readthedocs.io/en/latest/install/cplex.html>`_.

Quickstart
----------

A self-contained taste of the statistics layer (no solver or data files
required): the 5th/50th/95th percentile event counts for a Poisson process
with a log-normal prior on its rate.

.. code-block:: python

   import numpy as np
   from earthorbitplan.probability.rate import poisson_lognormal_rate_quantiles

   quantiles = poisson_lognormal_rate_quantiles(np.array([0.05, 0.5, 0.95]), mu=2.0, sigma=0.5)
   print(quantiles)  # -> [ 1.415...  6.826... 17.859... ]

For the full workflow, from a GW skymap to a scheduled UVEX/ULTRASAT
follow-up plan and its detection statistics, follow the
`M⁴OPT scheduler walkthrough <https://EarthOrbitPlan.readthedocs.io/en/latest/m4opt-scheduler/index.html>`_
in the documentation.

Key Features
------------

- **Educational Focus:** Tutorials, example notebooks, and documentation tailored for teaching multi-messenger astronomy, observation scheduling, and data analysis.
- **Realistic Simulations:** Recreates GW follow-up campaigns for various event types (BNS, NSBH, BBH) using state-of-the-art models.
- **Observatory Planning Tools:** Leverages `M⁴OPT <https://github.com/m4opt/m4opt>`_ to calculate telescope Field of Regard (FoR), optimal exposure times, and the most promising sky regions for follow-up.
- **Automation:** Streamlines workflows for event processing, data filtering, statistical analysis, and visualization, making it accessible for classroom or workshop use.
- **Impact Analysis:** Enables users to compare the performance of different observatories, assess the effect of scheduling strategies, and visualize outcomes.

EarthOrbitPlan demonstrates how to use `M⁴OPT <https://github.com/m4opt/m4opt>`_ to determine the optimal sky regions accessible to each telescope (Field of Regard, FoR), calculate the required exposure time to achieve a target signal-to-noise ratio, and identify the most probable sky coordinates for efficient follow-up observations.

Who is it for?
--------------

- **Students** who want to explore multi-messenger astrophysics through practical examples.
- **Educators** seeking ready-to-use tools and resources for teaching astronomy, physics, or data science.
- **Researchers** wishing to prototype observation strategies or demonstrate multi-messenger concepts.

Documentation and project files
-------------------------------

- Documentation: https://EarthOrbitPlan.readthedocs.io
- Contributing guide: `CONTRIBUTING.md <CONTRIBUTING.md>`_
- Changelog: `CHANGELOG.md <CHANGELOG.md>`_

How to cite
-----------

If you use EarthOrbitPlan in your research or teaching, please cite it using the
metadata in `CITATION.cff <CITATION.cff>`_ (GitHub's *Cite this repository*
button generates APA/BibTeX from it). The archived-release DOI will be added
here once available.
