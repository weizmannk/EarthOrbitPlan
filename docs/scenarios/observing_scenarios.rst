.. _observing_scenarios:

GWTC-3 Distribution and Observing Scenario Simulation Pipeline
=============================================================

Introduction
============

A robust understanding of the sensitivity, detection efficiency, and localization capabilities of the global gravitational wave detector network
is essential for optimizing observational strategies and guiding the development of future telescopes and missions.
To this end, **observing scenarios**, which simulate the detection and localization of gravitational wave (GW) events,
provide realistic forecasts of network performance during key science runs, notably O4 and O5 :footcite:`2023ApJ...958..158K,2014ApJ...795..105S,2018LRR....21....3A`.

Recent scenario studies :footcite:`2022ApJ...924...54P` have been carefully calibrated using public alerts from O3, improving the accuracy of localization by incorporating
realistic signal-to-noise ratio (SNR) thresholds and single-detector search strategies.
These efforts have advanced our ability to explore compact object populations, r-process nucleosynthesis, and cosmological measurements :footcite:`kiend2024,2013ApJ...767..124N,2017ApJ...848L..12A`.

The O4 observing run began with both LIGO Hanford (LHO) and LIGO Livingston (LLO) in operation, achieving a binary neutron star (BNS) range of 140–165 Mpc.
After a commissioning break, Virgo rejoined the network in March 2024 with a BNS range of 55 Mpc, followed later by KAGRA. O4 is scheduled to continue until October 7, 2025,
marking the first period with all four detectors operating together. This will significantly enhance both detection rates and localization accuracy for GW events.
For up-to-date information on detector status and sensitivity, see the `real-time detector status and range page <https://online.ligo.org>`_.

To model these capabilities, we simulate realistic astrophysical distributions of mass, spin, and sky locations for compact binary coalescences (CBCs).
The GWTC-3 distribution :footcite:`2022ApJ...931..108F,2023PhRvX..13a1048A` (Power Law + Dip + Break, PDB) serves as the foundation for population generation
in these simulations :footcite:`kiend2024,2023ApJ...958..158K`.

This section details the procedure for generating CBC populations from the GWTC-3 distribution and describes the simulation pipeline used to
reproduce observing scenarios for the O4, O5 and O6 campaigns of the LIGO-Virgo-KAGRA (LVK) network.


Population Modeling: GWTC-3 (PDB) Distribution
==============================================

The so-called GWTC-3 distribution  provides an empirical description of the compact object population based on the GWTC-3 catalog data.
This modeling has the following features:

- **Components**: All CBC systems (BNS, NSBH, BBH) are described by a continuous mass distribution, without prior separation.
- **Distribution shape**: The mass function follows a broken power law, with the possibility of a “dip” in the intermediate region to model a potential mass gap between neutron stars and black holes.
- **Tapering**: Cutoff functions are applied at low and high masses to reproduce the observed behavior at the ends of the distribution.
- **Pairing function**: Primary and secondary masses for each system are drawn from the distribution and paired via a law depending on the mass ratio, favoring realistic binary formation.
- **Spin distribution**: Spin magnitudes are drawn uniformly, and their directions are assumed isotropic. Ranges differ depending on the mass of the objects (see :footcite:`2016A&A...594A..13P` for details).

GWTC-3 Mass  and Spin distribution
==================================

The figure below shows the 1D Power Law + Dip + Break (PDB) mass distribution...

.. figure:: ../_static/supress_mass_gap.png
   :scale: 70 %
   :align: center
   :alt: PDB mass distribution
   :figclass: align-center

   The 1D PDB mass distribution :math:`p(m|\lambda)` in the range  :math:`[1, 100]\,\, {M}_\odot`, for a representative set of hyperparameters :math:`\lambda`.
   (See the  Table 7 in Appendix A.1 of :footcite:`2023ApJ...958..158K` for the full parameter values).

The model is based on :footcite:`2022ApJ...931..108F`, applied to the GWTC-3 distribution :footcite:`2023PhRvX..13a1048A`, and implemented
in our simulations as described in :footcite:`2023ApJ...958..158K`.

.. figure:: ../_static/Farah_gaussian_kde_distribution_of_masses_2.png
   :scale: 70 %
   :align: center
   :alt: Gaussian KDE Analysis of PDB/GWTC-3 distribution
   :figclass: align-center

   Gaussian kernel density estimator analysis of the PDB/GWTC-3 distribution, showing comparative mass and spin distributions across CBC categories.
   The left panel displays the 2D mass distributions for the components of each CBC category (all axes logarithmic).
   The right panel shows the spin distributions for the same CBC categories, providing insights into the spin characteristics of the population.
   The color density in both panels represents the number of CBC events per pixel, illustrating the density and variation within the distributions.

   Adapted from :footcite:`kiend2024,2023ApJ...958..158K`

PDB/GWTC-3 Masses and Spin Distribution
---------------------------------------

.. figure:: ../_static/Farah_gaussian_kde_distribution_of_masses_2.png
   :scale: 70 %
   :align: center
   :alt: Gaussian KDE analysis of PDB/GWTC-3 distribution

   Gaussian KDE analysis of the PDB/GWTC-3 distribution, showing comparative mass and spin distributions across CBC categories.
   The left panel displays the 2D mass distributions for the components of each CBC category (all axes logarithmic).
   The right panel shows the spin distributions for the same CBC categories, providing insights into the spin characteristics of the population.
   The color density in both panels represents the number of CBC events per pixel, illustrating the density and variation within the distributions.
   Reproduced from :footcite:`kiend2024,2023ApJ...958..158K`.

.. note::

   The figure above is based on a sample of **one million GWTC-3 events**.
   For efficiency, the full plot is **not generated during the documentation build**.
   You can reproduce it locally using the script below.
   The example uses only the first 1,000 events, but you may increase this number for a more detailed distribution (note: higher computation time).


.. plot::
   :caption: Mass distribution of the first 1000 events in PDB/GWTC-3 (demo)
   :include-source: False

   import os
   from astropy.table import Table
   import numpy as np
   from scipy.stats import gaussian_kde
   import matplotlib.pyplot as plt

    # Data path (edit if needed)
   data_dir = '../../scenarios/farah.h5'

   np.random.seed(42) # No need if running the all events
   Farah = Table.read(f"{data_dir}")[:1000] # Only first 1000 events

   Farah.sort('mass1')

   params = ['mass', 'spin']
   ns_max_mass = 3.0

   plt.clf()
   fig, axs = plt.subplots(nrows=1, ncols=2)

   for param in params:
       if param == 'mass':
           mass1 = np.log10(Farah['mass1'])
           mass2 = np.log10(Farah['mass2'])
           xy = np.vstack([mass1, mass2])
           z = gaussian_kde(xy)(xy)
           index = z.argsort()
           mass1, mass2, z = mass1[index], mass2[index], z[index]
           mass_density_Farah = axs[0].scatter(mass1, mass2, c=z, s=5)
           axs[0].set_xlabel(r'Primary mass, $ \log_{10}(m_1)$ [$M_\odot$]')
           axs[0].set_ylabel(r'Secondary mass, $\log_{10}(m_2)$ [$M_\odot$]')
       else:
           spin1z = Farah['spin1z']
           spin2z = Farah['spin2z']
           xy = np.vstack([spin1z, spin2z])
           z = gaussian_kde(xy)(xy)
           index = z.argsort()
           spin1z, spin2z, z = spin1z[index], spin2z[index], z[index]
           spin_density_Farah = axs[1].scatter(spin1z, spin2z, c=z, s=5)
           axs[1].set_xlabel(r'$\mathrm{spin}_1$')
           axs[1].set_ylabel(r'$\mathrm{spin}_2$')

   cbar1 = fig.colorbar(mass_density_Farah, ax=axs[0])
   cbar2 = fig.colorbar(spin_density_Farah, ax=axs[1])

   fig.text(0.5, 0.03, "PDB/GWTC-3 distribution", ha="center", va="center", fontsize=20, color='navy')
   plt.gcf().set_size_inches(14, 6)
   plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
   fig.tight_layout(rect=[0, 0.05, 1, 1])
   plt.savefig(f'{outdir}/Farah_gaussian_kde_distribution_of_masses_2.pdf')
   plt.savefig(f'{outdir}/Farah_gaussian_kde_distribution_of_masses_2.png')
   plt.close()

   # Print population counts
   print(f"number of BNS: {len(Farah[(Farah['mass1']<ns_max_mass)])}, "
         f"number of NSBH: {len(Farah[(Farah['mass1']>=ns_max_mass) & (Farah['mass2']<ns_max_mass)])}, "
         f"number of BBH: {len(Farah[(Farah['mass2']>=ns_max_mass)])}")

.. warning::

   For very large event samples, running the script above may require significant computing time and memory.



PDB/GWTC-3 Mass and Spin Distributions
======================================

.. figure:: ../_static/Farah_gaussian_kde_distribution_of_masses_2.png
   :scale: 70 %
   :align: center
   :alt: Mass and spin distributions in the PDB/GWTC-3 catalog

   **Gaussian KDE analysis** of the PDB/GWTC-3 catalog, comparing mass and spin distributions for different compact binary coalescence (CBC) categories.

   - **Left panel**: 2D distribution of primary and secondary masses (logarithmic scale).
   - **Right panel**: 2D spin distribution for the same events.

   Color density indicates the number of CBC events per pixel, highlighting the density and variation of the distributions.



Observing Scenario Simulation
===============================

Our simulations explore multiple detector configurations and signal-to-noise (SNR) thresholds to estimate GW detection rates under realistic observing conditions:

a. **SNR threshold of 8**

- The **HL configuration**, deployed during the O4a observing run, the simulations data are available at `Zenodo <https://doi.org/10.5281/zenodo.10078926>`_.
- The **HLVK configuration**, planned for O4 and O5,  results are in `Zenodo <https://doi.org/10.5281/zenodo.7026209>`_.
- We also simulate **HLV and HV configurations** for O5 to assess the effect of detector configurations, including scenarios where only one LIGO detector is operating, the simulation data are located in `zenodo <https://zenodo.org/records/15617982>`_.

b. **SNR threshold of 10**

- The **HLVK configuration** for the upcoming O5 and O6 runs is used to estimate detection rates based on a more conservative detectability threshold, reflecting planned pipeline improvements and noise rejection strategies.

The simulation pipeline follows these steps:

1. **Population sampling**
A large number of binary systems (e.g., :math:`10^6`) are generated by drawing their masses and spins according to the PDB distribution described above. Orientation parameters and comoving volume positions are also drawn uniformly and isotropically.

2. **Gravitational-wave detection simulation**:

- The generated signals are subject to a detectability threshold based on the signal-to-noise ratio (SNR) for each detector network, corresponding to the O4 or O5 configurations.
- Instrumental noise is simulated from the published `sensitivity curves (PSD) for each detector <https://dcc.ligo.org/T2200043-v3/public>`_.
- Detector duty cycles are realistically accounted for.

3. **Source localization**

- Events passing the SNR threshold are localized on the sky using the ``ligo.skymap`` toolchain (e.g., `bayestar-localize-coincs <https://lscsoft.docs.ligo.org/ligo.skymap/tool/bayestar_localize_coincs.html#offline-localization-bayestar-localize-coincs>`_), producing a sky probability map and distance estimate for each event.
- Credible regions (e.g., 90%) and the comoving distance distribution are extracted for each simulated event.

4. **Observing scenario preparation**

- The properties of the simulated events (localization, distance, etc.) serve as the basis for defining various electromagnetic (EM) observation scenarios, according to the capabilities of the planned follow-up instruments.
- This pipeline allows evaluation, for each instrumental configuration, of the probability of covering the EM counterpart of a given GW event.

The results of these simulations are used to update the :doc:`Observing Capabilities <userguide:capabilities>`


Tools and Resources
====================

- The simulation pipeline primarily relies on the `ligo.skymap <https://lscsoft.docs.ligo.org/ligo.skymap>`_ software suite.
- The scripts used to reproduce the entire population generation and simulation process are publicly available on GitHub (cf. https://github.com/lpsinger/observing-scenarios-simulations).
- Sensitivity curves and other configuration parameters are drawn from official IGWN consortium publications.

.. note::

    This page only describes the methodology for population generation and the simulation pipeline. For results and quantitative analysis,
    srefer to the corresponding section :footcite:`2023ApJ...958..158K`.



=================================
Zenodo GW Injection Data Unpacker
=================================

.. note::
    Here we show  how to easily unpacking, filtering, and conversion of injection datasets
    (e.g., GWTC-3) from Zenodo ZIP archives. It processes event tables and associated
    localization files for specific observing runs (e.g., O5, O6), and outputs
    filtered ECSV tables and organized FITS files.

Usage
-----

You can run this unpacker from the command line as follows:

.. code-block:: bash

   python scenarios/zenodo_unpacker.py --zip runs_SNR-10.zip --subdir runs_SNR-10 --runs O5 O6 --detectors HLVK --outdir ./data --mass-threshold 3


Or use a config file:

.. code-block:: bash

    python scenarios/zenodo_unpacker.py --config params.ini


Or import and call `process_zip()` in your Python code.

Config file example (`params.ini`)
----------------------------------
[params]
zip = runs_SNR-10.zip
subdir = runs_SNR-10
runs = O5 O6
detectors = HLVK
data_dir = data
skymap_dir = skymaps
mass_threshold = 3.0

Source
------

Zenodo Dataset: https://zenodo.org/records/14585837

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



==========
References
==========

.. footbibliography::
