.. _observing_scenarios:

Observing Scenario Simulation
=============================

Introduction
~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The so-called GWTC-3 distribution  provides an empirical description of the compact object population based on the GWTC-3 catalog data.
This modeling has the following features:

- **Components**: All CBC systems (BNS, NSBH, BBH) are described by a continuous mass distribution, without prior separation.
- **Distribution shape**: The mass function follows a broken power law, with the possibility of a “dip” in the intermediate region to model a potential mass gap between neutron stars and black holes.
- **Tapering**: Cutoff functions are applied at low and high masses to reproduce the observed behavior at the ends of the distribution.
- **Pairing function**: Primary and secondary masses for each system are drawn from the distribution and paired via a law depending on the mass ratio, favoring realistic binary formation.
- **Spin distribution**: Spin magnitudes are drawn uniformly, and their directions are assumed isotropic. Ranges differ depending on the mass of the objects (see :footcite:`2016A&A...594A..13P` for details).

Mass Distribution and the Mass Gap
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The figure below shows the 1D Power Law + Dip + Break (PDB) mass distribution...
The 1D PDB mass distribution :math:`p(m|\lambda)` in the range  :math:`[1, 100]\,\, {M}_\odot`,
for a representative set of hyperparameters :math:`\lambda`.(See the  Table 7 in Appendix A.1 of :footcite:`2023ApJ...958..158K` for the full parameter values).

.. plot::
   :caption: |
      Population model for the primary mass distribution.
      The blue shaded region and dashed lines highlight the "mass gap" between neutron stars and black holes.
   :include-source: False
   :scale: 0.43

   import numpy as np
   import matplotlib.pyplot as plt

   # Model parameters
   ALPHA_1 = -2.16
   ALPHA_2 = -1.46
   A = 0.97
   M_GAP_LO = 2.72
   M_GAP_HI = 6.13
   ETA_GAP_LO = 50
   ETA_GAP_HI = 50
   ETA_MIN = 50
   ETA_MAX = 4.91
   BETA = 1.89
   M_MIN = 1.16
   M_MAX = 54.38

   def lopass(m, m_0, eta):
       return 1 / (1 + (m / m_0) ** eta)

   def hipass(m, m_0, eta):
       return 1 - lopass(m, m_0, eta)

   def bandpass(m, m_lo, m_hi, eta_lo, eta_hi, A):
       return 1 - A * hipass(m, m_lo, eta_lo) * lopass(m, m_hi, eta_hi)

   def mass_distribution_1d(m):
       return (
           bandpass(m, M_GAP_LO, M_GAP_HI, ETA_GAP_LO, ETA_GAP_HI, A)
           * hipass(m, M_MIN, ETA_MIN)
           * lopass(m, M_MAX, ETA_MAX)
           * (m / M_GAP_HI) ** np.where(m < M_GAP_HI, ALPHA_1, ALPHA_2)
       )

   m = np.geomspace(1, 100, 2000)
   fig, ax = plt.subplots()
   ax.set_xscale("log")
   ax.set_yscale("log")

   # Violet: '#9400D3', Navy: '#001F75'
   ax.plot(m, m * mass_distribution_1d(m), color='navy', linewidth=2, label='Mass distribution')

   ax.set_xlim(1, 100)
   ax.set_ylim(0, 99.99)
   ax.set_xlabel(r"mass, $m$ [$M_\odot$]")
   ax.set_ylabel(r"$m\,p(m|\lambda)$")

   # Mass gap region in light blue, dashed lines
   y_min, y_max = ax.get_ylim()
   ax.fill_between([1.9, 2.9], y_min, y_max, color="blue", alpha=0.12, zorder=1)
   ax.axvline(x=1.9, color="blue", linestyle="--", alpha=0.7)
   ax.axvline(x=2.9, color="blue", linestyle="--", alpha=0.7)
   ax.text(
       2.4, y_min * 1.5, r"$2.4^{+0.5}_{-0.5}$",
       ha="center", va="bottom", fontsize=11, fontweight="bold", color="blue"
   )

   ax2 = ax.twiny()
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xscale(ax.get_xscale())
   ax2.set_xticks([M_MIN, M_GAP_LO, M_GAP_HI, M_MAX])
   ax2.set_xticklabels(
       [
           r"$M_\mathrm{min}$",
           r"$M^\mathrm{gap}_\mathrm{low}$",
           r"$M^\mathrm{gap}_\mathrm{high}$",
           r"$M_\mathrm{max}$",
       ]
   )
   ax2.grid(axis="x")
   fig.tight_layout()
   plt.show()
   plt.close()

The model is based on :footcite:`2022ApJ...931..108F`, applied to the GWTC-3 distribution :footcite:`2023PhRvX..13a1048A`, and implemented
in our simulations as described in :footcite:`2023ApJ...958..158K`.


.. plot::
   :caption: |
        Gaussian kernel density estimator analysis of the PDB/GWTC-3 distribution, showing comparative mass and spin distributions across CBC categories.
      **Left:** Logarithmic 2D distribution of primary vs. secondary masses for the first 1,000 PDB/GWTC-3 CBC events,
       based on Gaussian kernel density estimation.
      **Right:** Spin distribution of the same events, showing component spin correlations.
      Color scale indicates the event density per pixel.
   :include-source: False
   :scale: 0.43

   import os
   from astropy.table import Table
   import numpy as np
   from scipy.stats import gaussian_kde
   import matplotlib.pyplot as plt

   data_dir = '../../scenarios/farah.h5'
   outdir = '.'

   Farah = Table.read(data_dir)[:1000]
   Farah.sort('mass1')
   ns_max_mass = 3.0

   plt.clf()
   fig, axs = plt.subplots(1, 2, figsize=(14, 6))

   # Mass distribution (log scale)
   mass1 = np.log10(Farah['mass1'])
   mass2 = np.log10(Farah['mass2'])
   xy_mass = np.vstack([mass1, mass2])
   z_mass = gaussian_kde(xy_mass)(xy_mass)
   idx_mass = z_mass.argsort()
   mass1, mass2, z_mass = mass1[idx_mass], mass2[idx_mass], z_mass[idx_mass]
   msc = axs[0].scatter(mass1, mass2, c=z_mass, s=5)
   axs[0].set_xlabel(r'$\log_{10}(m_1)\ [M_\odot]$')
   axs[0].set_ylabel(r'$\log_{10}(m_2)\ [M_\odot]$')
   fig.colorbar(msc, ax=axs[0], label="Event density")

   # Spin distribution
   spin1z = Farah['spin1z']
   spin2z = Farah['spin2z']
   xy_spin = np.vstack([spin1z, spin2z])
   z_spin = gaussian_kde(xy_spin)(xy_spin)
   idx_spin = z_spin.argsort()
   spin1z, spin2z, z_spin = spin1z[idx_spin], spin2z[idx_spin], z_spin[idx_spin]
   ssc = axs[1].scatter(spin1z, spin2z, c=z_spin, s=5)
   axs[1].set_xlabel(r'$\mathrm{spin}_1$')
   axs[1].set_ylabel(r'$\mathrm{spin}_2$')
   fig.colorbar(ssc, ax=axs[1], label="Event density")

   plt.tight_layout()
   plt.show()
   plt.close()


.. note::

   This example uses only the first 1,000 events from the PDB/GWTC-3 catalog for clarity and fast plotting.
   For a full population analysis, you may increase this number (e.g., up to one million events), but this will require more time and memory.
   The documentation build does not run the full sample for efficiency—re-run locally for high-statistics plots.

.. warning::

   Using large samples (100,000+ events) may require significant computing resources.


Simulation process
~~~~~~~~~~~~~~~~~~

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

.. toctree::
   :maxdepth: 2

   ../notebooks/compute_GW_detection_rate


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
