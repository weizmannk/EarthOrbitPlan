.. _observing_scenarios:

===================
Observing Scenarios
===================

.. admonition:: Overview
   :class: tip

   Observing scenarios are used to generate quantitative forecasts of gravitational wave (gravitational-waves) network performance by simulating the detection
   and localization of gravitational-waves events during key science runs :footcite:`2023ApJ...958..158K,2014ApJ...795..105S,2018LRR....21....3A` (e.g., O4, O5).
   These simulations inform the optimization of observational strategies and the development of future instrumentation, based on population synthesis models and published detector sensitivity data.


.. dropdown:: Why Observing Scenarios Matter?

    .. admonition:: Main Purposes
       :class: important

        - Quantitative estimation of detection rates and sky localization accuracy for various detector network configurations and observing runs.

        - Optimization of electromagnetic follow-up strategies.

        - Assessment of requirements and trade-offs for future instrumentation and network design.

    Recent studies :footcite:`2022ApJ...924...54P,2023ApJ...958..158K` have calibrated these simulations using actual public alerts :footcite:`2023PhRvX..13a1048A`,
    thereby improving the reliability of forecasts through realistic Signal-to-noise ratio (:term:`SNR`) thresholds and advanced search methodologies.
    These developments have enabled population studies, :math:`r`-process nucleosynthesis analysis, and cosmological parameter inference :footcite:`kiendrebeogo:tel-04796327,2013ApJ...767..124N, 2023ApJ...958..158K`.


.. dropdown:: Current GW Network and Run Status

    .. admonition:: Current Network Status
       :class: info

       The O4 observing run began with both LIGO Hanford (:term:`LHO`) and LIGO Livingston (:term:`LLO`) in operation, achieving a :term:`BNS` range of 140–165 Mpc.
       After a commissioning break, Virgo rejoined the network in March 2024 with a BNS range of 55 Mpc, followed later by KAGRA. O4 is scheduled to continue until October 7, 2025,
       marking the first period with all four detectors operating together. This configuration significantly enhances both detection rates and localization accuracy for GW events.
       For up-to-date information on detector status and sensitivity, see the `real-time detector status and range page <https://online.ligo.org>`_.


.. dropdown:: Population Modeling

    All simulated populations are generated using the GWTC-3 "Power Law + Dip + Break" (:term:`PDB`) mass distribution :footcite:`2022ApJ...931..108F,2023PhRvX..13a1048A`, which empirically describes the properties of compact binaries detected by the LIGO-Virgo-KAGRA network.

    **Key characteristics:**

    - **Unified mass function:** Describes all CBC systems (:term:`BNS`, :term:`NSBH`, :term:`BBH`) continuously, without explicit subclass boundaries.

    - **Broken power law with dip:** Incorporates a broken power law and a "dip" to represent the observed mass gap between neutron stars and black holes.

    - **Tapering:** Cutoff functions at low and high masses reproduce observed behavior at distribution edges.

    - **Pairing law:** Component masses (Primary and secondary masses) are paired according to a mass-ratio-dependent prescription, producing physically plausible binaries.

    - **Spin properties:** Spin magnitudes are drawn from uniform distributions with isotropic orientations, and mass-dependent maximum values (see :footcite:`2016A&A...594A..13P`).

    This modeling framework is used to generate realistic populations of component masses, spins, and sky locations for all scenario simulations.



.. dropdown:: Mass Distribution and the Mass Gap

    The figure below shows the 1D Power Law + Dip + Break (PDB) mass distribution...
    The 1D PDB mass distribution :math:`p(m|\lambda)` in the range  :math:`[1, 100]\,\, {M}_\odot`.

    .. plot::
        :caption: Population model for the primary mass distribution.
                    The blue shaded region and dashed lines highlight the "mass gap" between neutron stars and black holes.
        :include-source: False

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
        ax.set_ylim(0, 100)
        ax.set_xlabel(r"mass, $m$ [$M_\odot$]")
        ax.set_ylabel(r"$m\,p(m|\lambda)$")

        # Mass gap region in light blue, dashed lines
        y_min, y_max = ax.get_ylim()
        ax.fill_between([1.9, 2.9], y_min, y_max, color="blue", alpha=0.12, zorder=1)
        ax.axvline(x=1.9, color="blue", linestyle="--", alpha=0.7)
        ax.axvline(x=2.9, color="blue", linestyle="--", alpha=0.7)
        ax.text(
            2.4, y_min-0.01, r"$2.4^{+0.5}_{-0.5}$",
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
        fig.show()


    The model is based on :footcite:`2022ApJ...931..108F`, applied to the GWTC-3 distribution :footcite:`2023PhRvX..13a1048A`, and implemented
    in our simulations as described in :footcite:`2023ApJ...958..158K`. The following table summarizes the full set of hyperparameters :math:`\lambda`.
    The first several entries describe the rate and mass distribution parameters,  and the last two entries describe the spin distribution parameters.


    .. tab-set::

        .. tab-item:: Hyperparameters

            .. table::  Hyperparameters used in the population model

                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | Parameter                            | Description                                                                       | Value                 |
                +======================================+===================================================================================+=======================+
                | :math:`\alpha_1`                     | Spectral index for the power law of the mass distribution at low mass             | -2.16                 |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\alpha_2`                     | Spectral index for the power law of the mass distribution at high mass            | -1.46                 |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\mathrm{A}`                   | Lower mass gap depth                                                              | 0.97                  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`M^\mathrm{gap}_\mathrm{low}`  | Location of lower end of the mass gap                                             | 2.72 :math:`M_\odot`  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`M^\mathrm{gap}_\mathrm{high}` | Location of upper end of the mass gap                                             | 6.13 :math:`M_\odot`  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\eta_\mathrm{low}`            | Parameter controlling how the rate tapers at the low end of the mass gap          | 50                    |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\eta_\mathrm{high}`           | Parameter controlling how the rate tapers at the low end of the mass gap          | 50                    |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\eta_\mathrm{min}`            | Parameter controlling tapering the power law at low mass                          | 50                    |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\eta_\mathrm{max}`            | Parameter controlling tapering the power law at high mass                         | 4.91                  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`\beta`                        | Spectral index for the power law-in-mass-ratio pairing function                   | 1.89                  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`M_{\rm min}`                  | Onset location of low-mass tapering                                               | 1.16 :math:`M_\odot`  |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`M_{\rm max}`                  | Onset location of high-mass tapering                                              | 54.38 :math:`M_\odot` |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`a_{\mathrm{max, NS}}`         | Maximum allowed component spin for objects with mass :math:`< 2.5\, M_\odot`      | 0.4                   |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+
                | :math:`a_{\mathrm{max, BH}}`         | Maximum allowed component spin for objects with mass :math:`\geq 2.5\, M_\odot`   | 1                     |
                +--------------------------------------+-----------------------------------------------------------------------------------+-----------------------+


            See :footcite:`2022ApJ...931..108F,2023ApJ...958..158K` for details, and :doc:`Observing Capabilities <userguide:capabilities>` for practical applications of the PDB distribution in network simulations.


        .. tab-item:: Gaussian kernel density estimaton

            .. plot::
                :caption: Gaussian kernel density estimator analysis of the PDB/GWTC-3 distribution, showing comparative mass and spin distributions across CBC categories.
                            **Left:** Logarithmic 2D distribution of primary vs. secondary masses for the first 10,000 PDB/GWTC-3 CBC events, based on Gaussian kernel density estimation.
                            **Right:** Spin distribution of the same events, showing component spin correlations. Color scale indicates the event density per pixel.
                :include-source: False

                    import os
                    from astropy.table import Table
                    import numpy as np
                    from scipy.stats import gaussian_kde
                    import matplotlib.pyplot as plt

                    # Load and process data
                    data_dir = '../../earthorbitplan/scenarios/farah.h5'
                    Farah = Table.read(data_dir)[:10000]
                    Farah.sort('mass1')

                    # Create subplots
                    fig, axs = plt.subplots(1, 2)

                    # increase the font size of the axes
                    for ax in axs:
                        for tick in ax.get_xticklabels() + ax.get_yticklabels():
                            tick.set_fontname("Times New Roman")
                            tick.set_fontsize(14)

                    # Mass distribution (log scale)
                    mass1 = np.log10(Farah['mass1'])
                    mass2 = np.log10(Farah['mass2'])
                    xy_mass = np.vstack([mass1, mass2])
                    z_mass = gaussian_kde(xy_mass)(xy_mass)
                    idx_mass = z_mass.argsort()
                    mass1, mass2, z_mass = mass1[idx_mass], mass2[idx_mass], z_mass[idx_mass]
                    msc = axs[0].scatter(mass1, mass2, c=z_mass, s=5)
                    axs[0].set_xlabel(r'$\log_{10}(m_1)\ [M_\odot]$',  fontname="Times New Roman", size=16, fontweight="bold")
                    axs[0].set_ylabel(r'$\log_{10}(m_2)\ [M_\odot]$', fontname="Times New Roman", size=16, fontweight="bold")
                    cbar1 = fig.colorbar(msc, ax=axs[0])
                    cbar1.set_label("Event density", fontname="Times New Roman", size=18)


                    # Spin distribution
                    spin1z = Farah['spin1z']
                    spin2z = Farah['spin2z']
                    xy_spin = np.vstack([spin1z, spin2z])
                    z_spin = gaussian_kde(xy_spin)(xy_spin)
                    idx_spin = z_spin.argsort()
                    spin1z, spin2z, z_spin = spin1z[idx_spin], spin2z[idx_spin], z_spin[idx_spin]
                    ssc = axs[1].scatter(spin1z, spin2z, c=z_spin, s=5)
                    axs[1].set_xlabel(r'$\mathrm{spin}_1$', fontname="Times New Roman", size=16, fontweight="bold")
                    axs[1].set_ylabel(r'$\mathrm{spin}_2$', fontname="Times New Roman", size=16, fontweight="bold")
                    cbar2 = fig.colorbar(ssc, ax=axs[1])
                    cbar2.set_label("Event density", fontname="Times New Roman", size=18)

                    # Adjust layout and figure size
                    fig.set_size_inches(14, 6)
                    plt.tight_layout()
                    plt.show()

            .. note::
                This example uses only the first 10,000 events from the PDB/GWTC-3 catalog for clarity and fast plotting.
                For a full population analysis, you may increase this number (e.g., up to one million events), but this will require more time and memory.
                The documentation build does not run the full sample for efficiency—re-run locally for high-statistics plots.

            .. warning::

                Using large samples (100,000+ events) may require significant computing resources.


.. dropdown:: Simulation Pipeline

   Our workflow consists of:

   1. **Population sampling**: Draw binaries from the PDB distribution, including mass, spin, orientation, and location.
   2. **Detection simulation**: Apply :term:`SNR` thresholds using each network’s published sensitivity curves and duty cycles.
   3. **Localization**: Use `ligo.skymap <https://lscsoft.docs.ligo.org/ligo.skymap>`_  tools to estimate sky position and distance for detected events.
   4. **Scenario preparation**: Characterize each event for electromagnetic` follow-up planning.

    .. seealso::

        Results are continually updated, see :doc:`Observing Capabilities <userguide:capabilities>` for the latest.


.. dropdown:: Data location on Zenodo

    Our simulations explore multiple detector configurations and :term: `SNR` thresholds to estimate GW detection rates under realistic observing conditions:

    .. tab-set::

        .. tab-item:: SNR threshold of 8

            - The **HL configuration**, deployed during the O4a observing run, the simulations data are available at `HL-config <https://doi.org/10.5281/zenodo.10078926>`_.
            - The **HLVK configuration**, planned for O4 and O5,  results are in `HLVK-config <https://doi.org/10.5281/zenodo.7026209>`_.
            - We also simulate **HLV and HV configurations** for O5 to assess the effect of detector configurations, including scenarios where only one LIGO detector is operating, the simulation data are located in `zenodo <https://zenodo.org/records/15617982>`_.

        .. tab-item:: SNR threshold of 10

            - The **HLVK configuration** for the upcoming O5 and O6 runs is used to estimate detection rates based on a more conservative detectability threshold, reflecting planned pipeline improvements and noise rejection strategies.


.. dropdown:: Notebooks & Resources

    .. tab-set::

        .. tab-item:: Population sampling

            A large number of binary systems (e.g., :math:`10^6`) are generated by drawing their masses and spins according to the PDB distribution described above.
            Orientation parameters and comoving volume positions are also drawn uniformly and isotropically.

        .. tab-item:: Gravitational-wave detection simulation

            - The generated signals are subjected to a detectability threshold based on the :term:`SNR` for each detector network, corresponding to the O4 or O5 configurations.
            - Instrumental noise is simulated using the published `sensitivity curves (PSD) for each detector <https://dcc.ligo.org/T2200043-v3/public>`_.
            - Detector duty cycles are realistically accounted for.

        .. tab-item:: Source localization

            - Events passing the SNR threshold are localized on the sky using the ``ligo.skymap`` toolchain
              (e.g., `bayestar-localize-coincs <https://lscsoft.docs.ligo.org/ligo.skymap/tool/bayestar_localize_coincs.html#offline-localization-bayestar-localize-coincs>`_),
              producing a sky probability map and distance estimate for each event.
            - Credible regions (e.g., 90%) and the comoving distance distribution are extracted for each simulated event.


    .. card:: Run the full analysis in interactive Jupyter notebooks

        You can explore the tutorials interactively

        ^^^
        - `Open in Binder <https://mybinder.org/v2/gh/weizmannk/EarthOrbitPlan/HEAD?urlpath=lab/tree/earthorbitplan/tutorials/observing_scenarios.ipynb>`__

        .. image:: https://mybinder.org/badge_logo.svg
            :target: https://mybinder.org/v2/gh/weizmannk/EarthOrbitPlan/HEAD?urlpath=lab/tree/earthorbitplan/tutorials/observing_scenarios.ipynb
        +++


.. dropdown:: Tools and Resources

    - The simulation pipeline primarily relies on the `ligo.skymap <https://lscsoft.docs.ligo.org/ligo.skymap>`_ software suite.
    - The scripts used to reproduce the entire population generation and simulation process are publicly available on GitHub (cf. https://github.com/lpsinger/observing-scenarios-simulations).
    - Sensitivity curves and other configuration parameters are drawn from official :term:`IGWN` consortium publications.


   .. note::
      This section covers only the simulation methodology.
      For results and quantitative comparisons, see :footcite:`2023ApJ...958..158K`.

.. dropdown:: Working With Zenodo Data

   .. note::
      This example demonstrates unpacking, filtering, and converting GW injection datasets (e.g., GWTC-3) from Zenodo archives.
      Outputs include ECSV tables and organized FITS files for O5/O6 runs.




References
==========

.. footbibliography::
