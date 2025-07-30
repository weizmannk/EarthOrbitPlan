.. _scheduler:

==========
Scheduling
==========


Scientific Rationale
====================

The emergence of all-sky surveys and the increasing number of real-time gravitational-waves alerts have created a pressing need for reliable scheduling frameworks.
These frameworks must efficiently coordinate follow-up observations across a range of telescopes and missions, both on the ground and in space.
To meet this need, we developed `|M⁴OPT| <https://m4opt.readthedocs.io/en/latest/>`_ : an open-source toolkit designed to optimize the scheduling of follow-up campaigns.

.. dropdown:: Why is scheduling so challenging?

  gravitational-waves alerts often have large localization uncertainties—sometimes hundreds of square degrees—so astronomers need to decide how to balance *sky coverage*
  (looking at as much area as possible) versus *depth* (spending enough time on each field to detect faint signals).
  Moreover, electromagnetic counterparts to gravitational-waves events can evolve rapidly, so observation plans must be generated quickly and efficiently.


  .. admonition:: Mixed Integer Linear Programming (MILP)
    :class: info

    :math:`\mathrm{M^4OPT}` addresses this by formulating the scheduling task as a **Mixed Integer Linear Programming (MILP)** problem.
    This approach dynamically allocates observation time across fields and optimizes the exposure time for each one, maximizing the overall
    probability of detection. Unlike fixed or manually tuned schedules, this method adapts to the conditions of each field, such as:

  - background noise and background light (including natural sky brightness from the galaxy, zodiacal light, and instrumental effects;
    for ground-based telescopes, this also includes the increased sky brightness during astronomical twilight—just before sunrise or after sunset),
  - distance uncertainty,
  - instrumental sensitivity.


.. important::

    In addition, :math:`\mathrm{M^4OPT}` takes into account practical telescope constraints, such as Sun and Moon exclusion zones, dynamic field-of-view, and slewing time between fields.

    This leads to a much more efficient use of telescope time and increases the chances of detecting faint or distant transients.


.. dropdown:: MILP Model and Constraints

    The MILP optimization considers several real-world constraints, including:

    1. Fixed latency between the alert and when observations can start;
    2. Field-of-Regard (FoR): the sky region accessible to the instrument, evolving over time due to the instrument’s motion and orientation;
    3. Slewing delays: time needed to move between different sky positions;
    4. Exposure time limits for each field (minimum and maximum).

    Detection probabilities are computed using synthetic photometry models, taking into account extinction, background light (zodiacal, galactic, etc.),
    instrument sensitivity, and the distribution of possible source magnitudes (for example, a Gaussian for kilonova absolute magnitudes).
    This modeling ensures realistic and robust predictions for observation planning.


Run process
===========

.. dropdown:: How to get your telescope observation plan ?

    Here is how you can run your first test with :math:`\mathrm{M^4OPT}` using the ULTRASAT mission.


    .. admonition:: With ULTRASAT telescope
        :class: info

        .. tab-set::

            .. tab-item:: 1. Running a Scheduling Simulation

                The simulation requires two inputs:
                ``./14.fits`` — the event skymap with localization probabilities
                ``./14.ecsv`` — the output file where the schedule will be saved

                Open your terminal and run the following command::

                    m4opt schedule ./data/14.fits ./data/14.ecsv \
                    --mission=ultrasat \
                    --skygrid=non-overlap \
                    --bandpass=NUV \
                    --absmag-mean=-16.0 \
                    --absmag-stdev=1.3 \
                    --exptime-min='300 s' \
                    --exptime-max='14400 s' \
                    --snr=10 \
                    --delay='15min' \
                    --deadline='24hour' \
                    --timelimit='2hour' \
                    --nside=128 \
                    --jobs 0 \
                    --cutoff=0.1

                This command launches a scheduling simulation for the ULTRASAT mission.
                You need to provide the main parameters, including the mission name, skygrid configuration, and observation settings.

                - The output file (e.g., ``14.ecsv``) will contain the observation schedule.
                - The simulation expects an event skymap file (usually a ``.fits`` file), which gives the localization probability of the event.

                .. note::
                    Missions like ULTRASAT support multiple skygrid models; use ``--skygrid`` to select (`non-overlap` and `allsky`).
                    Other missions (e.g., ZTF, UVEX, Rubin) support only a single skygrid and do not need this option.

                See the full list of parameters in the `CLI guide <https://m4opt.readthedocs.io/en/latest/guide/cli.html#m4opt-schedule>`_.



.. dropdown:: Output and Visualization

    .. admonition:: Understanding the Output
        :class: info

        The generated ECSV file (e.g. ``14.ecsv``) contains your observation plan, including:

        - Pointing coordinates,
        - Exposure times,
        - Slew (repositioning) times,
        - Visit (by default: two visits per field),
        - All relevant metadata.

        By default, the schedule includes **two visits per field**—so each coordinate may appear twice, corresponding to repeated observations.


    .. admonition:: Visualizing the Schedule
        :class: info

        .. tab-set::

            .. tab-item:: Visualizing the Schedule

                You can create an animation or a PDF showing the planned observations::

                    m4opt animate ./data/14.ecsv 14_MOVIE.gif --dpi 300 --still 14_MOVIE.pdf


            .. tab-item:: Animation

                The animation produces:

                - ``14_MOVIE.gif`` — an animation of the schedule
                - ``14_MOVIE.pdf`` — a static pdf,  of the observation sequence.

                .. image:: ../_static/14_MOVIE.gif
                    :alt: Example animation of the observation plan
                    :align: center

        .. tab-set::

            .. tab-item:: Explanation of the animation

                - The pink regions show the scheduled observation pointings the `footprints <https://m4opt.readthedocs.io/en/latest/api/m4opt.fov.footprint.html#footprint>`_.
                - The green outline marks the 90% credible region of the GW localization.
                - The deep blue areas are always outside the telescope’s Field of Regard; the light blue areas are temporarily out of view.
                - The lower panel shows how the detection probability and covered sky area accumulate over time, with different colors indicating
                    the number of times a region has been observed.
                - The symbol :math:`\oplus` shows the direction of the center of the Earth (sub-Earth point) projected onto the sky.
                - The symbol :math:`\odot` shows the direction of the Sun (sub-solar point) on the sky.

        .. seealso:: For more details of marker conventions

            For more details of marker conventions,
            see the `ligo.skymap plotting documentation <https://lscsoft.docs.ligo.org/ligo.skymap/plot/marker.html#module-ligo.skymap.plot.marker/>`_.


    .. note::

        This is a projection of the sky, **not a direct image of the Earth or the Moon**. The features shown correspond to sky coordinates,
        not to physical locations on Earth or lunar positions.


.. dropdown:: ECSV file inspection

    You can load and inspect a schedule file using Astropy:

    .. code-block:: python

        >>> from astropy.table import QTable
        >>> from earthorbitplan.utils.path import get_project_root
        >>> root = get_project_root()
        >>> output_file = root / "data" / "14.ecsv"
        >>> plan = QTable.read(output_file, format="ascii.ecsv")
        >>> obs = plan[plan["action"] == "observe"]
        >>> display = obs["start_time", "duration"]
        >>> display["ra"] = obs["target_coord"].ra
        >>> display["dec"] = obs["target_coord"].dec
        >>> display.round({'duration': 1, 'ra': 2, 'dec': 2})
        >>> print(display)


 .. dropdown:: ECSV Metadata Extraction

    Load a schedule, extract key metadata and visit counts:

    .. jupyter-execute::

        import sys, os
        sys.path.insert(0, os.path.abspath('../..'))
        from astropy.table import QTable
        from earthorbitplan.utils.path import get_project_root
        root = get_project_root()
        output_file = root / "data" / "14.ecsv"
        plan = QTable.read(output_file, format="ascii.ecsv")
        objective = plan.meta.get("objective_value")
        best_bound = plan.meta.get("best_bound")
        status = plan.meta.get("solution_status")
        time_used = plan.meta.get("solution_time")
        visits = plan.meta.get("args", {}).get("visits", 2)
        n_obs = len(plan[plan["action"] == "observe"])
        unique_fields = n_obs // visits
        print("Schedule metadata:")
        print(f" • Objective value: {objective:.4f}")
        print(f" • Best bound: {best_bound:.4f}")
        print(f" • Solver status: {status}")
        print(f" • Solution time: {time_used}")
        print(f" • Unique fields observed: {unique_fields}")


    .. list-table:: Schedule metadata summary
        :header-rows: 1
        :widths: 30 15

        * - Metric
          - Value
        * - Objective value
          - 0.9483
        * - Best bound
          - 0.9483
        * - Solver status
          - integer optimal solution
        * - Solution time (s)
          - 29.21
        * - Unique fields observed
          - 2


Statistics and predictions
==========================

.. dropdown:: Filtering from the :term:`CBC` events

    Here is how to filter :term:`BNS` and :term:`NSBH` events from the `Observing scenarios <https://m4opt.readthedocs.io/en/latest/guide/scenarios.html>`_.
    The following command will download the specified ZIP file, extract its contents, and filter the events based on your chosen criteria.

    .. admonition:: Zenodo API
        :class: info

        We have written a script for interacting with the Zenodo API, facilitating the download of files based on a DOI.
        This class provides functionality to retrieve the latest version DOI associated with a provided
        permanent DOI, and subsequently download the corresponding file from Zenodo.

        You can easily download another dataset from Zenodo by replacing the `permanent_doi`
        with a new one.

        Download the gravitational-waves simulation data from the `Zenodo database <https://zenodo.org/>`_

        .. tab-set::

            .. tab-item:: Using CLI

                .. code-block:: console

                    $ earthorbitplan.scenarios.zenodo_downloader --permanent-doi 14142969 --file-name runs_SNR-10.zip

            .. tab-item:: Using a config file

                .. code-block:: console

                    $ earthorbitplan.scenarios.zenodo_downloader --config ./earthorbitplan/config/params_ultrasat.ini


        .. note::
            For manual processing, see the source Zenodo dataset:
            `https://zenodo.org/records/14585837 <https://zenodo.org/records/14585837>`_


        .. admonition:: Unpack the zip file and filter the :term:`CBC` events
            :class: info

            This process automates the unpacking, filtering, and conversion of injection datasets
            (e.g., Farah / GWTC-3) from Zenodo ZIP archives. It processes event tables and associated
            localization files for specific observing runs (e.g., O5, O6), and outputs
            filtered ECSV tables and organized FITS files.

            The output will include an `.ecsv` file (`observing-scenarios.ecsv`) recording gravitational-waves` parameters such as mass, distance, and sky localization area.
            It will also copy the FITS files containing the gravitational-waves skymap probabilities into the directory specified by `--skymap-dir` (by default `./data/skymaps`) for each run.
            These outputs are useful for scheduling with :math:`\mathrm{M^4OPT}` and the statictric productions.


            .. tab-set::

                .. tab-item:: Using CLI

                    .. code-block:: console

                        $ earthorbitplan.workflow.unpacker --zip runs_SNR-10.zip --subdir runs_SNR-10 --runs O5 O6 --detectors HLVK --data-dir ./data --mass-threshold 3 --skymap-dir skymaps

                .. tab-item:: Using a config file

                    .. code-block:: console

                        $ earthorbitplan.workflow.unpacker --config ./earthorbitplan/config/params_ultrasat.ini


.. dropdown:: Submitting scheduling jobs in parallel or on a cluster


    .. admonition::  Why run scheduling jobs in parallel or on a cluster?
        :class: tip

        In gravitational-wave follow-up, researchers often need to process many sky maps or events quickly.
        Running scheduling jobs in parallel—using local multi-core processing, HTCondor, or SLURM—can significantly accelerate these computations.

        To select a specific execution backend, set the `backend` option in your configuration file.
        Available options include `"condor"`, `"parallel"`, `"slurm"`, or `"dask"`.


    .. admonition::  Submit jobs
        :class: tip

        The config file (e.g., :doc:`params_ultrasat.ini <../../config/params_ultrasat.ini>`)
        is used by default for ULTRASAT simulations.
        For other telescopes, use the relevant config files available in the
        :doc:`config <../../config>` directory.

        .. tab-set::

            .. tab-item::  Parallel execution

                **Local parallel execution** is ideal for small to medium workloads and can also be used on clusters.
                This approach distributes jobs across available CPU cores on a single machine or node.

                .. code-block:: console

                    python earthorbitplan.workflow.scheduler --config ../../config/params_ultrasat.ini --backend parallel

            .. tab-item::  SLURM

                Most of the clusters using  are cluster workload managers that handle large-scale job distribution across many compute nodes, making them ideal for processing many events efficiently.


                .. code-block:: console

                    python earthorbitplan.workflow.scheduler --config ./earthorbitplan/config/params_ultrasat.ini --backend slurm


            .. tab-item::  HTCondor

                **HTCondor** is a workload manager for high-throughput computing, suitable for running many independent jobs across a cluster or grid environment.
                We commonly used to run many independent jobs across clusters like the LIGO clusters at CIT, LHO, and LLO.

                This script is configured for use **only on LIGO clusters** (CIT, LHO, LLO).
                If you want to use HTCondor on a different cluster, you will need to update the backend implementation in
                :doc:`condor.py <../../backend/condor.py>`.

                .. code-block:: console

                    python earthorbitplan.workflow.scheduler --config ./earthorbitplan/config/params_ultrasat.ini --backend condor


            .. tab-item:: Dask

                **Dask** enables flexible parallel execution by dynamically distributing tasks across a cluster of workers managed by HTCondor.
                This backend is configured for use on LIGO clusters, but requires less cluster-specific configuration than the classic HTCondor mode.
                To adapt for another cluster running HTCondor, edit :doc:`dask.py <../../backend/dask.py>`.

                .. code-block:: console

                    python earthorbitplan.workflow.scheduler --config ./earthorbitplan/config/params_ultrasat.ini --backend dask


    .. note::

        The ini file :doc:`params_ultrasat.ini <../../config/params_ultrasat.ini>` contains all parameters and can be easily edited to adapt to another telescope's
        specifications, modify the observing campaign, or update output directories. The results are exactly the same as described in the :ref:`Run process` section,
        but here they are produced for multiple events for statistical analysis.
