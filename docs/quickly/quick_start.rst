.. _quick_start:

Scientific Rationale
====================

The emergence of all-sky surveys and the increasing number of real-time gravitational wave (GW) alerts have created a pressing need for reliable scheduling frameworks.
These frameworks must efficiently coordinate follow-up observations across a range of telescopes and missions, both on the ground and in space.
To meet this need, we developed `|M⁴OPT| <https://m4opt.readthedocs.io/en/latest/>`_ : an open-source toolkit designed to optimize the scheduling of follow-up campaigns.

Why is scheduling so challenging?
---------------------------------
GW alerts often have large localization uncertainties—sometimes hundreds of square degrees—so astronomers need to decide how to balance *sky coverage*
(looking at as much area as possible) versus *depth* (spending enough time on each field to detect faint signals).
Moreover, electromagnetic (EM) counterparts to GW events can evolve rapidly, so observation plans must be generated quickly and efficiently.

:math:`\mathrm{M^4OPT}` addresses this by formulating the scheduling task as a **Mixed Integer Linear Programming (MILP)** problem.
This approach dynamically allocates observation time across fields and optimizes the exposure time for each one, maximizing the overall
probability of detection. Unlike fixed or manually tuned schedules, this method adapts to the conditions of each field, such as:

- background noise and background light (including natural sky brightness from the galaxy, zodiacal light, and instrumental effects;
  for ground-based telescopes, this also includes the increased sky brightness during astronomical twilight—just before sunrise or after sunset),
- distance uncertainty,
- instrumental sensitivity.


In addition, :math:`\mathrm{M^4OPT}` takes into account practical telescope constraints, such as Sun and Moon exclusion zones, dynamic field-of-view, and slewing time between fields.

This leads to a much more efficient use of telescope time and increases the chances of detecting faint or distant transients.

MILP Model and Constraints
--------------------------
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

Here is how you can run your first test with :math:`\mathrm{M^4OPT}` using the ULTRASAT mission.

1. Running a Scheduling Simulation
---------------------------------

Open your terminal and run the following command::

    m4opt schedule 14.fits 14.ecsv \
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


The simulation requires two inputs:
``./14.fits`` — the event skymap with localization probabilities
``./14.ecsv`` — the output file where the schedule will be saved

This command launches a scheduling simulation for the ULTRASAT mission.
You need to provide the main parameters, including the mission name, skygrid configuration, and observation settings.

- The output file (e.g., ``14.ecsv``) will contain the observation schedule.
- The simulation expects an event skymap file (usually a ``.fits`` file), which gives the localization probability of the event.

.. note::

   Missions like ULTRASAT support multiple skygrid models; use ``--skygrid`` to select (non-overlap and allsky).
   Other missions (e.g., ZTF, UVEX, Rubin) support only a single skygrid and do not need this option.

See the full list of parameters in the :ref:`CLI guide <https://m4opt.readthedocs.io/en/latest/guide/cli.html#m4opt-schedule>`_.


2. Understanding the Output
---------------------------

The generated ECSV file (e.g. ``14.ecsv``) contains your observation plan, including:

- Pointing coordinates,
- Exposure times,
- Slew (repositioning) times,
- Visit (by default: two visits per field),
- All relevant metadata.

By default, the schedule includes **two visits per field**—so each coordinate may appear twice, corresponding to repeated observations.


3. Visualizing the Schedule
---------------------------

You can create an animation or a PDF showing the planned observations::

    m4opt animate 14.ecsv 14_MOVIE.gif --dpi 300 --still 14_MOVIE.pdf

This produces:

- ``14_MOVIE.gif`` — an animation of the schedule
- ``14_MOVIE.pdf`` — a static pdf,  of the observation sequence.

.. image:: ./14_MOVIE.gif
   :alt: Example animation of the observation plan
   :align: center

This workflow lets you quickly simulate and visualize follow-up plans for your favorite mission.
For more details and advanced options, check out the `full documentation <https://m4opt.readthedocs.io/en/latest/>`_.

**Explanation of the animation:**

- The pink regions show the scheduled observation pointings the :ref:`footprints <https://m4opt.readthedocs.io/en/latest/api/m4opt.fov.footprint.html#footprint>`_.
- The green outline marks the 90% credible region of the GW localization.
- The deep blue areas are always outside the telescope’s Field of Regard; the light blue areas are temporarily out of view.
- The lower panel shows how the detection probability and covered sky area accumulate over time, with different colors indicating
  the number of times a region has been observed.
- The symbol :math:`\oplus` shows the direction of the center of the Earth (sub-Earth point) projected onto the sky.
- The symbol :math:`\odot` shows the direction of the Sun (sub-solar point) on the sky.

For more details of marker conventions,
see the `ligo.skymap plotting documentation <https://lscsoft.docs.ligo.org/ligo.skymap/plot/marker.html#module-ligo.skymap.plot.marker/>`_.

.. note::

   This is a projection of the sky, **not a direct image of the Earth or the Moon**. The features shown correspond to sky coordinates,
   not to physical locations on Earth or lunar positions.


4. ECSV file inspection
-----------------------

You can load and inspect a schedule file using Astropy:

.. code-block:: console

   >>> from astropy.table import QTable
   >>> plan = QTable.read("14.ecsv", format="ascii.ecsv")
   >>> obs = plan[plan["action"] == "observe"]
   >>> display = obs["start_time", "duration"]
   >>> display["ra"] = obs["target_coord"].ra
   >>> display["dec"] = obs["target_coord"].dec
   >>> display.round({'duration': 1, 'ra': 2, 'dec': 2})
   >>> print(display)
        start_time          duration   ra    dec
                               s       deg   deg
    ----------------------- -------- ------ -----
    2012-07-14 16:04:59.480   1080.0 221.14 58.26
    2012-07-14 17:17:11.127   1080.0 221.14 58.26
    2012-07-14 17:35:59.480   3786.3 218.06 43.89
    2012-07-14 19:09:05.819   3786.3 218.06 43.89


5. ECSV Metadata Extraction
---------------------------

Load a schedule, extract key metadata and visit counts:

.. code-block:: console

   >>> from astropy.table import QTable
   >>> plan = QTable.read("14.ecsv", format="ascii.ecsv")
   >>> objective = plan.meta.get("objective_value")
   >>> best_bound = plan.meta.get("best_bound")
   >>> status = plan.meta.get("solution_status")
   >>> time_used = plan.meta.get("solution_time")
   >>> visits = plan.meta.get("args", {}).get("visits", 2)
   >>> n_obs = len(plan[plan["action"] == "observe"])
   >>> unique_fields = n_obs // visits
   >>> print("Schedule metadata:")
   >>> print(f" • Objective value: {objective:.4f}")
   >>> print(f" • Best bound: {best_bound:.4f}")
   >>> print(f" • Solver status: {status}")
   >>> print(f" • Solution time: {time_used}")
   >>> print(f" • Unique fields observed: {unique_fields}")
   Schedule metadata:
    • Objective value: 0.9483
    • Best bound: 0.9483
    • Solver status: integer optimal solution
    • Solution time: 29.206 s
    • Unique fields observed: 2

.. list-table:: Schedule metadata summary
   :header-rows: 1
   :widths: 30 15

   * - Metric                   - Value
   * - Objective value          - 0.9483
   * - Best bound               - 0.9483
   * - Solver status            - integer optimal solution
   * - Solution time (s)        - 29.21
   * - Unique fields observed   - 2



.. .. list-table:: Sample observations
..    :header-rows: 1
..    :widths: 20 12 8 8

..    * - start_time
..      - duration (s)
..      - ra (deg)
..      - dec (deg)
..    * - 2012‑07‑14 16:04:59.480
..      - 1080.0
..      - 221.14
..      - 58.26
..    * - 2012‑07‑14 17:17:11.127
..      - 1080.0
..      - 221.14
..      - 58.26
..    * - 2012‑07‑14 17:35:59.480
..      - 3786.3
..      - 218.06
..      - 43.89
..    * - 2012‑07‑14 19:09:05.819
..      - 3786.3
..      - 218.06
..      - 43.89
