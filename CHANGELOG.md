# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `CONTRIBUTING.md`, `CHANGELOG.md`, and `CITATION.cff`.
- `[dev]` optional-dependency group in `pyproject.toml` (`pre-commit` + test deps).
- Test suite under `tests/` (pure-helper unit tests plus a regression test
  pinning the Poisson–log-normal rate quantiles against a stored fixture).
- GitHub Actions CI: `pre-commit`, `mypy` (static type check), `pytest` on
  Python 3.12 / 3.13 / 3.14, and a Sphinx docs build.

### Changed
- **`src/` layout**: the package now lives at `src/earthorbitplan/` (matches
  `m4opt`). setuptools_scm writes `src/earthorbitplan/_version.py`.
- **Python floor raised to 3.12** (`requires-python = ">=3.12, <3.15"`), tracking
  `m4opt` v2.13; CI runs 3.12 and 3.13.
- `earthorbitplan.__version__` is now importable.
- `license` is declared as the SPDX string `BSD-3-Clause`; `pyproject.toml` is
  the single source of dependency truth.
- Non-package files moved out of the package: notebooks and helper scripts to
  `scripts/notebooks/`, paper figures/tables to `scripts/paper_plots/`,
  `build_uvex_followup_skymap.py` to `scripts/`.
- `get_project_root()` now keys off `pyproject.toml` + `data/` and uses the real
  `__file__` (was a quoted string literal, silently cwd-dependent).
- Observing-scenario downloader points at the GWTC-5.0 datasets (concept DOIs
  `22550047` FullPop / `22555948` PixelPop); `unpacker` gains `--pop
  {fullpop,pixelpop}`.
- **m4opt v2.13 API**: dropped the removed `m4opt.synphot.background.update_missions`
  calls in `detection_probability`, `mission_limmag` and `area_distance`. The
  Cerenkov / AE8 background is now part of the mission background model and is
  applied automatically from each `observing()` context. **Numerical outputs for
  ULTRASAT must be revalidated against the pre-v2.13 results.**

### Fixed
- Pre-commit `exclude` patterns written as `|` block scalars were silently
  disabled (trailing newline in the regex); the PSD data file is excluded again.
- Pre-commit `exclude` paths updated for the `src/` layout.
- `docs/tutorials/observing_scenarios.ipynb` reads `farah.h5` from
  `src/earthorbitplan/scenarios/` (was a stale in-tree path); the duplicate
  32 MB `docs/tutorials/farah.h5` copy is removed.
- sphinx-gallery output (`docs/auto_tutorials/`, `sg_execution_times.rst`) is
  no longer committed.

### Removed
- `old_requirements.txt` and `requirements.txt` (unused / superseded by
  `pyproject.toml`).

## [0.1.0] - 2026-08-28

Initial release. Reconstructed from the project history (first commit
2025-03-04).

### Added
- `earthorbitplan` package with the end-to-end follow-up workflow:
  - `workflow/`: scheduler wrapper, plan post-processing, table unpacking,
    limiting-magnitude computation, area–distance and skymap-volume analyses,
    LaTeX table generation, and visualization.
  - `probability/`: kilonova detection probability, and rate error propagation
    (beta-binomial, Poisson–log-normal quantiles) after Singer et al.
  - `backend/`: parallel execution backends (joblib, Dask, HTCondor, Slurm).
  - `scenarios/`: bundled PSDs, scenario parameters, and a resumable Zenodo
    downloader for the observing-scenarios data.
  - `config/`: mission parameter files for UVEX and ULTRASAT.
- Detection-probability path using the Cerenkov background for ULTRASAT.
- Educational tutorials and Jupyter notebooks (observing scenarios, kilonova
  detection rate, field of view, skymap volume) with Binder / sphinx-thebe
  execution.
- Sphinx documentation (Read the Docs) with a glossary, install guide,
  multi-messenger and scenario sections, an M4OPT scheduler walkthrough, and a
  tilepy integration section.

[Unreleased]: https://github.com/weizmannk/EarthOrbitPlan/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/weizmannk/EarthOrbitPlan/releases/tag/v0.1.0
