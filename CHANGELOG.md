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
- GitHub Actions CI: pre-commit, pytest (Python 3.11), and a warnings-as-errors
  docs build.

### Changed
- `earthorbitplan.__version__` is now importable; the setuptools_scm version
  file is written to `earthorbitplan/_version.py` instead of the repository root.
- `license` is declared as the SPDX string `BSD-3-Clause`; `pyproject.toml` is
  the single source of dependency truth.
- `build_uvex_followup_skymap.py` moved to `scripts/`.

### Fixed
- Pre-commit `exclude` patterns written as `|` block scalars were silently
  disabled (trailing newline in the regex); the PSD data file is excluded again.

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
