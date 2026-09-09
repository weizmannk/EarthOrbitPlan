# Contributing to EarthOrbitPlan

Thanks for taking the time to contribute. EarthOrbitPlan is an educational
framework built on top of [M4OPT](https://github.com/m4opt/m4opt); the notes
below are the minimum you need to get a working development setup and to get a
change merged.

## Development environment

EarthOrbitPlan targets **Python 3.12+** (floor tracks M4OPT).

```bash
git clone https://github.com/weizmannk/EarthOrbitPlan.git
cd EarthOrbitPlan
python3.12 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e ".[dev,docs]"
```

CPLEX is required for the M4OPT scheduler. See the
[M4OPT CPLEX install guide](https://m4opt.readthedocs.io/en/latest/install/cplex.html).

`[dev]` pulls in `pre-commit` and the test dependencies (`[test]`); add `[docs]`
only if you intend to build the documentation.

## Running the tests

```bash
pytest
```

Tests live in `tests/`. They cover the pure helpers (rate statistics, project
paths, table interleaving) plus one regression test that pins the output of the
Poisson–log-normal rate quantiles against a stored fixture in
`tests/data/`. Parts of the package that require CPLEX, network access, or large
skymaps are not exercised by the suite.

## Type checking

```bash
pip install mypy
mypy -p earthorbitplan
```

Configuration is in `[tool.mypy]` of `pyproject.toml`. Run it with the package
installed (`pip install -e ".[test]"`) so mypy resolves imports against the real
`m4opt` / `astropy` APIs — that is what catches a renamed or removed upstream
symbol before it reaches a user. CI runs the same check in the `mypy` job.

## Pre-commit

```bash
pre-commit install      # once, to run on every commit
pre-commit run --all-files
```

This runs ruff (lint + format), codespell, nbstripout, black on notebooks, and
the standard hygiene hooks. CI runs the exact same command and fails on any
change, so run it locally before pushing.

## Building the docs

```bash
pip install -e ".[docs]"
cd docs
make html            # output in docs/_build/html
```

CI builds the docs with `-W` (warnings are errors); keep the build clean.

## Branches, commits, and pull requests

- Branch off `main`. Use short prefixes: `feat/`, `fix/`, `docs/`, `chore/`,
  `ci/`, `test/`, `refactor/`.
- Write [Conventional Commit](https://www.conventionalcommits.org/) messages
  (`feat:`, `fix:`, `docs:`, `chore:`, `ci:`, `test:`, `refactor:`). Keep
  commits small and self-contained.
- Add a line under `## [Unreleased]` in `CHANGELOG.md` for any user-visible
  change.
- Never change scientific logic or numerical results in a hygiene/refactor PR
  without calling it out explicitly and updating the regression fixtures in the
  same PR.
- Open a pull request against `main`. CI (pre-commit, pytest, docs) must be
  green. At least one maintainer review is required before merge; the reviewer
  merges. Do not force-push shared branches.
