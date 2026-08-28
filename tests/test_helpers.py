"""Unit tests for the small pure helpers in the package."""

import sys
from pathlib import Path

import numpy as np

from earthorbitplan.utils.path import get_project_root
from earthorbitplan.workflow.mission_limmag import _masked_col, parse_arguments

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# get_project_root
# --------------------------------------------------------------------------
def test_get_project_root_finds_the_markers(monkeypatch):
    # get_project_root() walks up from the current working directory looking
    # for the data/, earthorbitplan/ and pyproject.toml markers.
    monkeypatch.chdir(REPO_ROOT)
    root = get_project_root()
    assert (root / "pyproject.toml").is_file()
    assert (root / "earthorbitplan").is_dir()
    assert (root / "data").is_dir()


# --------------------------------------------------------------------------
# _masked_col
# --------------------------------------------------------------------------
def test_masked_col_interleaves_and_masks_slew_rows():
    col = _masked_col([1.0, 2.0, 3.0], n_obs=3)
    # each value is duplicated (observe + slew) and the trailing slew is dropped
    np.testing.assert_array_equal(col.unmasked, [1.0, 1.0, 2.0, 2.0, 3.0])
    # odd (slew) positions are masked, even (observe) positions are not
    np.testing.assert_array_equal(col.mask, [False, True, False, True, False])


def test_masked_col_length_is_2n_minus_1():
    for n in (1, 2, 5):
        assert len(_masked_col(np.arange(n), n_obs=n)) == 2 * n - 1


# --------------------------------------------------------------------------
# parse_arguments
# --------------------------------------------------------------------------
def test_parse_arguments_reads_ini_params_section(tmp_path, monkeypatch):
    ini = tmp_path / "params.ini"
    ini.write_text(
        "[params]\n"
        "input_table = ultrasat_bns.ecsv\n"
        "plan_dir = data/O5\n"
        "output_dir = add_limmag\n"
    )
    monkeypatch.setattr(sys, "argv", ["mission_limmag.py", "--config", str(ini)])
    args = parse_arguments()
    assert args.input_table == "ultrasat_bns.ecsv"
    assert args.plan_dir == "data/O5"
    assert args.output_dir == "add_limmag"


def test_parse_arguments_falls_back_to_defaults(tmp_path, monkeypatch):
    ini = tmp_path / "empty.ini"
    ini.write_text("[params]\n")
    monkeypatch.setattr(sys, "argv", ["mission_limmag.py", "--config", str(ini)])
    args = parse_arguments()
    assert args.input_table == "bns.ecsv"
    assert args.plan_dir == "data/O5"
    assert args.output_dir == "add_limmag"
