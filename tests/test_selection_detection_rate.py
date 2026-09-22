"""Tests for the aggregated block and the P(N=0) column of the rate table.

The aggregate is the delicate part. The merger-rate prior is a single
astrophysical rate shared by every observing run, so its uncertainty is
*fully correlated* across them. Combining runs therefore means adding the
expected counts and re-deriving the quantiles once, at the same prior width.
Two tempting shortcuts are wrong and are guarded against below:

* summing the per-run quantiles -- the median is not additive, and the
  detected medians are all zero while their aggregate is not;
* convolving the per-run predictive distributions as if independent, which
  would treat one rate uncertainty as three.
"""

import re

import numpy as np
import pytest
from astropy.table import QTable
from scipy import stats

from earthorbitplan.probability.rate import (
    poisson_lognormal_rate_cdf,
    poisson_lognormal_rate_quantiles,
)
from earthorbitplan.utils.path import get_project_root
from earthorbitplan.workflow.selection_detection_rate import (
    summarize_selected_detected_events,
)

EVENTS = get_project_root() / "data" / "ultrasat" / "ultrasat_non-overlap.ecsv"
O5 = ("O5a", "O5b", "O5c")
PRIOR = dict(merger_rate_lo=50, merger_rate_mid=100, merger_rate_hi=203)

# Reference aggregate of O5a+O5b+O5c for the ULTRASAT LCS grid.
#   (row, class): (lambda, median, 5%, 95%, P(N=0) in percent)
REFERENCE = {
    ("Selected", "BNS"): (5.91, 5.3, 1.0, 13.1, 1.1),
    ("Selected", "NSBH"): (3.30, 2.7, 0.0, 7.8, 5.9),
    ("Selected", "All"): (9.21, 8.6, 2.5, 19.8, 0.2),
    ("Detected", "BNS"): (2.07, 1.4, 0.0, 5.2, 14.8),
    ("Detected", "NSBH"): (1.06, 0.4, 0.0, 2.9, 34.9),
    ("Detected", "All"): (3.14, 2.5, 0.0, 7.4, 6.6),
}
ROWS = ["Selected", "Detected"]
CLASSES = ["BNS", "NSBH", "All"]

pytestmark = pytest.mark.skipif(
    not EVENTS.exists(), reason=f"{EVENTS} not available in this checkout"
)


def summarize(**kwargs):
    return summarize_selected_detected_events(
        EVENTS,
        run_duration=1.0,
        poisson_lognormal_rate_quantiles=poisson_lognormal_rate_quantiles,
        verbose=False,
        **PRIOR,
        **kwargs,
    )


def aggregate_arrays():
    """Recompute the aggregate independently of the table renderer."""
    table = QTable.read(EVENTS)
    table = table[table["objective_value"] >= table["cutoff"][0]]
    effective = {
        run: value.to_value("Gpc-3 yr-1")
        for run, value in table.meta["effective_rate"].items()
    }
    (width_90,) = np.diff(stats.norm.interval(0.9))
    sigma = np.log(PRIOR["merger_rate_hi"] / PRIOR["merger_rate_lo"]) / width_90

    is_class = {
        "BNS": lambda t: t["source_class"] == "BNS",
        "NSBH": lambda t: t["source_class"] == "NSBH",
        "All": lambda t: np.ones(len(t), dtype=bool),
    }
    lam = np.zeros((2, len(CLASSES)))
    for run in O5:
        events = table[table["run"] == run]
        scale = PRIOR["merger_rate_mid"] / effective[run]
        for i_cls, cls in enumerate(CLASSES):
            subset = events[is_class[cls](events)]
            lam[0, i_cls] += scale * len(subset)
            lam[1, i_cls] += scale * np.sum(
                subset["detection_probability_known_position"]
            )
    return lam, np.log(lam), sigma


# --------------------------------------------------------------------------
# Default behaviour is untouched
# --------------------------------------------------------------------------
def test_defaults_add_no_aggregate_row_and_no_zero_column():
    latex = summarize()
    assert "O5 total" not in latex
    assert "P(N{=}0)" not in latex
    # two sub-columns per class: lambda and the interval
    assert r"\begin{tabular}{ll" + "cc" * len(CLASSES) + "}" in latex


def test_aggregate_leaves_the_per_run_rows_untouched():
    without = summarize()
    with_aggregate = summarize(aggregate_runs=O5)
    per_run = [line for line in without.splitlines() if line.endswith(r"\\")]
    still_there = [line for line in with_aggregate.splitlines() if line.endswith(r"\\")]
    for line in per_run:
        assert line in still_there


# --------------------------------------------------------------------------
# Reference values
# --------------------------------------------------------------------------
@pytest.mark.parametrize("row", ROWS)
@pytest.mark.parametrize("cls", CLASSES)
def test_aggregate_matches_reference(row, cls):
    lam, mu, sigma = aggregate_arrays()
    i_row, i_cls = ROWS.index(row), CLASSES.index(cls)
    want_lam, want_med, want_lo, want_hi, want_zero = REFERENCE[row, cls]

    assert lam[i_row, i_cls] == pytest.approx(want_lam, abs=0.1)

    median, lo, hi = poisson_lognormal_rate_quantiles(
        np.array([0.5, 0.05, 0.95]), mu[i_row, i_cls], sigma
    )
    assert median == pytest.approx(want_med, abs=0.1)
    assert lo == pytest.approx(want_lo, abs=0.1)
    assert hi == pytest.approx(want_hi, abs=0.1)

    zero = 100 * poisson_lognormal_rate_cdf(0, mu[i_row, i_cls], sigma)
    assert zero == pytest.approx(want_zero, abs=0.2)


def test_prior_width_is_the_expected_sigma():
    _, _, sigma = aggregate_arrays()
    assert sigma == pytest.approx(0.4259, abs=1e-4)


# --------------------------------------------------------------------------
# The aggregate must not be built by summing quantiles
# --------------------------------------------------------------------------
def test_median_is_not_additive():
    """Guard against re-implementing the aggregate as a sum of quantiles.

    Every per-run Detected median rounds to zero, so a quantile-summing
    implementation would report zero for the aggregate too. The correct
    aggregate is about 2.5.
    """
    table = QTable.read(EVENTS)
    table = table[table["objective_value"] >= table["cutoff"][0]]
    effective = {
        run: value.to_value("Gpc-3 yr-1")
        for run, value in table.meta["effective_rate"].items()
    }
    _, mu_aggregate, sigma = aggregate_arrays()

    per_run_medians = []
    for run in O5:
        events = table[table["run"] == run]
        lam_run = (PRIOR["merger_rate_mid"] / effective[run]) * np.sum(
            events["detection_probability_known_position"]
        )
        per_run_medians.append(
            float(poisson_lognormal_rate_quantiles(0.5, np.log(lam_run), sigma))
        )

    correct = float(
        poisson_lognormal_rate_quantiles(
            0.5, mu_aggregate[1, CLASSES.index("All")], sigma
        )
    )

    # Every printed per-run median is zero, so reading them off the table and
    # adding them gives exactly zero.
    assert all(np.rint(m) == 0 for m in per_run_medians), per_run_medians
    assert sum(np.rint(m) for m in per_run_medians) == 0

    # Even adding the unrounded medians falls far short: the median of a sum
    # is not the sum of the medians.
    assert sum(per_run_medians) == pytest.approx(1.16, abs=0.1)

    assert correct == pytest.approx(2.5, abs=0.1)
    assert correct > 2 * sum(per_run_medians), "aggregate must not sum medians"


def aggregate_rows(latex):
    """The two body rows of the aggregate block, from the tabular only."""
    tabular = latex.split(r"\begin{tabular}")[1]
    body = [line for line in tabular.splitlines() if line.rstrip().endswith(r"\\")]
    selected = next(line for line in body if "O5 total" in line)
    return selected, body[body.index(selected) + 1]


def test_aggregate_detected_row_is_not_zero_like_the_per_run_rows():
    """The per-run Detected medians all print as 0; the aggregate must not."""
    latex = summarize(aggregate_runs=O5)
    _, detected = aggregate_rows(latex)

    per_run_detected = [
        line
        for line in latex.split(r"\begin{tabular}")[1].splitlines()
        if "& Detected &" in line and "O5 total" not in line
    ]
    assert per_run_detected
    for line in per_run_detected:
        assert "$0^" in line or "$<1$" in line

    # The All column is where the contrast shows. Its per-run Detected median
    # prints as 0 in every sub-period, while the aggregate prints 3. NSBH
    # still rounds to 0 even aggregated (its median is 0.4), which is why the
    # assertion targets the last class rather than all of them.
    all_column = detected.split("&")[-1]
    assert "$<1$" not in all_column, detected
    median = int(re.search(r"\$(\d+)\^", all_column).group(1))
    assert median >= 2, detected


# --------------------------------------------------------------------------
# P(N = 0)
# --------------------------------------------------------------------------
def test_zero_probability_adds_one_column_per_class():
    latex = summarize(include_zero_probability=True)
    assert r"$P(N{=}0)$" in latex
    assert r"\begin{tabular}{ll" + "ccc" * len(CLASSES) + "}" in latex
    assert re.search(r"\d+\.\d\\%", latex)


def test_zero_probability_is_a_percentage_with_one_decimal():
    latex = summarize(include_zero_probability=True)
    for value in re.findall(r"(\d+\.\d)\\%", latex):
        assert 0.0 <= float(value) <= 100.0


def test_zero_probability_decreases_with_lambda():
    """P(N=0) must fall as the expected count rises."""
    _, mu, sigma = aggregate_arrays()
    detected = poisson_lognormal_rate_cdf(0, mu[1], sigma)
    selected = poisson_lognormal_rate_cdf(0, mu[0], sigma)
    assert np.all(selected < detected)


# --------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------
def test_unknown_aggregate_run_raises():
    with pytest.raises(ValueError, match="absent"):
        summarize(aggregate_runs=("O5a", "O9z"))
