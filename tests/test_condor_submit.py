"""Checks on the HTCondor submit description generated for m4opt jobs.

The memory and disk requests are ClassAd expressions rather than constants, so
that a job held for exceeding a limit is released and retried with a larger
request. Getting the units or the escalation ceiling wrong fails silently --
the job is simply never matched -- so both are pinned here.
"""

import subprocess

import pytest

from earthorbitplan.backend import condor


@pytest.fixture
def submit_description(monkeypatch):
    """Return the submit file that ``submit_condor_job`` feeds to condor_submit."""
    captured = {}

    class FakeProc:
        returncode = 0

        def communicate(self, input=None):
            captured["script"] = input
            return "", ""

    monkeypatch.setattr(subprocess, "Popen", lambda *args, **kwargs: FakeProc())

    def render(**kwargs):
        condor.submit_condor_job("O5a", "1234", "/logs", "/wrapper.sh", **kwargs)
        return dict(
            line.split(" = ", 1)
            for line in captured["script"].splitlines()
            if " = " in line
        )

    return render


@pytest.mark.parametrize(
    "value,expected_kb",
    [
        ("4096", 4096),
        ("8000 MB", 8192000),
        ("8 GB", 8388608),
        ("80000MB", 81920000),
    ],
)
def test_to_kb(value, expected_kb):
    assert condor._to_kb(value) == expected_kb


def test_to_kb_rejects_garbage():
    with pytest.raises(ValueError):
        condor._to_kb("plenty")


def test_default_requests_escalate_within_a_ceiling(submit_description):
    description = submit_description()
    # request_memory is in MiB, request_disk in KiB: 8000 MB -> 8192000 KiB.
    assert description["request_memory"] == (
        "ifthenelse(MemoryUsage =!= undefined, "
        "min({200000, max({40000, MemoryUsage * 2.0})}), 40000)"
    )
    assert description["request_disk"] == (
        "ifthenelse(DiskUsage =!= undefined, max({8192000, DiskUsage * 2.0}), 8192000)"
    )
    assert description["periodic_release"] == (
        "(HoldReasonCode == 34 || HoldReasonCode == 104) && (NumJobStarts <= 3)"
    )


def test_requests_are_configurable(submit_description):
    description = submit_description(
        request_memory="16 GB",
        request_disk="2 GB",
        max_request_memory="64 GB",
        growth_factor=1.5,
        max_retries=5,
    )
    assert description["request_memory"] == (
        "ifthenelse(MemoryUsage =!= undefined, "
        "min({65536, max({16384, MemoryUsage * 1.5})}), 16384)"
    )
    assert description["request_disk"] == (
        "ifthenelse(DiskUsage =!= undefined, max({2097152, DiskUsage * 1.5}), 2097152)"
    )
    assert "NumJobStarts <= 5" in description["periodic_release"]


def test_ceiling_never_falls_below_the_floor(submit_description):
    """A ceiling below the floor would make the job unmatchable from the start."""
    description = submit_description(
        request_memory="100 GB", max_request_memory="10 GB"
    )
    assert "min({102400, max({102400," in description["request_memory"]
