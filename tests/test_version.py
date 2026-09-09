"""Import / version smoke test."""

import earthorbitplan


def test_version_is_exposed():
    assert isinstance(earthorbitplan.__version__, str)
    assert earthorbitplan.__version__
    # setuptools_scm always produces at least "MAJOR.MINOR..."
    assert earthorbitplan.__version__[0].isdigit()
