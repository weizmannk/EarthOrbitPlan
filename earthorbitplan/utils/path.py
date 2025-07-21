from pathlib import Path


def get_project_root() -> Path:
    """
    Return the project root directory (containing 'data/', 'earthorbitplan/', and 'pyproject.toml').

    This function works from anywhere in the project—script, notebook, etc.
    It checks for the presence of these marker directories/files to identify the root.

    Returns
    -------
    Path
        Path to the project root directory.

    Raises
    ------
    RuntimeError
        If the project root cannot be found.
    """
    try:
        path = Path("__file__").resolve()
    except NameError:
        path = Path.cwd().resolve()
    for parent in [path] + list(path.parents):
        if (
            (parent / "data").exists()
            and (parent / "earthorbitplan").exists()
            and (parent / "pyproject.toml").is_file()
        ):
            return parent
    raise RuntimeError(
        "Could not find the project root (looked for 'data/', 'earthorbitplan/', and 'pyproject.toml')."
    )
