from pathlib import Path


def get_project_root() -> Path:
    """
    Return the project root directory (the one holding ``pyproject.toml`` and ``data/``).

    Walks up from this file. Falls back to the current working directory when
    ``__file__`` is unavailable (e.g. an interactive session).

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
        path = Path(__file__).resolve()
    except NameError:
        path = Path.cwd().resolve()
    for parent in [path, *path.parents]:
        if (parent / "pyproject.toml").is_file() and (parent / "data").is_dir():
            return parent
    raise RuntimeError(
        "Could not find the project root (looked for 'pyproject.toml' and 'data/')."
    )
