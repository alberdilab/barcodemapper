"""Installed package version, sourced from pyproject.toml metadata."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("barcodemapper")
except PackageNotFoundError:
    # A source checkout can be imported before it has been installed.
    __version__ = "0+uninstalled"
