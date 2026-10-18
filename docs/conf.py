import re
import tomllib
from pathlib import Path

root = Path(__file__).resolve().parent.parent
release = tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"]

project = "BarcodeMapper"
author = "Antton Alberdi"
copyright = "2026, Antton Alberdi"
version = release

extensions = ["myst_parser"]
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build", "requirements.txt", "validation.md", "releasing.md"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "../images/barcodemapper.png"
html_theme_options = {"logo_only": False}
