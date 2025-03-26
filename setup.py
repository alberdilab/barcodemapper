import os
import re
from setuptools import setup, find_packages

def get_version():
    with open(os.path.join("barcodemapper", "__version__.py")) as f:
        content = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name="barcodemapper",
    version=get_version(),
    author="Antton Alberdi",
    author_email="antton.alberdi@sund.ku.dk",
    description="Metabarcoding marker gene mapper for metagenomic data",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "argparse",
        "PyYAML",
        "requests",
        "plotly",
        "biopython"
    ],
    entry_points={
        "console_scripts": [
            "barcodemapper=barcodemapper.cli:main"
        ],
    },
    python_requires=">=3.12",
)
