from setuptools import setup, find_packages

setup(
    name="dietscan",
    version="1.0.0",
    author="Antton Alberdi",
    author_email="antton.alberdi@sund.ku.dk",
    description="Dietary profiling from metagenomic data",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "argparse",
        "PyYAML",
        "requests",
        "snakemake"
    ],
    entry_points={
        "console_scripts": [
            "dietscan=dietscan.cli:main"
        ],
    },
    python_requires=">=3.6",
)
