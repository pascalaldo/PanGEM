from setuptools import find_packages, setup

setup(
    name="PanGEM",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas == 2.2.2",
        "numpy==1.26.4",
        "biopython == 1.83",
    ],
    entry_points={"console_scripts": ["pangem=PanGEM.cli:main"]},
)
