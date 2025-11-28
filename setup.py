from setuptools import setup, find_packages

setup(
    name="hcv-vaccine",
    description="A simulation model for an HCV vaccine",
    version="1.0",
    url="https://github.com/Burnet-Modelling/hcv-vaccine",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "sciris",
        "pandas",
        "scipy",
        "pathlib",
        "openpyxl",
        ],
    python_requires='~=3.11',
)
