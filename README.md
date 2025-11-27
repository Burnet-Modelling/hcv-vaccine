# hcv-vaccine
Compartmental model and analysis pipeline for hepatitis C (HCV) vaccine evaluations.

## What this project does

This repository contains a simulation framework to evaluate the epidemiological and economic impact of hypothetical HCV vaccines. It includes code to:

- Generate country-level databooks used to parameterise models
- Run automatic calibration workflows (YAML-driven) to fit model parameters
- Run vaccine scenario ensembles and sensitivity analyses
- Perform economic evaluation and produce aggregated outputs and plots
- Post-process and plot calibration, scenario and sensitivity results

The primary entrypoints are `run_main.py` for running simulations per country and `run_postprocess.py` for aggregating and plotting results.

## Key features and benefits

- Reproducible, scripted pipeline for calibration → scenario → sensitivity → post-processing
- Country-level databook generation and YAML-driven calibration for scalable runs
- Uses standard scientific Python stack for numerical work and plotting
- Designed for batch execution (see `job.conf` example for running many countries)

## Quick start

Prerequisites

- Python 3.11 

Install (recommended: create a virtual environment first)

```powershell
python -m pip install -e .
```

This installs the package in editable mode and pulls the runtime dependencies listed in `setup.py` (numpy, pandas, scipy, sciris, atomica, openpyxl, etc.).

## Running the model

The primary script for running the end-to-end pipeline for a single country is `run_main.py`.

Usage (example for Afghanistan `AFG`):

```powershell
python run_main.py AFG
```

What it does:

- Calls `generate_databook(country)` to create the country databook used by the model
- Runs calibration and plots calibration diagnostics
- Runs vaccine scenario sampling and economic evaluation
- Runs sensitivity and genotype coverage analyses

For batch runs, `job.conf` contains an example list of commands that call a wrapper (`srun_main.py` in your cluster environment) across many countries.

## Post-processing and plotting

Use `run_postprocess.py` to aggregate scenario outputs, produce regional and global summary tables and generate publication-ready figures.

```powershell
python run_postprocess.py
```

The post-processing script expects results to be under the `results/` folder with the structure produced by `run_main.py` (calibration, scenarios, sensitivity_analyses, figures).

## Project structure

- `hcv/` — package with model code, databook generation, utilities and plotting helpers
- `data/` — input spreadsheets and shapefiles used to build databooks and for mapping
- `databooks/` — generated per-country databooks 
- `results/` — output folder used by scripts for calibration, scenarios and figures
- `run_main.py` — per-country pipeline entrypoint
- `run_postprocess.py` — aggregation and plotting 
- `job.conf` — example batch commands for running many countries
- `setup.py` — package and dependency information

## Dependencies

Declared in `setup.py`. Key packages:

- numpy, pandas, scipy: numerical and data handling
- sciris: utilities
- atomica: used for compartmental modelling
- openpyxl: read/write Excel files

Install via pip as shown above.

## Maintainers and contribution

Maintainer: Farah Houdroge

Contributors: Phillip Luong, Jessica Zuk, Chris Seaman, Kelly Maynard

If you want to contribute:

1. Fork the repository and create a feature branch
2. Run locally 
3. Open a pull request with a clear description of changes


## Acknowledgements

This project was developed by the Burnet Modelling group in collaboration with 1DaySooner.

