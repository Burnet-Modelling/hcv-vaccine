# -*- coding: utf-8 -*-
"""
Created on Tue May 27 13:50:59 2025

@author: farah.houdroge

Script to:
    - Generate databook
    - Run YAML auto-calibration
    - Plot calibration
    - Run vaccine scenarios
    - Run sensitivity analyses
"""
from hcv.generate_databooks import generate_databook
from hcv import utils as ut
from hcv import utils_plotting as ut_plt
import sys

dir = ut.get_project_root()

country = sys.argv[1]
rand_seed = 250711
n_samples = 100
results_folder = dir / "results"
cal_folder = results_folder / "calibration" / "y-factors"
savedir_calib = results_folder / "calibration" / "plots"
savedir_calib.mkdir(parents=True, exist_ok=True)
savedir_scens = results_folder / "scenarios" / country
savedir_scens.mkdir(parents=True, exist_ok=True)
sens_folder = results_folder / "sensitivity_analyses"
sens_folder.mkdir(parents=True, exist_ok=True)
# %%
# Generate databook
generate_databook(country)

# Run point estimate calibration
ut.run_calibration(country, savedir=cal_folder)
ut_plt.plot_calibration(country, cal_folder=cal_folder, savedir=savedir_calib) # plot

# Run vaccine scenarios
ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens)
ut.econ_eval(country, savedir_scens, results_folder, rand_seed=rand_seed, n_samples=n_samples)

# Run sensitivity analyses
ut.run_sensitivity_analyses(country, cal_folder, sens_folder, results_folder)
ut.run_genotype_analyses(country, cal_folder, sens_folder, results_folder)
ut.run_coverage_analyses(country, cal_folder, sens_folder)
ut.econ_eval_central(country, sens_folder)
