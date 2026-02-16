# -*- coding: utf-8 -*-
"""
Created on Tue May 27 13:50:59 2025

@author: farah.houdroge

Script to:
    - Generate databook
    - Run YAML auto-calibration
    - Plot calibration
    - Run vaccine scenarios
    - Post-process & plot results
"""
from hcv.generate_databooks import generate_databook
from gitlab_utils import _get_sharepoint_folder
from hcv import utils as ut
from hcv import utils_plotting as ut_plt
from pathlib import Path
import sys
# dir = pathlib.Path(_get_sharepoint_folder())
dir = ut.get_project_root()

country = sys.argv[1]
rand_seed=250711
n_samples=100
results_folder = dir/'results'
results_folder.mkdir(parents=True, exist_ok=True)
cal_folder_parent = results_folder / 'calibration'
cal_folder_parent.mkdir(parents=True, exist_ok=True)
cal_folder = cal_folder_parent / 'y-factors'
cal_folder.mkdir(parents=True, exist_ok=True)
savedir_calib = cal_folder_parent / 'plots'
savedir_calib.mkdir(parents=True, exist_ok=True)
savedir_scens = results_folder / 'scenarios' / country
savedir_scens.mkdir(parents=True, exist_ok=True)
sens_folder = results_folder / 'sensitivity_analyses'
sens_folder.mkdir(parents=True, exist_ok=True)

#%%
# Generate databook
# generate_databook(country)
#%%
# Run calibration
ut.run_calibration(country, savedir=cal_folder)
ut_plt.plot_calibration(country, cal_folder=cal_folder, savedir=savedir_calib,cal_version='v1') # plot

# Run vaccine scenarios
# ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens)
# ut.econ_eval(country, savedir_scens, results_folder, rand_seed=rand_seed, n_samples=n_samples)

# # Run sensitivity analyses
# ut.run_sensitivity_analyses(country, cal_folder, sens_folder, results_folder)
# ut.run_genotype_analyses(country, cal_folder, sens_folder, results_folder)
# ut.run_coverage_analyses(country, cal_folder, sens_folder)
# ut.econ_eval_central(country, sens_folder)
