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
P = ut.project(country, load_calibration=True, cal_folder=cal_folder, load_programs=False,cal_version='v1')
_, D = ut.return_fw_db(country)
parset = P.make_parset()
parset.load_calibration(cal_folder/f'{country}_calibration.xlsx')   
result = P.run_sim(parset=parset) 
ut_plt.plot_outcomes2(P, result, country, start_year=2000, end_year=2050)

cal = P.parsets.copy('default','Calibrated')

par_adj = ['18-64_males','PWID_males']
pop_adj = ['Prisoners_males','Prisoners_PWID_males']
# measurable
par_meas = ['alive','alive']
pop_meas = ['Prisoners_males','Prisoners_PWID_males']
bounds = [(0, None),(0, None)]

ut.optimise_y_factor(country, cal_folder, parset, P, D, par_adj, pop_adj, par_meas, pop_meas, bounds,transfer='inc') 
P = ut.project(country, load_calibration=True, cal_folder=cal_folder, load_programs=False)
cal = P.parsets.copy('default','Calibrated')

par_adj = ['18-64_males']
pop_adj = ['PWID_males']
# measurable
par_meas = ['alive']
pop_meas = ['PWID_males']
bounds = [(0, None)]

ut.optimise_y_factor(country, cal_folder, parset, P, D, par_adj, pop_adj, par_meas, pop_meas, bounds,transfer='idu_status') 
#%%
# Generate databook
# generate_databook(country)
#%%
# Run calibration
# ut.run_calibration(country, savedir=cal_folder)
ut_plt.plot_calibration(country, cal_folder=cal_folder, savedir=savedir_calib) # plot

# Run vaccine scenarios
# ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens)
# ut.econ_eval(country, savedir_scens, results_folder, rand_seed=rand_seed, n_samples=n_samples)

# Run sensitivity analyses
# ut.run_sensitivity_analyses(country, cal_folder, sens_folder, results_folder)
# ut.run_genotype_analyses(country, cal_folder, sens_folder, results_folder)
# ut.run_coverage_analyses(country, cal_folder, sens_folder)
# ut.econ_eval_central(country, sens_folder)
