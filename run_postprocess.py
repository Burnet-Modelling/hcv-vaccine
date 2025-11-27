# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 11:37:29 2025

@author: farah.houdroge

Post-process results from calibration, scenarios and sensitivity analyses.
"""

from hcv import utils as ut
from hcv import utils_plotting as ut_plt

dir = ut.get_project_root()
#%%
# Simulation parameters
n_samples = 100
results_folder = dir / "results"
scens_folder = results_folder / "scenarios"
cal_folder = results_folder / "calibration"
sens_folder = results_folder / "sensitivity_analyses"

fig_folder = results_folder / "figures"
fig_folder.mkdir(parents=True, exist_ok=True)
# %% pkl files

## Aggregate epi and econ pkl outcomes
# requires "scens_folder/agg_outputs" from econ_eval
ut.aggregate_ensembles(scens_folder, n_samples)

# ## Outcomes by region
# requires "scens_folder/econ_agg.pkl" from aggregate_ensembles
ut.econ_analysis(scens_folder, n_samples)
ut.write_print_table(scens_folder)

#%% calibration validation
ut.calibration_outputs(cal_folder,scens_folder)

#%% Sensitivity analyses
ut.aggregate_ensembles_central(sens_folder,results_folder)
ut.econ_analysis_central(sens_folder)
ut.write_sensitivity_table(scens_folder,sens_folder)

#%% BCR correlation
ut.bcr_correlation()

#%% plots
## Impact bar panel plot
# requires "scens_folder/epi_agg.pkl" from aggregate_ensembles
ut_plt.plot_bars_impact(scens_folder) # plot all regions on one page (for manuscript)

for regions in [['global','top10'],['global'],['top10'],['AFR','EMR'],['EUR','AMR'],['SEAR','WPR']]:
    ut_plt.plot_bars_impact(scens_folder,regions=regions) # plot select region combinations (for slides)
for regions in [['top10','AFR','EMR','EUR','AMR','SEAR','WPR']]:
    ut_plt.plot_bars_impact(scens_folder,regions=regions)

## BCR maps
# requires "scens_folder/bcr_map.csv" from econ_analysis
ut_plt.plot_bcr_map(scens_folder)
#%%
## Time-series panel plot of outcomes
for scenarios in [['scenario_0'],['scenario_1'],['scenario_2','scenario_3','scenario_4']]: # bundle scenarios with same counterfactual
    ut_plt.plot_outcomes_timeseries(scens_folder,scenarios=scenarios)
#%%
ut_plt.plot_calibration_forest(cal_folder)

# %%
ut_plt.plot_calibration_panel(scens_folder, cal_folder)
