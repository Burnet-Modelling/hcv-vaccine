# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:53:37 2026

@author: farah.houdroge

Script to run convergence analysis on sample size (number of simulations).
"""
from hcv import utils as ut
from hcv.convergence_checks import convergence_check


dir = ut.get_project_root()

country = "VNM"
rand_seed = 250711
n_samples = 1000
results_folder = dir / "results"
cal_folder = results_folder / "calibration" / "y-factors"
savedir_scens = results_folder / "convergence" / country / f"{n_samples}_samples"
savedir_scens.mkdir(parents=True, exist_ok=True)
savedir_pp = results_folder / "convergence" / "results" / country / f"{n_samples}_samples"
savedir_pp.mkdir(parents=True, exist_ok=True)

scens_folder = results_folder / "convergence"
#%%
ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens,scenarios_to_run=["counter_0", "scenario_0"])
ut.econ_eval(country, savedir_scens, savedir_pp, rand_seed=rand_seed, n_samples=n_samples,scenarios_to_run=["counter_0", "scenario_0"])

#%%
pkl_file = scens_folder / f"{country}_econ_eval.pkl"
# save_file = scens_folder / f"{country}_convergence_check.xlsx"

convergence_check(pkl_file, scens_folder, ["counter_0"])