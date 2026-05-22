# -*- coding: utf-8 -*-
"""
Created on Fri May 22 13:53:37 2026

@author: farah.houdroge

Script to run convergence analysis on sample size (number of simulations).
"""
from hcv.generate_databooks import generate_databook
from hcv import utils as ut
from hcv import utils_plotting as ut_plt
import sys

dir = ut.get_project_root()

country = "VNM"
rand_seed = 250711
n_samples = int(sys.argv[1])
results_folder = dir / "results"
cal_folder = results_folder / "calibration" / "y-factors"
savedir_scens = results_folder / "convergence" / country / n_samples
savedir_scens.mkdir(parents=True, exist_ok=True)
savedir_pp = results_folder / "convergence" / "results" / country / n_samples
savedir_pp.mkdir(parents=True, exist_ok=True)

ut.run_scenario_sampling(country, cal_folder, rand_seed=rand_seed, n_samples=n_samples, savedir=savedir_scens,nb_scenarios=1)
ut.econ_eval(country, savedir_scens, savedir_pp, rand_seed=rand_seed, n_samples=n_samples,nb_scenarios=1)