import time
import os

import sys

in_path = '../'
sys.path.append(in_path)

loc = 'C:/Users/Phil/Documents/GitHub/atomica-tools/'
sys.path.append(loc)

import atomica as at
import at_tools as att
import matplotlib.pyplot as plt
from os import sep
import numpy as np
import sciris as sc
import yaml

import datetime
import time
import shutil
import winsound
out_path = '../autocalibration_tool'
sys.path.append(out_path)

#run settings
testing = False
dis_model = 'hepc'
db_loc = 'databooks_revised'

ctry_list = ['Australia', 'Georgia'] # ['Armenia', 'Australia' , 'Georgia', 'Myanmar', 'Pakistan', 'Tanzania']

for country in ctry_list:
    print(f'\n\nCurrent country: {country} --------------------------------')
    F = at.ProjectFramework(f"{in_path}/framework/hcv_framework.xlsx") # Edit Framework
    D = at.ProjectData.from_spreadsheet(f"{in_path}/{db_loc}/databook_{country}.xlsx", framework=F) # Edit databook
    # D.rename_pop('18–64', '18-64', 'General Population 18-64')
    P = at.Project(framework=F, databook=D, do_run=False)
    P.settings.update_time_vector(start=2000, end=2031, dt=1 / 10)
    cal = P.make_parset()

    config_file_path = in_path + '/yaml_calibrations/model_config_HepC_v2_1_1.yaml' # Edit YAML to specific location on drive
    config = att.ProjectConfig.from_yaml(config_file_path)

    # make new folder to put calibrations in
    date = time.strftime("%Y-%m-%d_%H")
    calspath = f'{out_path}/{dis_model}_YAML_autocalibrations_nocal_{date}'
    savecals = calspath
    if testing:
        calspath += '_TESTING'
        savecals = False
    else:
        #save everything in folder
        os.makedirs(calspath)
        # os.chdir(calspath)

        # save yaml file for reference
        fname = 'model_config_HepC_v2_1.yaml'
        new_fname = f'{country}_{dis_model}_YAML-config_{date}.yaml'
        shutil.copy(config_file_path, calspath)
        #save logfile
        at.start_logging(f'{calspath}/logfile_{date}.txt')

    #run calibration w yaml instructions
    newcal = config.run_calibration(P,
                                    parset=cal,
                                    save_intermediate_cals=savecals)

    #save final calibration
    date = time.strftime("%Y-%m-%d_%H%M")
    newcal.save_calibration(f'{calspath}/{country}_{dis_model}_YAML-autocalibration')
    at.stop_logging()

winsound.Beep(frequency = 2500, duration = 200)
winsound.Beep(frequency = 2750, duration = 200)
winsound.Beep(frequency = 3050, duration = 200)