# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 10:35:13 2024

@author: kelly.maynard
"""
import os
import socket
from platform import system

def _get_sharepoint_folder():

    
    user = socket.gethostname()
    platform = system() 


    if user in ['KellyM-Mar22']:
        folder = 'C:/Users/kelly.maynard/Burnet Institute/WG-Modelling-Hepatitis C - Documents/Applications/Multi-country/'
    elif user in ['JessicaZ-Nov23']:
        folder = 'C:/Users/jessica.zuk/Burnet Institute/WG-Modelling-Hepatitis C - Multi-country/'
    elif user in ['PhilL-Desktop-May23']:
        folder = 'C:/Users/Phil/Burnet Institute/WG-Modelling - HCV vaccine model - Documents/HCV vaccine model/HCV vaccine model/'
    elif user in ['phillipluong']:
        folder = 'C:/Users/iamph/Burnet Institute/WG-Modelling - HCV vaccine model - Documents/HCV vaccine model/HCV vaccine model/'
    elif user in ['FarahH-May20']:
        folder = 'C:/Users/farah.houdroge/Burnet Institute/WG-Modelling - HCV vaccine model - HCV vaccine model/HCV vaccine model'
    elif user in ['NickS-Sep24']:
        folder = 'C:/Users/nick.scott/Burnet Institute/WG-Modelling - HCV vaccine model - HCV vaccine model/'
    elif user in ['Chris-Desktop']:
        folder = 'C:/Users/seama/Documents/Burnet Institute/WG-Modelling - HCV vaccine model - Documents/HCV vaccine model/HCV vaccine model/'
    elif user in ['CSeaman-Oct18']:
        folder = 'C:/Users/chris.seaman/Burnet Institute/WG-Modelling - HCV vaccine model - Documents/HCV vaccine model/HCV vaccine model/'
    elif user in ['Your username from socket.gethostname()']:
        folder = 'Your hcv-multicountry file path' 
    else:
        if platform == 'Linux':
            folder = '/home/farah.houdroge/HCV-vaccine-model/'
        else:
            raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')

    return os.path.join(os.path.abspath(folder),'')


def _get_gitlab_folder():

    user = socket.gethostname()
    platform = system() 

    if user in ['FarahH-May20']:
        folder = 'C:/Users/farah.houdroge/Documents/GitHub/hcv-vaccine/'
    elif user in ['TharinduW-Oct20']:
        folder = 'D:/Covid-India-UP/hcv-vaccine/'
    elif user in ['KellyM-Mar22']:
        folder = 'C:/Users/kelly.maynard/Documents/GitHub/hcv-vaccine/'
    elif user in ['NickS-Sep24']:
        folder = 'C:/Users/nick.scott/Desktop/GitHub/hcv-vaccine/'
    elif user in ['TomTidhar-Jun18']:
        folder = 'C:/Users/tom.tidhar/Documents/GitHub/hcv-vaccine'
    elif user in ['Your username from socket.gethostname()']:
        folder = 'Your covasim-australia file path'
    elif user in ['JessicaZ-Nov23']:
        folder = 'C:/Users/jessica.zuk/Documents/GitHub/hcv-vaccine/'
    elif user in ['phillipluong']:
        folder = 'C:/Users/iamph/Documents/GitHub/hcv-vaccine/'
    elif user in ['PhilL-Desktop-May23']:
        folder = 'C:/Users/Phil/Documents/GitHub/hcv-vaccine/'
    elif user in ['Chris-Desktop']:
        folder = 'C:/Users/seama/Documents/GitHub/hcv-vaccine'
    elif user in ['CSeaman-Oct18']:
        folder = 'C:/Users/chris.seaman/Documents/hcv-vaccine'
    else:
        if platform == 'Linux':
            folder = '/home/farah.houdroge/hcv-vaccine/'
        else:
            raise Exception(f'Error: unknown user "{user}", please add user information for future convenience!')

    return os.path.join(os.path.abspath(folder),'')
