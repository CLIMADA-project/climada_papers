#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of CLIMADA.

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---
This script was used to generate results as displayed in the paper by
Lüthi et al. (2021). It depends on the open-source and -access model
CLIMADA available at https://github.com/CLIMADA-project/climada_python

Please refer to details on installation and model set-up to
https://climada-python.readthedocs.io/en/v0.1.0/install.html

@author: Samuel Luethi - samuel.luethi@usys.ethz.ch
"""
import os

from climada.hazard.wildfire import WildFire

from climada.hazard import Centroids
from climada.util.constants import ONE_LAT_KM
from climada.util.dates_times import date_to_str
from climada.entity.exposures.base import Exposures
from climada.entity.exposures.litpop import LitPop
from climada.entity.exposures.litpop import exposure_set_admin1

import numpy as np
import pandas as pd
import pickle

''' ------------- LOAD DATA -----------------------------------------------'''
# Data folder
DP = "/path/to/data/"

# Read events file
events = pd.read_csv(os.path.join(DP, 'EMDAT_WF_MATCHED_COMPARE_RESOLUTION.csv'))
events['MODIS_ID'] = events['MODIS_ID'].apply(str)

# function to check if province names are correctly written (utf-8)
def check_provinces(list_prov_event, list_prov_all):
    result =  all(elem in list_prov_all  for elem in list_prov_event)
    if result:
        print('prov ok')
    else:
        print('issue with provinces')

''' ------------- Prepare exposure ----------------------------------------'''
# Prepare exposure on Province level 300 arcsec
exp_prov_level = []
for i, ctry in enumerate(events.ISO.values):
    
    if events.ISO.values[i-1]!=ctry:
        print('Setting up - ' + ctry)
        exp_ctry = LitPop()
        exp_ctry.set_country(ctry, res_arcsec=300, reference_year=2019)
        print('Setting admin1')
        exposure_set_admin1(exp_ctry, 300) # function is contained within LitPop
    try:
        provinces = np.char.split(events.Admin1.values[i],sep=",").item()
        check_provinces(provinces, exp_ctry.gdf.admin1.unique().tolist())
        exp_prov = exp_ctry.gdf[exp_ctry.gdf['admin1'].str.contains('|'.join(provinces), na=False)]
        exp = Exposures(exp_prov)
        exp_prov_level.append(exp)
    except TypeError:
        print('ISSUE')
        exp_prov_level.append(exp_ctry)
        
    print('-------------- done:')
    print(i)
    
with open(os.path.join(DP, 'all_exp_LitPop_admin1_300arc_compare_new_exp.pkl'), 'wb') as f:
    pickle.dump(exp_prov_level, f)

# Prepare exposure on Province level 120 arcsec
exp_prov_level = []
for i, ctry in enumerate(events.ISO.values):
    
    if events.ISO.values[i-1]!=ctry:
        print('Setting up - ' + ctry)
        exp_ctry = LitPop()
        exp_ctry.set_country(ctry, res_arcsec=120, reference_year=2019)
        print('Setting admin1')
        exposure_set_admin1(exp_ctry, 120) # function is contained within LitPop
    try:
        provinces = np.char.split(events.Admin1.values[i],sep=",").item()
        check_provinces(provinces, exp_ctry.admin1.unique().tolist())
        exp_prov = exp_ctry[exp_ctry['admin1'].str.contains('|'.join(provinces), na=False)]
        exp = Exposures(exp_prov)
        exp_prov_level.append(exp)
    except TypeError:
        print('ISSUE')
        exp_prov_level.append(exp_ctry)
        
    print('-------------- done:')
    print(i)
    
with open(os.path.join(DP, 'all_exp_LitPop_admin1_120arc_compare.pkl'), 'wb') as f:
    pickle.dump(exp_prov_level, f)

# Prepare exposure on Province level 30 arcsec
# you might require cluster computing for this part
exp_prov_level = []
for i, ctry in enumerate(events.ISO.values):
    
    if events.ISO.values[i-1]!=ctry:
        print('Setting up - ' + ctry)
        exp_ctry = LitPop()
        exp_ctry.set_country(ctry, res_arcsec=30, reference_year=2019)
        print('Setting admin1')
        exposure_set_admin1(exp_ctry, 30)
    try:
        provinces = np.char.split(events.Admin1.values[i],sep=",").item()
        check_provinces(provinces, exp_ctry.gdf.admin1.unique().tolist())
        exp_prov = exp_ctry.gdf[exp_ctry.gdf['admin1'].str.contains('|'.join(provinces), na=False)]
        exp = Exposures(exp_prov)
        exp_prov_level.append(exp)
    except TypeError:
        print('ISSUE')
        exp_prov_level.append(exp_ctry)
        
    print('-------------- done:')
    print(i)
    
with open(os.path.join(DP, 'all_exp_LitPop_admin1_30arc_compare_new.pkl'), 'wb') as f:
    pickle.dump(exp_prov_level, f)

''' ------------- Prepare hazard ----------------------------------------'''
# prepare hazard
# calcs for 10 km hazard (300 arcsec)
all_wf = []
for i in range(events.shape[0]):
    
    # file locations
    folder_name_haz = 'MODIS/DL_FIRE_M6_' + events.MODIS_ID.values[i] + '_' + events.identifier.values[i]
    file_name_haz = 'fire_archive_M6_' + events.MODIS_ID.values[i] + '.csv'
    
    HAZ_DIR = os.path.join(DP, folder_name_haz, file_name_haz)
    
    # load and cut FIRMS to save computing time
    exp = exp_prov_level[i]
    FIRMS = pd.read_csv(HAZ_DIR)
    FIRMS = FIRMS[FIRMS['latitude']<exp.latitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['latitude']>exp.latitude.min()-1.0]
    FIRMS = FIRMS[FIRMS['longitude']<exp.longitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['longitude']>exp.longitude.min()-1.0]
    # Hazard calcs
    wf = WildFire()
    wf.set_hist_fire_FIRMS(FIRMS, centr_res_factor=1/10.)
    wf.combine_fires()
    
    all_wf.append(wf)
    print('-------------- done:')
    print(i)

with open(os.path.join(DP, 'all_wf_10km_res_compare.pkl'), 'wb') as f:
    pickle.dump(all_wf, f)

# calcs for 4 km hazard (120 arcsec)
all_wf = []
for i in range(events.shape[0]):
    
    # file locations
    folder_name_haz = 'MODIS/DL_FIRE_M6_' + events.MODIS_ID.values[i] + '_' + events.identifier.values[i]
    file_name_haz = 'fire_archive_M6_' + events.MODIS_ID.values[i] + '.csv'
    
    HAZ_DIR = os.path.join(DP, folder_name_haz, file_name_haz)
    
    # load and cut FIRMS to save computing time
    exp = exp_prov_level[i]
    FIRMS = pd.read_csv(HAZ_DIR)
    FIRMS = FIRMS[FIRMS['latitude']<exp.latitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['latitude']>exp.latitude.min()-1.0]
    FIRMS = FIRMS[FIRMS['longitude']<exp.longitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['longitude']>exp.longitude.min()-1.0]
    # Hazard calcs
    wf = WildFire()
    wf.set_hist_fire_FIRMS(FIRMS, centr_res_factor=1/4.)
    wf.combine_fires()
    
    all_wf.append(wf)
    print('-------------- done:')
    print(i)

with open(os.path.join(DP, 'all_wf_4km_res_compare.pkl'), 'wb') as f:
    pickle.dump(all_wf, f)


# calcs for 1 km hazard (30 arcsec)  
all_wf = []
for i in range(events.shape[0]):
    # file locations
    folder_name_haz = 'MODIS/DL_FIRE_M6_' + events.MODIS_ID.values[i] + '_' + events.identifier.values[i]
    file_name_haz = 'fire_archive_M6_' + events.MODIS_ID.values[i] + '.csv'
    
    HAZ_DIR = os.path.join(DP, folder_name_haz, file_name_haz)
    
    # load and cut FIRMS to save computing time
    exp = exp_prov_level[i]
    FIRMS = pd.read_csv(HAZ_DIR)
    FIRMS = FIRMS[FIRMS['latitude']<exp.latitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['latitude']>exp.latitude.min()-1.0]
    FIRMS = FIRMS[FIRMS['longitude']<exp.longitude.max()+1.0]
    FIRMS = FIRMS[FIRMS['longitude']>exp.longitude.min()-1.0]
    # Hazard calcs
    wf = WildFire()
    wf.set_hist_fire_FIRMS(FIRMS, centr_res_factor=1/1.)
    wf.combine_fires()
    
    all_wf.append(wf)
    print('-------------- done:')
    print(i)

with open(os.path.join(DP, 'all_wf_1km_res_compare.pkl'), 'wb') as f:
    pickle.dump(all_wf, f)
