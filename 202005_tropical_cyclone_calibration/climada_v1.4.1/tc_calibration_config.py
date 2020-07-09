#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the CLIMADA papers repository:
    https://github.com/CLIMADA-project/climada_papers

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Created on Wed Jul  1 16:21:10 2020

description: Configuration script for tropical cyclone (TC) impact function calibration.
    You need so set some paths and parameters here before running tc_calibration_main.py.
    Please refer to the paper (submitted) and README for more information.
    Not required for the Jupyer notebook.

@author: Samuel Eberenz
"""

import os
import numpy as np

from climada.util.constants import SYSTEM_DIR

### Configuration of calibration routine:
CALIB = [1, 2, 3, 4, 5, 6] # List of calibration steps to be run (default: [1, 2, 3, 4, 5, 6])
HAZ = 'TC' # hazard type abbreviation (CLIMADA convention), default: 'TC'
v_step = .1 # step length of velocity parameters in m/s. set larger steps to improve performance. Default=.1
v0 = [25.7] # IMPACT FUNCTION PARAMETER SPACE FOR CALIBRATION: V_threshold in m/s
scale = [1.] # IMPACT FUNCTION PARAMETER SPACE FOR CALIBRATION: linear scaling along y-axis
v_half_range = np.arange(v0[0]+.1, v0[0]+300.1, v_step) # IMPACT FUNCTION PARAMETER SPACE FOR CALIBRATION: V_half in m/s

""" CALIB:
    1:  Loading or initiating HAZARD and EXPOSURE sets (required for 2, 3, 5, 6, 7)
    2:  With 1 global impact functions, calculate EDR. Extract NRD from EM-DAT.
        Compute ratio EDR=SED/NRD. Save to CSV file.
        Required for CALIB 3.
    3:  Core calibration. Loop over v_half_range to compute impact and ratio for each matched event/country.
        Required for CALIB 4.
    4:  Optimization: Compute cost functions and find optimized v_half for each cost function and region,
        Save calibration results to CSV.
    5:  Compute annual average damage (AAD) per country for comparison with GAR 2013 
        (calibrated and uncalibrated)
    6:  Compute damage time series, trend and significance, as well as
        standard deviation of annual damage for EM-DAT and CLIMADA
        per country and region.
"""

## Set directories (different on cluster):
if 'cluster' in SYSTEM_DIR:
    on_cluster = True
    DATA_DIR = '/<SET PATH>' # set data directory on cluster
    HAZARD_DIR = DATA_DIR + '/tc_tracks'
    TRACK_DIR = HAZARD_DIR
    ENTITY_DIR = DATA_DIR + '/litpop'
    RES_DIR = DATA_DIR + '/calibration/20191203'
else:
    on_cluster = False
    DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../data")
    HAZARD_DIR = os.path.join(SYSTEM_DIR, 'hazard') # hazard data directory
    TRACK_DIR = DATA_DIR # os.path.join(SYSTEM_DIR, 'data_master_sl_short') #set tc track data directory
    ENTITY_DIR = os.path.join(SYSTEM_DIR, 'litpop') # exposure data directory
    RES_DIR = os.path.join(DATA_DIR, 'results') #!!! results. Change according to your file structure.

TRACK_FOLDER = 'ibtracs_calibration' #!!! subfolder in TRACK_DIR with subset of tracks used. IBTrACS need to be downloaded seperately from data provider.
EMDAT_CSV = os.path.join(DATA_DIR, '201810_emdat_TC_1980-2017.csv') #!!! Change accoridngly. EM-DAT data needs to be downloaded seperately from emdat.be
EMDAT_MAP = os.path.join(DATA_DIR, '202002_TC_event_assignment.csv')

### Basic configuration of model setup:
version_date = '2019-04-05' # set reference date for versioning
INIT_HAZARD = True # set to TRUE for first init of hazard set, later to FALSE  to reuse saved hazard set
REF_YEAR = 2014 # Reference year used to initiate exposure data and for damage normalization, default: 2014
YEAR_RANGE = [1980, 2017] # Range of years to be considered for calibration and analysis (default: [1980, 2017])
RES_ARCSEC = 300 # Resoultion of hazard and exposure in arc seconds, default: 300 arcsec
dist_cst_lim = 1500000 # distance to coast threshold in meteres: data further from coast is discarded for performance reasons
keep_variables = True # keep exposure and hazard variables once initiated?

### Define regions
basins_short = ['NA', 'NI', 'OC', 'SI', 'WP']
regions_short = ['NA1', 'NA2', 'NI', 'OC', 'SI', 'WP1', 'WP2', 'WP3', 'WP4']
regions_short = np.sort(regions_short).tolist()
regions_long = dict()
regions_long[regions_short[0]] = 'North Atlantic 1' # NA1; 'Caribbean and Mexico'  # CAM
regions_long[regions_short[1]] = 'North Atlantic 2' # NA2; 'USA and Canada' # NAM
regions_long[regions_short[2]] = 'North Indian' # NI
regions_long[regions_short[3]] = 'Oceania' # OC, 'Southern Hemisphere' # SHE
regions_long[regions_short[4]] = 'South Indian' # SI, 'Southern Hemisphere' # SHE
regions_long[regions_short[5]] = 'West Pacific 1' # WP1; 'South East Asia' # SEA
regions_long[regions_short[6]] = 'West Pacific 2' # WP2; 'Philippines' # PHL
regions_long[regions_short[7]] = 'West Pacific 3' # WP3; 'China' # CHN
regions_long[regions_short[8]] = 'West Pacific 4' # WP4; 'North West Pacific' # NWP
regions_long['all'] = 'Global'
regions_long['GLB'] = 'Global'
# countries by region:
region_ids_cal = {'NA1': ['AIA', 'ATG', 'ARG', 'ABW', 'BHS', 'BRB', 'BLZ', 'BMU', 'BOL', 'CPV', 'CYM', 'CHL', 'COL', 'CRI', 'CUB', 'DMA', 'DOM', 'ECU', 'SLV', 'FLK', 'GUF', 'GRD', 'GLP', 'GTM', 'GUY', 'HTI', 'HND', 'JAM', 'MTQ', 'MEX', 'MSR', 'NIC', 'PAN', 'PRY', 'PER', 'PRI', 'SHN', 'KNA', 'LCA', 'VCT', 'SXM', 'SUR', 'TTO', 'TCA', 'URY', 'VEN', 'VGB', 'VIR'], \
                  'NA2': ['CAN', 'USA'], \
                  'NI': ['AFG', 'ARM', 'AZE', 'BHR', 'BGD', 'BTN', 'DJI', 'ERI', 'ETH', 'GEO', 'IND', 'IRN', 'IRQ', 'ISR', 'JOR', 'KAZ', 'KWT', 'KGZ', 'LBN', 'MDV', 'MNG', 'MMR', 'NPL', 'OMN', 'PAK', 'QAT', 'SAU', 'SOM', 'LKA', 'SYR', 'TJK', 'TKM', 'UGA', 'ARE', 'UZB', 'YEM'], \
                  'OC': ['ASM', 'AUS', 'COK', 'FJI', 'PYF', 'GUM', 'KIR', 'MHL', 'FSM', 'NRU', 'NCL', 'NZL', 'NIU', 'NFK', 'MNP', 'PLW', 'PNG', 'PCN', 'WSM', 'SLB', 'TLS', 'TKL', 'TON', 'TUV', 'VUT', 'WLF'], \
                  'SI': ['COM', 'COD', 'SWZ', 'MDG', 'MWI', 'MLI', 'MUS', 'MOZ', 'ZAF', 'TZA', 'ZWE'], \
                  'WP1': ['KHM', 'IDN', 'LAO', 'MYS', 'THA', 'VNM'], \
                  'WP2': ['PHL'], \
                  'WP3': ['CHN'], \
                  'WP4': ['HKG', 'JPN', 'KOR', 'MAC', 'TWN'], 
                  'ROW': ['ALB', 'DZA', 'AND', 'AGO', 'ATA', 'AUT', 'BLR', 'BEL', 'BEN', 'BES', 'BIH', 'BWA', 'BVT', 'BRA', 'IOT', 'BRN', 'BGR', 'BFA', 'BDI', 'CMR', 'CAF', 'TCD', 'CXR', 'CCK', 'COG', 'HRV', 'CUW', 'CYP', 'CZE', 'CIV', 'DNK', 'EGY', 'GNQ', 'EST', 'FRO', 'FIN', 'FRA', 'ATF', 'GAB', 'GMB', 'DEU', 'GHA', 'GIB', 'GRC', 'GRL', 'GGY', 'GIN', 'GNB', 'HMD', 'VAT', 'HUN', 'ISL', 'IRL', 'IMN', 'ITA', 'JEY', 'KEN', 'PRK', 'XKX', 'LVA', 'LSO', 'LBR', 'LBY', 'LIE', 'LTU', 'LUX', 'MLT', 'MRT', 'MYT', 'MDA', 'MCO', 'MNE', 'MAR', 'NAM', 'NLD', 'NER', 'NGA', 'MKD', 'NOR', 'PSE', 'POL', 'PRT', 'ROU', 'RUS', 'RWA', 'REU', 'BLM', 'MAF', 'SPM', 'SMR', 'STP', 'SEN', 'SRB', 'SYC', 'SLE', 'SGP', 'SVK', 'SVN', 'SGS', 'SSD', 'ESP', 'SDN', 'SJM', 'SWE', 'CHE', 'TGO', 'TUN', 'TUR', 'UKR', 'GBR', 'UMI', 'ESH', 'ZMB', 'ALA']}

### Criteria for automatic and manual exclusion of countries / events:
def_log_ratio_threshold = np.log(1000) # threshold of log(EDR) to exclude events from calibration
exclude_cntries = [] # countries to be excluded manually from calibration
drop_excluded_cntries = True # drop excluded countries from calibration? Default=True
mismatched_emdat_events = ['1985-0124', '1985-0375', '1995-0247', '2000-0788', '2013-0515'] # mismatched events to be excluded manually from calibration
new_reg_names = True # for backward compatibility: update regions in dataframe?
trunc_rus = True # remove Russian exposure west of 100E for better performance

### Plotting settings:
make_plots = True
if on_cluster:
    make_plots = False
colors8 = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf']
colors8b = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','black','#a65628','#f781bf','#ffff33']
colors9b = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf', '#252525']
colors9 = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#999999','#a65628','#f781bf','#ffff33']

### Strings for file names when written to hard drive ashdf5
ENTITY_STR = "litpop_%04das_%04d_%s_iter_coast_" + version_date + ".hdf5" # exposure files
    # (RES_ARCSEC, REF_YEAR, HAZ + '_emdat')
HAZARD_STR = "TC_GLB_%04das_%03d.hdf5" # (RES_ARCSEC, 0) # hazard files
RES_STR = "Cal_r2_%04das_%04d_%s_%s_v" + version_date +"_%s.%s" # Result files
track_bin_size = 0 # amount of tracks processed at once, set to 0 for all.
