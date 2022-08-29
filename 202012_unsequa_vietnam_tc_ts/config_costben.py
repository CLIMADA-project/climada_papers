#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:31:14 2021

@author: ckropf
"""

from pathlib import Path
from entity_VNM import (measures_20, measures_50, disc_rates)

#
run_name = Path(__file__).with_suffix('').name

#Define exposures, hazard
haz_file = "HAZ_ts_VNM_2020.hdf5"
haz_fut_file = "HAZ_ts_VNM_2050_cc85.hdf5"
impf_id = 3

#Future year
horizon = 2050

#Define the measures and disc rates
meas_pres = measures_20()
meas_fut = measures_50()
disc = disc_rates()['PE']

assign_centr_kwargs = {}
population = 97330000.0
litpop_kwargs = {
    'countries' : ['VNM'],
    'res_arcsec' : 30,
    'reference_year' : 2020,
    'fin_mode' : 'norm',
    'total_values' : [population]
    }
value_unit = 'people'

#Exposures unc
bounds_totval = [0.9, 1.1]
mean_std = (1, 0.05)
choice_mn = [
    [0, 0.75],
    [0, 1],
    [0, 1.25],
    [0.5, 0.75],
    [0.5, 1],
    [0.5, 1.25],
    [1, 0.75],
    [1, 1],
    [1, 1.25],
    ]
bounds_eg = [1.1, 1.16]

#Impact function
bounds_mdd = None
bounds_impfi = [-0.5, 2]

bounds_int = [0.9, 1.1] #for cc
bounds_freq = [0.5, 2] #for cc
n_ev = 1 #percent of total number of events

#Adaptation bounds
bounds_disk = None
bounds_cost = [0.5, 2]


#Define uncertainty Salib configuration
#Sampling
N_samples = 2**10
sampling_method = "saltelli"
sampling_kwargs = {'skip_values': 2**11}
#Sensitivity
sensitivity_method = "sobol"
sensitivity_kwargs = None

