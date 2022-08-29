#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:31:14 2021

@author: ckropf
"""

from pathlib import Path

#
run_name = Path(__file__).with_suffix('').name

#Define exposures, hazard
haz_file = "HAZ_ts_VNM_2020.hdf5"

fun_id = 3
population = 97330000.0
litpop_kwargs = {
    'countries' : ['VNM'],
    'res_arcsec' : 30,
    'reference_year' : 2020,
    'fin_mode' : 'norm',
    'total_values' : [population]
    }
value_unit = 'people'
assign_centr_kwargs = {}

#define uncertainty bounds
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

bounds_mdd = None
bounds_impfi = [-0.5, 2]

bounds_int = None #for cc
bounds_freq = None #for cc
n_ev = 1 #percent of total number of events

#Define uncertainty Salib configuration
#Sampling
N_samples = 2**10
sampling_method = "saltelli"
sampling_kwargs = {'skip_values': 2**11}
#Sensitivity
sensitivity_method = "sobol"
sensitivity_kwargs = None
#Dsitribution
rp = None
calc_eai_exp = True
calc_at_event = True
