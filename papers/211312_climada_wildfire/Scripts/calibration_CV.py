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
import logging


from climada.hazard.wildfire import WildFire

from climada.hazard import Centroids
from climada.entity.exposures.base import Exposures
from climada.entity.exposures.litpop import LitPop

from climada.entity.impact_funcs.wildfire import IFWildfire
from climada.entity.impact_funcs import ImpactFuncSet
from climada.engine import Impact

from matplotlib import colors
import numpy as np
import pandas as pd
import geopandas as gp
import pickle
from itertools import compress
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, KFold
import skopt
from skopt import dummy_minimize, gbrt_minimize, gp_minimize, forest_minimize
from skopt.plots import plot_convergence, plot_evaluations, plot_objective

LOGGER = logging.getLogger(__name__)

from sklearn.base import BaseEstimator, RegressorMixin

%matplotlib inline

''' ------------- LOAD DATA -----------------------------------------------'''
# calibration was performed for all combinations of hazard on exposure.
# Data folder
DP = "/path/to/data/"

# Read events file
events = pd.read_csv(os.path.join(DP, 'EMDAT_WF_MATCHED_COMPARE_RESOLUTION.csv'))
events['MODIS_ID'] = events['MODIS_ID'].apply(str)

# Load list of exposure
with open(os.path.join(DP, 'all_exp_LitPop_admin1_30arc_compare_new.pkl'), 'rb') as f:
    exp_prov_level = pickle.load(f)


# Load list of hazards
with open(os.path.join(DP, 'all_wf_1km_res_compare.pkl'), 'rb') as f:
    haz = pickle.load(f)
    
''' ------------- Functions & class for calibration -----------------------'''
class CalibratorWildFire(BaseEstimator, RegressorMixin):
    
    def __init__(self, haz_list, exp_list, i_half=310):
        #include all default parameters here    
        self.haz_type = 'WFsingle'
        self.i_half = i_half
        self.haz_list = haz_list
        self.exp_list = exp_list

    
    def fit(self):
        # MDD Emanuel type function with one degree of freedom
        
        self.if_wf = IFWildfire()
        self.if_wf.set_default_FIRMS(self.i_half)
        self.if_set = ImpactFuncSet()
        self.if_set.append(self.if_wf)
        
        # calc all damages
        estimated_damage = np.zeros(len(self.haz_list))
        # calc impact
        for idx, (haz, exp) in enumerate(zip(self.haz_list, self.exp_list)):

            exp.check()
            
            imp_wf = Impact()
            imp_wf.calc(exp, self.if_set, haz)
            estimated_damage[idx] = np.sum(imp_wf.at_event)
    
        self.estimated_damage = estimated_damage
        self.fit = estimated_damage

    def predict(self):
        # MDD Emanuel type function with one degree of freedom
        self.if_wf = IFWildfire()
        self.if_wf.set_default_FIRMS(self.i_half)
        self.if_set = ImpactFuncSet()
        self.if_set.append(self.if_ex)
        
        # calc all damages
        estimated_damage = np.zeros(len(self.haz_list))
        # calc impact
        for idx, (haz, exp) in enumerate(zip(self.haz_list, self.exp_list)):

            exp.check()
            
            imp_wf = Impact()
            imp_wf.calc(exp, self.if_set, haz)
            estimated_damage[idx] = np.sum(imp_wf.at_event)
    
        self.estimated_damage = estimated_damage
        self.fit = estimated_damage
        
    
def CF_RMSF(y_true, y_est):
    assert(len(y_est)==len(y_true))
    a = 0
    # add 1 to make sure that there are no 0 (leads to -inf)
    y_est=y_est+1.
    
    for i in range(0,len(y_true)):
        a = a + np.log(y_est[i]/y_true[i])**2

    error = np.exp(np.sqrt(a/len(y_true)))
    return(error)



# function for skopt input
HPO_PARAMS = {'n_calls':    500,
              'n_random_starts':25,
              'base_estimator':'ET',
              'acq_func':'EI',
              'xi':0.02,
              'kappa':1.96,
              'n_points':10000,
             }

SPACE = [skopt.space.Real(295, 1000, name='i_half', prior='uniform')]

@skopt.utils.use_named_args(SPACE)
def opt_climada_RMSF(i_half):
    
    # init calibrator
    calibrator = CalibratorWildFire(haz_list=haz_l, exp_list=exp_l,
                                      i_half=i_half)

    # run climada
    calibrator.fit()
    
    y_est = calibrator.estimated_damage
    # calc error score
    err = CF_RMSF(y, y_est)
    
    return err

''' ------------- Start script -------------------------------------------'''
# clean data -> exclude events that don't cause damage in CLIMADA
d = CalibratorWildFire(haz, exp_prov_level,i_half=523.8)
d.fit()

exclude = np.array(d.estimated_damage!=0)

events_clean = events[exclude]
haz_clean = list(compress(haz, exclude))
exp_clean = list(compress(exp_prov_level, exclude))

haz_l = haz_clean
exp_l = exp_clean
y = events_clean.tot_damage_inflated.values

# set up calibration and cross validation

# lists for saving params
i_half_CV = []
score_training = []
score_testing = []

# in order for opt_climada_RMSF to work, hazard (haz_l) and exposure (exp_l) 
# list as well as inflated damages (y) need to be named accordingly

# run optimisation algorithm with all data
res = gbrt_minimize(opt_climada_RMSF, SPACE, **HPO_PARAMS)

# save parameter and scores (same for test and train, as not differentiated)
i_half_CV.append(res.x[0])
score_training.append(res.fun)
score_testing.append(res.fun)

# initiated CV
n_splits = 10
cv = KFold(n_splits=n_splits, shuffle=True)
for train_index, test_index in cv.split(events_clean.tot_damage_inflated):
    
    # train set
    haz_l = list(map(haz_clean.__getitem__, train_index))
    exp_l = list(map(exp_clean.__getitem__, train_index))
    
    y = events_clean.tot_damage_inflated.values[train_index]
     
    # run optimisation with training data           
    res = gbrt_minimize(opt_climada_RMSF, SPACE, **HPO_PARAMS)
    i_half_CV.append(res.x[0])
    score_training.append(res.fun)
    
    # test set
    haz_test = list(map(haz_clean.__getitem__, test_index))
    exp_test = list(map(exp_clean.__getitem__, test_index))
    y_test = events_clean.tot_damage_inflated.values[test_index]
    
    # calc score for test data 
    calc_test = CalibratorWildFire(haz_test, exp_test, i_half=res.x[0])
    calc_test.fit()
    score_testing.append(CF_RMSF(y_test, calc_test.estimated_damage))

# one more calc with best parameter of each CV
estimated_damages_CV = np.zeros((len(events.tot_damage_inflated), n_splits+1))

for i in range(n_splits+1):
    d = CalibratorWildFire(haz, exp_prov_level,i_half=i_half_CV[i])
    d.fit()
    estimated_damages_CV[:,i] = d.estimated_damage
    
# order an save results
res_CV = np.vstack((np.array(i_half_CV), np.array(score_training),
                 np.array(score_testing)))
results_CV = pd.DataFrame(res_CV.transpose(),
                          columns= ['i_half', 'train_RMSF', 'test_RMSF'])
results_CV.to_csv(os.path.join(DP, 'results_CV.csv'))

estimated_damages_CV = pd.DataFrame(estimated_damages_CV)
estimated_damages_CV.rename(columns={estimated_damages_CV.columns[0]:
                                     "estimated_damage" }, inplace = True)
estimated_damages_CV.to_csv(os.path.join(DP, 'estimated_damage.csv'))
