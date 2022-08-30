#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2021

@author: C. M. Kropf
"""

import importlib
import datetime as dt

import pandas as pd
import scipy as sp

from pathlib import Path
from pathos.pools import ProcessPool as Pool


from climada.hazard import Hazard
from climada.engine.unsequa import CalcImpact, InputVar #CLIMADA >= 3.1.0
#from climada.engine.uncertainty_quantification import CalcCostBenefit, InputVar #CLIMADA < 3.1.0

from entity_VNM import generate_litpop_base
from impf_VNM import impf


#Load configuration of uncertainty run

DATA_DIR = Path("VNM_Data/")
OUT_DIR = Path("./")

CONFIG_FILE  = "config_impact"


def main(config_file):

    config = importlib.import_module(config_file)

    uncertainty_filename = (OUT_DIR /
        ("uncertainty_" + str(config.run_name) + '_' +
         dt.datetime.now().strftime("%Y-%m-%d") + '.hdf5')
    )
    if uncertainty_filename.is_file():
        print("Data files for this run already exist. Abort")
        return None
    sensitivity_filename = (OUT_DIR /
        ("sensitivity_" + str(config.run_name) + '_' +
         dt.datetime.now().strftime("%Y-%m-%d") + '.hdf5')
    )

    haz = Hazard()
    haz.read_hdf5(DATA_DIR / config.haz_file)
    haz_type = haz.tag.haz_type

    impf_set = impf()

    n_ev = config.n_ev if config.n_ev is None else int(config.n_ev * haz.size)
    litpop_list = generate_litpop_base(
        config.fun_id, config.value_unit, haz, config.assign_centr_kwargs,
        config.choice_mn, **config.litpop_kwargs
        )
    exp_iv = InputVar.exp(litpop_list, bounds_totval=config.bounds_totval
        )

    myclip_a, myclip_b = config.bounds_totval
    my_mean, my_std = config.mean_std
    a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std
    exp_iv.distr_dict['ET'] = sp.stats.truncnorm(a=a, b=b, loc=my_mean, scale=my_std)

    haz_iv = InputVar.haz(haz, n_ev=n_ev, bounds_int=config.bounds_int,
                          bounds_freq=config.bounds_freq)

    impfset_iv = InputVar.impfset(
        impf_set, bounds_impfi=config.bounds_impfi,
        bounds_mdd=config.bounds_mdd, haz_id_dict={haz_type : [config.fun_id]})
    unc_calc = CalcImpact(exp_iv, impfset_iv, haz_iv)

    #Make sample
    unc_data = unc_calc.make_sample(
        N=config.N_samples,
        sampling_method=config.sampling_method,
        sampling_kwargs=config.sampling_kwargs
        )

    pool = Pool()
    #Compute distribution
    unc_data = unc_calc.uncertainty(
        unc_data,
        rp = config.rp,
        calc_eai_exp=config.calc_eai_exp,
        calc_at_event=config.calc_at_event,
        pool=pool)

    #Save uncertainty data
    unc_data.to_hdf5(uncertainty_filename)

    pool.close()
    pool.join()
    pool.clear()

    if config.calc_at_event:
        unc_data.at_event_unc_df = pd.DataFrame([])

    #Compute sensitivity
    unc_data = unc_calc.sensitivity(
        unc_data,
        sensitivity_method=config.sensitivity_method,
        sensitivity_kwargs=config.sensitivity_kwargs,
        )

    if config.calc_eai_exp:
        unc_data.eai_exp_unc_df = pd.DataFrame([])

    #Save data
    unc_data.to_hdf5(sensitivity_filename)


main(CONFIG_FILE)
