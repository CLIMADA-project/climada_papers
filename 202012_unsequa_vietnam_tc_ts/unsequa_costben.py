#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 14:31:14 2021

@author: ckropf
"""

import importlib
import datetime as dt

import scipy as sp

from pathlib import Path
from pathos.pools import ProcessPool as Pool

from climada.hazard import Hazard
from climada.engine.unsequa import CalcCostBenefit, InputVar #CLIMADA >= 3.1.0
#from climada.engine.uncertainty_quantification import CalcCostBenefit, InputVar #CLIMADA < 3.1.0

from impf_VNM import impf
from entity_VNM import generate_litpop_base


DATA_DIR = Path("VNM_Data/")
OUT_DIR = Path("./")

CONFIG_FILE = "config_costben"

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

    haz_pres = Hazard()
    haz_pres.read_hdf5(DATA_DIR / config.haz_file)

    haz_fut = Hazard()
    haz_fut.read_hdf5(DATA_DIR / config.haz_fut_file)

    impf_set = impf()

    n_ev = config.n_ev if config.n_ev is None else int(config.n_ev * haz_pres.size)
    #Create uncertainty object
    haz_iv = InputVar.haz(haz_pres, n_ev=n_ev)

    haz_fut_iv = InputVar.haz(haz_fut, n_ev=n_ev,
                                bounds_int=config.bounds_int,
                                bounds_freq=config.bounds_freq)

    haz_id_dict = {haz_pres.tag.haz_type: [config.impf_id]}

    litpop_list = generate_litpop_base(
        config.impf_id, config.value_unit, haz_pres, config.assign_centr_kwargs,
        config.choice_mn, **config.litpop_kwargs
        )

    ent_iv = InputVar.ent(
        bounds_disc=config.bounds_disk,
        bounds_cost=config.bounds_cost,
        bounds_totval=config.bounds_totval,
        bounds_mdd=config.bounds_mdd,
        bounds_impfi=config.bounds_impfi,
        exp_list=litpop_list,
        impf_set=impf_set,
        disc_rate=config.disc,
        meas_set=config.meas_pres,
        haz_id_dict=haz_id_dict
        )

    myclip_a, myclip_b = config.bounds_totval
    my_mean, my_std = config.mean_std
    a, b = (myclip_a - my_mean) / my_std, (myclip_b - my_mean) / my_std
    ent_iv.distr_dict['ET'] = sp.stats.truncnorm(a=a, b=b, loc=my_mean, scale=my_std)

    litpop_fut_list = generate_litpop_base(
        config.impf_id, config.value_unit, haz_pres, config.assign_centr_kwargs,
        config.choice_mn, fut_year=config.horizon, **config.litpop_kwargs
        )

    ent_fut_iv = InputVar.entfut(
        bounds_cost=config.bounds_cost,
        bounds_eg=config.bounds_eg,
        bounds_mdd=config.bounds_mdd,
        bounds_impfi=config.bounds_impfi,
        exp_list=litpop_fut_list,
        impf_set=impf_set,
        meas_set=config.meas_fut,
        haz_id_dict=haz_id_dict
        )

    unc_calc = CalcCostBenefit(haz_iv, ent_iv,
                             haz_fut_iv, ent_fut_iv)

    #Make sample
    unc_output = unc_calc.make_sample(
        N=config.N_samples,
        sampling_method=config.sampling_method,
        sampling_kwargs=config.sampling_kwargs
        )

    pool = Pool()
    unc_output = unc_calc.uncertainty(unc_output, pool=pool)

    #Save data
    unc_output.to_hdf5(uncertainty_filename)

    pool.close()
    pool.join()
    pool.clear()

    #Compute sensitivity
    unc_output = unc_calc.sensitivity(
        unc_output,
        sensitivity_method=config.sensitivity_method,
        sensitivity_kwargs=config.sensitivity_kwargs,
        )

    #Save data
    unc_output.to_hdf5(sensitivity_filename)


main(CONFIG_FILE)
