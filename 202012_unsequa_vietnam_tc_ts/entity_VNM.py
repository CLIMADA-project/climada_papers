#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 2021

@author: C. M. Kropf
"""

import numpy as np

from climada.entity import Measure
from climada.entity import MeasureSet
from climada.entity import DiscRates

BASE_YEAR = 2020
HORIZON = 2050

#Define a generic method to make litpop instances with different exponent pairs.
def generate_litpop_base(impf_id, value_unit, haz, assign_centr_kwargs,
                          choice_mn, fut_year=None, **litpop_kwargs):
    from climada.entity import LitPop
    litpop_base = []
    for [m, n] in choice_mn:
        print('\n Computing litpop for m=%f, n=%f \n' %(m, n))
        litpop_kwargs['exponents'] = (m, n)
        exp = LitPop.from_countries(**litpop_kwargs)
        exp.gdf['impf_' + haz.tag.haz_type] = impf_id
        exp.gdf.drop('impf_', axis=1, inplace=True)
        if value_unit is not None:
            exp.value_unit = value_unit
        exp.assign_centroids(haz, **assign_centr_kwargs)
        exp.gdf.region_id = 0
        exp.gdf.loc[(exp.gdf.latitude < 11.029167) & (exp.gdf.longitude < 106.820833), 'region_id'] = 1 #define north/south regions
        exp.gdf.loc[exp.gdf['latitude'] > 19.912500, 'region_id'] = 2
        if fut_year is not None:
            exp.ref_year = fut_year
        litpop_base.append(exp)
    return litpop_base

def measures_50():
    #measures in 2050

    n_years = HORIZON - BASE_YEAR

    TS_meas1 = Measure()
    TS_meas1.name = 'Mangrove'
    TS_meas1.haz_type = 'TS'
    TS_meas1.exp_region_id = [1, 2]
    TS_meas1.color_rgb = np.array([0.16, 0.62, 0.56])
    mangrove_maintain = 4000000
    TS_meas1.cost = 172 * 1000000 * 0.000043 * 1160 * 1000 * 150 / 10000 + mangrove_maintain * n_years  # 361,664,400$
    TS_meas1.mdd_impact = (1, 0)
    TS_meas1.paa_impact = (1, 0)
    TS_meas1.hazard_inten_imp = (1, -0.5)

    TS_meas2 = Measure()
    TS_meas2.name = 'Seadykes'
    TS_meas2.haz_type = 'TS'
    TS_meas2.color_rgb = np.array([0.91, 0.77, 0.42])
    TS_meas2.exp_region_id = [1]
    seadyke_maintain = 200
    TS_meas2.cost = (110 * 1000000 * 0.000043 + seadyke_maintain * n_years) * 150 * 1000
    TS_meas2.hazard_inten_imp = (1, -2)

    TS_meas3 = Measure()
    TS_meas3.name = 'Gabions'
    TS_meas3.haz_type = 'TS'
    TS_meas3.color_rgb = np.array([0.65, 0.65, 0.55])
    TS_meas3.exp_region_id = [1]
    gabion_maintain = 130
    TS_meas3.cost = (1300 + gabion_maintain) * 150 * 1000
    TS_meas2.cost = (110 * 1000000 * 0.000043 + 0.7 * seadyke_maintain * n_years) * 150 * 1000  # protect against event once every 20 years
    TS_meas3.hazard_inten_imp = (1, -0.5)

    TS_meas_set = MeasureSet()
    TS_meas_set.append(TS_meas1)
    TS_meas_set.append(TS_meas2)
    TS_meas_set.append(TS_meas3)
    TS_meas_set.check()

    MeasSet = MeasureSet()
    MeasSet.append(TS_meas1)
    MeasSet.append(TS_meas2)
    MeasSet.append(TS_meas3)

    return  MeasSet

def measures_20():
    #Measures in 2020


    TS_meas1 = Measure()
    TS_meas1.name = 'Mangrove'
    TS_meas1.haz_type = 'TS'
    TS_meas1.exp_region_id = [1, 2]
    TS_meas1.color_rgb = np.array([0.16, 0.62, 0.56])
    TS_meas1.cost = 172 * 1000000 * 0.000043 * 1160 * 1000 * 150 / 10000
    TS_meas1.mdd_impact = (1, 0)
    TS_meas1.paa_impact = (1, 0)
    TS_meas1.hazard_inten_imp = (1, -0.5)

    TS_meas2 = Measure()
    TS_meas2.name = 'Seadykes'
    TS_meas2.haz_type = 'TS'
    TS_meas2.color_rgb = np.array([0.91, 0.77, 0.42])
    TS_meas2.exp_region_id = [1]
    TS_meas2.cost = (110 * 1000000 * 0.000043) * 150 * 1000
    TS_meas2.hazard_inten_imp = (1, -2)

    TS_meas3 = Measure()
    TS_meas3.name = 'Gabions'
    TS_meas3.haz_type = 'TS'
    TS_meas3.color_rgb = np.array([0.65, 0.65, 0.55])
    TS_meas3.exp_region_id = [1]
    gabion_maintain = 130
    TS_meas3.cost = (1300 + gabion_maintain) * 150 * 1000
    TS_meas2.cost = (110 * 1000000 * 0.000043) * 150 * 1000  # protect against event once every 20 years
    TS_meas3.hazard_inten_imp = (1, -0.5)

    MeasSet = MeasureSet()
    MeasSet.append(TS_meas1)
    MeasSet.append(TS_meas2)
    MeasSet.append(TS_meas3)

    return  MeasSet


def disc_rates():

    # number of people does not discount
    people_disc = DiscRates()
    people_disc.years = np.arange(2000, 2100)
    people_disc.rates = np.zeros(people_disc.years.size)

    return {
        'PE': people_disc
        }

