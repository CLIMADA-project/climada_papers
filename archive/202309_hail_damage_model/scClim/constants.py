# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 14:16:50 2022
by Timo Schmid

Constants to be used for hail calculations
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy

VARNAME_DICT = {'meshs':'MZC',
                'MESHS':'MZC',
                'POH':'BZC',
                'poh':'BZC',
                }
UNIT_DICT = {'meshs':'mm',
            'MESHS':'mm',
            'POH':'%',
            'poh':'%',
            'crowd':'mm',
            'crowdFiltered':'mm',
            'E_kin':'J/m2',
            'E_kinCC':'J/m2',
            'MESHSdBZ':'mm',
            'MESHSdBZ_p3':'mm',
            'MESHS_4km':'mm',
            'dBZ':'dBZ',
            'VIL':'g/m2',
            }
CH_EXTENT = [5.8, 10.6, 45.7, 47.9]#x0 x1 y0 y1
CH_EXTENT_EPSG2056 = (2485014.011782823, 2837016.892254034, 1075214.1686203221, 1299782.7670088513)
ZRH_EXTENT = [8.35, 9, 47.15, 47.7]
ZRH_EXTENT_EPSG2056 = (2668000, 2721000, 1222000, 1284000)
SUB_CH_EXTENT_2056 = (2.55e6, 2.75e6, 1.125e6,1.29e6)
CMAP_VIR = copy.copy(plt.cm.get_cmap('viridis'))
CMAP_VIR.set_under('white',alpha=0)

# Variable_dictionaries
INT_RANGE_DICT = {
    'MESHS':    np.concatenate(([0],np.arange(20,100))),
    'MESHS_4km':np.concatenate(([0],np.arange(20,100))),
    'HKE':      np.arange(0,3001,40),
    'dBZ':      np.arange(45, 75, 0.5),
    'crowd':    np.arange(0, 70, 1),
    'crowdFiltered': np.arange(0, 70, 1),
    'POH':      np.arange(0, 100, 1),
    'E_kin':    np.concatenate(([0],np.arange(60,1000,10))),
    'E_kinCC':  np.arange(0,2000,50),
    'MESHSdBZ': np.arange(0, 90, 2),
    'MESHSdBZ_p3': np.arange(0, 90, 2),
    'VIL':      np.concatenate(([0],np.arange(10,70))),
}

INT_LABEL_DICT = {
    'MESHS_4km': 'Intensity: MESHS [mm]',
    'MESHSdBZ_p3': 'Intensity: MESHS (dBZ-scaled) [mm]',
    'MESHSdBZ': 'Intensity: MESHS (dBZ-scaled) [mm]',
    'MESHS': 'Intensity: MESHS [mm]',
    'HKE': 'Intensity: HKE (MESHS-derived) [J m$^{-2}$]',
    'dBZ': 'Intensity: Reflectivity [dBZ]',
    'crowdFiltered': 'Intensity: Crowd-sourced hail size [mm]',
    'crowd': 'Intensity: Crowd-sourced hail size [mm]',
    'POH': 'Intensity: Probability of hail [%]',
    'E_kinCC': 'Intensity: E$_{kin}$ [J m$^{-2}$]',
    'E_kin': 'Intensity: E$_{kin}$ [J m$^{-2}$]',
    'VIL': 'Intensity: VIL [g m$^{-2}$]',
}

CUT_OFF_DICT = {
    'MESHS':    60,
    'MESHS_4km':60,
    'HKE':      2000,
    'dBZ':      65,
    'crowd':    45,
    'crowdFiltered': 45,
    'POH':      100,
    'E_kin':    800,
    'E_kinCC':  800,
    'MESHSdBZ': 60,
    'MESHSdBZ_p3': 60,
    'VIL':      55,
}

ID_COL_DICT = {
    'GVL':          'Kantonale Versicherungsnummer',
    'AGV':          'VertragsNr',
    'GVB':          'Vertragsnummer',
    'GVZ':          'VersicherungsID',
    'MFZrandom_':   'POLNR',
    'KGV':          'id_col'
}

#Plot extent: should be roughly square
PLOT_EXTENT_DICT = { #x0 x1 y0 y1
    'GVZ':          [8.3, 9, 47.1, 47.8],
    '':             [8.3, 9, 47.1, 47.8],
    'MFZrandom_':   [6, 10, 44.5, 48.5],
    'GVL':          [7.85, 8.5,46.7,47.35],
    'AGV':          [7.75, 8.45,47.05,47.65],
    'GVB':          [6.9, 8.4,46.3,47.4],
    'KGV':          [6.9, 9,46.2,48],
}

CANTON_DICT = {
    'GVZ':          'Zürich',
    '':             'Zürich', #GVZ is also called as '' only
    'MFZrandom_':   None,
    'GVL':          'Luzern',
    'AGV':          'Aargau',
    'GVB':          'Bern',
    'KGV':          ['Zürich','Bern','Luzern','Aargau']
}

#Dictionary of windowsize for rolling window (in #steps)
W_SIZE_DICT = {
    'MESHS':    11,
    'MESHS_4km':11,
    'HKE':      15,
    'dBZ':      3,
    'crowd':    5,
    'crowdFiltered': 5,
    'POH':      7,
    'E_kin':    3,
    'E_kinCC':  3,
    'MESHSdBZ': 5,
    'MESHSdBZ_p3':5,
    'VIL':      7,
}

DMG_BIN_DICT = {
    'MESHS': 5,
    'MESHS_4km':5,
    'HKE': 10,
    'dBZ': 2,
    'crowd': 2,
    'crowdFiltered': 2,
    'POH': 2,
    'E_kin': 50,
    'E_kinCC': 50,
    'MESHSdBZ': 5,
    'MESHSdBZ_p3':5,
    'VIL': 5,
}

BAUINDEX = pd.DataFrame(
    index=np.arange(2000,2021+1),
    data = {
    #See PDF in ../data/AGV
    'AGV':[402,417,436,436,422,422,436,436,464,482,482,482,498,498,498,498,498,486,486,486,486,486],
    #https://gvb.ch/de/versicherungen/baukostenindex.html
    'GVB': [120,125,127,123,124,127,130,134,140,139,138,141,141,141,141,141,141,140,141,143,144,147],
    #source: https://www.gvl.ch/versicherung/versicherungswert/
    'GVL':[810,810,815,808,810,820,846,880,917,902,912,926,927,921,918,916,912,906,911,907,901,928],
    #GVZ source: GVZ hail loss data (.csv)
    'GVZ':[840,  900,  900,  900,  900,  900,  900,  900,  970, 1025, 1025,
       1025, 1025, 1025, 1025, 1025, 1025, 1025, 1025, 1025, 1025, 1025],
    #KGV damages are already indexed! see scClim/grid_cantonal_data.py
    'KGV':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], 
    }
)

DATA_RANGE_DICT = {
    'MESHS':    np.arange(2002,2021+1),
    'MESHS_4km':np.arange(2002,2021+1),
    'HKE':      np.arange(2002,2021+1),
    'dBZ':      np.arange(2013,2021+1),
    'crowd':    np.arange(2017, 2021+1),
    'crowdFiltered': np.arange(2017, 2021+1),
    'POH':      np.arange(2002,2021+1),
    'E_kin':    np.arange(2013,2021+1),
    'E_kinCC':  np.arange(2013,2021+1),
    'MESHSdBZ': np.arange(2013,2021+1),
    'MESHSdBZ_p3': np.arange(2013,2021+1),
    'VIL':      np.arange(2013,2021+1),
    #damage data
    '': np.arange(2002,2021+1),
    'GVZ': np.arange(2002,2021+1),
    'GVL': np.arange(2002,2021+1),
    'AGV': np.arange(2002,2021+1),
    'GVB': np.arange(2002,2021+1),
    'KGV': np.arange(2002,2021+1),
    'KGV_nonExtreme': np.arange(2002,2021+1),
    'KGV_1e5': np.arange(2002,2021+1),
    'scaledKGV_': np.arange(2002,2021+1),
    'gridKGV_': np.arange(2002,2021+1),
    'MFZrandom_': np.arange(2017,2021+1),

}

DROP_DATES_DICT={'Weizen': ['ev_2017-08-01'],
           'Gerste': ['ev_2017-07-08', 'ev_2017-08-01'],
           'Raps': ['ev_2017-08-01'],
           'Mais': [],
           'Reben': [],
           'Aepfel': []}

FRACTION_INSURED_DICT ={'Weizen': 0.69,
           'Gerste': 0.69,
           'Raps': 0.69,
           'Mais': 0.69,
           'Reben': 0.43,
           'Aepfel': 0.32}

PRE_PROC_PARAM_DICT = {
    7: {
        'version_id' : 7,
        'min_day' : -2, # maximum -2
        'max_day' : 2, # maximum 2
        'min_POH_diff' : 50, #min POH difference to change a date
        'poh_level' : 10, #POH level for plausible hail damage
        'buffer_km' : 5, #in km. buffer around POH level above
        'delete_nonPOH_dmgs' : True, #wether to delete non-plausible damages
        'use_likelihood' : False, #use boolean
    },

    #Version 8: to be used for MFZrandom ONLY
    8: {
        'version_id' : 8,
        'min_day' : -2, # maximum -2
        'max_day' : 2, # maximum 2
        'min_POH_diff' : None, #min POH difference to change a date: so date will never be changed based on local POH difference 
        'poh_level' : 10, #POH level for plausible hail damage
        'buffer_km' : 50, #in km. buffer around POH level above
        'delete_nonPOH_dmgs' : True, #wether to delete non-plausible damages
        'use_likelihood' : False, #use boolean
    },

    #Version 9: to be used for MFZrandom ONLY, does not change any dates!
    9: {
        'version_id' : 9,
        'min_day' : 0, # maximum -2
        'max_day' : 0, # maximum 2
        'min_POH_diff' : None, #min POH difference to change a date: so date will never be changed based on local POH difference 
        'poh_level' : 10, #POH level for plausible hail damage
        'buffer_km' : 50, #in km. buffer around POH level above
        'delete_nonPOH_dmgs' : False, #wether to delete non-plausible damages
        'use_likelihood' : False, #use boolean
    },

    #Version 1: to be used for MFZrandom ONLY
    #assigns likelyhood for each date, instead of boolean (possible hail = yes/no)
    1: {
        'version_id' : 1,
        'min_day' : -2, # maximum -2
        'max_day' : 2, # maximum 2
        'min_POH_diff' : None, #min POH difference to change a date: so date will never be changed based on local POH difference 
        'poh_level' : 10, #POH level for plausible hail damage
        'buffer_km' : 50, #in km. buffer around POH level above
        'delete_nonPOH_dmgs' : True, #wether to delete non-plausible damages
        'use_likelihood' : True, #wether to use likelihood instead of boolean
    },
}

BAUJAHR_DICT = {
    '' : (0,2022),
    'before1960' : (0,1959), 
    '1960-2002' : (1960,2002), 
    'after2002' : (2003,2022)
}