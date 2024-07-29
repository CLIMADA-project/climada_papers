"""
"""
import sys
import os
import pickle
import numpy as np
import pandas as pd
from iso3166 import countries as iso_cntry

from climada.entity import BlackMarble, ImpactFuncSet, IFTropCyclone
from climada.hazard import TropCyclone, Centroids
from climada.engine import Impact
from climada.util.save import save

from constants import CNTRIES, CNTRIES_ISO, RESOL, YEAR, GDP_NLD_ISL, INC_GRP, GDP, POLY_VAL
from plt_exposure_irma import fig03_fig04
from plt_analysis import fig06


DATA_DIR = os.path.abspath(os.path.dirname(__file__))
""" Input/output data folder relative path """

IBTRACS_DIR = os.path.join(DATA_DIR, 'tracks')
""" Tracks data in DATA_DIR """

FIG_DIR = DATA_DIR
""" Folder where the images are written """

def calc_tracks(data_dir):
    """ Compute tracks from ibtracs data, if not contained in data_dir.
    This functions is the longest one to execute."""
    try:
        abs_path = os.path.join(data_dir, 'sel_hist_syn_1h.p')
        with open(abs_path, 'rb') as f:
            sel_ibtracs = pickle.load(f)
        print('Loaded sel_hist_syn_1h:', sel_ibtracs.size)
    except FileNotFoundError:
        abs_path = os.path.join(data_dir, 'sel_hist.p')
        with open(abs_path, 'rb') as f:
            sel_ibtracs = pickle.load(f)
        print('Loaded sel_hist:', sel_ibtracs.size)

        sel_ibtracs.calc_random_walk(49)
        print('num tracks hist+syn:', sel_ibtracs.size)
        save(os.path.join(data_dir, 'sel_hist_syn.p'), sel_ibtracs)

        sel_ibtracs.equal_timestep(1)
        save(os.path.join(data_dir, 'sel_hist_syn_1h.p'), sel_ibtracs)

    return sel_ibtracs

def calc_exposure(data_dir):
    """ Compute exposure in every island group, if not contained in data_dir."""
    try:
        abs_path = os.path.join(data_dir, 'exp_irma.p')
        with open(abs_path, 'rb') as f:
            expo_dict = pickle.load(f)
        print('Loaded exp_irma:', len(expo_dict))
    except FileNotFoundError:
        expo_dict = dict()
        for cntry, cntry_iso in zip(CNTRIES, CNTRIES_ISO):
            if cntry == 'Netherlands':
                ent = BlackMarble()
                ent.set_countries({cntry: ['St. Eustatius', 'Saba']}, YEAR,
                                   res_km=RESOL, poly_val=POLY_VAL)
                ent.value = ent.value/ent.value.sum()*GDP_NLD_ISL*(INC_GRP['NLD']+1)
            else:
                ent = BlackMarble()
                ent.set_countries([cntry], YEAR, res_km=RESOL, poly_val=POLY_VAL,
                                  **{'gdp': GDP, 'inc_grp': INC_GRP})

            # set impact functions for TCs from if_ column (ones)
            expo_dict[cntry_iso] = ent.rename(columns={"if_": "if_TC"})

        save(os.path.join(data_dir, 'exp_irma.p'), expo_dict)

    return expo_dict

def calc_tc(expo_dict, tracks, data_dir):
    """ Compute tropical cyclone events from tracks at every island group,
    if not contained in data_dir. """
    try:
        abs_path = os.path.join(data_dir, 'tc_isl.p')
        with open(abs_path, 'rb') as f:
            tc_dict = pickle.load(f)
        print('Loaded tc_isl:', len(tc_dict))
    except FileNotFoundError:
        all_isl = BlackMarble(pd.concat(list(expo_dict.values())))

        centr = Centroids()
        centr.coord = np.zeros((all_isl.latitude.size, 2))
        centr.coord[:, 0] = all_isl.latitude
        centr.coord[:, 1] = all_isl.longitude
        centr.id = np.arange(centr.lat.size) + 1
        centr.region_id = all_isl.region_id

        tc = TropCyclone()
        tc.set_from_tracks(tracks, centr)

        tc_dict = dict()
        for ent_iso, ent_val in expo_dict.items():
            reg_id = np.unique(ent_val.region_id)[0]
            tc_dict[ent_iso] = tc.select(reg_id=reg_id)

        save(os.path.join(data_dir, 'tc_isl.p'), tc_dict)

    return tc_dict

def calc_imp(expo_dict, tc_dict, data_dir):
    """ Compute impacts of TCs in every island group. """
    try:
        abs_path = os.path.join(data_dir, 'imp_isl.p')
        with open(abs_path, 'rb') as f:
            imp_dict = pickle.load(f)
        print('Loaded imp_isl:', len(imp_dict))
    except FileNotFoundError:
        if_exp = ImpactFuncSet()
        if_em = IFTropCyclone()
        if_em.set_emanuel_usa()
        if_exp.add_func(if_em)

        imp_dict = dict()
        for isl_iso in expo_dict:
            imp = Impact()
            imp.calc(expo_dict[isl_iso], if_exp, tc_dict[isl_iso])
            imp_dict[isl_iso] = imp

        save(os.path.join(data_dir, 'imp_isl.p'), imp_dict)

    return imp_dict

def get_irma_damage(imp_dict):
    """ Compute Irma's damage per island group """
    msg = "Damage Irma on {0:<30s}: {1:1.3e}"
    for imp_isl, imp_val in imp_dict.items():
        id_irma = imp_val.event_name.index('2017242N16333')
        print(msg.format(iso_cntry.get(imp_isl).name,
              imp_val.at_event[id_irma]))

def aai_isl(imp_dict):
    """ Compute average annual impact"""
    msg = "Average annual impact on {0:<30s}: {1:1.3e}"
    for imp_isl, imp_val in imp_dict.items():
        print(msg.format(iso_cntry.get(imp_isl).name, imp_val.aai_agg))

def get_efc_isl(imp_dict):
    """ Compute impact exceedance level for different return periods"""
    msg = "Impact level for return period [{0:d} {1:d} {2:d} {3:d} {4:d}] on " \
          +"{5:<30s}: [{6:1.3e} {7:1.3e} {8:1.3e} {9:1.3e} {10:1.3e}]"
    for imp_isl, imp_val in imp_dict.items():
        imp_rp = imp_val.calc_freq_curve([25, 100, 250, 1000, 3000]).impact
        print(msg.format(25, 100, 250, 1000, 3000, iso_cntry.get(imp_isl).name,
              imp_rp[0], imp_rp[1], imp_rp[2], imp_rp[3], imp_rp[4]))

def main(argv):
    print('Input/Output data folder: ', DATA_DIR)

    # exposures
    expo_dict = calc_exposure(DATA_DIR)

    # tracks
    sel_tr = calc_tracks(DATA_DIR)

    # dictionary of tc per island
    tc_dict = calc_tc(expo_dict, sel_tr, DATA_DIR)

    # damage per isl
    imp_dict = calc_imp(expo_dict, tc_dict, DATA_DIR)

    # damage irma
    get_irma_damage(imp_dict)

    # average annual impact
    aai_isl(imp_dict)

    # compute impact exceedance frequency
    get_efc_isl(imp_dict)

    # FIG03 and FIG04
    fig03_fig04(IBTRACS_DIR, DATA_DIR, FIG_DIR) # 5min

    # FIG 06
    fig06(DATA_DIR, FIG_DIR)


if __name__ == "__main__":
   main(sys.argv[1:])
