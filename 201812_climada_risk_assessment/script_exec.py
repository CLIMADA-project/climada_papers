"""
"""
import sys
import os
import pickle
import numpy as np
from iso3166 import countries as iso_cntry

from climada.entity import BlackMarble, ImpactFuncSet, IFTropCyclone
from climada.hazard import TCTracks, TropCyclone, Centroids
from climada.engine import Impact
from climada.util.save import save
from climada.util.files_handler import get_file_names

DATA_DIR = os.path.abspath(os.path.dirname(__file__))
""" Input/output data folder relative path """

IBTRACS_DIR = os.path.join(DATA_DIR, 'tracks')
""" Tracks data in DATA_DIR """

RESOL = 0.5
""" Approx. resolution in km """

YEAR = 2016
""" Year of exposure data """

CNTRIES = ['Saint Barthelemy', 'Saint Martin', 'Sint Maarten', 'Anguilla',
           'British Virgin Islands', 'United States Virgin Islands',
           'Turks And Caicos Islands', 'Saint Kitts And Nevis',
           'Antigua And Barbuda', 'Netherlands']
""" Country (island groups) names """

CNTRIES_ISO = ['BLM', 'MAF', 'SXM', 'AIA', 'VGB', 'VIR', 'TCA', 'KNA', 'ATG', 'NLD']
""" Country (island groups) ISO3 codes """

GDP = {'BLM': 414710000, 'MAF': 614258169, 'SXM': 1081577185, \
       'AIA': 337201995, 'VGB': 971237110, \
       'VIR': 3765000000, 'TCA': 917550492, \
       'KNA': 909854630, 'ATG': 1460144703, 'NLD': ''}
""" GDP at YEAR per island group """

GDP_NLD_ISL = 48.0e6 + 100.0e6
""" GDP Saba and St. Eustatius """

INC_GRP_DEF = 4
INC_GRP = {'BLM': INC_GRP_DEF, 'MAF': INC_GRP_DEF, 'SXM': INC_GRP_DEF,
           'AIA': INC_GRP_DEF, 'VGB': INC_GRP_DEF, 'VIR': INC_GRP_DEF,
           'TCA': INC_GRP_DEF, 'KNA': INC_GRP_DEF, 'ATG': INC_GRP_DEF,
           'NLD': INC_GRP_DEF}
""" income group level at YEAR per island group """

POLY_VAL = [0, 0, 1]
""" Polygonal transformation in night lights """

def calc_tracks(data_dir, ibtracs_dir):
    """ Compute tracks from ibtracs data, if not contained in data_dir.
    This functions is the longest one to execute."""
    try:
        abs_path = os.path.join(data_dir, 'sel_hist_syn_1h.p')
        with open(abs_path, 'rb') as f:
            sel_ibtracs = pickle.load(f)
        print('Loaded sel_hist_syn_1h:', sel_ibtracs.size)
    except FileNotFoundError:
        sel_ibtracs = TCTracks()
        for file in get_file_names(ibtracs_dir):
            tmp_tr = TCTracks()
            tmp_tr.read_ibtracs_csv(file)
            sel_ibtracs.append(tmp_tr.data)
        print('num tracks hist:', sel_ibtracs.size)
        save(os.path.join(data_dir, 'sel_hist.p'), sel_ibtracs)

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

            expo_dict[cntry_iso] = ent

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
        all_isl = BlackMarble()
        for ent_iso, ent_val in expo_dict.items():
            all_isl.append(ent_val)

        centr = Centroids()
        centr.coord = all_isl.coord
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
    sel_tr = calc_tracks(DATA_DIR, IBTRACS_DIR)

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

if __name__ == "__main__":
   main(sys.argv[1:])
