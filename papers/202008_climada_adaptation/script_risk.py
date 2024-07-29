"""
"""
import sys
import os
import numpy as np
from rasterio import Affine
from pathos.pools import ProcessPool as Pool

from climada.entity import BlackMarble, ImpactFuncSet, IFTropCyclone, ImpactFunc
from climada.hazard import TCTracks, TropCyclone, TCSurge, Centroids
from climada.engine import Impact
from climada.util.constants import ENT_TEMPLATE_XLS, DEF_CRS
from constants import RESOL, YEAR, GDP_NLD_ISL, INC_GRP, \
GDP, POLY_VAL, CNTRIES, CNTRIES_ISO
from climada.util.save import save, load
import climada.hazard.centroids.centr as cr
cr.MAX_DEM_TILES_DOWN = 800

CURR_DIR = os.path.abspath(os.path.dirname(__file__))
""" Current directory """

OUT_FOLDER = os.path.join(CURR_DIR, 'results_eca')
""" Folder where outputs will be written """

try:
    os.mkdir(OUT_FOLDER)
except FileExistsError:
    pass

def calc_exposure():
    """ Compute or load existing exposure islands using BlackMarble """
    expo_dict = dict()
    for cntry, cntry_iso in zip(CNTRIES, CNTRIES_ISO):
        file_name = os.path.join(OUT_FOLDER, 'exp_'+cntry_iso+'.h5')
        ent = BlackMarble()
        if os.path.isfile(file_name):
            ent.read_hdf5(file_name)
            expo_dict[cntry_iso] = ent
            print('Loaded:', cntry, expo_dict[cntry_iso].shape[0])
            continue
        if cntry == 'Netherlands':
            ent.set_countries({cntry: ['St. Eustatius', 'Saba']}, YEAR,
                               res_km=RESOL, poly_val=POLY_VAL)
            ent.value = ent.value/ent.value.sum()*GDP_NLD_ISL*(INC_GRP['NLD']+1)
        else:
            ent.set_countries([cntry], YEAR, res_km=RESOL, poly_val=POLY_VAL,
                              **{'gdp': GDP, 'inc_grp': INC_GRP})
        ent.rename(columns={"if_": "if_TC"}, inplace=True)
        ent.check()
        ent.write_hdf5(file_name)
        expo_dict[cntry_iso] = ent
    return expo_dict

def calc_tracks(pool=None):
    """ Compute or load existing tracks from ibtracs data."""
    # 13:05 - 14:13 --> 1 hour
    file_name = os.path.join(OUT_FOLDER, 'tracks_na_1950_2018_05.p')
    if os.path.isfile(file_name):
        sel_ibtracs = load(file_name)
        print('Loaded tracks:', sel_ibtracs.size)
    else:
        sel_ibtracs = TCTracks(pool)
        sel_ibtracs.read_ibtracs_netcdf(provider='usa', basin='NA',
                                        year_range=(1950, 2018), correct_pres=True)
        print('num tracks hist:', sel_ibtracs.size)
        sel_ibtracs.calc_random_walk(49)
        print('num tracks hist+syn:', sel_ibtracs.size)
        sel_ibtracs.equal_timestep(0.5)
        save(file_name, sel_ibtracs)
    return sel_ibtracs

def _add_raster_sea(meta, new_cells):
    """ Add sea raster surrounding each island """
    new_trans = [0] * 6
    new_trans[0] = meta['transform'][0]
    new_trans[4] = meta['transform'][4]
    new_trans[2] = meta['transform'][2] - meta['transform'][0]*new_cells
    new_trans[5] = meta['transform'][5] - meta['transform'][4]*new_cells
    new_trans = Affine(new_trans[0], new_trans[1], new_trans[2], new_trans[3], new_trans[4], new_trans[5])

    new_meta = meta
    new_meta['transform'] = new_trans
    new_meta['height'] += new_cells*2
    new_meta['width'] += new_cells*2

    return new_meta

def calc_tc(expo_dict, tracks, pool=None):
    """ Compute tropical cyclone wind gust from tracks at every island group,
    if not contained in OUT_FOLDER. """
    # 23:09 - 02:01 --> 3 horas
    if os.path.isfile(os.path.join(OUT_FOLDER, 'tc_'+'BLM'+'.h5')):
        tc_dict = dict()
        for ent_iso in expo_dict.keys():
            tc_dict[ent_iso] = TropCyclone()
            tc_dict[ent_iso].read_hdf5(os.path.join(OUT_FOLDER, 'tc_'+ent_iso+'.h5'))
            tc_dict[ent_iso].check()
        print('Loaded tc_isl:', len(tc_dict))
    else:
        new_cells = 10 # sea cells to add as envelope
        centr_all = Centroids()
        centr_all.geometry.crs = DEF_CRS
        for key, value in expo_dict.items():
            centr = Centroids()
            centr.meta = _add_raster_sea(value.meta, new_cells)
            centr.set_meta_to_lat_lon()
            centr.region_id = np.ones(centr.size)*value.region_id[0]
            centr.meta = dict()
            centr_all.append(centr)
        centr_all.check()

        tc = TropCyclone(pool)
        tc.set_from_tracks(tracks, centr_all)
        tc.check()

        tc_dict = dict()
        for ent_iso, ent_val in expo_dict.items():
            reg_id = np.unique(ent_val.region_id)[0]
            tc_dict[ent_iso] = tc.select(reg_id=reg_id)
            tc_dict[ent_iso].centroids.set_lat_lon_to_meta()
            tc_dict[ent_iso].centroids.lat = np.array([])
            tc_dict[ent_iso].centroids.lon = np.array([])
            tc_dict[ent_iso].write_hdf5(os.path.join(OUT_FOLDER, 'tc_'+ent_iso+'.h5'))
    return tc_dict

def calc_impact(expo_dict, tc_dict):
    """ Compute or load resulting impact statistics """
    # impact function TC
    ifs = ImpactFuncSet()
    if_tc = IFTropCyclone()
    if_tc.set_emanuel_usa()
    ifs.append(if_tc)

    imp_tc_dict = dict()
    for isl_iso in expo_dict.keys():
        imp_tc = Impact()
        file_name = os.path.join(OUT_FOLDER, 'imp_tc_'+isl_iso+'.csv')
        if os.path.isfile(file_name):
            imp_tc.read_csv(file_name)
            print('Loaded imp_tc:', imp_tc.at_event.size)
        else:
            imp_tc.calc(expo_dict[isl_iso], ifs, tc_dict[isl_iso])
            imp_tc.write_csv(file_name)
        imp_tc_dict[isl_iso] = imp_tc
    return imp_tc_dict, ifs

def main(argv):
    print('Input/Output data folder: ', OUT_FOLDER)

    # set parallel processing
    pool = Pool()

    # exposures
    expo_dict = calc_exposure()

    # tracks
    sel_tr = calc_tracks(pool)

    # dictionary of tc per island
    tc_dict = calc_tc(expo_dict, sel_tr, pool)

    # damage per isl
    calc_impact(expo_dict, tc_dict)

    # close processes
    pool.close()
    pool.join()

if __name__ == "__main__":
   main(sys.argv[1:])
