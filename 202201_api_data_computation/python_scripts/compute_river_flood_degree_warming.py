import os
import sys
import glob

from climada.hazard import Centroids
from climada_petals.hazard.river_flood import RiverFlood
import numpy as np
from config import DATA_DIR
from create_log_file import log_msg

YEARS_WARMING = {'1': {
    'rcp26': {'gfdl-esm2m': [2006, 2024], 'miroc5': [2006, 2025], 'hadgem2-es': [2006, 2022], 'ipsl-cm5a-lr': None},
    'rcp60': {'gfdl-esm2m': [2006, 2026], 'miroc5': [2013, 2033], 'hadgem2-es': [2006, 2024], 'ipsl-cm5a-lr': None},
    'rcp85': {'gfdl-esm2m': [2006, 2024], 'miroc5': [2015, 2024], 'hadgem2-es': [2006, 2022], 'ipsl-cm5a-lr': None}},
    '1.5': {
    'rcp26': {'gfdl-esm2m': None, 'miroc5': [2038, 2058], 'hadgem2-es': [2016, 2036], 'ipsl-cm5a-lr': [2006, 2019]},
    'rcp60': {'gfdl-esm2m': [2046, 2066], 'miroc5': [2042, 2062], 'hadgem2-es': [2022, 2042], 'ipsl-cm5a-lr': [2006, 2020]},
    'rcp85': {'gfdl-esm2m': [2026, 2046], 'miroc5': [2023, 2043], 'hadgem2-es': [2015,2035], 'ipsl-cm5a-lr': [2006, 2019]}},
    '2': {'rcp26': {'gfdl-esm2m': None, 'miroc5': None, 'hadgem2-es': None, 'ipsl-cm5a-lr': [2019, 2039]},
          'rcp60': {'gfdl-esm2m': [2066, 2086], 'miroc5': [2061, 2081], 'hadgem2-es': [2040, 2060], 'ipsl-cm5a-lr': [2019, 2039]},
          'rcp85': {'gfdl-esm2m': [2043, 2063], 'miroc5': [2038, 2058], 'hadgem2-es': [2027, 2047], 'ipsl-cm5a-lr': [2014, 2024]}},
    '3': {'rcp26': {'gfdl-esm2m': None, 'miroc5': None, 'hadgem2-es': None, 'ipsl-cm5a-lr': [2019, 2039]},
    'rcp60': {'gfdl-esm2m': [2066, 2086], 'miroc5': [2061, 2081], 'hadgem2-es': [2030, 2050],
              'ipsl-cm5a-lr': [2019, 2039]},
    'rcp85': {'gfdl-esm2m': [2043, 2063], 'miroc5': [2028, 2048], 'hadgem2-es': [2027, 2047],
              'ipsl-cm5a-lr': [2014, 2024]}}}


CENT_FILE_PATH = os.path.join(DATA_DIR, "centroids/earth_centroids_150asland_1800asoceans_distcoast_region.hdf5")

OUT_FILE = 'river_flood_150arcsec_{warming}.hdf5'


def main(warming='1'):
    LOG_FILE = "progress_make_river_flood_degrees_warming.txt"
    log_msg(f"Started computing floods for {warming} degrees warming level\n", LOG_FILE)
    centroids = Centroids.from_hdf5(CENT_FILE_PATH)
    filename = OUT_FILE.format(warming=warming)
    path = os.path.join(DATA_DIR, 'flood_warming_level', 'global', warming)
    file = os.path.join(path, filename)
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path, exist_ok=True)
    rf_list = []
    for scenario in ['rcp26', 'rcp60', 'rcp85']:
        flddph_data_dir = os.path.join(DATA_DIR, "".join(['flood/flood_flddph/', scenario]))
        files_list = glob.glob(os.path.join(flddph_data_dir, '*.nc*'))
        for flddph_file in files_list:
            gcm = [gcm for gcm in YEARS_WARMING[warming][scenario].keys() if gcm in flddph_file][0]
            years = YEARS_WARMING[warming][scenario][gcm]
            if years is None:
                continue
            flddph_file = flddph_file.split('/')[-1]
            flddph_file_path = os.path.join(flddph_data_dir, flddph_file)
            fldfrc_file_path = flddph_file_path.replace('flddph', 'fldfrc')
            log_msg(f"Initating river flood hazard for {scenario}.\n", LOG_FILE)
            rf = RiverFlood()
            rf.set_from_nc(years=list(np.arange(int(years[0]), int(years[1]))), dph_path=flddph_file_path, frc_path=fldfrc_file_path, centroids=centroids)
            rf.event_name = ["_".join([str(y), flddph_file.split('_')[3], scenario, flddph_file.split('_')[2]]) for y in np.arange(int(years[0]),int(years[1]))]
            rf.tag.haz_type = 'RF'
            rf_list.append(rf)
    log_msg(f"Concatenating hazards\n", LOG_FILE)

    rf_concat = RiverFlood.concat(rf_list)
    rf_concat.frequency = np.ones(len(rf_concat.event_id))/len(rf_concat.event_id)
    rf_concat.write_hdf5(file)


if __name__ == "__main__":
    main(warming=sys.argv[1])
