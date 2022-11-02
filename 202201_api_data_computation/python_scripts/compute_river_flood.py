import os
import sys
import glob

from climada.hazard import Centroids
from climada_petals.hazard.river_flood import RiverFlood
from config import DATA_DIR
from create_log_file import log_msg

CENT_FILE_PATH = os.path.join(DATA_DIR, "centroids/earth_centroids_150asland_1800asoceans_distcoast_region.hdf5")
OUT_FILE_PATH = "".join(['river_flood_150arcsec_{scenario}_{years_str}.hdf5'])


def main(years=None, scenario='hist'):
    LOG_FILE = "progress_make_river_flood_global.txt"
    if years is None:
        years = [1980, 2000]
    if scenario == 'hist':
        flddph_data_dir = os.path.join(DATA_DIR, 'flood/flood_flddph_hist')
    else:
        flddph_data_dir = os.path.join(DATA_DIR, "".join(['flood/flood_flddph/', scenario]))
    log_msg(f"Started computing floods for scenario {scenario}  and years {years}\n", LOG_FILE)

    centroids = Centroids.from_hdf5(CENT_FILE_PATH)
    years_str = "_".join([str(years[0]), str(years[1])])
    years = range(int(years[0]), int(years[1]))
    filename = OUT_FILE_PATH.format(scenario=scenario, years_str=years_str)
    path = os.path.join(DATA_DIR, 'flood_v2', 'global', scenario, years_str)
    file = os.path.join(path, filename)
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    rf_list = []
    files_list = glob.glob(os.path.join(flddph_data_dir, '*.nc*'))
    for flddph_file in files_list:
        flddph_file = flddph_file.split('/')[-1]
        flddph_file_path = os.path.join(flddph_data_dir, flddph_file)
        fldfrc_file_path = flddph_file_path.replace('flddph', 'fldfrc')
        rf = RiverFlood()
        rf.set_from_nc(years=years, dph_path=flddph_file_path, frc_path=fldfrc_file_path, centroids=centroids)
        rf.event_name = ["_".join([str(y), flddph_file.split('_')[2], flddph_file.split('_')[3]]) for y in years]
        rf.tag.haz_type = 'RF'
        rf_list.append(rf)
    rf_concat = rf.concat(rf_list)
    rf_concat.frequency = rf_concat.frequency / len(files_list)
    rf_concat.write_hdf5(file)
    log_msg(f"Computing of river floods for scenario {scenario}  and years {years} done\n", LOG_FILE)


if __name__ == "__main__":
    main(years=[sys.argv[1], sys.argv[2]], scenario=sys.argv[3])