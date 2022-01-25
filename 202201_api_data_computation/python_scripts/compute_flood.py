import os
import sys
import glob

from climada.hazard import Centroids
from climada_petals.hazard.river_flood import RiverFlood

from config import IN_DATA_DIR
from config import OUT_DATA_DIR


def main(years=[1980, 2000], scenario='hist'):
    if scenario == 'hist':
        flddph_data_dir = os.path.join(IN_DATA_DIR, 'flood/flood_flddph_hist')
    else:
        flddph_data_dir = os.path.join(IN_DATA_DIR, 'flood/flood_flddph')
    centroids = Centroids.from_hdf5(
        os.path.join(OUT_DATA_DIR, "centroids/earth_centroids_150asland_1800oceans_distcoast_region_nopoles.hdf5"))
    years_str  = "_".join([str(years[0]), str(years[1])])
    years = range(int(years[0]), int(years[1]))
    filename = "".join(['river_flood_150arcsec', '_', scenario, '_', years_str, '.hdf5'])
    path = os.path.join(OUT_DATA_DIR, 'flood','global', scenario, years_str)
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


if __name__ == "__main__":
    main(years=[sys.argv[1], sys.argv[2]], scenario=sys.argv[3])