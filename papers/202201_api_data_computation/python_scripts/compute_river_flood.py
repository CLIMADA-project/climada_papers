import os
import sys
import glob
import datetime
from climada.hazard import Centroids
from climada_petals.hazard.river_flood import RiverFlood
from config import DATA_DIR
from create_log_file import log_msg
from climada.util.api_client import Client
OUT_FILE_PATH = "".join(['river_flood_150arcsec_{scenario}_{years_str}.hdf5'])

# Data sources for river flood script. Ensure these datasets are downloaded and available:
DATA_LINKS = {
    'ISIMIP2a': 'https://zenodo.org/record/4446364#.YbCOjbso_mE', # Observed climate
    'ISIMIP2b': 'https://zenodo.org/record/4627841#.Ysb_iHhBw5k'
}
# The hazards on the CLIMADA API were computed with the files that end with flopros.nc,
# indicating the consideration of protection measures.
# Reference: https://nhess.copernicus.org/articles/16/1049/2016/
# Note: Calculations can also be done without protection measures or with 'flopros_modeled',
# which means protection measures are modeled where observations are not available.


def main(years=None, scenario='hist'):
    """
        Compute river floods hazard for given years and scenario and write results to an HDF5 file.

        Parameters:
            years (list of int, optional): List of start and end year. Default is [1980, 2010].
            scenario (str, optional): Scenario to consider. Default is 'hist'.

        """
    LOG_FILE = "progress_make_river_flood_global.txt"
    if years is None:
        years = [1980, 2010]
    if scenario == 'hist':
        flddph_data_dir = os.path.join(DATA_DIR, 'river_flood/flood_flddph_hist')
    else:
        flddph_data_dir = os.path.join(DATA_DIR, "".join(['river_flood/flood_flddph/', scenario]))
    log_msg(f"Started computing floods for scenario {scenario}  and years {years}\n", LOG_FILE)
    client = Client()
    centroids = client.get_centroids(res_arcsec_land=150,
    res_arcsec_ocean=1800,
    extent=(-180, 180, -90, 90),)
    years_str = "_".join([str(years[0]), str(years[1])])
    years = range(int(years[0]), int(years[1]))
    filename = OUT_FILE_PATH.format(scenario=scenario, years_str=years_str)
    path = os.path.join(DATA_DIR, 'river_flood', 'global', scenario, years_str)
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
        rf = RiverFlood.from_nc(years=years, dph_path=flddph_file_path, frc_path=fldfrc_file_path, centroids=centroids)
        date_year = [datetime.date.fromordinal(ordinal).year for ordinal in rf.date]
        rf.event_name = ["_".join([str(y), flddph_file.split('_')[3], flddph_file.split('_')[2], scenario]) for y in date_year]
        rf.tag.haz_type = 'RF'
        rf_list.append(rf)
        log_msg(f"computation for one GCM done \n", LOG_FILE)
    rf_concat = rf.concat(rf_list)
    rf_concat.frequency = rf_concat.frequency / len(files_list)
    rf_concat.write_hdf5(file)
    log_msg(f"Computing of river floods for scenario {scenario}  and years {years} done\n", LOG_FILE)


if __name__ == "__main__":
    main(years=[sys.argv[1], sys.argv[2]], scenario=sys.argv[3])