from climada.hazard import TropCyclone
import os
import numpy as np
from config import DATA_DIR
from create_log_file import log_msg


BASINS = ["SI", "NA", "SP", "WP", "SA", "EP"]
LOG_FILE = "concatenate_basins_tc.txt"
FILE_NAME = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_{scenario}_genesis_{basin}_{year}.hdf5'
FILE_NAME_HIST = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_genesis_{basin}_{year}.hdf5'
FILE_NAME_GLOBAL = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_{scenario}_global_{year}.hdf5'
FILE_NAME_GLOBAL_HIST = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_global_{year}.hdf5'


def main(climate_scenarios=None, n_tracks=10, years=None):
    if climate_scenarios is None:
        climate_scenarios = ['historical', 'rcp26', 'rcp60', 'rcp45', 'rcp85']
    tracks_str = "".join([str(n_tracks), 'synth_tracks'])
    for scenario in climate_scenarios:
        if years is None:
            if scenario == 'historical':
                years = ['1980_2020']
            else:
                years = ['2040', '2060', '2080']

        for year in years:
            log_msg(f"Starting concatenating basins for year {year} and scenario {scenario}\n", LOG_FILE)

            tc = TropCyclone()
            path0 = os.path.join(DATA_DIR, 'tropical_cyclones/genesis_basin/',
                                 tracks_str)

            path = os.path.join(path0, 'NI')
            if scenario == 'historical':
                tc_file = FILE_NAME_HIST.format(n_tracks=n_tracks, basin='NI', year=year)
                tc_file_path = os.path.join(path, scenario, tc_file)
            else:
                tc_file = FILE_NAME.format(n_tracks=n_tracks, basin='NI', year=year, scenario=scenario)
                tc_file_path = os.path.join(path, scenario, year, tc_file)
            tc.read_hdf5(tc_file_path)
            max_event_id = np.max(tc.event_id)
            for basin in BASINS:
                if scenario == 'historical':
                    tc2_file = FILE_NAME_HIST.format(n_tracks=n_tracks, basin=basin, year=year)
                    tc2_file_path = os.path.join(path0, basin, scenario, year, tc2_file)
                else:
                    tc2_file = FILE_NAME.format(n_tracks=n_tracks, basin=basin, year=year, scenario=scenario)
                    tc2_file_path = os.path.join(path0, basin, scenario, year, tc2_file)
                tc2 = TropCyclone()
                tc2.read_hdf5(tc2_file_path)
                tc2.event_id = tc2.event_id+max_event_id #keep same event ids as in the basin files
                tc2.write_hdf5(tc2_file_path)
                max_event_id = np.max(tc2.event_id)
                tc.append(tc2)
                tc_file = FILE_NAME_GLOBAL.format(n_tracks=n_tracks, scenario=scenario, year=year)
                if scenario == 'historical':
                    tc_file = FILE_NAME_GLOBAL_HIST.format(n_tracks=n_tracks, scenario=scenario, year=year)
                path = os.path.join(DATA_DIR, 'tropical_cyclones/global/',
                                    tracks_str, scenario, year)

                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                tc.write_hdf5(os.path.join(path, tc_file))
                log_msg(f"Finished concatenating basins for year {year} and scenario {scenario}\n", LOG_FILE)


if __name__ == "__main__":
    main(n_tracks=10, years=['2060', '2080'], climate_scenarios=['rcp85'])
