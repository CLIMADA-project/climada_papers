import os
from config import DATA_DIR
import sys
from climada.hazard import TropCyclone
from create_log_file import log_msg

FILE_NAME = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_rcp{scenario}_genesis_{basin}_{year}.hdf5'
FILE_NAME_HIST = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_genesis_{basin}_{year}.hdf5'
LOG_FILE = "concatenate_basins_tc.txt"


def main(basin='EP', climate_scenarios=None, future_years=None, n_tracks=10):
    if future_years is None:
        future_years = [2040, 2060, 2080]
    if climate_scenarios is None:
        climate_scenarios = [85]
    path0 = os.path.join(DATA_DIR, 'tropical_cyclones/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
    path = os.path.join(path0, "historical")
    file_name = FILE_NAME_HIST.format(n_tracks=n_tracks,
    basin=basin, year='1980_2020')
    tc_haz = TropCyclone.from_hdf5(os.path.join(path, file_name))
    for climate_scenario in climate_scenarios:
        for year in future_years:
            log_msg(f"Started computing climate change for scenario {climate_scenario} and "
                    f"year {year}.\n", LOG_FILE)

            rcp_str = 'rcp' + str(climate_scenario)
            path = os.path.join(path0, rcp_str, str(year))
            file_name = FILE_NAME.format(n_tracks=n_tracks, basin=basin, year=year,scenario=climate_scenario)
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            tc_haz_future = tc_haz.set_climate_scenario_knu(ref_year=year, rcp_scenario=climate_scenario)
            tc_haz_future.write_hdf5(os.path.join(path, file_name))

            log_msg(f"Finished computing climate change for scenario {climate_scenario} and "
                    f"year {year}.\n", LOG_FILE)


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])