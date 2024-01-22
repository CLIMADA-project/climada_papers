import os
from config import DATA_DIR
import sys
from climada.hazard import TropCyclone, Centroids, TCTracks
from create_log_file import log_msg
import numpy as np

FILE_NAME = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_rcp{scenario}_genesis_{basin}_{year}.hdf5'
FILE_NAME_HIST = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_genesis_{basin}_{year}.hdf5'
LOG_FILE = "compute_tc_cc.txt"


def main(basin='EP', climate_scenarios=None, future_years=None, n_tracks=10):
    """
    Compute climate change scenarios for tropical cyclones for specified basins, scenarios, and years.

    The function reads the historical tropical cyclone hazard data, applies the climate change scenarios,
    and saves the modified hazard data.

    Parameters:
        basin (str, optional): Name of the tropical cyclone genesis basin. Default is 'EP'.
        climate_scenarios (list of int, optional): List of climate scenarios (RCPs) to consider. Default is [85].
        future_years (list of int, optional): List of years to consider for the climate change scenarios. Default is [2040, 2060, 2080].
        n_tracks (int, optional): Number of synthetic tracks. Default is 10.
    """
    if future_years is None:
        future_years = [2040, 2060, 2080]
    if climate_scenarios is None:
        climate_scenarios = [85]
    path0 = os.path.join(DATA_DIR, 'tropical_cyclones_v3/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
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
            file_name = FILE_NAME.format(n_tracks=str(n_tracks), basin=basin, year=str(year), scenario=str(climate_scenario))
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            tc_haz_future = tc_haz.apply_climate_scenario_knu(ref_year=year, rcp_scenario=climate_scenario)
            tc_haz_future.fraction.data = np.zeros(len(tc_haz_future.fraction.data))
            tc_haz_future.fraction.eliminate_zeros()
            path = os.path.join(path0, "".join(["rcp",str(climate_scenario)]))
            tc_haz_future.write_hdf5(os.path.join(path, file_name))

            log_msg(f"Finished computing climate change for scenario {str(climate_scenario)} and "
                    f"year {str(year)}.\n", LOG_FILE)


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])