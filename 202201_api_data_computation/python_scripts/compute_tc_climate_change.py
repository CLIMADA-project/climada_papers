import os
from config import OUT_DATA_DIR
import sys
from climada.hazard import TropCyclone


def main(basin='EP', climate_scenarios=None, future_years=None, n_tracks=10):
    if future_years is None:
        future_years = [2040, 2060, 2080]
    if climate_scenarios is None:
        climate_scenarios = [26, 60, 45, 85]
    path0 = os.path.join(OUT_DATA_DIR, 'tropical_cyclones/genesis_basin', str(n_tracks) + 'synth_tracks', basin)
    tracks_str = str(n_tracks) + 'synth_tracks'
    path = os.path.join(path0, "historical")
    file_name = "_".join(['tropical_cyclone', tracks_str, '150arcsec_genesis', basin, '1980_2020.hdf5'])
    tc_haz = TropCyclone.from_hdf5(os.path.join(path, file_name))
    for climate_scenario in climate_scenarios:
        for year in future_years:
            rcp_str = 'rcp' + str(climate_scenario)
            path = os.path.join(path0, rcp_str, str(year))
            file_name = "_".join(['tropical_cyclone', tracks_str, '150arcsec', rcp_str, 'genesis', basin, str(year)])
            file_name = ".".join([file_name, 'hdf5'])
            isExist = os.path.exists(path)
            if not isExist:
                os.makedirs(path)
            tc_haz_future = tc_haz.set_climate_scenario_knu(ref_year=year, rcp_scenario=climate_scenario)
            tc_haz_future.write_hdf5(os.path.join(path, file_name))


if __name__ == "__main__":
    print(sys.argv)
    main(basin=sys.argv[1], n_tracks=sys.argv[2])