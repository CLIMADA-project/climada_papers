from climada.hazard import TropCyclone
import os
import numpy as np
from config import OUT_DATA_DIR

BASINS = ["SI", "NA", "SP", "WP", "SA", "EP"]


def main(scenarios=None, n_tracks=10):
    if scenarios is None:
        scenarios = ['historical', 'rcp26', 'rcp60', 'rcp45', 'rcp85']
    tracks_str = "".join([str(n_tracks), 'synth_tracks'])
    for scenario in scenarios:
        if scenario == 'historical':
            years = ['1980_2020']
        else:
            years = ['2040', '2060', '2080']

        for year in years:
            tc = TropCyclone()
            path0 = os.path.join(OUT_DATA_DIR, 'tropical_cyclones/genesis_basin/',
                                 tracks_str)
            tc_file = "".join(['tropical_cyclone_'+str(n_tracks)+'synth_tracks_150arcsec_', scenario, '_genesis_NI_', year, '.hdf5'])
            if scenario == 'historical':
                tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_genesis_NI_', year, '.hdf5'])
            path = os.path.join(path0, 'NI')
            if scenario == 'historical':
                tc_file = os.path.join(path, scenario, tc_file)
            else:
                tc_file = os.path.join(path, scenario, year, tc_file)
            tc.read_hdf5(tc_file)
            max_event_id = np.max(tc.event_id)
            for basin in BASINS:
                tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_',scenario,'_genesis_', basin, '_', year,'.hdf5'])
                if scenario == 'historical':
                    tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_genesis_', basin, '_',year,'.hdf5'])

                tc2 = TropCyclone()
                if scenario == 'historical':
                    tc2_file = os.path.join(path0, basin, scenario, tc_file)
                else:
                    tc2_file = os.path.join(path0, basin, scenario, year, tc_file)
                tc2.read_hdf5(tc2_file)
                tc2.event_id = tc2.event_id+max_event_id
                tc2.write_hdf5(tc2_file)
                max_event_id = np.max(tc2.event_id)
                tc.append(tc2)
                tc_file = "".join(['tropical_cyclone_'+str(n_tracks)+'synth_tracks_150arcsec_', scenario, '_global_', year, '.hdf5'])
                if scenario == 'historical':
                    tc_file = "".join(
                        ['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_global_', year, '.hdf5'])
                path = os.path.join(OUT_DATA_DIR,'tropical_cyclones/global/',
                                    tracks_str, scenario, year)

                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                tc.write_hdf5(os.path.join(path, tc_file))


if __name__ == "__main__":
    main(n_tracks=10)
#    main(n_tracks=50)