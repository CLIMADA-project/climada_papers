from climada.hazard import TropCyclone
import os

from python_scripts.config import OUT_DATA_DIR


def main(years=['2040', '2060', '2080'], scenarios=['historical','rcp26', 'rcp60', 'rcp45', 'rcp85'], n_tracks=10):
    tracks_str = "".join([str(n_tracks), 'synth_tracks'])
    for scenario in scenarios:
        if scenario == 'historical':
            years = ['1980_2020']
        for year in years:
            tc = TropCyclone()
            path0 = os.path.join(OUT_DATA_DIR, '/tropical_cyclones/genesis_basin/',
                                 tracks_str)
            tc_file = "".join(['tropical_cyclone_'+str(n_tracks)+'synth_tracks_150arcsec_', scenario, '_genesis_NI_', year, '.hdf5'])
            if scenario == 'historical':
                tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_genesis_NI_', year, '.hdf5'])
            path = os.path.join(path0, 'NI')
            tc.read_hdf5(os.path.join(path, scenario, year, tc_file))

            for basin in ["SI", "NA", "SP", "WP", "SA", "EP"]:

                tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_',scenario,'_genesis_', basin, '_', year,'.hdf5'])
                if scenario == 'historical':
                    tc_file = "".join(['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_genesis_', basin, '_',year,'.hdf5'])

                tc2 = TropCyclone()
                tc2.read_hdf5(os.path.join(path0, basin, scenario, year, tc_file))
                tc.append(tc2)
                tc_file = "".join(['tropical_cyclone_'+str(n_tracks)+'synth_tracks_150arcsec_', scenario, '_global_', year, '.hdf5'])
                if scenario == 'historical':
                    tc_file = "".join(
                        ['tropical_cyclone_', str(n_tracks), 'synth_tracks_150arcsec_genesis_global_', year, '.hdf5'])
                path = os.path.join('DATA_DIR','tropical_cyclones/global/',
                                    tracks_str, scenario, year)

                isExist = os.path.exists(path)
                if not isExist:
                    os.makedirs(path)
                tc.write_hdf5(os.path.join(path, tc_file))


if __name__ == "__main__":
    main(n_tracks=50)
    main(n_tracks=10)