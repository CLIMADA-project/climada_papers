import os
from pathlib import Path
from climada.hazard import TropCyclone
from pycountry import countries
from config import DATA_DIR
from create_log_file import log_msg


FILE_NAME = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_rcp{scenario}_{country}_{year}.hdf5'
FILE_NAME_HIST = 'tropical_cyclone_{n_tracks}synth_tracks_150arcsec_{scenario}_{country}_{year}.hdf5'


def main(years_list=None, scenarios=None, n_tracks=10, replace=True):
    if years_list is None:
        years_list = [2040, 2060, 2080]
    if scenarios is None:
        scenarios = [26, 60, 45, 85]
    for scenario in scenarios:
        str_scenario = "".join(['rcp',str(scenario)])
        if scenario == 'historical':
            years_list = ['1980_2020']
            str_scenario = scenario
        for year in years_list:
            tracks_str = "".join([str(n_tracks), 'synth_tracks'])
            path0 = os.path.join(DATA_DIR, 'tropical_cyclones')
            path = os.path.join(path0, 'global', tracks_str, str_scenario, str(year))
            for file in os.listdir(path):
                file = os.path.join(path, file)
                tc = TropCyclone()
                tc.read_hdf5(file)
                path_country = os.path.join(path0, 'countries', tracks_str, str_scenario, str(year))
                isExist = os.path.exists(path_country)
                if not isExist:
                    os.makedirs(path_country)
                for country in countries:
                    if scenario != 'historical':
                        file_country = FILE_NAME.format(scenario=scenario, year=year, country=country.alpha_3,
                                                        n_tracks=n_tracks)

                    else:
                        file_country = FILE_NAME_HIST.format(scenario=scenario, year=year, country=country.alpha_3,
                                                             n_tracks=n_tracks)
                    file_country = os.path.join(path_country, file_country)

                    if Path(file_country).exists() and replace is False:
                        continue
                    tc_country = tc.select(reg_id=int(country.numeric))
                    if tc_country is None:
                        continue
                    tc_country.write_hdf5(file_country)


if __name__ == "__main__":
    main(n_tracks=10)
