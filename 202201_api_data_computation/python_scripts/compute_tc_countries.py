import os
from pathlib import Path
from climada.hazard import TropCyclone
from pycountry import countries
from config import OUT_DATA_DIR


def main(years_list=None, scenarios=None, n_tracks=10, replace=True):
    if years_list is None:
        years = [2040, 2060, 2080]
    if scenarios is None:
        scenarios = [26, 60, 45, 85]
    for scenario in scenarios:
        str_scenario = "".join(['rcp',str(scenario)])
        if scenario == 'historical':
            years = ['1980_2020']
            str_scenario = scenario
        for year in years:
            tracks_str = "".join([str(n_tracks), 'synth_tracks'])
            path0 = os.path.join(OUT_DATA_DIR, 'tropical_cyclones_v2')
            path = os.path.join(path0, 'global', tracks_str, str_scenario, str(year))
            for file in os.listdir(path):
                f = file.split('_', 7)
                file = os.path.join(path, file)
                tc = TropCyclone()
                tc.read_hdf5(file)
                path_country = os.path.join(path0, 'countries', tracks_str, str_scenario, str(year))
                isExist = os.path.exists(path_country)
                if not isExist:
                    os.makedirs(path_country)
                for country in countries:
                    if scenario!='historical':
                        file_country = "".join((f[0], '_', f[1], '_', f[2], '_', f[3], '_', f[4], '_', f[5], '_',
                                        country.alpha_3, '_', f[7]))

                    else:
                        file_country = "".join(
                            [f[0], '_', f[1], '_', f[2], '_', f[3], '_', f[4], '_', country.alpha_3, '_',
                             f[6], '_', f[7]])
                    file_country = os.path.join(path_country, file_country)

                    if Path(file_country).exists() and replace is False:
                        continue
                    tc_country = tc.select(reg_id=int(country.numeric))
                    if tc_country is None:
                        continue
                    tc_country.write_hdf5(file_country)


if __name__ == "__main__":
    main(n_tracks=10, scenarios=['historical'])
#    main(n_tracks=50)

