import os
from pycountry import countries
from climada.entity import LitPop
from climada.entity import Exposures

from config import OUT_DATA_DIR

data_dir_countries = os.path.join(OUT_DATA_DIR, 'litpop', 'countries', 'default')
data_dir = os.path.join(OUT_DATA_DIR, 'litpop', 'global', 'default')
missing_country = []
exposures_list = []
for country in countries:
    try:
        exposures_list.append(LitPop.from_hdf5("".join([data_dir_countries,'/LitPop_150arcsec_',country.alpha_3,'.hdf5'])))
    except KeyError as e:
        missing_country.append(country)
exposures_concat = Exposures()
exposures_concat = exposures_concat.concat(exposures_list)
exposures_concat.tag.description = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pc, exp: (1, 1), admin1_calc: False'
exposures_concat.write_hdf5(os.path.join(data_dir,'LitPop_150arcsec.hdf5'))


data_dir_countries = os.path.join(OUT_DATA_DIR, 'litpop', 'countries', 'pop')
data_dir = os.path.join(OUT_DATA_DIR, 'litpop', 'pop')
exposures_list = []
for country in countries:
    exposures_list.append(LitPop.from_hdf5(data_dir_countries,'/LitPop_pop_150arcsec_',country.alpha_3,'.hdf5'))
exposures_concat = Exposures()
exposures_concat = exposures_concat.concat(exposures_list)
exposures_concat.tag = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pop, exp: (0, 1), admin1_calc: False'
exposures_concat.write_hdf5(os.path.join(data_dir,'LitPop_pop_150arcsec.hdf5'))


data_dir_countries = os.path.join(OUT_DATA_DIR, 'litpop', 'countries', 'nightlight')
data_dir = os.path.join(OUT_DATA_DIR, 'litpop', 'global', 'nightlight')
exposures_list = []
for country in countries:
    exposures_list.append(LitPop.from_hdf5("".join([data_dir_countries,'/LitPop_pop_150arcsec_',country.alpha_3,'.hdf5'])))
exposures_concat = Exposures()
exposures_concat = exposures_concat.concat(exposures_list)
exposures_concat.tag = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pc, exp: (1, 0), admin1_calc: False'
exposures_concat.write_hdf5(os.path.join(data_dir, 'LitPop_nightlight_150arcsec.hdf5'))
