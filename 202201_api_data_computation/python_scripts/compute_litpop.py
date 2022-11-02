import os
from pycountry import countries
from climada.entity import LitPop
from climada.hazard import Centroids
from config import DATA_DIR
from create_log_file import log_msg

missing_country = []

OUT_DIR_COUNTRIES = os.path.join(DATA_DIR, 'litpop', 'countries')
OUT_DIR = os.path.join(DATA_DIR, 'litpop', 'global', 'pop')
EXPONENTS = {'pop': '(0,1)', 'default': '(1,1)', 'assets': '(3,0)'}
FIN_MODE = {'pop': 'pop', 'default': 'pc', 'assets': 'pc'}
EXP_STR = {'pop': 'pop', 'default': '', 'assets': 'assets_pc'}

OUT_FILE_COUNTRY = 'LitPop_{exposure}_150arcsec_{country}.hdf5'
OUT_FILE = 'LitPop_{exposure}_150arcsec.hdf5'
LOG_FILE = "progress_make_litpop.txt"


def make_litpop(exposure):
    exposures_list = []
    for country in countries:
        try:
            litpop = LitPop.from_countries(country.alpha_3, res_arcsec=150, exponents=EXPONENTS[exposure],
                                           fin_mode=FIN_MODE[exposure])
            litpop.check()
            notnull = litpop.gdf['values'].notnull()
            litpop.gdf = litpop.gdf[notnull]
            exposures_list.append(litpop)
            litpop.write_hdf5(os.path.join(OUT_DIR_COUNTRIES, OUT_FILE_COUNTRY.
                                           format(exposure=EXP_STR[exposure],country=country)))
            log_msg(f"Country {country} done.\n", LOG_FILE)
        except ValueError:
            missing_country.append(country.alpha_3)
            log_msg(f"Country {country} failed.\n", LOG_FILE)
    log_msg(f"Start creating global exposure.\n", LOG_FILE)
    exposures_concat = LitPop.concat(exposures_list)
    exposures_concat.gdf = exposures_concat.gdf.set_geometry(exposures_concat.gdf.geometry)
    exposures_concat.set_lat_lon()
    cent = Centroids.from_lat_lon(exposures_concat.gdf.latitude, exposures_concat.gdf.longitude)
    cent.set_lat_lon_to_meta()
    exposures_concat.meta = cent.meta
    exposures_concat.tag.description = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pop, exp: (0, 1), ' \
                                       'admin1_calc: False'
    exposures_concat.write_hdf5(os.path.join(OUT_DIR, exposure, OUT_FILE.
                                           format(exposure=EXP_STR[exposure])))

    log_msg(f"Global file saved.\n", LOG_FILE)
    log_msg(f"The following countries were not successfull: {missing_country}.\n", LOG_FILE)



