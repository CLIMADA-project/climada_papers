import os
from pycountry import countries
from climada.entity import LitPop
from climada.entity import Exposures
from climada.hazard import Centroids

from config import OUT_DATA_DIR




missing_country = []

data_dir_countries = os.path.join(OUT_DATA_DIR, 'litpop', 'countries', 'pop')
data_dir = os.path.join(OUT_DATA_DIR, 'litpop', 'global','pop')
exposures_list = []
for country in countries:
   try:
       litpop = LitPop.from_countries(country.alpha_3, res_arcsec=150, exponents=(0, 1), fin_mode='pop')
       litpop.check()
       if litpop.gdf.value.isna().sum()<1:
           exposures_list.append(litpop)
           litpop.write_hdf5("".join([data_dir_countries,'/LitPop_pop_150arcsec_',country.alpha_3,'.hdf5']))
   except:
       missing_country.append(country.alpha_3)
exposures_concat = LitPop.from_hdf5("/Users/szelie/OneDrive - ETH Zurich/data/climada_api/litpop/global/pop/LitPop_pop_150arcsec_global.hdf5")
exposures_concat.gdf = exposures_concat.gdf.dropna()
exposures_concat.gdf = exposures_concat.gdf.set_geometry(exposures_concat.gdf.geometry)
exposures_concat.set_lat_lon()
cent = Centroids.from_lat_lon(exposures_concat.gdf.latitude, exposures_concat.gdf.longitude)
cent.set_lat_lon_to_meta()
exposures_concat.meta = cent.meta
exposures_concat.tag.description = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pop, exp: (0, 1), admin1_calc: False'
exposures_concat.write_hdf5(os.path.join(data_dir, 'LitPop_pop_150arcsec.hdf5'))

data_dir_countries = os.path.join(OUT_DATA_DIR, 'litpop', 'countries', 'assets')
data_dir = os.path.join(OUT_DATA_DIR, 'litpop', 'global', 'assets')
exposures_list = []
for country in countries:
   if country in ['ATA', 'IOT', 'ATF', 'LBY', 'SGS', 'SYR']:
       continue
   try:
       litpop = LitPop.from_countries(country.alpha_3, res_arcsec=150, exponents=(3, 0), fin_mode='pc')
       litpop.check()
       if litpop.gdf.value.isna().sum()<1:
           exposures_list.append(litpop)
           litpop.write_hdf5("".join([data_dir_countries,'/LitPop_assets_pc_150arcsec_',country.alpha_3,'.hdf5']))
   except:
       missing_country.append(country.alpha_3)

exposures_concat = LitPop.from_hdf5("/Users/szelie/OneDrive - ETH Zurich/data/climada_api/litpop/global/default/LitPop_150arcsec_global.hdf5")
exposures_concat.gdf = exposures_concat.gdf.dropna()

#exposures_concat = LitPop.concat(exposures_list)
exposures_concat.gdf = exposures_concat.gdf.set_geometry(exposures_concat.gdf.geometry)
exposures_concat.set_lat_lon()
cent = Centroids.from_lat_lon(exposures_concat.gdf.latitude, exposures_concat.gdf.longitude)
cent.set_lat_lon_to_meta()
exposures_concat.meta = cent.meta
exposures_concat.tag.description = 'Global LitPop Exposure at 150 as, year: 2018, financial mode: pc, exp: (3, 0), admin1_calc: False'
exposures_concat.write_hdf5(os.path.join(data_dir, 'LitPop_assets_pc_150arcsec.hdf5'))

