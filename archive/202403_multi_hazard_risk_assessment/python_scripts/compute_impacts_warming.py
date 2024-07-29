import os
import pandas as pd
import pycountry
import numpy as np
import scipy as sp
from climada.util.api_client import Client
from climada_petals.entity.impact_funcs.river_flood import ImpfRiverFlood, flood_imp_func_set, RIVER_FLOOD_REGIONS_CSV
from climada.engine.impact import Impact
from climada.hazard import Hazard
from define_impf import impact_fun_set_pop, impact_function_set_assets, get_impf_id
from compute_tc_mit import make_tc_hazard
from config import BASE_DIR, WARMING_LEVELS, EXPOSURE_TYPES, HAZARD_TYPES
# 
OUTPUT_DIR = BASE_DIR / "impacts_multi_risk"
PATH_HAZARD = {
    "TC": str(BASE_DIR / "hazards/hazard_TC_MIT/hazard_tc_{warming}.hdf5"),
    "RF": str(BASE_DIR / 'hazards/flood_warming_level/global/{warming}/river_flood_150arcsec_{warming}.hdf5')
}

def get_exposure(client, exposure_str, country="global"):
    properties = {
        "res_arcsec": "150",
        "exponents": "(0,1)" if exposure_str == "pop" else "(1,1)",
        "spatial_coverage": "global"
    }

    exposure = client.get_exposures("litpop", properties=properties)
    exposure.gdf = exposure.gdf.dropna().reset_index(drop=True)
    return exposure

def process_impf(exposure, exposure_str):
    regions_df = pd.read_csv(RIVER_FLOOD_REGIONS_CSV)
    imp_fun_set = impact_fun_set_pop() if exposure_str == 'pop' else impact_function_set_assets()
    exposure.gdf['impf_TC'] = 1
    exposure.gdf['impf_RF'] = 1
    if exposure_str == "assets":
        for cnt in regions_df['ISO']:
            try:
                country_numeric = int(pycountry.countries.get(alpha_3=cnt).numeric)
                exposure.gdf.loc[exposure.gdf['region_id'] == country_numeric, 'impf_RF'] =\
                    int(regions_df[regions_df['ISO'] == cnt]['impf_RF'])
            except:
                continue
        for cnt in np.unique(exposure.gdf.region_id):
            exposure.gdf.loc[exposure.gdf['region_id'] == cnt, 'impf_TC'] = get_impf_id(int(cnt))[1]
    return imp_fun_set, exposure

def compute_impacts_warming_level(hazard_type, exposure_str, warming_level):
    client = Client()
    exposure = get_exposure(client, exposure_str)
    imp_fun_set, exposure = process_impf(exposure, exposure_str)

    hazard = Hazard.from_hdf5(PATH_HAZARD[hazard_type].format(warming=warming_level))
    if hazard_type == 'RF':
        events_rf = [event_name for event_name in hazard.event_name if 'rcp85' not in event_name] # remove all rcp85
        hazard.select(event_names=events_rf)
    
    impact = Impact()
    impact.calc(exposure, imp_fun_set, hazard, save_mat=True)

    for ext in ['csv', 'npz']:
        file_output = str(BASE_DIR / f"impacts/global/{ext}/{hazard_type}_impact_{exposure_str}_150arcsec_{warming_level}_global.{ext}")

        if ext == 'csv':
            impact.write_csv(file_output)
        else:
            impact.write_sparse_csr(file_output)

if __name__ == "__main__":
    for warming_level in WARMING_LEVELS:
        for exposure_str in EXPOSURE_TYPES:
            for hazard_type in HAZARD_TYPES:
                compute_impacts_warming_level(hazard_type, exposure_str, warming_level)
