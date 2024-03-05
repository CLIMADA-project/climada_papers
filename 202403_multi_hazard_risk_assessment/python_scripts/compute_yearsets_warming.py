import logging
import numpy as np
import sys

from climada.util.api_client import Client
from climada.engine import Impact
from create_log_msg import log_msg
from compute_tc_mit import make_tc_hazard
from config import BASE_DIR, WARMING_LEVELS, EXPOSURE_TYPES, HAZARD_TYPES
from multi_risk import *

OUTPUT_DIR = BASE_DIR / "impact_yearsets/global"
INPUT_DIR = BASE_DIR / "impacts/global"
LOG_FILE = "../logs/yearsets_warming.txt"

correction_factor_assets = {'TC': 1, 'RF': 48.3 / 481}


def load_exposure_data(client, exposure_type):
    properties = {'res_arcsec': '150', 'spatial_coverage': 'global', 
                  'exponents': '(1,1)' if exposure_type == 'assets' else '(0,1)',
                  'fin_mode': 'pc' if exposure_type == 'assets' else None}
    exposure = client.get_exposures('litpop', properties=properties, version='v2')
    exposure.gdf = exposure.gdf.dropna().reset_index(drop=True)
    return exposure

def load_and_prepare_impact(hazard_type, exposure_str, warming):
    file_path = f"{INPUT_DIR}/csv/{exposure_str}/{hazard_type}_impact_{exposure_str}_150arcsec_{warming}_global.csv"
    impact = Impact.from_csv(file_path)
    impact.imp_mat = Impact.read_sparse_csr(file_path.replace('csv', 'npz'))
    return impact

def apply_corrections_and_mask(impacts_per_year):
    for hazard_type, correction_factor in correction_factor_assets.items():
        impacts_per_year[hazard_type].imp_mat *= correction_factor
        mask = np.abs(impacts_per_year[hazard_type].imp_mat.data) < 100
        impacts_per_year[hazard_type].imp_mat.data[mask] = 0
        impacts_per_year[hazard_type].imp_mat.eliminate_zeros()

def common_elements(list1, list2):
    return [element for element in list1 if element in list2]

def compute_yearsets(exposure_str, warming):
    impacts_per_year = {}
    impacts = {}
    for hazard_type in HAZARD_TYPES:
        impacts[hazard_type] = load_and_prepare_impact(hazard_type, exposure_str, warming)
    impacts_per_year['TC'] = aggregate_impact_from_event_name(
        impacts['TC'])
    impacts_per_year['RF'] = impacts['RF']
    logging.info("Impact objects loaded")

    if exposure_str == 'assets':
        apply_corrections_and_mask(impacts_per_year)

    # Additional processing such as event name manipulation and impact ordering can go here
    event_names = [["_".join(event.split("_")[0:3]) for event in impacts_per_year[hazard_type].event_name] for
                   hazard_type
                   in HAZARD_TYPES]

    event_names = np.unique(common_elements(event_names[0], event_names[1]))
    log_msg(f"len common event names is {len(event_names)} \n", LOG_FILE)
    log_msg(f"shape tc imp_mat {impacts_per_year['TC'].imp_mat.shape} \n", LOG_FILE)


    # order the impacts by event name
    # write the impacts to CSV and NPZ files
    log_msg(f"starting combining events \n", LOG_FILE)

    impacts_yearsets = order_by_event_name(impacts_per_year, 2, list_string=event_names)


    for hazard_type in HAZARD_TYPES:
        csv_path = OUTPUT_DIR / f"csv/{exposure_str}/{hazard_type}_{exposure_str}_impacts_yearsets_150arcsec_{warming}_global.csv"
        npz_path = OUTPUT_DIR / f"npz/{exposure_str}/{hazard_type}_{exposure_str}_impacts_yearsets_150arcsec_{warming}_global.npz"
        impacts_yearsets[hazard_type].write_csv(csv_path)
        impacts_yearsets[hazard_type].write_sparse_csr(npz_path)
        log_msg(f"Processed and saved data for {hazard_type}, warming level {warming}, exposure {exposure_str}.\n", LOG_FILE)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    for warming in WARMING_LEVELS:
        for exposure_str in EXPOSURE_TYPES:
            compute_yearsets(exposure_str, warming)
