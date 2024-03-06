import itertools
import logging
import numpy as np
import copy
from climada.engine import Impact
from create_log_msg import log_msg
from sklearn.utils import shuffle
from climada.util.api_client import Client
from multi_risk import *
from config import BASE_DIR


INPUT_DIR = BASE_DIR / "impact_yearsets"
LOG_FILE = '../logs/compute_combined_impact_logs.txt'
OCCUR_TOGETHER = False

def get_exposure_data(client, exposure_type, cap_at_exp):
    """
    Fetch and prepare exposure data based on the specified type and capping option.
    
    Parameters:
    - client (Client): An instance of the API client to fetch exposure data.
    - exposure_type (str): Type of exposure ('pop' or 'assets').
    - cap_at_exp (bool): Indicates whether the impact data should be capped at exposure levels.
    
    Returns:
    - dict: A dictionary containing the exposure data.
    """
    properties = {'res_arcsec': '150', 'spatial_coverage': 'global'}
    if exposure_type == 'assets':
        properties.update({'exponents': '(1,1)', 'fin_mode': 'pc'})
    else:  # 'pop'
        properties.update({'exponents': '(0,1)'})
    
    exposure_data = client.get_exposures('litpop', properties=properties, version='v2')
    exposure_data.gdf = exposure_data.gdf.dropna().reset_index(drop=True)
    return {'pop': exposure_data if exposure_type == 'pop' else None,
            'assets': exposure_data if exposure_type == 'assets' else None}

def load_impact_data(hazard, exposure, warming, cap_str):
    """
    Load impact data for a specific hazard, exposure, and warming scenario.
    
    Parameters:
    - hazard (str): Type of hazard ('TC' or 'RF').
    - exposure (str): Type of exposure ('pop' or 'assets').
    - warming (str): Warming scenario ('1' or '2').
    - cap_str (str): String indicating whether impacts are capped ('_caped' or '').
    
    Returns:
    - Impact: An instance of the Impact class loaded with data.
    """
    if hazard == "RF":
        cap_str = "" #RF yearsets do not need to be caped, as they anyways are yearly impacts to start with
    csv_path = INPUT_DIR /  f"global/csv/{exposure}/{hazard}_{exposure}_impacts_yearsets_150arcsec_{warming}_global{cap_str}.csv"
    npz_path = INPUT_DIR / f"global/npz/{exposure}/{hazard}_{exposure}_impacts_yearsets_150arcsec_{warming}_global{cap_str}.npz"
    
    impact = Impact.from_csv(csv_path)
    impact.imp_mat = Impact.read_sparse_csr(npz_path)
    return impact

def save_combined_impact(combined_impact, exposure, combi, warming, cap_str):
    """
    Save the combined impact data to CSV and NPZ files.
    
    Parameters:
    - combined_impact (Impact): The combined impact object to be saved.
    - exposure (str): Type of exposure.
    - combi (tuple): Tuple of combined hazards.
    - warming (str): Warming scenario.
    - cap_str (str): String indicating whether impacts are capped.
    """
    combi_str = "_".join(combi)
    csv_filename = BASE_DIR / f"impact_aggr/global/csv/{exposure}/{exposure}_combined_impact_{combi_str}_150arcsec_{warming}_global{cap_str}.csv"
    npz_filename = BASE_DIR / f"impact_aggr/global/npz/{exposure}/{exposure}_combined_impact_{combi_str}_150arcsec_{warming}_global{cap_str}.npz"
    
    combined_impact.write_csv(csv_filename)
    combined_impact.write_sparse_csr(npz_filename)
    log_msg(f"Saved combined impact for {combi_str} under exposure {exposure}, warming {warming}.\n", LOG_FILE)

def combine_impacts(impacts_yearsets_ordered, exposure_dict, exposure, warming, cap_str):
    """
    Combine impacts for different hazards, log messages, and save the results.
    
    Parameters:
    - impacts_yearsets_ordered (dict): Dictionary of ordered Impact objects.
    - exposure_dict (dict): Dictionary containing exposure data.
    - exposure (str): Type of exposure.
    - warming (str): Warming scenario.
    - cap_str (str): String indicating whether impacts are capped.
    """
    hazards = ['TC', 'RF']
    combinations = list(itertools.combinations(hazards, 2))
    
    for combi in combinations:
        log_msg(f"Start combining {combi} for exposure {exposure}.\n", LOG_FILE)
        
        combined_impact_ordered = combine_yearsets(
            [impacts_yearsets_ordered[c] for c in combi], how='sum',
            occur_together=OCCUR_TOGETHER, exposures=exposure_dict[exposure]
        )
        
        combined_impact_ordered.unit = impacts_yearsets_ordered[combi[0]].unit
        save_combined_impact(combined_impact_ordered, exposure, combi, warming, cap_str)

def compute_aggr_impact(exposure, warming, cap_at_exp=True):
    """
    Adjusted main function to compute combined impacts for specific exposure and warming settings.
    
    Parameters:
    - exposure (str): The exposure setting to compute impacts for.
    - warming (str): The warming scenario to compute impacts for.
    - cap_at_exp (bool): Flag to indicate whether to cap impacts at exposure levels each year. Defaults to True.
    """
    client = Client()
    
    cap_str = "_caped" if cap_at_exp else ""
    exposure_dict = get_exposure_data(client, exposure, cap_at_exp)
    
    impacts_yearsets_ordered = {
        hazard: load_impact_data(hazard, exposure, warming, cap_str) for hazard in ['TC', 'RF']
    }
    
    # Your existing logic here to combine impacts and log
    combine_impacts(impacts_yearsets_ordered, exposure_dict, exposure, warming, cap_str)

if __name__ == "__main__":
    for exposure in ['assets', 'pop']:  # Example, add other exposures as needed
        for warming in ['1', '2']:  # Example warming levels
            compute_aggr_impact(exposure, warming, cap_at_exp=True)

