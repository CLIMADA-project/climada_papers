import itertools
import logging
import numpy as np
import copy
from climada.util.api_client import Client
from climada.engine import Impact
from create_log_msg import log_msg
from multi_risk import *
from sklearn.utils import shuffle

from config import BASE_DIR, WARMING_LEVELS, EXPOSURE_TYPES, HAZARD_TYPES

INPUT_DIR = BASE_DIR / "impact_yearsets"

OUTPUT_DIR = BASE_DIR / "impact_compound"
LOG_FILE = 'compute_combined_impact_logs.txt'
OCCUR_TOGETHER = {'pop': True, 'assets': True}



def generate_file_paths(hazard, exposure, warming, file_type='csv'):
    """
    Generate file paths for CSV and NPZ files based on hazard, exposure, and warming.
    
    Parameters:
    - hazard (str): The type of hazard (e.g., 'TC', 'RF').
    - exposure (str): The type of exposure (e.g., 'pop', 'assets').
    - warming (str): The warming scenario (e.g., '1', '2').
    - file_type (str): Type of the file ('csv' or 'npz'). Defaults to 'csv'.
    
    Returns:
    - A string representing the file path.
    """
    base_path = INPUT_DIR / f"global/{file_type}/{exposure}/"
    file_name = f"{hazard}_{exposure}_impacts_yearsets_150arcsec_{warming}_global_caped.{file_type}"
    return base_path / file_name

def load_impacts(hazard, exposure, warming):
    """
    Load impacts from CSV and NPZ files for a given hazard, exposure, and warming scenario.
    
    Parameters:
    - hazard (str): The type of hazard.
    - exposure (str): The type of exposure.
    - warming (str): The warming scenario.
    
    Returns:
    - An Impact object loaded with data from the files.
    """
    csv_path = generate_file_paths(hazard, exposure, warming, 'csv')
    npz_path = generate_file_paths(hazard, exposure, warming, 'npz')

    impact = Impact.from_csv(csv_path)
    impact.imp_mat = Impact.read_sparse_csr(npz_path)
    return impact


def combine_and_save_impacts(impacts, exposure, warming, combination, occur_together, order_str = 'ordered_', how='sum'):
    """
    Combine impacts for given hazards and save the combined impact to CSV and NPZ files.
    
    Parameters:
    - impacts (dict): Dictionary of Impact objects with exposure types as keys.
    - exposure (str): The exposure type.
    - warming (str): The warming level.
    - combination (tuple): Tuple of hazard types to combine.
    - occur_together (bool): Indicates whether events from different hazards occur together.
    - how (str): Method to combine impacts, default is 'sum'.
    """
    combined_impact = combine_yearsets([impacts[hazard] for hazard in combination], how=how, occur_together=occur_together)
    combined_impact.unit = impacts[list(combination)[0]].unit

    # Construct file paths for saving combined impacts
    combi_str = "_".join(combination)
    file_identifiers = f"{exposure}_combined_impact_{order_str}{combi_str}_150arcsec_{warming}_global"

    combined_impact.write_csv(OUTPUT_DIR / f"global/csv/{exposure}/{file_identifiers}.csv")
    combined_impact.write_sparse_csr(OUTPUT_DIR / f"global/npz/{exposure}/{file_identifiers}.npz")

    
def compute_compound(exposure, warming):
    impacts_yearsets_ordered = {hazard: load_impacts(hazard, exposure, warming) for hazard in ['TC', 'RF']}
    impacts_yearsets = copy.deepcopy(impacts_yearsets_ordered)
    # Randomly shuffle RF events
    impacts_yearsets['RF'] = order_events_by_indices(
        impacts_yearsets_ordered['RF'], shuffle(np.arange(impacts_yearsets_ordered['RF'].imp_mat.shape[0])))


    combinations = list(itertools.combinations(['TC', 'RF'], 2))
    for combi in combinations:
        log_msg(f"Start combining {combi} for exposure {exposure} and warming {warming}.\n", LOG_FILE)
        #combine_and_save_impacts(impacts_yearsets_ordered, exposure, warming, combi, OCCUR_TOGETHER[exposure])
        combine_and_save_impacts(impacts_yearsets, exposure, warming, combi, OCCUR_TOGETHER[exposure], order_str='')


if __name__ == "__main__":
    exposures = ['pop', 'assets']
    warmings = ['1', '2']
    for exposure, warming in itertools.product(EXPOSURE_TYPES, WARMING_LEVELS):
        compute_compound(exposure, warming)
