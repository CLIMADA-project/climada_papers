import itertools
import logging
import numpy as np
import scipy as sp
from sklearn.utils import shuffle
from climada.util.api_client import Client
from climada.engine import Impact
from create_log_msg import log_msg
import copy
from config import BASE_DIR, WARMING_LEVELS, EXPOSURE_TYPES


LOG_FILE = "log/yearsets_Cc.txt"  # Ensure this is correctly defined

def cap_yearsets(exposures_str, warming):
    client = Client()
    if exposures_str=='pop':
        exposures = client.get_exposures('litpop', properties={'res_arcsec': '150', 'exponents': '(0,1)', 'spatial_coverage': 'global'}, version='v2')
    if exposures_str == 'assets':
        exposures = client.get_exposures('litpop', properties={'res_arcsec': '150', 'exponents': '(1,1)', 'fin_mode': 'pc', 'spatial_coverage': 'global'}, version='v2')

    exposures.gdf = exposures.gdf.dropna().reset_index(drop=True)
    impacts_yearsets_ordered = {}

    for hazard in ['TC']:
        csv_file = BASE_DIR / f"impact_yearsets/global/csv/{exposures_str}/{hazard}_{exposures_str}_impacts_yearsets_150arcsec_{warming}_global.csv"
        npz_file = BASE_DIR / f"impact_yearsets/global/npz/{exposures_str}/{hazard}_{exposures_str}_impacts_yearsets_150arcsec_{warming}_global.npz"

        impacts_yearsets_ordered[hazard] = Impact.from_csv(csv_file)
        impacts_yearsets_ordered[hazard].imp_mat = Impact.read_sparse_csr(npz_file)

        imp_mat = impacts_yearsets_ordered[hazard].imp_mat

        m1 = imp_mat.data
        m2 = exposures.gdf.value[imp_mat.nonzero()[1]]
        imp_mat = sp.sparse.csr_matrix((np.minimum(m1, m2), imp_mat.indices, imp_mat.indptr), shape=imp_mat.shape)

        imp_mat.eliminate_zeros()

        imp = copy.deepcopy(impacts_yearsets_ordered[hazard])
        # Assuming set_imp_mat is a function to update imp's imp_mat, need to define it or update accordingly
        imp.imp_mat = imp_mat
        imp.frequency = np.ones(imp_mat.shape[0]) / imp_mat.shape[0]
        imp.date = np.arange(1, len(imp.at_event) + 1)
        imp.event_id = np.arange(1, len(imp.at_event) + 1)

        csv_output = BASE_DIR / f"impact_yearsets/global/csv/{exposures_str}/{hazard}_{exposures_str}_impacts_yearsets_150arcsec_{warming}_global_caped.csv"
        npz_output = BASE_DIR / f"impact_yearsets/global/npz/{exposures_str}/{hazard}_{exposures_str}_impacts_yearsets_150arcsec_{warming}_global_caped.npz"
        
        imp.write_csv(csv_output)
        imp.write_sparse_csr(npz_output)

if __name__ == "__main__":

    for exposures, warming in itertools.product(EXPOSURE_TYPES, WARMING_LEVELS):
        cap_yearsets(exposures, warming)
