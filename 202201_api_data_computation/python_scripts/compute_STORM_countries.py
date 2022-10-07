import os
from pathlib import Path

import numpy as np

from climada.hazard import TropCyclone
from pycountry import countries
from config import OUT_DATA_DIR


def main(replace=True):
    for scenario in ['historical', 'RCP85']:
        path0 = os.path.join('/nfs/n2o/wcr/szelie/CLIMADA_api_data/STORM')
        path = os.path.join(path0, 'global', scenario)
        files = os.listdir(path)
        for file in files:
            tc = TropCyclone.from_hdf5(os.path.join(path, file))
            tc.centroids.set_region_id()
            path_country = os.path.join(path0, 'countries', scenario)
            isExist = os.path.exists(path_country)
            if not isExist:
                os.makedirs(path_country)
            for country in countries:
                file_country = file.replace('global', country.alpha_3)
                file_country = os.path.join(path_country, file_country)
                if Path(file_country).exists() and replace is False:
                    continue
                tc_country = tc.select(reg_id=int(country.numeric))
                if tc_country is None:
                    continue
                tc_country.write_hdf5(file_country)


if __name__ == "__main__":
    main(replace=False)
#    main(n_tracks=50)

