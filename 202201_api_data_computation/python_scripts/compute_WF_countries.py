import os
from pathlib import Path

import numpy as np

from climada.hazard import TropCyclone, Hazard, Centroids
from pycountry import countries
from config import OUT_DATA_DIR


def main(replace=True):
    for scenario in ['historical']:
        path0 = os.path.join('/nfs/n2o/wcr/szelie/CLIMADA_api_data/wildfire')
        path = os.path.join(path0, 'global', scenario)
        files = os.listdir(path)
        for file in files:
            wf = Hazard.from_hdf5(os.path.join(path, file))
            wf.centroids = Centroids.from_hdf5("/nfs/n2o/wcr/szelie/CLIMADA_api_data/"
                                               "centroids/earth_centroids_150asland_1800asoceans_distcoast_region.hdf5")
            path_country = os.path.join(path0, 'countries', scenario)
            isExist = os.path.exists(path_country)
            if not isExist:
                os.makedirs(path_country)
            for country in countries:
                file_country = file.replace('global',country.alpha_3)
                file_country = os.path.join(path_country, file_country)
                if Path(file_country).exists() and replace is False:
                    continue
                wf_country = wf.select(reg_id=int(country.numeric))
                if wf_country is None:
                    print("country not found:"+ country.numeric)
                    continue
                wf_country.write_hdf5(file_country)


if __name__ == "__main__":
    main()
#    main(n_tracks=50)

