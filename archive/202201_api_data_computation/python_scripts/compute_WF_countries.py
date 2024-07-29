import os
from pathlib import Path

from climada.hazard import Hazard, Centroids
from pycountry import countries
from config import DATA_DIR
from create_log_file import log_msg
from climada.util.api_client import Client

def main(replace=True):
    """
        Process wildfire hazard data from a global scale down to individual countries.

        The function reads the global wildfire hazard data and saves the hazard data for each individual country.
    """
    client = Client()
    for scenario in ['historical']:
        path0 = os.path.join('/nfs/n2o/wcr/szelie/CLIMADA_api_data/wildfire')
        path = os.path.join(path0, 'global', scenario)
        files = os.listdir(path)
        for file in files:
            wf = Hazard.from_hdf5(os.path.join(path, file))
            wf.centroids = client.get_centroids(res_arcsec_land=150,
            res_arcsec_ocean=1800,
            extent=(-180, 180, -90, 90),)
            path_country = os.path.join(path0, 'countries', scenario)
            isExist = os.path.exists(path_country)
            if not isExist:
                os.makedirs(path_country)
            for country in countries:
                file_country = file.replace('global', country.alpha_3)
                file_country = os.path.join(path_country, file_country)
                if Path(file_country).exists() and replace is False:
                    continue
                wf_country = wf.select(reg_id=int(country.numeric))
                if wf_country is None:
                    print("country not found:" + country.numeric)
                    continue
                wf_country.write_hdf5(file_country)


if __name__ == "__main__":
    main()
#    main(n_tracks=50)

