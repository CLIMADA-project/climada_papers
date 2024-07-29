import os
from pathlib import Path

import numpy as np

from climada.hazard import TropCyclone
from pycountry import countries
from config import DATA_DIR
from create_log_file import log_msg

PATH_STORM_global = os.path.join(DATA_DIR, 'STORM_v2')
PATH_STORM_countries = os.path.join(DATA_DIR, 'STORM_v3')

LOG_FILE = "progress_make_STORM_countries.txt"


def main(replace=True):
    """
            Process STORM TC hazard data from a global scale down to individual countries.

            The function reads the global hazard data and saves the hazard data for each individual country.
        """
    for scenario in ['historical', 'RCP85']:
        log_msg(f"Started computing floods for scenario {scenario}\n", LOG_FILE)
        path = os.path.join(PATH_STORM_global, 'global', scenario)
        files = os.listdir(path)
        for file in files: #going through the different GCMs
            tc = TropCyclone.from_hdf5(os.path.join(path, file))
            tc.centroids.set_region_id()
            path_country = os.path.join(PATH_STORM_countries, 'countries', scenario)
            isExist = os.path.exists(path_country)
            if not isExist:
                os.makedirs(path_country)
            for country in countries:
                file_country = file.replace('global', country.alpha_3)
                file_country = os.path.join(path_country, file_country)
                if Path(file_country).exists() and replace is False:
                    log_msg(f"{file_country} already exist, and replace is False.\n", LOG_FILE)
                    continue
                tc_country = tc.select(reg_id=int(country.numeric))
                if tc_country is None:
                    continue
                tc_country.write_hdf5(file_country)
            log_msg(f"Countries where computed for the file {file}\n", LOG_FILE)


if __name__ == "__main__":
    main(replace=False)
#    main(n_tracks=50)

