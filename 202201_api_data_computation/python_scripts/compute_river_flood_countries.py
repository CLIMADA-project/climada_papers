import os
from pathlib import Path
from climada_petals.hazard.river_flood import RiverFlood
from pycountry import countries
import sys
from config import DATA_DIR
from create_log_file import log_msg


def main(years=None, scenario='rcp26', replace=True):
    """
     Process river flood hazard data from a global scale down to individual countries.

     The function reads the global river flood hazard data for a specified scenario and years,
     then filters and saves the hazard data for each individual country.

     Parameters:
         years (list of str, optional): List of start and end year. Default is ['2010', '2030'].
         scenario (str, optional): Scenario to consider (e.g., 'rcp26'). Default is 'rcp26'.
         replace (bool, optional): If True, existing country-specific river flood hazard
                                   files will be overwritten. Default is True.

     Returns:
         None
     """

    LOG_FILE = "progress_make_river_flood_countries.txt"

    if years is None:
        years = ['2010', '2030']
    years_str = "_".join([str(years[0]), str(years[1])])
    path0 = os.path.join(DATA_DIR, 'flood_v2')
    path = os.path.join(path0, 'global', scenario, years_str)
    log_msg(f"Reading in global flood file for scenario {scenario}  and years {years}\n", LOG_FILE)

    for file in os.listdir(path):
        file_path = os.path.join(path, file)
        path_country = os.path.join(path0, 'country', scenario, years_str)
        isExist = os.path.exists(path_country)
        if not isExist:
            os.makedirs(path_country)
        f = file.split('_', 4)
        rf = RiverFlood()
        rf.read_hdf5(file_path)
        for country in countries:
            log_msg(f"Begin computing river floods for {country}. \n", LOG_FILE)

            file_country = "".join(
                [path_country, '/', f[0], '_', f[1], '_', f[2], '_', f[3], '_', country.alpha_3, '_',
                 f[4]])
            if Path(file_country).exists() and replace is False:
                continue
            rf2 = rf.select(reg_id=int(country.numeric))
            if rf2 is None:
                continue
            rf2.write_hdf5(file_country)
            log_msg(f"Country {country} saved. \n", LOG_FILE)


if __name__ == "__main__":
    main(years=[sys.argv[1], sys.argv[2]], scenario=sys.argv[3])