import os
from pathlib import Path
from climada_petals.hazard.river_flood import RiverFlood
from pycountry import countries
import sys
from config import OUT_DATA_DIR


def main(years=['2010', '2030'], scenario='rcp26', replace=False):
    years_str = "_".join([str(years[0]), str(years[1])])
    path0 = os.path.join(OUT_DATA_DIR, 'flood')
    path = os.path.join(path0, 'global', scenario, years_str)
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
            file_country = "".join(
                [path_country, '/', f[0], '_', f[1], '_', f[2], '_', f[3], '_', country.alpha_3, '_',
                 f[4]])
            if Path(file_country).exists() and replace is False:
                continue
            rf2 = rf.select(reg_id=int(country.numeric))
            if rf2 is None:
                continue
            rf2.write_hdf5(file_country)


if __name__ == "__main__":
    main(years=[sys.argv[1], sys.argv[2]], scenario=sys.argv[3])