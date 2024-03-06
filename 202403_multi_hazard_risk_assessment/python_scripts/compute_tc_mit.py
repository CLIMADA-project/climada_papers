####################### import basic modules #################################
import os
import copy
import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd
import scipy.sparse as sp
from datetime import datetime as dt
from collections import Counter
import datetime
from climada.hazard import Hazard
import climada.util.coordinates as u_coord
from climada.hazard import TropCyclone
from config import BASE_DIR

MIT_DIR = BASE_DIR / "MIT_windfields/windfields"

MIT_STR = "%s_rcp%02d_%s_100.csv"

GCM = ['GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5']  # MIT tracks

MAX_VAL = 1e11  # max. value for cotbour plot
MIN_VAL = 2e2  # min. value
BUFFER_DEG_LON = .5  # map buffer
BUFFER_DEG_LAT = .5

RCP = [26, 60]
RCP_STR = ['rcp26', 'rcp60']

BASIN_BOUNDS = {
    # Eastern Pacific Basin
    'EP': [-180.0, -75.0, 0.0, 60.0],

    # North Atlantic Basin
    'NA': [-105.0, -30.0, 0.0, 60.0],

    # Northern Indian Basin
    'NI': [37.0, 99.0, 0.0, 35.0],

    # Southern Indian Basin
    'SIW': [20.0, 75.0, -50.0, 0.0],
    'SIE': [75.0, 135.0, -50.0, 0.0],

    # Southern Pacific Basin
    'SPW': [135.0, 172.0, -50.0, 0.0],
    'SPE': [172.0, -60.0, -50.0, 0.0],

    # Western Pacific Basin
    'WPN': [99.0, 180.0, 20.0, 60.0],
    'WPS': [99.0, 180.0, 0.0, 20.0],
}
years_warming_level = {'1': {
    'rcp26': {'gfdl-esm2m': [2006, 2023], 'miroc5': [2006, 2025], 'hadgem2-es': [2006, 2019], 'ipsl-cm5a-lr': None},
    'rcp60': {'gfdl-esm2m': [2006, 2027], 'miroc5': [2013, 2034], 'hadgem2-es': [2006, 2023], 'ipsl-cm5a-lr': None}},

    '2': {'rcp26': {'gfdl-esm2m': None, 'miroc5': None, 'hadgem2-es': None, 'ipsl-cm5a-lr': [2019, 2040]},
          'rcp60': {'gfdl-esm2m': [2066, 2087], 'miroc5': [2061, 2082], 'hadgem2-es': [2040, 2061], 'ipsl-cm5a-lr': [2019, 2040]}}}

def haz_select_bounds(haz, bounds):
    """
    Selects hazard events within specified geographical bounds.

    Parameters:
    - haz (Hazard): A CLIMADA Hazard object.
    - bounds (tuple): A tuple of (lonmin, lonmax, latmin, latmax) defining the geographical bounds.

    Returns:
    - haz (Hazard): A modified CLIMADA Hazard object with events outside the specified bounds removed.
    """
    
    lonmin, lonmax, latmin, latmax = bounds

    lat, lon = haz.centroids.lat, haz.centroids.lon
    if lonmin < lonmax:
        haz.centroids.region_id = (latmin <= lat) & (lat <= latmax) \
                                  & (lonmin <= lon) & (lon <= lonmax)
    else:
        # basin crossing the -180/180 longitude threshold
        haz.centroids.region_id = (latmin <= lat) & (lat <= latmax) \
                                  & ((lonmin <= lon) | (lon <= lonmax))
    haz.select(reg_id=1)
    return haz


def load_tc_realization(path, year_start, year_end, basin=None, draws_subset=None):
    """
    Loads tropical cyclone hazard realizations from CSV files and filters them by year and basin.

    Parameters:
    - path (str): Path to the CSV file containing tropical cyclone data.
    - year_start (int): Start year for filtering.
    - year_end (int): End year for filtering.
    - basin (str): Specific basin to filter the cyclones (default is None).
    - draws_subset (array-like): Specific subset of realization IDs to include (default is None).

    Returns:
    - hazard (Hazard): A CLIMADA Hazard object containing filtered tropical cyclone events.
    """
    
    draws = pd.read_csv(path)
    draws = draws[(year_start <= draws['year']) & (draws['year'] <= year_end)]
    if draws_subset is not None:
        draws = draws[draws['real_id'].isin(draws_subset)].reset_index(drop=True)

    # load master hazard object containing all relevant Emanuel tracks
    filenames = np.unique(draws['name'].str.replace("-[0-9]+$", "", regex=True))
    hazards = []
    for fname in filenames:
        path = os.path.join(MIT_DIR, fname)
        haz = TropCyclone.from_mat(path)
        haz.event_name = [f"{fname}-{n}" for n in haz.event_name]
        if basin is not None:
            haz = haz_select_bounds(haz, BASIN_BOUNDS[basin])
            haz.basin = [basin] * len(haz.event_name)
            mask = np.abs(haz.intensity.data) < 17
            haz.intensity.data[mask] = 0
            haz.intensity.eliminate_zeros()
        hazards.append(haz)
    master_hazard = TropCyclone.concat(hazards)
    event_names = list(draws['name'].values)
    hazard = master_hazard.select(event_names=event_names)
    hazard.event_name = [f"{n}-{i}" for n, i in zip(hazard.event_name, draws['real_id'].values)]
    years = [dt.fromordinal(d).year for d in hazard.date]
    date_diff = [dt(y_dst, 1, 1).toordinal() - dt(y, 1, 1).toordinal() \
                 for y_dst, y in zip(draws['year'].values, years)]
    hazard.date += date_diff
    return hazard


def make_tc_hazard(n_draws=25):
    """
    Generates and stores tropical cyclone hazards for specified RCP scenarios and warming levels.
    
    This function iterates over different RCP scenarios, GCM models, and warming levels to generate
    a set of tropical cyclone hazards. The hazards are then concatenated and stored in a specified
    HDF5 file.
    """
    basins = list(BASIN_BOUNDS.keys())

  # Renamed from 'rcps' to 'RCP_STR'
    for warming_level in years_warming_level:
        haz_list = []
        for r, rcp in enumerate(RCP):
            for g, gcm in enumerate(GCM):
                gcm = gcm.lower()
                if years_warming_level[warming_level][RCP_STR[r]][gcm] is None:
                    continue
                year_start, year_end = years_warming_level[warming_level][RCP_STR[r]][gcm]
                draws_subset = np.random.choice(np.arange(1, 100), n_draws)
                for basin in basins:
                    path = os.path.join(MIT_DIR, 'draws', MIT_STR % (GCM[g], RCP[r], basin))
                    haz = load_tc_realization(path, basin=basin, year_start=year_start, year_end=year_end, draws_subset=draws_subset)
                    haz_list.append(haz)
                    draw = [e.split('-')[-1] for e in haz.event_name]
                    event_name = [f"{str(datetime.date.fromordinal(date).year)}_{gcm.lower()}_{RCP_STR[r]}" for date in haz.date]
                    haz.event_name = ["_".join([name, draw[n], basin, str(n)]) for n, name in enumerate(event_name)]
        hazard = Hazard.concat(haz_list)
        hazard.event_name = ["_".join(event_name.split("_")[0:4]) for event_name in hazard.event_name]
        hazard.centroids.lon = u_coord.lon_normalize(hazard.centroids.lon)
        hazard_path = BASE_DIR / f"MIT_windfields/hazard_TC_MIT/hazard_tc_{warming_level}.hdf5"
        hazard.write_hdf5(hazard_path)

if __name__ == "__main__":
    make_tc_hazard()