"""
This file is part of CLIMADA.

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Script to create population exposure and match with habitat protective capacity data
Author: Sarah Hülsen
"""

import geopandas as gpd
from climada.entity import Exposures
import pandas as pd

# This script creates an exposure file from vectorised WorldPop population data (per TC basin)
# It also merges habitat data with the population data (4 permutations of two habitat datasets and two pop datasets)
# The merged files can later be used to analyse impacts per habitat category

# variables
geo_epsg = 4326  # WGS84
proj_eck4 = '+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'  # World Eckert IV
prot_dist = 2000
pop_years = ['2000', '2020']
hab_years = ['1992', '2020']
basins = ['AP', 'WP', 'IO', 'SH']
input = '/Users/sarah/Documents/CLIMADA/coastal_ecosystems/data/' 
output = '/Users/sarah/Documents/CLIMADA/coastal_ecosystems/results/intermediate/'

# create exposure files for both years per basin
for basin in basins:
    for pop_year in pop_years:
        # Load data
        pop_fn = f'wp_{pop_year}_{basin}.shp'
        pop = gpd.read_file(f'{input}{pop_fn})

        # prepare dataframe
        pop['longitude'] = pop['geometry'].x
        pop['latitude'] = pop['geometry'].y

        # write exposures
        exp_pop = Exposures(pop)

        # save exposure data
        exp_pop.write_hdf5(f'{output}exp_wp_{basin}_{pop_year}.hdf5')

# create merged habitat and population files per basin
for basin in basins:
    for pop_year in pop_years:
        for hab_year in hab_years:
            # Load data
            rhab_fn = f'rhab_{hab_year}_{basin}.shp'    # ecosystem rank data filename
            pop_fn = f'wp_{pop_year}_{basin}.shp'   # population data filename
            pop = gpd.read_file(f'{input}{pop_fn})
            rhab = gpd.read_file(f'{input}{rhab_fn})

            # prepare dataframes
            rhab['rhab'] = rhab['Rhab_all'].round().astype(int)   # discrete habitat categories
            rhab['rhab_cont'] = rhab['Rhab_all']   # continuous habitat values
            pop['longitude'] = pop['geometry'].x
            pop['latitude'] = pop['geometry'].y

            # buffer rhab
            rhab_buffer = rhab.to_crs(proj_eck4)
            rhab_buffer['buffer_geometry'] = rhab_buffer.geometry.buffer(prot_dist)
            rhab_buffer = rhab_buffer.rename(columns={'geometry': 'point', 'buffer_geometry': 'geometry'})
            rhab_buffer = rhab_buffer.to_crs(epsg=geo_epsg)

            col_select = ['geometry', 'rhab', 'rhab_cont']
            rhab_sel = rhab_buffer[col_select]
            exp_join = gpd.sjoin(pop, rhab_sel, how="inner", op="within")
            exp_join['exp_index'] = exp_join.index
            join_filter = exp_join.sort_values(by='rhab', ascending=True)
            join_filter = join_filter.drop_duplicates(subset='exp_index', keep='first')    # drop duplicate matches between population and habitat

            # save data
            join_filter.to_csv(f'{output}exp_rhab_{hab_year}_wp_{pop_year}_{basin}.csv')

# create global merged habitat and population files
for pop_year in pop_years:
    for hab_year in hab_years:
        # load files
        AP = pd.read_csv(f'{output}exp_rhab_{hab_year}_wp_{pop_year}_AP.csv')
        WP = pd.read_csv(f'{output}exp_rhab_{hab_year}_wp_{pop_year}_WP.csv')
        IO = pd.read_csv(f'{output}exp_rhab_{hab_year}_wp_{pop_year}_IO.csv')
        SH = pd.read_csv(f'{output}exp_rhab_{hab_year}_wp_{pop_year}_SH.csv')

        # concatenate basins
        basins_all = pd.concat([AP, WP, IO, SH])
        basins_all.to_csv(f'{output}wp{pop_year}_rhab{hab_year}.csv')
