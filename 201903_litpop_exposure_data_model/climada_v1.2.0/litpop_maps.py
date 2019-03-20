"""
This file is part of CLIMADA-papers.

Global LitPop: An Exposure Data Model for Disaster Risk Assessment based on Nightlight and Population Data

Section 3.1
Figures 2, A1
Compute LitPop, Lit, and Pop for metropolitan areas and plot maps

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0
or later

@author: ThomasRoosli
"""

# Import required packages:
import os
import numpy as np
# import pandas as pd
from matplotlib import colors
import matplotlib.pyplot as plt
# from iso3166 import countries as iso_cntry
import cartopy.crs as ccrs
#import pickle

from climada.entity.exposures.litpop import LitPop
from climada.util.constants import DATA_DIR
from climada.util.config import CONFIG

# set output directory
output_path = os.path.join(DATA_DIR, 'results')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
CONFIG['local_data']['save_dir'] = output_path


def get_entity(country,exponents,description_str,date_str):
    
    file_name = (country + '_litpop_' + description_str 
                 + '_' + date_str + '.hdf5')
    exposure_file = os.path.abspath(os.path.join( \
                CONFIG['local_data']['save_dir'], file_name))
    if os.path.isfile(exposure_file):
        exposure = LitPop()
        exposure.read_hdf5(exposure_file)
#        with open(exposure_file, 'rb') as f:
#            exposure = pickle.load(f)
    else:
        exposure = LitPop()
        # The arguments 'exponent' is used to set the power with which Lit and Pop go into LitPop:
        exposure.set_country(country, exponents=exponents, reference_year=2014)
        exposure.set_geometry_points()
        exposure.check()
        exposure.tag.description = (exposure.tag.description + '_' 
                                    + description_str) 
        exposure.set_lat_lon()
        exposure.write_hdf5(exposure_file)
        #pickle.dump(exposure,open(exposure_file, 'wb'))
    return exposure


date_str = '20190308'
exponents = [[1,0], [0,1], [1,1]]
country = ['USA','GBR','ZAF','IND','MEX']
city_name =['New York','London','Cape Town','Mumbai','Mexico City']
extent = [(-74.6, -73, 40, 41),(-0.6,0.4,51,52),
          (17.7,19.1,-34.65,-33.2),(72,73.35,18.8,19.4),
          (-99.8,-98.6,18.9,20)]
markersize = [1,6,1,4,1]
col_max_value = [1.2*10**10,3*10**9,5*10**8,6.5*10**8,1.1*10**9]
col_min_value = np.true_divide(col_max_value,10**3)
plot_method = 1

for country_i,city_name_i,extent_i,markersize_i,col_max_value_i,col_min_value_i \
            in zip(country,city_name,extent,markersize,col_max_value,col_min_value):
                   
    for exponents_i in exponents:
        if exponents_i == [1,0]:
            description_str = 'nightlight'
        elif exponents_i == [0,1]:
            description_str = 'population'
        elif exponents_i == [1,1]:
            description_str = 'litpop'
            
        exposure = get_entity(country_i,exponents_i,description_str,date_str)

        index_plot = (exposure.latitude > extent_i[2]) & (exposure.latitude < extent_i[3]) & (exposure.longitude > extent_i[0]) & (exposure.longitude < extent_i[1])

        
        plt.figure()
        if plot_method == 0:
            scatter_params = {
                    'cmap':'plasma',
                    'marker':'s',
                    's':markersize_i,
                    'norm':colors.LogNorm(10**3,col_max_value_i)
                    }
            value_log = exposure.value[index_plot].replace(0,0.0001)
            plt.scatter(exposure.longitude[index_plot], exposure.latitude[index_plot], c=value_log , **scatter_params)
    #        plt.scatter(exposure.longitude[index_plot], exposure.latitude[index_plot], c=exposure.value[index_plot], **scatter_params)
            plt.colorbar()
        elif plot_method == 1:
            lon_vec = np.sort(exposure.longitude[index_plot].unique())
            lat_vec = np.sort(exposure.latitude[index_plot].unique())
            value_mesh = np.empty([lat_vec.size,lon_vec.size])
            value_mesh.fill(np.nan)
            for value_i,lon_i,lat_i in zip(exposure.value[index_plot],exposure.longitude[index_plot], exposure.latitude[index_plot]):
                value_mesh[lat_vec==lat_i,lon_vec==lon_i]=value_i
            ax = plt.axes(projection=ccrs.PlateCarree())
            plt.pcolormesh(lon_vec,lat_vec,value_mesh,
                           transform=ccrs.PlateCarree(),
                           cmap='plasma',vmin=col_min_value_i,vmax=col_max_value_i,
                           norm=colors.LogNorm(col_min_value_i,col_max_value_i))
            plt.colorbar()
        elif plot_method == 2:
            exposure_crop = exposure.cx[extent_i[2]:extent_i[3],extent_i[0]:extent_i[1]]
        plt.title(country_i + ', ' + city_name_i + ' - ' + description_str + ' distributed')
        plot_name = os.path.abspath(os.path.join( \
                CONFIG['local_data']['save_dir'], country_i + '_' + city_name_i.replace(' ', '') + '_' + description_str +'.pdf'))
        plt.savefig(plot_name, bbox_inches='tight',dpi=300)
