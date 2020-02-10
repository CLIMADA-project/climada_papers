"""
This file is part of CLIMADA-papers.

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019. 

Compute Lit^1Pop^1, Lit^1, and Pop^1 for four metropolitan areas and plot maps

Section 3.3
Figures 4, A1

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0
or later

The required gridded population data GPWv4.10 is available from SEDAC's Beta site, please see
https://beta.sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse

For more guidance on the LitPop module please refer to the CLIMADA tutorial:
https://climada-python.readthedocs.io/en/latest/tutorial/climada_entity_LitPop.html

@authors: Thomas Roosli, Samuel Eberenz
"""

# Import required packages:
import os
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from climada.entity import Exposures
from climada.entity.exposures.litpop import LitPop
from climada.util.constants import DATA_DIR
from climada.util.config import CONFIG
from climada.util.constants import SYSTEM_DIR


# set output directory
save_figures = True
output_path = os.path.join(DATA_DIR, 'results')
if not os.path.isdir(output_path):
    os.mkdir(output_path)
CONFIG['local_data']['save_dir'] = output_path

DPI = 600 # Graph export resolution


### WORLD MAP:

version_date = 'v1_1'
REF_YEAR = 2014
RES_ARCSEC = 30
ENTITY_DIR = os.path.join(SYSTEM_DIR, 'litpop_%i' % (REF_YEAR))
filename_start = 'LitPop_pc_%iarcsec_' % (RES_ARCSEC)
RES_DIR = output_path # '/Users/eberenzs/Documents/Projects/climada_papers/201903_litpop_exposure_data_model/climada_v1.3.1/results'

write_to_hdf5 = True
try_read_from_hdf5 = True

if not os.path.exists(RES_DIR):
    os.makedirs(RES_DIR)

files = [i for i in os.listdir(ENTITY_DIR) if os.path.isfile(os.path.join(ENTITY_DIR,i)) and \
         filename_start in i]
files = np.unique(files)
print('Number of country exposure files: %i' %(len(files)))

# 
exposure_data = LitPop()

if try_read_from_hdf5 and os.path.exists(os.path.join(ENTITY_DIR, 'LitPop_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC))):
    exposure_data.read_hdf5(os.path.join(ENTITY_DIR, 'LitPop_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC)))
else:
    grid_stats = pd.DataFrame(index=np.arange(0, len(files)), columns=['country', 'grid_count', 'sum', 'max', 'mean', 'median'])
    for idx, fi in enumerate(files):
        print('Loading: %s ...' %(fi))
        exposure_tmp = LitPop()
        exposure_tmp = exposure_tmp.from_csv(os.path.join(ENTITY_DIR, fi), index_col=None)
        exposure_data = exposure_data.append(exposure_tmp)
        # print('Max. grid cell value: USD %1.0f' %(exposure_tmp.value.max()))
        grid_stats.loc[idx, 'country'] = fi[-7:-4]
        grid_stats.loc[idx, 'grid_count'] = exposure_tmp.value[exposure_tmp.value > 0].count()
        grid_stats.loc[idx, 'sum'] = exposure_tmp.value.sum()
        grid_stats.loc[idx, 'max'] = exposure_tmp.value.max()
        grid_stats.loc[idx, 'mean'] = exposure_tmp.value.mean()
        grid_stats.loc[idx, 'median'] = exposure_tmp.value.median()

    del exposure_tmp
    grid_stats.to_csv(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_grid_stats.csv' % (RES_ARCSEC, REF_YEAR)))
    grid_stats.to_csv(os.path.join(ENTITY_DIR, 'LitPop_pc_%iarcsec_%i_grid_stats.csv' % (RES_ARCSEC, REF_YEAR)))
    exposure_data = Exposures(exposure_data)
    print('\n' + '\x1b[1;03;30;30m' + 'exposure_data is now an Exposures:', str(type(exposure_data)) + '\x1b[0m')
    exposure_data.set_geometry_points() # set geometry attribute (shapely Points) from GeoDataFrame from latitude and longitude
    print('\n' + '\x1b[1;03;30;30m' + 'check method logs:' + '\x1b[0m')
    exposure_data.check() # puts metadata that has not been assigned
    
    print('\n' + '\x1b[1;03;30;30m'  + 'exposure_data looks like:' + '\x1b[0m')
    print(exposure_data.head())
    print('\n' + '\x1b[1;03;30;30m'  + 'plotting global map...' + '\x1b[0m')
    print('Global max. grid cell value: USD %1.0f' %(exposure_data.value.max()))
    print('Global mean grid cell value: USD %1.0f' %(exposure_data.value.mean()))
    print('Global median grid cell value: USD %1.0f' %(exposure_data.value.median()))
    if write_to_hdf5:
        exposure_data.write_hdf5(os.path.join(ENTITY_DIR, 'LitPop_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC)))

ax_exp = exposure_data.plot_hexbin(pop_name=False, cmap='plasma', norm=LogNorm(vmin=10, vmax=0.1*exposure_data.value.max()))
if save_figures:
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map.pdf' % (RES_ARCSEC, REF_YEAR)), \
                  dpi=DPI, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='pdf', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map.png' % (RES_ARCSEC, REF_YEAR)), \
                  dpi=DPI, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='png', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)




### DETAIL MAPS:

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
country = ['GBR','ZAF','IND','MEX','USA']
city_name =['London', 'Cape Town', 'Mumbai', 'Mexico City', 'New York']
extent = [(-0.6,0.4,51,52),
          (17.7,19.1,-34.65,-33.2), (72,73.35,18.8,19.4),
          (-99.8,-98.6,18.9,20), (-74.6, -73, 40, 41)]
markersize = [6, 1, 4, 1, 1]
col_max_value = [3e9, 5e8, 6.5e8, 1.1e9, 1.2e10]
col_min_value = np.true_divide(col_max_value,1000)
plot_method = 1

for country_i,city_name_i,extent_i,markersize_i,col_max_value_i,col_min_value_i \
            in zip(country,city_name,extent,markersize,col_max_value,col_min_value):
    print('++++++++++++ %s ++++++++++++' %(city_name_i))               
    for exponents_i in exponents:
        if exponents_i == [1,0]:
            description_str = r'$Lit^1$'
        elif exponents_i == [0,1]:
            description_str = r'$Pop^1$'
        elif exponents_i == [1,1]:
            description_str = r'$Lit^1Pop^1$'
            
        exposure = get_entity(country_i,exponents_i,description_str,date_str)

        index_plot = (exposure.latitude > extent_i[2]) & (exposure.latitude < extent_i[3]) & (exposure.longitude > extent_i[0]) & (exposure.longitude < extent_i[1])

        
        plt.figure()
        if plot_method == 0:
            scatter_params = {
                    'cmap':'plasma',
                    'marker':'s',
                    's':markersize_i,
                    'norm':colors.LogNorm(1000,col_max_value_i)
                    }
            value_log = exposure.value[index_plot].replace(0,0.0001)
            plt.scatter(exposure.longitude[index_plot], exposure.latitude[index_plot], c=value_log , **scatter_params)
    #        plt.scatter(exposure.longitude[index_plot], exposure.latitude[index_plot], c=exposure.value[index_plot], **scatter_params)
            cbar = plt.colorbar()
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
            cbar = plt.colorbar()
            cbar.ax.set_ylabel('Exposure value [USD]')
        elif plot_method == 2:
            exposure_crop = exposure.cx[extent_i[2]:extent_i[3],extent_i[0]:extent_i[1]]
        plt.title(city_name_i + ' (' + country_i  + ') - ' + description_str)
        plot_name = os.path.abspath(os.path.join( \
                CONFIG['local_data']['save_dir'], country_i + '_' + city_name_i.replace(' ', '') + '_' + description_str +'.pdf'))
        if save_figures:
            plt.savefig(plot_name, bbox_inches='tight',dpi=DPI)

