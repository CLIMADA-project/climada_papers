#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of CLIMADA-papers.

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019. 


Plot global asset exposure data for 224 countries
Sections 3.1
Figure 3

Requires asset exposure data in the form of CSV as available from the ETH research repository:
https://doi.org/10.3929/ethz-b-000331316
Save CSV data in local folder ENTITY_DIR or set path in variable ENTITY_DIR before executing this script.

Specify result folder path in variable RES_DIR before running this script.

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.3.1
or later

The required gridded population data GPWv4.10 is available from SEDAC's Beta site, please see
https://beta.sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse

For more guidance on the LitPop module please refer to the CLIMADA tutorial:
https://climada-python.readthedocs.io/en/latest/tutorial/climada_entity_LitPop.html

@author: Samuel Eberenz
"""


import os


import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cartopy.crs as ccrs
from cartopy.io import shapereader

from climada.util.constants import SYSTEM_DIR
from climada.entity import Exposures
# from climada.entity.exposures.litpop import LitPop

version_date = 'v1_1'
REF_YEAR = 2014
RES_ARCSEC = 30

# List of target resolutions in arcsec
res_targets = [600] 
ENTITY_DIR = os.path.join(SYSTEM_DIR, 'litpop_%i' % (REF_YEAR))
ENTITY_DIR_HDF5 = os.path.join(SYSTEM_DIR, 'litpop_%i_hdf5' % (REF_YEAR))
filename_start = 'LitPop_pc_%iarcsec_' % (RES_ARCSEC)
RES_DIR = '' # SPECIFY RES_DIR BEFORE USE

# Import and export options:
write_to_tiff = True
try_read_from_tiff = True
save_plots = True
write_to_hdf5 = False
try_read_from_hdf5 = False

fix_zeros = True # CHANGES TOTAL VALUES OF EXPOSURE; USE FOR PLOTTING ONLY

if fix_zeros:
    fadd = '_fixzeros_'
else:
    fadd = ''


plot_minimum = 100

if not os.path.exists(RES_DIR):
    os.makedirs(RES_DIR)

files = [i for i in os.listdir(ENTITY_DIR) if os.path.isfile(os.path.join(ENTITY_DIR,i)) and \
         filename_start in i]
files = np.unique(files)
print('\n' + '\x1b[1;03;30;30m' + 'Number of country exposure files: %i' %(len(files)) + '\x1b[0m')

"""LOADING DATA FROM CSV AND REGRIDDING TO TARGET RESOLUTION PER COUNTRY:"""
print('\n' + '\x1b[1;03;30;30m' + 'REGRIDDING TO TARGET RESOLUTION PER COUNTRY' + '\x1b[0m')
for res_target in res_targets:
    for idx, fi in enumerate(files):
        exposure_tmp = Exposures()
        if os.path.exists(os.path.join(RES_DIR, '%s_%ias.tiff' %(fi[0:-4]+fadd, res_target))):
            print('\n' + '\x1b[1;03;30;30m' + 'TIFF exists already, skipping: %s_%ias.tiff' %(fi[0:-4], res_target) + '\x1b[0m')

            continue
        else:
            print('\n' + '\x1b[1;03;30;30m' + 'Loading: %s ...' %(fi) + '\x1b[0m')
            exposure_tmp = exposure_tmp.from_csv(os.path.join(ENTITY_DIR, fi), index_col=None)

            if np.isnan(exposure_tmp.value.max()):
                continue
            exposure_tmp = Exposures(exposure_tmp)
            exposure_tmp.set_geometry_points() # set geometry attribute (shapely Points) from GeoDataFrame from latitude and longitude
            exposure_tmp.check() # puts metadata that has not been assigned        
            if write_to_hdf5:
                exposure_tmp.write_hdf5(os.path.join(ENTITY_DIR_HDF5, '%s.hdf5' %(fi[0:-4]+fadd)))
            if write_to_tiff:
                if fix_zeros:
                    exposure_tmp.value[exposure_tmp.value<1] = 1
                # exposure_tmp.plot_raster(res=RES_ARCSEC/3600, save_tiff=\
                #                          os.path.join(RES_DIR, '%s_%ias.tiff' %(fi[0:-4], RES_ARCSEC)))
                exposure_tmp.plot_raster(res=RES_ARCSEC/3600, raster_res=res_target/3600, save_tiff=\
                                         os.path.join(RES_DIR, '%s_%ias.tiff' %(fi[0:-4]+fadd, res_target)))

"""COMBINE AND PLOT EXPOSURE AT TARGET RESOLUTION:"""
print('\n' + '\x1b[1;03;30;30m' + 'COMBINE AND PLOT EXPOSURE AT TARGET RESOLUTION' + '\x1b[0m')
for res_target in res_targets:
    exposure_data = Exposures()
    for idx, fi in enumerate(files):
        exposure_tmp = Exposures()
        if try_read_from_tiff and os.path.exists(os.path.join(RES_DIR, '%s_%ias.tiff' %(fi[0:-4]+fadd, res_target))):
            print('\n' + '\x1b[1;03;30;30m' + 'Loading: %s_%ias.tiff' %(fi[0:-4]+fadd, res_target) + '\x1b[0m')
            exposure_tmp.set_from_raster(os.path.join(RES_DIR, '%s_%ias.tiff' %(fi[0:-4]+fadd, res_target)))
            exposure_tmp = Exposures(exposure_tmp)
            exposure_tmp.set_geometry_points() # set geometry attribute (shapely Points) from GeoDataFrame from latitude and longitude
            exposure_tmp.check() # puts metadata that has not been assigned
            exposure_data = exposure_data.append(exposure_tmp)
            if fix_zeros:
                exposure_tmp.value[exposure_tmp.value<plot_minimum] = plot_minimum
        else:
            print('\n' + '\x1b[1;03;30;30m' + 'ERROR Loading: %s_%ias.tiff' %(fi[0:-4]+fadd, res_target) + '\x1b[0m')
        
    print('\n' + '\x1b[1;03;30;30m' + 'Checking combined data...' + '\x1b[0m')
    exposure_data.check()
    if write_to_tiff:
        print('\n' + '\x1b[1;03;30;30m' + 'Writing combined data to TIFF...' + '\x1b[0m')
        exposure_data.plot_raster(res=res_target/3600, raster_res=res_target/3600, save_tiff=\
            os.path.join(RES_DIR, '%s_000_%ias.tiff' %(files[0][0:-8]+fadd, res_target)))
    ax_exp = exposure_data.plot_hexbin(pop_name=False, cmap='plasma', norm=LogNorm(vmin=plot_minimum, vmax=0.1*exposure_data.value.max()))
    if save_plots:
        plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_%s_world_map.png' % (res_target, REF_YEAR, fadd)), \
                      dpi=300, facecolor='w', edgecolor='w', \
                      orientation='portrait', papertype=None, format='png', \
                      transparent=False, bbox_inches=None, pad_inches=0.1, \
                      frameon=None, metadata=None)
        plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_%s_world_map.pdf' % (res_target, REF_YEAR, fadd)), \
                      dpi=300, facecolor='w', edgecolor='w', \
                      orientation='portrait', papertype=None, format='pdf', \
                      transparent=False, bbox_inches=None, pad_inches=0.1, \
                          frameon=None, metadata=None)
    
    
    ax_exp_scatter = exposure_data.plot_scatter(pop_name=False, cmap='plasma', s=.01, shapes=False, \
                                                norm=LogNorm(vmin=plot_minimum, vmax=np.max([0.1*exposure_data.value.max(), 10**9])))

    shp_file = shapereader.natural_earth(resolution='10m', \
                    category='cultural', name='admin_0_countries')
    shp = shapereader.Reader(shp_file)
    for geometry in shp.geometries():
        ax_exp_scatter.add_geometries([geometry], crs=ccrs.PlateCarree(), facecolor='', \
                                edgecolor='black', linewidth=.5)
    
    if save_plots:
        plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map_scatter001pt_lw3%s.png' % (res_target, REF_YEAR, fadd)), \
                      dpi=300, facecolor='w', edgecolor='w', \
                      orientation='portrait', papertype=None, format='png', \
                      transparent=False, bbox_inches=None, pad_inches=0.1, \
                      frameon=None, metadata=None)
        plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map_scatter001pt_lw3%s.pdf' % (res_target, REF_YEAR, fadd)), \
                      dpi=300, facecolor='w', edgecolor='w', \
                      orientation='portrait', papertype=None, format='pdf', \
                      transparent=False, bbox_inches=None, pad_inches=0.1, \
                          frameon=None, metadata=None)