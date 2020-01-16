#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:36:02 2019

@author: eberenzs
"""
import os


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from climada.util.constants import SYSTEM_DIR
from climada.entity import Exposures
from climada.entity.exposures.litpop import LitPop


version_date = 'v1_1'
REF_YEAR = 2014
RES_ARCSEC = 30
ENTITY_DIR = os.path.join(SYSTEM_DIR, 'litpop_%i' % (REF_YEAR))
filename_start = 'LitPop_pc_%iarcsec_' % (RES_ARCSEC)
RES_DIR = '/Users/eberenzs/Documents/Projects/climada_papers/201903_litpop_exposure_data_model/climada_v1.3.1/results'

write_to_hdf5 = False
try_read_from_hdf5 = True

if not os.path.exists(RES_DIR):
    os.makedirs(RES_DIR)

files = [i for i in os.listdir(ENTITY_DIR) if os.path.isfile(os.path.join(ENTITY_DIR,i)) and \
         filename_start in i]
files = np.unique(files)
print('Number of country exposure files: %i' %(len(files)))

# 
exposure_data = LitPop()

if try_read_from_hdf5 and os.path.exists(os.path.join(ENTITY_DIR, 'LitPop_hdf5_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC))):
    exposure_data.read_hdf5(os.path.join(ENTITY_DIR, 'LitPop_hdf5_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC)))
else:
    grid_stats = pd.DataFrame(index=np.arange(0, len(files)), columns=['country', 'grid_count', 'sum', 'max', 'mean', 'median'])
    for idx, fi in enumerate(files[0:3]):
        print('Loading: %s ...' %(fi))
        exposure_tmp = LitPop()
        exposure_tmp = exposure_tmp.from_csv(os.path.join(ENTITY_DIR, fi), index_col=None)
        exposure_data = exposure_data.append(exposure_tmp)

        grid_stats.loc[idx, 'country'] = fi[-7:-4]
        grid_stats.loc[idx, 'grid_count'] = exposure_tmp.value[exposure_tmp.value > 0].count()
        grid_stats.loc[idx, 'sum'] = exposure_tmp.value.sum()
        grid_stats.loc[idx, 'max'] = exposure_tmp.value.max()
        grid_stats.loc[idx, 'mean'] = exposure_tmp.value.mean()
        grid_stats.loc[idx, 'median'] = exposure_tmp.value.median()

    del exposure_tmp
    grid_stats.to_csv(os.path.join(RES_DIR, 'LitPop_meta_%iarcsec_%i_grid_stats.csv' % (RES_ARCSEC, REF_YEAR)))
    grid_stats.to_csv(os.path.join(ENTITY_DIR, 'LitPop_meta_%iarcsec_%i_grid_stats.csv' % (RES_ARCSEC, REF_YEAR)))
    exposure_data = Exposures(exposure_data)
    print('\n' + '\x1b[1;03;30;30m' + 'exposure_data is now an Exposures:', str(type(exposure_data)) + '\x1b[0m')
    exposure_data.set_geometry_points() # set geometry attribute (shapely Points) from GeoDataFrame from latitude and longitude
    print('\n' + '\x1b[1;03;30;30m' + 'check method logs:' + '\x1b[0m')
    exposure_data.check() # puts metadata that has not been assigned
    
    print('\n' + '\x1b[1;03;30;30m'  + 'exposure_data looks like:' + '\x1b[0m')
    print(exposure_data.head())

    print('Global number of grid cells with value > USD 0: %i' %(exposure_data.value[exposure_data.value > 0].count()))
    print('Global max. grid cell value: USD %1.0f' %(exposure_data.value.max()))
    print('Global mean grid cell value: USD %1.0f' %(exposure_data.value.mean()))
    print('Global median grid cell value: USD %1.0f' %(exposure_data.value.median()))

    if write_to_hdf5:
        print('\n' + '\x1b[1;03;30;30m'  + 'Write to hdf5...' + '\x1b[0m')
        exposure_data.write_hdf5(os.path.join(ENTITY_DIR, 'LitPop_hdf5_pc_%iarcsec_000_all.hdf5' % (RES_ARCSEC)))

print('\n' + '\x1b[1;03;30;30m'  + 'regridding raster & plotting global maps...' + '\x1b[0m')

plot_res = [1800, 300]
for res_target in plot_res:
    print('\n' + '\x1b[1;03;30;30m'  + 'raserize to %ias...' %(res_target) + '\x1b[0m')
    ax_exp_raster1 = exposure_data.plot_raster(res=RES_ARCSEC/3600, raster_res=res_target/3600, save_tiff=\
                os.path.join(RES_DIR, 'LitPop_pc_1800as_%i_world_map.tiff' % (REF_YEAR)))
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map_raster_%ias.png' % (RES_ARCSEC, REF_YEAR, res_target)), \
                  dpi=600, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='png', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map_raster_%ias.pdf' % (RES_ARCSEC, REF_YEAR, res_target)), \
                  dpi=600, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='pdf', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)
    
    print('\n' + '\x1b[1;03;30;30m'  + 'plot hexbin %ias...' %(res_target) + '\x1b[0m')
    exp_plot = Exposures()    
    exp_plot.set_from_raster(os.path.join(RES_DIR, 'LitPop_pc_300as_%i_world_map.tiff' % (REF_YEAR)))
    
    ax_exp = exp_plot.plot_hexbin(pop_name=False, cmap='plasma', norm=LogNorm(vmin=10, vmax=0.1*exp_plot.value.max()))
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map.png' % (res_target, REF_YEAR)), \
                  dpi=600, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='png', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)
    plt.savefig(os.path.join(RES_DIR, 'LitPop_pc_%iarcsec_%i_world_map.pdf' % (res_target, REF_YEAR)), \
                  dpi=600, facecolor='w', edgecolor='w', \
                  orientation='portrait', papertype=None, format='pdf', \
                  transparent=False, bbox_inches=None, pad_inches=0.1, \
                  frameon=None, metadata=None)
