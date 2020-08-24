#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the CLIMADA papers repository:
    https://github.com/CLIMADA-project/climada_papers

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Created on Mon Aug 24 09:56:35 2020

Description:
    This script is to compute the ABSOLUTE number of displaced people due to river floods at 5km resolution globally.
    The output raster files represent the number of displaced people at every 10 years from 1970-2000 using the baseline population at year 2000 and historical simulation of flooding hazards.
    The flow of the script is similar to fl_cc_displacement_raster.py
    Please refer to the paper (submitted) and README for more information.
    For the baseline calculation please refer to flood_historical_displacement_raster.py.
    JUPYTER NOTEBOOK documented the calculation using Vietnam as a case study.

@author: manniekam
"""
import logging
import numpy as np
from rasterio.warp import Resampling
from rasterio import Affine

from climada.hazard import Hazard
from climada.entity import Exposures, ImpactFuncSet, ImpactFunc
from climada.entity.exposures.base import INDICATOR_CENTR, INDICATOR_IF
from climada.engine import Impact
from climada.util.coordinates import write_raster
from climada.util.files_handler import get_file_names

logging.basicConfig(format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                    level=logging.INFO)
HAZ_TYPE = 'FL' # hazard type abbreviation (CLIMADA convention), default: 'FL'

HAZARD_PTH = '/<SET PATH>' # The directory which contains hazard files

HAZ_FILES = get_file_names(HAZARD_PTH)

EXP_POP_PTH = '/<SET PATH>' # The directory which contains exposure files

SAVE_PTH =   '/<SET PATH>' # Destinated path for the output raster files

RCP = 'historical' # historical simulation

HYDRO_MODEL = ['H08', 'LPJmL', 'MPI-HM', 'ORCHIDEE', 'PCR-GLOBWB', 'WaterGAP2']     # 6 * hydrological model

GCM_MODEL = ['gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc5']      # 4 * climate model

YEAR = np.arange(1970,2010,10) # Years for the computation

DST_META = {'width': 8640, 'height': 3347, 'crs': {'init':'epsg:4326'}, \
            'transform': Affine(0.04166666666666666, 0.0, -180.0,
                                0.0, -0.041666666666666664, 83.64166666610012)}

class IF_FL_1m(ImpactFunc):
    """ Define impact function which is 0 for flood of less than 1m height and
    1 otherwise. """

    def __init__(self):
        ImpactFunc.__init__(self)
        self.id = 1
        self.name = 'flood_1m'
        self.intensity_unit = 'm'
        self.haz_type = HAZ_TYPE

        self.intensity = np.arange(0, 1000, 0.5)
        self.paa = np.ones(self.intensity.size)
        self.mdd = np.zeros(self.intensity.size)
        for idx in range(self.mdd.size):
            if self.intensity[idx] >= 1:
                self.mdd[idx] = 1.0

    def calc_mdr(self, inten):
        if not isinstance(inten, np.ndarray):
            inten = np.array([inten])
        mdr = np.zeros(inten.shape)
        mdr[inten>=1] = 1.0
        return mdr

def iffl():
    """ Define impact functions """
    if_1m = IF_FL_1m()
    if_fl = ImpactFuncSet()
    if_fl.tag.description = '1m step function'
    if_fl.append(if_1m)

    return if_fl

ifs_step = iffl()

yy_start = 106
yy_end = 116

idx_band = 1

exp = Exposures()

ssp_file = EXP_POP_PTH +'baseYr_total_2000.tif'
exp.set_from_raster(ssp_file, transform=DST_META['transform'], height=DST_META['height'], 
                width=DST_META['width'], resampling=Resampling.average)

exp.value *= 25     # sum of the grids after upscaling
if np.any(exp.value<0) == True:
    raise ValueError
exp.value_unit = 'N people per pixel'
exp.ref_year = 2000
exp[INDICATOR_CENTR+HAZ_TYPE] = np.arange(len(exp), dtype=int)
exp[INDICATOR_IF+HAZ_TYPE] = np.ones(len(exp), dtype=int)
exp.check()

for year in YEAR:    
    
    imp_model = []  # initite a list fro all models in one year
    imp_save = None

    for hydro, gcm in [(hydro, gcm) for hydro in HYDRO_MODEL for gcm in GCM_MODEL]:
        
        haz_file = HAZARD_PTH +'flddph_' +hydro +'_' +gcm +'_' +RCP +'_flopros_gev_picontrol_2006_2300_0.1.nc'
        haz_frac = HAZARD_PTH +'fldfrc_' +hydro +'_' +gcm +'_' +RCP +'_flopros_gev_picontrol_2006_2300_0.1.nc'
    
        if haz_file not in HAZ_FILES:
            logging.error('no file: flddph_' +hydro +'_' +gcm +'_' +RCP +'_flopros_gev_picontrol_2006_2300_0.1.nc')
            continue
        
        haz = Hazard('FL')
        haz.set_raster([haz_file], [haz_frac], band=range(yy_start,yy_end), 
                          transform=DST_META['transform'], height=DST_META['height'],
                          width=DST_META['width'], resampling=Resampling.bilinear,
                          attrs={'frequency':np.ones(10)/10}) 
        
        imp_tmp = Impact()
        imp_tmp.calc(exp, ifs_step, haz, save_mat=True)
        imp_model.append(imp_tmp)
        
        if imp_save is None:
            imp_save = np.reshape(imp_tmp.eai_exp, [1,-1])
        else:
            imp_save_tmp = np.reshape(imp_tmp.eai_exp, [1,-1])
            imp_save = np.append(imp_save, imp_save_tmp, axis=0)

        
        logging.info('%s %s %s %s %s %s', year, hydro, gcm, imp_tmp.eai_exp.sum(), imp_tmp.aai_agg)
        
    save_file = SAVE_PTH +'imp_' +RCP +'_BaseYearPOP_' +str(year) +'.tif'
    write_raster(save_file, imp_save, DST_META)
    
    yy_start += 10
    yy_end += 10