# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 09:05:37 2022

@author: timo_
"""
# %%
import numpy as np
import pandas as pd 
import xarray as xr
import sys, os
import datetime as dt
import geopandas as gpd
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
from scClim.constants import ID_COL_DICT,PRE_PROC_PARAM_DICT

data_dir = str(CONFIG.local_data.data_dir)
out_dir = str(CONFIG.local_data.out_dir)
read_extent = [5.8, 10.6, 45.7, 47.9] #all of switzerland
########################### Definitions ######################################

######################## Selection parameters ################################
version_id = 1
if not version_id ==7:
    print("Warning: Version id 7 is the standard pre-processing. Only use others for MFZ")
    raise ValueError("avoid accidental cacluation")
param_dict = PRE_PROC_PARAM_DICT[version_id].copy()
years = np.arange(2002,2021+1) 
if version_id==1:
    #Version_id=1 is only for cars where data is available 2017+
    years = np.arange(2017,2021+1) 
##############################################################################

    
#read exposure and hazard (POH)
pathsPOH = sc.E.get_hazard_files_TS('POH',years,version_id,data_dir)
haz_poh = sc.hazard_from_radar(pathsPOH,varname='POH',
                                extent = read_extent)
poh = xr.open_mfdataset(pathsPOH, concat_dim = 'time',combine='nested',
                          coords = 'minimal').BZC
pathsMESHS = sc.E.get_hazard_files_TS('MESHS',years,version_id,data_dir)
haz_meshs = sc.hazard_from_radar(pathsMESHS,varname='MESHS',
                                extent = read_extent)

#create / read PH (possible hail) netcdfs
if os.path.exists(data_dir+'/V5/possible_hail/ph_v%d.nc'%version_id):
    ds_PH = xr.open_dataset(data_dir+'/V5/possible_hail/ph_v%d.nc'%version_id)
else:
    ds_PH= sc.get_possible_hail('no date', poh, haz_poh, extent=read_extent, 
                                poh_thresh=param_dict['poh_level'], 
                                buffer_km=param_dict['buffer_km'],
                                get_likelihood=param_dict['use_likelihood'],
                                return_type='all')
    ds_PH.to_netcdf(data_dir+'/V5/possible_hail/ph_v%d.nc'%version_id,
                    encoding={"possible_hail":{'zlib':True,'complevel':9}})
    
haz_PH = sc.hazard_from_radar(ds_PH,varname='possible_hail', extent=read_extent)
np.testing.assert_array_equal(haz_poh.centroids.lat,haz_PH.centroids.lat)
np.testing.assert_array_equal(haz_poh.date,haz_PH.date)

# %% GVZ 
#use GVZ_damages .csv file as template for the processed file
gvz_template = pd.read_csv(data_dir+'/GVZ/GVZ_Hail_Loss_200001_to_202203_nonZero.csv',
                           sep=';',nrows=100)#read first 100 rows to infer data types
gvz_template['date_dt'] = pd.to_datetime('1800-01-01')
gvz_template=gvz_template.iloc[:0,:]

gvz_gdf = sc.read_gvz_dmg(data_dir+'/GVZ/GVZ_Hail_Loss_200001_to_202203_nonZero.csv',
                           'GVZ', w_rel_dmg=True,return_type='gdf', index_dmgs =False)

#correct date using the prescribed POH-based filter
corrected_dmg = sc.correct_impact_dates(gvz_gdf,haz_PH,haz_poh,haz_meshs,poh,
                                        gvz_template,param_dict,
                                        log_file=data_dir +'/GVZ/readMe.txt')  

#save new damage file to csv
corrected_dmg.to_csv(data_dir+'/GVZ/GVZ_Hail_Loss_date_corrected%d.csv'%version_id,
                     sep=';',index=False)    

# %% GVL 
years = (2002,2021)
#use GVL_damages .csv file as template for the processed file
gvl_template = pd.read_csv(data_dir+'/GVL/GVL_damages_geocoded.csv',
                           sep=';',nrows=100)
gvl_template['date_dt'] = pd.to_datetime('1800-01-01')
gvl_template=gvl_template.iloc[:0,:]

gvl_gdf = sc.read_gvz_dmg(data_dir+'/GVL/GVL_damages_geocoded.csv',
                           cant_insurance='GVL', crs='EPSG:4326',
                           w_rel_dmg=True,return_type='gdf',years=years,
                           id_col=ID_COL_DICT['GVL'],index_dmgs=False)

#use only April to September
sel = (gvl_gdf.date_dt.dt.month>=4) & (gvl_gdf.date_dt.dt.month<=9)
print(f'{len(gvl_gdf)-sum(sel)} dmg entries outside of April-Sept are ignored')
gvl_gdf = gvl_gdf.loc[sel,:]

#correct date using the prescribed POH-based filter
corrected_dmg_GVL = sc.correct_impact_dates(gvl_gdf,haz_PH,haz_poh,haz_meshs,
                                        poh,gvl_template,param_dict,
                                        log_file=data_dir +'/GVL/_ReadMe.txt')  

#save new damage file to csv
corrected_dmg_GVL.to_csv(data_dir+'/GVL/GVL_Hail_Loss_date_corrected%d.csv'%version_id,
                         sep=';',index=False)    

# %% AGV 
years = (2002,2021)
#use AGV_damages .csv file as template for the processed file
agv_template = pd.read_csv(data_dir+'/AGV/AGV_damages_geocoded.csv',
                           sep=';',nrows=100)
agv_template['date_dt'] = pd.to_datetime('1800-01-01')
agv_template=agv_template.iloc[:0,:]

agv_gdf = sc.read_gvz_dmg(data_dir+'/AGV/AGV_damages_geocoded.csv','AGV', crs='EPSG:4326',
                           w_rel_dmg=True,return_type='gdf',years=years,
                           id_col=ID_COL_DICT['AGV'],index_dmgs =False)
#use only April to September
sel = (agv_gdf.date_dt.dt.month>=4) & (agv_gdf.date_dt.dt.month<=9)
print(f'{len(agv_gdf)-sum(sel)} dmg entries outside of April-Sept are ignored')
agv_gdf = agv_gdf.loc[sel,:]

#correct date using the prescribed POH-based filter
corrected_dmg_AGV = sc.correct_impact_dates(agv_gdf,haz_PH,haz_poh,haz_meshs,poh,
                                        agv_template,param_dict,
                                        log_file=data_dir +'/AGV/_ReadMe.txt')  

#save new damage file to csv
corrected_dmg_AGV.to_csv(data_dir+'/AGV/AGV_Hail_Loss_date_corrected%d.csv'%version_id,
                         sep=';',index=False) 
# %% GVB
years = (2002,2021)
#use GVB_damages .csv file as template for the processed file
gvb_template = pd.read_csv(data_dir+'/GVB/GVB_damages_geocoded.csv',
                           sep=';',nrows=100)
gvb_template['date_dt'] = pd.to_datetime('1800-01-01')
gvb_template=gvb_template.iloc[:0,:]

gvb_gdf = sc.read_gvz_dmg(data_dir+'/GVB/GVB_damages_geocoded.csv','GVB', crs='EPSG:4326',
                           w_rel_dmg=True,return_type='gdf',years=years,
                           id_col=ID_COL_DICT['GVB'],index_dmgs =False)
#use only April to September
sel = (gvb_gdf.date_dt.dt.month>=4) & (gvb_gdf.date_dt.dt.month<=9)
print(f'{len(gvb_gdf)-sum(sel)} dmg entries outside of April-Sept are ignored') #929
gvb_gdf = gvb_gdf.loc[sel,:]

#correct date using the prescribed POH-based filter
corrected_dmg_GVB = sc.correct_impact_dates(gvb_gdf,haz_PH,haz_poh,haz_meshs,poh,
                                        gvb_template,param_dict,
                                        log_file=data_dir +'/GVB/_ReadMe.txt')  

#save new damage file to csv
corrected_dmg_GVB.to_csv(data_dir+'/GVB/GVB_Hail_Loss_date_corrected%d.csv'%version_id,
                         sep=';',index=False)    


# %% Car damage data 

#Note: Insurance name is replaced by *** here
log_file_MFZ = data_dir +'/***/readMe.txt'
out_file_MFZ = data_dir+'/***/MFZ_Hail_Loss_date_corrected%d.csv'%version_id

MFZ_gdf = gpd.read_file(data_dir+'/***/MFZ_claims_random.gpkg') #random points

if 'CLAIMS_DATE' in MFZ_gdf.columns:
    MFZ_gdf['date_dt'] = pd.to_datetime(MFZ_gdf['CLAIMS_DATE'])
    MFZ_gdf['Schadendatum'] =  [dt.datetime.strftime(d, "%d%m%Y")  for d in MFZ_gdf['date_dt']]

#Use only claims until 2021 (for now)
years = np.array([ d.year  for d in MFZ_gdf['date_dt']])
MFZ_gdf = MFZ_gdf.loc[years<=2021,:]

MFZ_gdf=MFZ_gdf.to_crs(epsg=4326)
MFZ_gdf['latitude'] = MFZ_gdf.geometry.y
MFZ_gdf['longitude'] = MFZ_gdf.geometry.x
if not 'value' in MFZ_gdf.columns:
    MFZ_gdf = MFZ_gdf.rename(columns={'PAID':'value'})

corrected_dmg_MFZ = sc.correct_impact_dates(MFZ_gdf,haz_PH,haz_poh,haz_meshs,poh,
                                        MFZ_gdf.iloc[:0,],param_dict,log_file_MFZ)  

#rename columns 'PAID' to 'Schadensumme'
corrected_dmg_MFZ = corrected_dmg_MFZ.rename(columns={'value':'Schadensumme',
                                                      'VEHVALUE':'Versicherungssumme'})

#save new damage file to csv
corrected_dmg_MFZ.to_csv(out_file_MFZ,sep=';',index=False)  