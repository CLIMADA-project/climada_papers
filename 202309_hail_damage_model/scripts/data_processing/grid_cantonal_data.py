# -*- coding: utf-8 -*-
"""
Created on Tue May 24 22:03:54 2022

Script to create gridded damage data from cantonal data
"""

# %% Imports
import geopandas as gpd
import pandas as pd
import numpy as np
import sys
import shapely.geometry as sg
import xarray as xr
from scipy import optimize
import cartopy.crs as ccrs
import datetime as dt
from climada import util
from climada.entity import LitPop, Exposures
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim.hail_climada as fct
import scClim as sc
from scClim.constants import SUB_CH_EXTENT_2056,CH_EXTENT_EPSG2056,ID_COL_DICT
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import ticker, cm

data_dir = str(CONFIG.local_data.data_dir)
out_dir = str(CONFIG.local_data.out_dir)
poh_ex = xr.open_dataset(data_dir +"/V5/BZC/BZC_X1d66_2021.nc")
event_def_version = 7

# %% Merge damage and exposure data from all cantons

#Damage
years = (2002,2021)

#loop over damage data from 4 cantons (GVZ separately) and add them to a dictionary
dmg_dict = {}
for exposure in ['GVL','GVB','AGV']:
    f_path = (f'{data_dir}/{exposure}/{exposure}_Hail_Loss_date_corrected'
              f'{event_def_version}.csv')
    dmg_gdf = fct.read_gvz_dmg(f_path, exposure,years=years,return_type='gdf',
                               index_dmgs=True,crs='EPSG:4326',
                               id_col=ID_COL_DICT[exposure])
    dmg_gdf['id_col']=exposure+dmg_gdf[ID_COL_DICT[exposure]].astype(str)
    dmg_dict.update({exposure:dmg_gdf})
    
dmg_gvz = fct.read_gvz_dmg(data_dir+'/GVZ/GVZ_Hail_Loss_date_corrected%d.csv'%event_def_version, 'GVZ',
                                return_type='gdf',years=years,index_dmgs=True) 
dmg_gvz['id_col']='GVZ'+dmg_gvz[ID_COL_DICT['GVZ']].astype(str)
dmg_dict.update({'GVZ':dmg_gvz})

#merge all cantons
col_select=['date_dt','value','Versicherungssumme','geometry','id_col','Baujahr']
dmg_gdf = pd.concat([dmg_dict['GVL'][col_select],dmg_dict['AGV'][col_select],
           dmg_dict['GVB'][col_select],dmg_dict['GVZ'][col_select]],axis=0)

#save as dataframe in EPSG 4236
dmg_gdf_save = dmg_gdf.copy(deep=True)
dmg_gdf_save['longitude'] = dmg_gdf_save.geometry.x
dmg_gdf_save['latitude'] = dmg_gdf_save.geometry.y
dmg_gdf_save=dmg_gdf_save.rename(columns={'value':'Schadensumme'})
dmg_gdf_save.Versicherungssumme[dmg_gdf_save.Versicherungssumme==0]=np.nan
assert(dmg_gdf_save.crs=='EPSG:4326')
f_path = f'{data_dir}/KGV/KGV_Hail_Loss_date_corrected{event_def_version}.csv'
dmg_gdf_save.drop(columns=['geometry']).to_csv(f_path,sep=';')

#transform to epsg 2056 for gridded data
dmg_gdf = dmg_gdf.to_crs('epsg:2056')

# %% Exposure
exp_dict={}
for exposure in ['GVL','GVB','AGV']:
    exp = fct.read_gvz_exposure(f'{data_dir}/{exposure}/Gebäudewert_geocoded_v3.csv',
                                crs = 'EPSG:4326')
    exp.gdf['id_col']=exposure+exp.gdf[ID_COL_DICT[exposure]].astype(str)
    exp_dict.update({exposure:exp})
    if exposure == 'GVL':
        exp.gdf = exp.gdf.rename(columns = {'Volumen *':'Volumen'})
exp_zrh = fct.read_gvz_exposure(data_dir+'/GVZ/GVZ_Exposure_202201.csv')
exp_zrh.gdf['id_col']='GVZ'+exp_zrh.gdf[ID_COL_DICT['GVZ']].astype(str)
exp_dict.update({'GVZ':exp_zrh})

# merge all cantons
col_select=['value','geometry','Baujahr','id_col']
exp_gdf = pd.concat([exp_dict['GVL'].gdf[col_select],exp_dict['AGV'].gdf[col_select],
           exp_dict['GVB'].gdf[col_select],exp_dict['GVZ'].gdf[col_select]],axis=0)

# ensure no duplicates in ID col
if exp_gdf.id_col.duplicated().sum()<0.001*exp_gdf.shape[0]:
    exp_gdf = exp_gdf.loc[~exp_gdf.id_col.duplicated(),:] #remove duplicates (<0.1%)
else:
    raise UserWarning('more than 0.1% are duplicates. please check  data')
    
#save as dataframe in EPSG 4236
exp_gdf_save = exp_gdf.copy(deep=True)
exp_gdf_save['longitude'] = exp_gdf_save.geometry.x
exp_gdf_save['latitude'] = exp_gdf_save.geometry.y
assert(exp_gdf_save.crs=='EPSG:4326')
exp_gdf_save.drop(columns=['geometry']).to_csv(f'{data_dir}/KGV/Gebäudewert_geocoded_v3.csv',sep=';')
#transform to epsg 2056 for gridded data
exp_gdf = exp_gdf.to_crs('epsg:2056')

#save scale factors
year_range = np.arange(years[0],years[1]+1)
tot_val = exp_gdf.value.sum()
scale_factor = np.array([1-exp_gdf.loc[exp_gdf.Baujahr>year,:].value.sum()/tot_val 
                         for year in year_range])
tot_num_buildings = exp_gdf.value.count()
scale_factorPAA = np.array([1-exp_gdf.loc[exp_gdf.Baujahr>year,:].value.count()/
                            tot_num_buildings for year in year_range])

#save scale factor as dataframe
df_data = {'year':year_range, 'scale_factor':scale_factor,
           'scale_factorPAA':scale_factorPAA}
pd.DataFrame(data=df_data).to_csv(data_dir+'/out_files/constants/KGV_scale_factor.csv',
                                  index=False)

# %% grid damage data 
cell_size =1000 #in m
cell_CH = sc.regrid.create_cell(SUB_CH_EXTENT_2056,cell_size = cell_size,epsg=2056)

#Expsoure
cell_out = sc.regrid.group_to_cell(cell_CH, exp_gdf,
                                   aggfunc={'value':"sum",'Baujahr':'mean'},
                                   count_var='id_col',plot=True)
cell_multiIndex = cell_out.drop(columns='geometry').set_index(['KoordinateNord','KoordinateOst'])
ds_exp=cell_multiIndex.to_xarray().set_coords(('lat','lon'))

#damage
aggfunc={'value':"sum",'Versicherungssumme':'sum'}
for datum in dmg_gdf.date_dt.unique():
    print(datum)
    gdfNow = dmg_gdf.loc[dmg_gdf.date_dt == datum,:]
    cell_now = sc.regrid.group_to_cell(cell_CH,gdfNow,aggfunc=aggfunc,
                                       count_var='id_col',plot=False)
    cell_multiIndex = cell_now[['KoordinateOst','KoordinateNord','value','Versicherungssumme','n_count']].set_index(['KoordinateNord','KoordinateOst'])
    ds_now=cell_multiIndex.to_xarray().expand_dims(dim={'date': [datum]})
    ds_dmg = ds_now if datum==dmg_gdf.date_dt.unique()[0] else xr.concat([ds_dmg, ds_now], dim='date')

ds_dmg[['value','Versicherungssumme','n_count']] = ds_dmg[['value','Versicherungssumme','n_count']].fillna(0) #fill na
ds_dmg = ds_dmg.where(ds_exp['n_count']>0) # make sure it's NA outside of Exposure area

#calculate PAA, MDR by compareing exposure to damages (FIXED EXPOSURE FOR NOW)
ds_dmg['PAA']=ds_dmg['n_count']/ds_exp['n_count']
ds_dmg['MDR']=ds_dmg['value']/ds_exp['value']
for var in ['n_count','Baujahr','value']: #,'lat','lon'
    newVar = var if var in ['lat','lon'] else var+'_exposure'
    ds_dmg[newVar]=ds_exp[var]

ds_dmg = ds_dmg.rename({'KoordinateOst':'chx','KoordinateNord':'chy'})
#sort by dates
ds_dmg = ds_dmg.sortby('date')

#save as netcdf
ds_dmg.to_netcdf(f'{data_dir}/KGV/ds_building_dmg_v{event_def_version}_{cell_size}m.nc',encoding={
    "value":{'zlib':True,'complevel':5},
    "Versicherungssumme":{'zlib':True,'complevel':5},
    "n_count":{'zlib':True,'complevel':5},
    "PAA":{'zlib':True,'complevel':5},
    "MDR":{'zlib':True,'complevel':5},
    "Baujahr_exposure":{'zlib':True,'complevel':5},
    "value_exposure":{'zlib':True,'complevel':5},
    "n_count_exposure":{'zlib':True,'complevel':5},
})

sc.plot_nc(ds_dmg.sel(date='2021-06-28').value,crs=ccrs.epsg(2056),
           canton=['Bern','Zürich','Luzern','Aargau'],vmin=1,vmax=5e6)
extr_dates =ds_dmg.date.values[np.argsort(ds_dmg.value.sum(dim=['chx','chy'])).values[-10:]]
for date in extr_dates:
    sc.plot_nc(ds_dmg.sel(date=date).PAA,crs=ccrs.epsg(2056),canton=['Bern','Zürich','Luzern','Aargau'],
               title= f"{str(date)[:10]}, {ds_dmg.sel(date=date).value.sum(dim=['chx','chy']).values:.1e}CHF, ")
    plt.savefig(out_dir+f'/explorative/KGV/{str(date)[:10]}_PAA.png',bbox_inches='tight',dpi=220)


#confirming that the dmg_gdf and the ds_dmg are equal (dmg/counts per day)
assert(max(ds_dmg.value.sum(dim=['chx','chy'])-
           dmg_gdf.groupby('date_dt').value.sum())<1e-5)
assert(max(ds_dmg.n_count.sum(dim=['chx','chy'])-
           dmg_gdf.groupby('date_dt').value.count()<1e-5))