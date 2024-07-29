# -*- coding: utf-8 -*-
"""
Main script for hail damage calculations with climada
"""
# %%
import numpy as np #1.20.3 is default
import pandas as pd 
import xarray as xr
import sys
import os
import copy
import matplotlib.pyplot as plt
import datetime as dt
import warnings
import cartopy.crs as ccrs
from climada.entity import Exposures
from climada.engine import Impact, ImpactCalc
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
from scClim.constants import ZRH_EXTENT,CANTON_DICT,PLOT_EXTENT_DICT,ID_COL_DICT,DATA_RANGE_DICT

data_dir = str(CONFIG.local_data.data_dir)
out_dir = str(CONFIG.local_data.out_dir)
mode = str(CONFIG.mode)  # plot / save_hdf5

event_def_version = 7
haz_var = 'MESHS'
impf = 'emp_emanuel'
exposure = 'gridKGV_'

# select year range
years = np.intersect1d(DATA_RANGE_DICT[exposure],DATA_RANGE_DICT[haz_var])
n_years = len(years)
exp_str = exposure

# Hazard
paths = sc.E.get_hazard_files_TS(haz_var, years, event_def_version, data_dir)
haz = sc.hazard_from_radar(
    paths, extent=[5.8, 10.6, 45.7, 48], varname=haz_var)

# Exposure 
if exposure == '' :
    exp = sc.read_gvz_exposure(data_dir+'/GVZ/GVZ_Exposure_202201.csv')
elif exposure == 'scaledGVZ_':
    exp = sc.read_gvz_exposure(data_dir+'/GVZ/GVZ_Exposure_202201.csv',
                                exterior_only=True)
elif exposure == 'MFZrandom_':
    exp_dict = {}
    for year in years:
        exp = Exposures.from_hdf5(f'{data_dir}/***/hdf5/MFZ_exp_random_{year}.hdf5')
        exp.value_unit = 'CHF'
        exp_dict.update({year:exp})
elif exposure in ['GVL','AGV','GVB','KGV','scaledKGV_','KGV_nonExtreme','KGV_1e5']:
    if 'scaled' in exposure:
        exterior_only = True; exp_str = exposure[6:9]
    elif 'nonExtreme' in exposure or '1e5' in exposure:
        exp_str = exposure[0:3]; exterior_only=False
    else:
        exterior_only=False; exp_str = exposure
        
    exp = sc.read_gvz_exposure(f'{data_dir}/{exp_str}/Gebäudewert_geocoded_v3.csv',
                                crs = 'EPSG:4326',exterior_only = exterior_only)
    if '1e5' in exposure:
        exp.gdf.value = 1e5
elif exposure == 'gridKGV_' and event_def_version == 7:
    exp = sc.read_xr_exposure(data_dir+'/KGV/ds_building_dmg_v7_1000m.nc','value_exposure')
    imp_measured = sc.read_xr_impact(data_dir+'/KGV/ds_building_dmg_v7_1000m.nc',
                                      'value',years=(years[0],years[-1]))
    exp_str = 'KGV'
    
#plotting values
pl_extent = PLOT_EXTENT_DICT[exp_str]; pl_canton = CANTON_DICT[exp_str]
intensity_label = f'{haz_var} [{sc.constants.UNIT_DICT[haz_var]}]'

# %% Vulnerability (impact functions)

impf_exp_str = 'KGV' if exposure == 'gridKGV_' else exposure 
impf_path = (f"{data_dir}/out_files/paa_mdd_smooth_{impf_exp_str}{haz_var}"
             f"_v{event_def_version}.csv")
if impf == '':
    imp_fun_set = sc.impf_from_csv(impf_path,smooth=False,plot=True)
elif impf == 'smooth':
    imp_fun_set = sc.impf_from_csv(impf_path,smooth=True,plot=True)
elif impf == 'emp_emanuel':
    imp_fun_set = sc.impf_from_csv(impf_path, smooth=False,
                                    emanuel_fit=True, plot=True)
elif impf == 'emanuel':
    impf_path = (f"{data_dir}/out_files/calib_opt/impf_opt_R2_{impf_exp_str}"
                 f"{haz_var}_v{event_def_version}.csv")
    imp_fun_set = sc.impf_from_csv(impf_path, smooth=False,
                                    emanuel_fit=True, plot=True)

# %% Imact

# 1) observed: read in damages from file as Impact object
if exposure in ['','scaledGVZ_']:
    imp_measured = sc.read_gvz_dmg(data_dir+'/GVZ/GVZ_Hail_Loss_date_corrected%d.csv'%event_def_version, 'GVZ',
                                    return_type='imp',years=(years[0],years[-1])) 
    exp_id = 'VersicherungsID'
elif exposure == 'MFZrandom_':
    imp_path = (f"{data_dir}/***/***_MFZ_Hail_Loss_date_corrected"
                f"{event_def_version}.csv")
    imp_measured=sc.read_gvz_dmg(imp_path,exposure,years=(years[0],years[-1]),
                                  index_dmgs=False,crs='EPSG:4326',
                                  id_col=ID_COL_DICT[exposure])    
    exp_id = 'POLNR'
elif exposure in ['GVL','AGV','GVB','KGV','scaledKGV_','KGV_nonExtreme','KGV_1e5']:
    index_dmgs = False if 'KGV' in exposure else True #KGV exposure is already indexed
    add_string = '_nonExtreme_' if 'nonExtreme' in exposure else ''
    imp_path = (f"{data_dir}/{exp_str}/{exp_str}_Hail_Loss_date_corrected"
                f"{add_string}{event_def_version}.csv")
    imp_measured = sc.read_gvz_dmg(imp_path, exp_str, return_type='imp',
                                    years=(years[0],years[-1]),index_dmgs=index_dmgs,
                                    crs='EPSG:4326',id_col=ID_COL_DICT[exp_str])

#extend hazard in case it doesn't cover all dates with reported impact
if any([imp_date not in haz.date for imp_date in imp_measured.date]):
    warnings.warn('Some dates with measured impact, are not contained in the \
                  Hazard file. Intensity=0 is assumed for these events')
    #identify dates with damages, but no hazard info
    dates = imp_measured.date[np.array([imp_date not in haz.date for imp_date in imp_measured.date])]
    #set hazard intensity to zero at these dates
    haz2 = sc.calibration.extend_haz(haz,dates)
    haz=haz2
    
    
    
    
# 2) calculated: 
if exposure == 'MFZrandom_':
    imp_dict={}
    for year in years:
        haz_now = haz.select(date=[f'{year}-01-01',f'{year}-12-31'])
        imp_now = ImpactCalc(exp_dict[year], imp_fun_set, haz_now).impact(save_mat=True)
        imp_dict.update({year:imp_now})
else:
    imp = ImpactCalc(exp, imp_fun_set, haz).impact(save_mat=True)
    if mode == 'plot': imp.plot_hexbin_eai_exposure(ignore_zero=False, gridsize=40)


if mode == 'save_hdf5':
    imp.write_sparse_csr(str(CONFIG.local_data.save_dir)+'/imp_test.npz')

# %% scatter plots of true and modelled damages
dur_path = data_dir + '/V0/hailgrids_duration/XXX/XXX_DURd66_V2_YEAR.nc'
imp_var = 'PAA' #'PAA' or 'MDR'
dmg_thresh = 1e5 #threshold for damages to calculate skill scores

if imp_var == 'PAA':
    dmg_thresh = 1e2

    #impact function
    if impf == 'smooth':
        imp_fun_set_PAA = sc.impf_from_csv(impf_path,PAA_only=True,smooth=True,plot=True)
    elif impf == 'emp_emanuel':
        imp_fun_set_PAA = sc.impf_from_csv(impf_path,PAA_only=True,smooth=False,
                                            emanuel_fit=True,plot=True)
    
    #Damage data
    if exposure in ['','scaledGVZ_']:
        imp_measured2 = sc.read_gvz_dmg(data_dir+'/GVZ/GVZ_Hail_Loss_date_corrected%d.csv'%event_def_version, 'GVZ',
                                        return_type='imp',years=(years[0],years[-1]),exp_n_assets=True) 
    elif exposure in ['GVL','AGV','GVB','KGV','scaledKGV_']:
        imp_measured2 = sc.read_gvz_dmg(imp_path,exp_str, return_type='imp',
                                         years=(years[0],years[-1]),
                                         exp_n_assets=True,index_dmgs=True,
                                         crs='EPSG:4326',id_col=ID_COL_DICT[exp_str])  
    elif exposure == 'gridKGV_':
        exp_PAA = sc.read_xr_exposure(data_dir+'/KGV/ds_building_dmg_v7_1000m.nc','n_count_exposure',val_unit='')
        imp_measured2 = sc.read_xr_impact(data_dir+'/KGV/ds_building_dmg_v7_1000m.nc',
                                           'n_count',unit='',years=(years[0],years[-1]))
        imp_PAA = ImpactCalc(exp_PAA, imp_fun_set_PAA, haz).impact(save_mat=True)
    elif exposure == 'MFZrandom_':
        imp_measured2=sc.read_gvz_dmg(f'{data_dir}/***/***_MFZ_Hail_Loss_date_corrected{event_def_version}.csv',
                                       exposure,years=(years[0],years[-1]),index_dmgs=False,
                                       crs='EPSG:4326',id_col=ID_COL_DICT[exposure],exp_n_assets=True)    
    
        
    #Exposure & impact calculation 
    if exposure in ['GVZ','GVL','AGV','GVB','KGV','scaledKGV_']:
        exp_PAA=copy.deepcopy(exp)
        exp_PAA.gdf.value = 1 #set values of all exposure points to 1 (only for non-gridded expousure)
        exp_PAA.value_unit = ''
        
        imp_PAA = ImpactCalc(exp_PAA, imp_fun_set_PAA, haz).impact(save_mat=True)
        
    elif exposure == 'MFZrandom_':
        exp_dict_PAA = copy.deepcopy(exp_dict)
        imp_PAA_dict={}
        for year in years:
            exp_dict_PAA[year].gdf.value = 1
            exp_dict_PAA[year].value_unit=''
            #impact calculation
            haz_now = haz.select(date=[f'{year}-01-01',f'{year}-12-31'])
            imp_PAA_now = ImpactCalc(exp_dict_PAA[year], imp_fun_set_PAA, haz_now).impact(save_mat=True)
            imp_PAA_dict.update({year:imp_PAA_now})        
    
    
#choose correct impact objects to compare for MDR or PAA comparisons
if imp_var=='PAA':
    if exposure == 'MFZrandom_':
        imp_now = Impact()
        imp_now.date = np.concatenate([imp_PAA_dict[key].date for key in imp_PAA_dict.keys()])
        imp_now.at_event = np.concatenate([imp_PAA_dict[key].at_event for key in imp_PAA_dict.keys()])
        imp_now.unit = ''
    else:
        imp_now = imp_PAA
    xmin = 0.5 #threshold for per-event damages to be completely disregarded
    imp_obs_now = imp_measured2
elif imp_var == 'MDR':
    if exposure == 'MFZrandom_':
        imp_now = Impact()
        imp_now.date = np.concatenate([imp_dict[key].date for key in imp_dict.keys()])
        imp_now.at_event = np.concatenate([imp_dict[key].at_event for key in imp_dict.keys()])
        imp_now.unit = 'CHF'
    else:
        imp_now=imp
    imp_obs_now = imp_measured
    xmin=1e2
    

#create impact dataframe with modelled and reported impacts per event (day)
#get all dates with modelled OR reported damages above xmin
ord_dates_nonZero = np.sort(np.unique(np.concatenate((imp_now.date[imp_now.at_event>xmin],
                                                      imp_obs_now.date[imp_obs_now.at_event>xmin]))))
imp_df = pd.DataFrame(index=ord_dates_nonZero,data={
    'date':[dt.datetime.fromordinal(d) for d in ord_dates_nonZero],
    "prePost2012":[dt.datetime.fromordinal(d).year>2012 for d in ord_dates_nonZero]})

imp_dfMod = pd.DataFrame(data={"imp_modelled":imp_now.at_event},index = imp_now.date)
imp_dfObs = pd.DataFrame(data={"imp_obs":imp_obs_now.at_event},index = imp_obs_now.date)

imp_df= imp_df.join(imp_dfMod,how='left') #join on index
imp_df= imp_df.join(imp_dfObs,how='left') #join on index
imp_df = imp_df.fillna(0) #fill NaN with zeros (for days with no reported or modelled damages)


#scale modelled impact to reflect that before some buildings were not built yet in previous years (exposure is from 2021)
#Note: explicitly changing the exposure each year only has a negligible difference to this (faster) approach
if any([name in exposure for name in ['GVL','AGV','GVB','GVZ','KGV']]):
    scale_factor_year = pd.read_csv(data_dir+'/out_files/constants/KGV_scale_factor.csv',index_col=0)
    scale_factor_str = 'scale_factor' if imp_var=='MDR' else 'scale_factorPAA'
    imp_df['imp_modelled_raw'] = imp_df['imp_modelled']
    imp_df['imp_modelled'] = scale_factor_year.loc[imp_df.date.dt.year,scale_factor_str].values*imp_df['imp_modelled'].values
elif exposure != 'MFZrandom_':
    raise Warning("scaling for changing exposure not implemented \
                  and not explicitly considered (as with MFZrandom)")
    
#calculate skill measures and save as csv 
expHazImp= f'{exposure}_{haz_var}_{impf}'
if not event_def_version == 7: expHazImp = f'{expHazImp}_v{event_def_version}'
rmse,rmsf,rmsf_weighted,FAR,POD,p_within_OOM,n_ev = sc.E.calc_and_save_skill_scores(
    imp_df,dmg_thresh,imp_var,expHazImp)

#create dictionary with evaluation metrics to pass to plotting function
eval_dict = {var_name: globals()[var_name] for var_name in ["rmse","rmsf","rmsf_weighted",
            "FAR","POD","p_within_OOM","n_ev","haz_var","exposure","impf"]}

#plot scatter plot
fig = sc.plot_funcs.scatter_from_imp_df(imp_df,imp_now.unit,xmin,dmg_thresh,eval_dict)
fig.savefig(out_dir+'/skill_plots/%s/%s_scatter_v%d_%s_%s.png'%(exposure,haz_var, event_def_version,impf,imp_var),
            dpi=250,bbox_inches = 'tight')
    
# %% plot Example for paper
sel_dates = [dt.datetime.strptime(str_date, "%Y-%m-%d") for str_date in ["2020-06-26","2021-06-21","2021-06-28"]]
fig,ax = sc.plot_event(sel_dates[0], imp_now, haz, imp_obs_now, exp,dur_path, haz_var = intensity_label,#draw_grid=True,
                      extent=pl_extent,pl_type='sel_dates',canton=pl_canton,vmax=0,pop_name=False,
                      vmax_haz=80,cmap_log=False,sel_dates=sel_dates)#,,vmax=1
for ax in [fig.axes[i] for i in (0,3,6)]:
    ax.set(title=ax.get_title()[-10:]) #remove Event_id from title
fig.savefig(out_dir+'/paper1/fig10.png',dpi=300,bbox_inches='tight')