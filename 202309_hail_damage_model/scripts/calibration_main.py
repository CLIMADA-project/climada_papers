# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 10:37:47 2022

Main script for calibration
"""

# Load packages and set base_dir
import sys
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import datetime as dt

from climada.entity import ImpactFunc, ImpactFuncSet, Exposures
from climada.engine import ImpactCalc
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
from scClim.constants import CUT_OFF_DICT, INT_RANGE_DICT,ID_COL_DICT,W_SIZE_DICT,DATA_RANGE_DICT,DMG_BIN_DICT,UNIT_DICT,INT_LABEL_DICT

############################Definitions#######################################
event_def_version = 7 #7, (1,8,81,9 for MFZ)
haz_var = 'MESHS' # 'MESHS', 'dBZ' ,'E_kinCC'
exposure = 'KGV' # 'MFZrandom_','GVL','AGV','GVB','KGV','scaledKGV_',
get_per_gridcell = False  # if True: calculates additional per-grid_cell dataset
baujahr= '' #'', 'before1960', '1960-2002', 'after2002'
##############################################################################

# set directories
data_dir = str(CONFIG.local_data.data_dir)
out_dir = str(CONFIG.local_data.out_dir)

# load data into CLIMADA objects
years = np.intersect1d(DATA_RANGE_DICT[exposure],DATA_RANGE_DICT[haz_var])
n_years = len(years)

# configuration_path
c_path = '%s%s_v%d%s' % (exposure,haz_var, event_def_version,baujahr)
if baujahr=='':
    fig_path = '%s%s/v%d' %(exposure,haz_var, event_def_version)
else:
    fig_path = '%s%s/%s/v%d' %(exposure,haz_var, baujahr, event_def_version)
    
#defaults:
dmg_bin_size=DMG_BIN_DICT[haz_var]
#give warning if not whole year range is selected
if not (years[0]==2002 and years[-1]==2021):
    print(f'INFO: Only subrange selected: {years[0]}-{years[-1]}')

#Hazaard
paths = sc.E.get_hazard_files_TS(haz_var, years, event_def_version, data_dir)
    
#set cut_off and intensity_range
unit = UNIT_DICT[haz_var]
cut_off = CUT_OFF_DICT[haz_var]
intensity_range=INT_RANGE_DICT[haz_var]
intensity_label = INT_LABEL_DICT[haz_var]
#windowsize (rolling) in hazard units:
assert(len(np.unique(np.diff(intensity_range[1:])))==1)
stepsize = np.unique(np.diff(intensity_range[1:]))[0]
window_size = W_SIZE_DICT[haz_var]
unit_windowsize = stepsize*window_size
#xlim for plots
xlim = [min(intensity_range)-2, max(intensity_range)]
if haz_var=='HKE': xlim=[-2,160] #plot HKE only until 160j/m2 for better visibility

haz = sc.hazard_from_radar(
    paths, extent=[5.8, 10.6, 45.7, 47.9], varname=haz_var)  # country_code=756)

if exposure == '' :
    exp = sc.read_gvz_exposure(data_dir+'/GVZ/GVZ_Exposure_202201.csv')
    exp.assign_centroids(haz)
elif exposure == 'scaledGVZ_':
    exp = sc.read_gvz_exposure(data_dir+'/GVZ/GVZ_Exposure_202201.csv',
                               exterior_only=True)
    exp.assign_centroids(haz)
elif exposure in ['GVL','AGV','GVB','KGV','scaledKGV_','KGV_nonExtreme','KGV_1e5']:
    if 'scaled' in exposure:
        exterior_only = True; exp_str = exposure[6:9]
    elif 'nonExtreme' in exposure or '1e5' in exposure:
        exp_str = exposure[0:3]; exterior_only=False
    else:
        exterior_only=False; exp_str = exposure
        
    exp = sc.read_gvz_exposure(f'{data_dir}/{exp_str}/Gebäudewert_geocoded_v3.csv',
                               crs = 'EPSG:4326',exterior_only = exterior_only)
    exp.assign_centroids(haz)
    exp.gdf=exp.gdf.rename(columns={ID_COL_DICT[exp_str]:'VersicherungsID'})
    if '1e5' in exposure:
        exp.gdf.value = 1e5
elif exposure == 'MFZrandom_':
    exp_dict = {}
    for year in years:
        exp = Exposures.from_hdf5(f'{data_dir}/***/hdf5/MFZ_exp_random_{year}.hdf5')
        exp.assign_centroids(haz)
        exp_dict.update({year:exp})

if exposure in ['','scaledGVZ_']:
    imp_measured = sc.read_gvz_dmg(data_dir+'/GVZ/GVZ_Hail_Loss_date_corrected%d.csv'%event_def_version, 'GVZ',
                                    return_type='imp',years=(years[0],years[-1]),index_dmgs=True) 
    exp_id = 'VersicherungsID'
elif exposure == 'MFZrandom_':
    imp_path = (f"{data_dir}/***/***_MFZ_Hail_Loss_date_corrected"
                f"{event_def_version}.csv")
    imp_measured=sc.read_gvz_dmg(imp_path,exposure,years=(years[0],years[-1]),
                                 index_dmgs=False,crs='EPSG:4326',
                                 id_col=ID_COL_DICT[exposure])    
    
    #rename ID_col to match VersicherungsID
    exp.gdf=exp.gdf.rename(columns={ID_COL_DICT[exposure]:'VersicherungsID'})
    
elif exposure in ['GVL','AGV','GVB','KGV','scaledKGV_','KGV_nonExtreme','KGV_1e5']:
    index_dmgs = False if 'KGV' in exposure else True #KGV exposure is already indexed
    add_string = '_nonExtreme_' if 'nonExtreme' in exposure else ''
    imp_path = (f"{data_dir}/{exp_str}/{exp_str}_Hail_Loss_date_corrected"
                f"{add_string}{event_def_version}.csv")
    imp_measured = sc.read_gvz_dmg(imp_path, exp_str,return_type='imp',
                                   years=(years[0],years[-1]), 
                                   baujahr_filter=baujahr,index_dmgs=index_dmgs,
                                   crs='EPSG:4326',id_col=ID_COL_DICT[exp_str])     
   
    month= np.array([dt.datetime.fromordinal(d).month for d in imp_measured.date])
    event_ids= imp_measured.event_id[(month>=4) & (month<=9)]
    imp_measured = imp_measured.select(event_ids=event_ids)
    print(f'Warning: Dmgs from Oct-March are discarded: {sum((month<4) | (month>9))} dates')

# create idendity impact function (only for selcetion of dates with nonzero 
# hazard intensity over an exposed asset)
imp_fun_set = ImpactFuncSet()
imp_fun_identity = ImpactFunc.from_step_impf((0, 1, max(intensity_range)*2),haz_type='HL')
imp_fun_set.append(imp_fun_identity)
# calculate impact (with identitiy impact function)
imp = ImpactCalc(exp, imp_fun_set, haz).impact(save_mat=True)

dates_modeled_imp = np.array([dt.datetime.fromordinal(d) for d in imp.date[imp.at_event > 0]])

# %%

if any([imp_date not in haz.date for imp_date in imp_measured.date]):
    warnings.warn('Some dates with measured impact, are not contained in the \
                  Hazard file. Intensity=0 is assumed for these events')
    
    #identify dates with damages, but no hazard info
    dates = imp_measured.date[np.array([imp_date not in haz.date for imp_date 
                                        in imp_measured.date])]
    #set hazard intensity to zero at these dates
    haz2 = sc.calibration.extend_haz(haz,dates)
    haz=haz2
    
if 'Baujahr' in exp.gdf.columns:
    filter_year = sc.constants.BAUJAHR_DICT[baujahr] #(0,2022)
    if any(exp.gdf.Baujahr.isna()):
        if sum(exp.gdf.Baujahr.isna())/len(exp.gdf.Baujahr)>0.05:
            raise ValueError('Over 5% of exposure points without valid Baujahr')
        # Missing entries are mainly from Bern, when buildings were renovated!
        # i.e. they were already standing before 2002, when the radar data begins!
        exp.gdf.Baujahr[exp.gdf.Baujahr.isna()] = 2000
else:
    print('Warning: No "baujahr" in Exposure')
    filter_year = None

calib_tuple = sc.empirical_calibration_per_exposure(hazard_object = haz,
    exposure_object = exp, damages = imp_measured, exposure_type = 'GVZ', 
    variable = haz_var,filter_year=filter_year,dates_modeled_imp=dates_modeled_imp,
    roll_window_size=window_size,get_PVA=True)
ds, df_all, ds_roll, ds_roll_cut, values_at_centroid_all, intensity_range = calib_tuple

#make sure zero values are included (not NaN)
values_at_centroid_all['PAA'] = (values_at_centroid_all['n_dmgs']/
                                 values_at_centroid_all['n_exp'])
values_at_centroid_all['MDR'] = (values_at_centroid_all['dmg_val']/
                                 values_at_centroid_all['exp_val'])
#save values_at_centroids to csv
values_at_centroid_all.to_csv(f'{data_dir}/out_files/calib_emp/at_centr_{c_path}.csv',index=False)

n_samples=1000
out_tuple=sc.bootstrapping(ds, ds_roll,haz_var,n_samples,intensity_range,
                           log_fit=False,cut_off=cut_off,keep_raw_values=True)
ds_boot,ds_boot_roll,ds_boot_roll_cut, fit_data = out_tuple

# df, df_roll, df_roll_cut, df_roll_cut_fit, n_dmgs = sc.compute_empirical_damage_functions(ds,ds_roll,ds_roll_cut)
df, df_roll, df_roll_cut, n_dmgs = sc.compute_empirical_damage_functions(ds,ds_roll,
                                                                         ds_roll_cut,
                                                                         get_monotonic_fit=False)

# fit Sigmoidal function
y_bounds = [min(intensity_range),max(intensity_range)]
v_tresh_bounds = y_bounds
if haz_var=='MESHS' or haz_var == 'MESHS_4km':
    v_tresh_bounds = (0,20) #visually detemined (MDR/PAA is >0 at MESHS=20)
elif haz_var == 'E_kin' or haz_var == 'E_kinCC':
    v_tresh_bounds = (0,200) #visually determined
pbounds={'v_thresh': v_tresh_bounds, 'v_half': y_bounds, 
         'scale': [df.MDR.max()/10, min(1,df.MDR.max()*10)],
          'power': (3,3)
          }
p,res,impf_emanuel = sc.calib_opt.fit_emanuel_impf_to_emp_data(df,pbounds)

######### Emanuel fit for every bootstrapped sample ##########################
df_boot_emanuel = pd.DataFrame(index = intensity_range, dtype=float,
                               columns = [f'b_{i}' for i in range(n_samples)])
for i in range(n_samples):
    print(i)
    df_now = ds_boot.isel(b_sample=i).to_dataframe()
    pNow,resNow,impf_emanuelNow = sc.calib_opt.fit_emanuel_impf_to_emp_data(df_now,pbounds,plot=False)
    df_boot_emanuel.loc[impf_emanuelNow.intensity,f'b_{i}']=impf_emanuelNow.mdd

fig,ax=plt.subplots()
ax.plot(df.index,impf_emanuel.mdd*100)
sc.plot_funcs.fill_df_quantile(df_boot_emanuel,0.1,0.9,ax)
ax.legend()
##############################################################################

pboundsPAA={'v_thresh': v_tresh_bounds, 'v_half': y_bounds,
            'scale': [df_roll.PAA.max()/10, min(1,df.PAA.max()*10)]}
p_PAA,resPAA,impf_emanuelPAA = sc.calib_opt.fit_emanuel_impf_to_emp_data(df,pboundsPAA,opt_var='PAA')

######### Emanuel fit for every bootstrapped sample PAA ######################
df_boot_emanuel_PAA = pd.DataFrame(index = intensity_range,dtype=float,
                                   columns = [f'b_{i}' for i in range(n_samples)])
for i in range(n_samples):
    print(i)
    df_now = ds_boot.isel(b_sample=i).to_dataframe()
    pNow,resNow,impf_emanuelNow = sc.calib_opt.fit_emanuel_impf_to_emp_data(df_now,pboundsPAA,
                                                                            opt_var='PAA',plot=False)
    df_boot_emanuel_PAA.loc[impf_emanuelNow.intensity,f'b_{i}']=impf_emanuelNow.mdd

fig,ax=plt.subplots()
ax.plot(df.index,impf_emanuelPAA.mdd*100)
sc.plot_funcs.fill_df_quantile(df_boot_emanuel_PAA,0.1,0.9,ax)
ax.legend()
##############################################################################

df_roll.loc[df_roll.MDD==np.inf,'MDD']=np.nan
pboundsMDD={'v_thresh': y_bounds, 'v_half': y_bounds,
            'scale': [df.MDD.max()/10, df.MDD.max()*5]}
p_MDD,resMDD,impf_emanuelMDD = sc.calib_opt.fit_emanuel_impf_to_emp_data(df_roll,pboundsMDD,opt_var='MDD')


# # #convert to rolling xarray dataset
if get_per_gridcell:
    if not haz_var == 'MESHS':
        raise Warning("windowsize for gc_df is only done for MESHS")
    for i, start in enumerate(intensity_range[1:]):
        print(start)
        now_cond = np.logical_and(values_at_centroid_all[haz_var] > 
                                  start-np.floor(unit_windowsize/2), 
                                  values_at_centroid_all[haz_var] < 
                                  start+np.floor(unit_windowsize/2))
        gc_now = values_at_centroid_all.where(now_cond).dropna(axis=0, how='all')
        ds_temp = gc_now.to_xarray().drop_vars(
            haz_var).expand_dims(dim={haz_var: [start]})
        ds_temp['index'] = np.arange(0, len(gc_now.index))
        if i == 0:
            ds_roll_gc = ds_temp.copy(deep=True)
        else:
            ds_roll_gc = xr.concat([ds_roll_gc, ds_temp], dim=haz_var)


# %%#Plots
relative_bins=True if exposure=='MFZrandom_' else False
# plot rolling PAA
title=(f'Empirical Percent of Assets Affected ({exposure});\n '
       f'{unit_windowsize} {unit} moving average, {baujahr}')
fig,ax=sc.impf_plot(df_all,df_roll,df_roll_cut,ds_boot_roll,ds_boot_roll_cut,
                    haz_var,impf_emanuelPAA,cut_off,'PAA',title,dmg_bin_size,
                    intensity_label)
fig.savefig(out_dir + '/impf/%s_PAA_roll_.png' %(fig_path), dpi=300, 
            bbox_inches='tight')


#plot combined figure
fig = plt.figure(figsize=(11.5,3.5))
gs = GridSpec(1,2, figure=fig,wspace=0.4)
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
fig,_=sc.impf_plot2(df_all,df_roll,df_roll_cut,ds_boot_roll,ds_boot_roll_cut,
                    haz_var,impf_emanuelPAA,cut_off,'PAA',"",dmg_bin_size,
                    intensity_label,df_boot_emanuel=df_boot_emanuel_PAA,
                    relative_bins=relative_bins,figAx=(fig,ax0))
fig,_=sc.impf_plot2(df_all,df_roll,df_roll_cut,ds_boot_roll,ds_boot_roll_cut,
                    haz_var,impf_emanuel,cut_off,'MDR',"",dmg_bin_size,
                    intensity_label,color='blue',df_boot_emanuel=df_boot_emanuel,
                    relative_bins=relative_bins,figAx=(fig,ax1))
fig.savefig(f"{out_dir}/impf/{fig_path}_combined.png", dpi=300, bbox_inches='tight')


# empirical MDR impact function
if haz_var == 'MESHS':
    intensity = np.arange(0, 100, 1)
elif haz_var == 'MESHSweigh':
    intensity = np.arange(0, max(intensity_range)+1, 1)
else:
    intensity = intensity_range
df_impf = pd.DataFrame(index=intensity, columns=['PAA', 'MDD', 'MDR'])
for var in ['PAA', 'MDD', 'MDR']:
    df_impf.loc[1:, var] = df_roll[var].loc[intensity_range[1:]]

    # add smoothed function
    # df_roll_cut.apply(sc.)
    smooth = sc.smooth_monotonic(
        df_roll_cut.index[1:], df_roll_cut[var].loc[intensity_range[1:]])
    # plt.plot(df_roll_cut.index[1:],smooth)

    df_impf.loc[intensity_range[1]:, var+'_smooth'] = smooth

#add Sigmoidal function fit
np.testing.assert_array_equal(impf_emanuelPAA.intensity,intensity_range)
df_impf.loc[intensity_range[1:], 'PAA_emanuel'] = impf_emanuelPAA.mdd[1:]
df_impf.loc[intensity_range[1:], 'MDR_emanuel'] = impf_emanuel.mdd[1:]
df_impf.loc[intensity_range[1:], 'MDD_emanuel'] = impf_emanuelMDD.mdd[1:]

# set values above cut_off to values at cut_off
cut_vars = ['PAA', 'MDD', 'MDR']
if haz_var == 'MESHS':  # cut off at 59 instead of 60 for realistic values
    df_impf.loc[cut_off-1:,
                cut_vars] = df_impf.loc[cut_off-1, cut_vars].values
elif haz_var != 'POH': # for POH there is no cut off!
    df_impf.loc[cut_off:, cut_vars] = df_impf.loc[cut_off, cut_vars].values
    
for var in ['PAA','MDR','MDD']:
    df_impf[[var,f'{var}_smooth',f'{var}_emanuel']].plot()

#Save impact function data as .csv file
df_impf.fillna(0).to_csv(
    data_dir+'/out_files/paa_mdd_smooth_%s%s_v%d%s.csv' % 
    (exposure,haz_var, event_def_version,baujahr))