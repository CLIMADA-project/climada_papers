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
Created on Monday, June 29 2020
    
description: Main for tropical cyclone (TC) impact function calibration

Requires:
    CLIMADA repository version 1.4.1+:
        https://wcr.ethz.ch/research/climada.html
        https://github.com/CLIMADA-project/climada_python
    TC track data from IBTrACS v4, 1980-2017 (get data online or ask the author):
        https://www.ncdc.noaa.gov/ibtracs/
    EM-DAT data for tropical cyclone 1980-2017 (get data from EM-DAT or ask the author):
        https://www.emdat.be
        https://public.emdat.be/

Please refer to the paper (submitted) and README for more information.

@author: Samuel Eberenz, samuel.eberenz@usys.ethz.ch
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from iso3166 import countries as iso_cntry

from climada.engine import Impact
from climada.entity import IFTropCyclone, ImpactFuncSet
#from climada.entity.impact_funcs.trop_cyclone import IFSTropCyclone
from if_trop_cyclone_stable_202006 import IFSTropCyclone, IFTropCyclone # stable version
# from climada.engine.impact_data import emdat_countries_by_hazard, emdat_to_impact
from impact_data_stable_202006 import emdat_countries_by_hazard, emdat_to_impact # stable version
import tc_calibration_functions as tc_cal

""" Variable initiation """
from tc_calibration_config import CALIB, DATA_DIR, EMDAT_CSV, EMDAT_MAP, \
    ENTITY_DIR, ENTITY_STR, HAZ, HAZARD_DIR, HAZARD_STR, REF_YEAR, \
    RES_ARCSEC, RES_DIR, RES_STR, TRACK_DIR, TRACK_FOLDER, YEAR_RANGE, \
    basins_short, def_log_ratio_threshold, \
    drop_excluded_cntries, exclude_cntries, keep_variables, \
    make_plots, mismatched_emdat_events, new_reg_names, region_ids_cal, \
    regions_short, scale, v0, v_half_range

# init dict with countries by basin:
basins_ids_cal = dict()
for basin in basins_short:
    basins_ids_cal[basin] = list()
    for region in regions_short:
        if basin in region:
            basins_ids_cal[basin] += region_ids_cal[region]
    basins_ids_cal[basin].sort()
# exclude countries
if drop_excluded_cntries:
    for reg in region_ids_cal:
        for cntry_ex in exclude_cntries:
            if cntry_ex in region_ids_cal[reg]:
                region_ids_cal[reg].remove(cntry_ex)


if (1 in CALIB or not (keep_variables and 'hazard' in locals() and 'exp_coast' in locals())) and not -1 in CALIB:
    """ CALIB 1: Loading or initiating HAZARD and EXPOSURE sets: """
    print('\n...........................................................\n')
    print('CALIB 1: Loading or initiating HAZARD and EXPOSURE sets')
    print('.............................................................')
    hazard, exp_coast, _ = tc_cal.prepare_calib_data(EMDAT_CSV, \
                                os.path.join(TRACK_DIR, TRACK_FOLDER), ENTITY_DIR, \
                                HAZARD_DIR, ENTITY_STR, HAZARD_STR, \
                                make_plots=make_plots, init_hazard=True, \
                                hazard_type=HAZ, res_arcsec=RES_ARCSEC, ref_year=REF_YEAR)

else:
    print("--------Keeping existing variables: hazard, exp_coast-------------")

if 2 in CALIB:
    """ CALIB 2: With one global impact function (V_0=25.7 m/s; V_half=74.7m/s, scale=1),
    calculate simulated event damage (SED) per TC event and get normalized 
    reported damage (NRD) from EM-DAT.
    compute the event damage ratio EDR=SED/NRD.
    Save to CSV file.
    Required as a starting point for calibration, i.e. CALIB 3 to 5"""
    print('\n...........................................................\n')
    print('Calibration step CALIB 2')
    print('.............................................................')
    # init data:
    calib2_file_str = 'CALIB2_event_analysis'
    if os.path.isfile(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, \
                                   HAZ, calib2_file_str, 'ALL', 'csv'))):
        if (not 'results_CVA' in locals()) or not keep_variables:
            print("----------------------Loading Country Vulnerability Analysis----------------------")
            results_CVA = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, \
                                          HAZ, calib2_file_str, 'ALL', 'csv')), \
                                          encoding="ISO-8859-1", header=0)
            # reset basin to get read of import error ('NA'-->nan):
            for basin in basins_short:
                for idx in results_CVA.index:
                    if results_CVA.loc[idx, 'country'] in basins_ids_cal[basin]:
                        results_CVA.loc[idx, 'cal_region'] = basin

    else:
        results_CVA, fail_list = tc_cal.country_vulnerability_analysis(exp_coast, \
                     hazard, [], basins_ids_cal, \
                     EMDAT_CSV, EMDAT_MAP, RES_DIR, RES_STR, RES_ARCSEC, \
                     parameters = (25.7, 49, 1), \
                     yearly_impact=False, hazard_type=HAZ, \
                     year_range = YEAR_RANGE, ref_year=REF_YEAR)
        print('Fail list: ')
        print(fail_list)
        # results = results[results.EM_ID != 'NONE']
        results_CVA = results_CVA.reset_index(drop=True)
        results_CVA.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, calib2_file_str, 'ALL', 'csv')), index=False)

    # compute EDR (ratio between CLIMADA and EMDAT):
    edr = np.asarray(results_CVA['climada_impact']/results_CVA['emdat_impact_scaled']).astype(np.float64)
    log_ratio = np.log(edr)
    results_CVA_ratio = results_CVA
    results_CVA_ratio['ratio'] = edr
    results_CVA_ratio['log_ratio'] = log_ratio
    # clean up invalid numbers:
    results_CVA_ratio.replace([np.inf, -np.inf], np.nan)
    results_CVA_ratio = results_CVA_ratio[results_CVA_ratio['log_ratio'] > -np.inf]
    results_CVA_ratio = results_CVA_ratio[results_CVA_ratio['log_ratio'] < np.inf] # 534 data points
    
    # Compute number of data points (N) with ratio and median ratio EDR per calibration region:
    regions_N = dict()
    regions_Median = dict()
    regions_N['ALL'] = results_CVA_ratio.shape[0]
    regions_Median['ALL'] = results_CVA_ratio.median(axis=0)['ratio']
    
    # N and Median per country and region, introducing cal_region2:
    countries_N_median = pd.DataFrame(index=np.arange(0, np.unique(results_CVA_ratio['country']).size), \
                               columns=['country', 'region_id', 'cal_region', \
                                        'cal_region2', \
                                        'median_ratio', 'N'])
    for cal_reg in np.unique(results_CVA_ratio['cal_region']):
        results_CVA_ratio_reg = results_CVA_ratio[results_CVA_ratio['cal_region']==cal_reg]
        regions_N[cal_reg] = results_CVA_ratio_reg.shape[0]
        regions_Median[cal_reg] = results_CVA_ratio_reg.median(axis=0)['ratio']
    for ind, cntry in enumerate(np.unique(results_CVA_ratio['country'])):
        # print(str(ind) + ': ' + cntry)
        countries_N_median.loc[ind,"country"] = cntry
        countries_N_median.loc[ind,"region_id"] = iso_cntry.get(cntry).numeric
        countries_N_median.loc[ind,"cal_region"] = results_CVA_ratio.loc\
            [results_CVA_ratio['country']==cntry,"cal_region"].values[0]
        countries_N_median.loc[ind,"N"] = \
            results_CVA_ratio.loc[results_CVA_ratio['country']==cntry,"ratio"].shape[0]
        countries_N_median.loc[ind,"median_ratio"] = \
            results_CVA_ratio.loc[results_CVA_ratio['country']==cntry,"ratio"].median()

    # Attribute events / countries to 2nd set of calibration region (further
    # splitting up of basins NA and WP into sub-regions:
    for ind, cntry in enumerate(countries_N_median['country']):
        for r2 in regions_short:
            if cntry in region_ids_cal[r2]:
                countries_N_median.loc[ind,"cal_region2"] = r2
                continue
            
    results_CVA_ratio = results_CVA_ratio.reset_index(drop=True)
    for ind in np.arange(results_CVA_ratio.shape[0]):
        for r2 in regions_short:
            if results_CVA_ratio.loc[ind, "country"] in region_ids_cal[r2]:
                results_CVA_ratio.loc[ind,"cal_region2"] = r2
                continue
    for cntry_ex in exclude_cntries: # change cal_region2 of excluded countries to 'None'
        results_CVA_ratio.loc[results_CVA_ratio['country']==cntry_ex,'cal_region2'] = 'None'
    
    if new_reg_names:
        results_CVA_ratio = tc_cal.update_regions(results_CVA_ratio, region_ids_cal)[0]
    # save results to CSV:
    results_CVA_ratio.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                             REF_YEAR, HAZ, 'event_analysis_ratio', 'ALL', 'csv')), index=False)

    regions2_N = dict()
    regions2_Median = dict()
    regions2_N['ALL'] = results_CVA_ratio.shape[0]
    regions2_Median['ALL'] = results_CVA_ratio.median(axis=0)['ratio']
    for r2 in region_ids_cal.keys():
        countries_N_median_tmp = countries_N_median[countries_N_median['cal_region2']==r2]
        regions2_N[r2] = countries_N_median_tmp['N'].sum()
        regions2_Median[r2] = results_CVA_ratio_reg.median(axis=0)['ratio']
    
    for cal_reg in np.unique(results_CVA_ratio['cal_region']):
        countries_region = countries_N_median.loc[countries_N_median['cal_region'] == cal_reg, :]
        countries_region = countries_region.sort_values('median_ratio')

if 3 in CALIB:
    """ Core calibration (slow!)
    Loop over v_half_range to compute event damage (SED) and ratio EDR for each 
    matched event."""
    print('\n...........................................................\n')
    print('Calibration step CALIB 3')
    print('.............................................................')
    if (not 'results_CVA_ratio' in locals()) or not keep_variables:
        try:
            results_CVA_ratio = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                              REF_YEAR, HAZ, 'event_analysis_ratio', 'ALL', 'csv')),\
                                          encoding="ISO-8859-1", header=0)
            results_CVA_ratio.shape
        except FileNotFoundError:
            print('ERROR: File not found: ' + os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                              REF_YEAR, HAZ, 'event_analysis_ratio', 'ALL', 'csv')))

    results_CALIB3 = pd.DataFrame() # init results DataFRame
    # loop over regions:
    for cal_reg in regions_short:
        regions_dict_tmp  = dict()
        regions_dict_tmp[cal_reg] = region_ids_cal[cal_reg]

        ### CORE CALIBRATION ENGINE: Computing SED (simulated event damage) for all events and complete
        ### impact function parameter set (brute force takes long. run on cluster):
        results_CALIB3_tmp, fail_list = tc_cal.country_calib_3(exp_coast, hazard, \
                                            results_CVA_ratio, regions_dict_tmp, \
                                            parameter_space = \
                                            [np.array(v0), \
                                            v_half_range-v0[0], \
                                            np.array(scale)], \
                                            yearly_impact=False, \
                                            year_range=YEAR_RANGE)
        # save results at each step (in case script is interrupted)
        results_CALIB3_tmp.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, \
                                HAZ, 'event_CALIB3_v_half', cal_reg, 'csv')), index=False)
        # append results DataFrame an clean up:
        results_CALIB3 = results_CALIB3.append(results_CALIB3_tmp)
        results_CALIB3 = results_CALIB3.reset_index(drop=True)
        for cntry_ex in exclude_cntries: # change cal_region2 of excluded countries to 'None'
            results_CALIB3.loc[results_CALIB3['country']==cntry_ex, 'cal_region2'] = 'None'
        # ensure backward compatibility:
        if new_reg_names:
            results_CALIB3 = tc_cal.update_regions(results_CALIB3, region_ids_cal)[0]

        results_CALIB3 = tc_cal.fix_data(results_CALIB3, mismatches=mismatched_emdat_events, \
                                          year_range=YEAR_RANGE, \
                                          save_path=None, rm_mismatches=True)
        # save combined results of CALIB 3:
        results_CALIB3.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, \
                                    HAZ, 'event_CALIB3_v_half', 'all', 'csv')), index=False)


if 4 in CALIB:
    """Optimization: Compute cost functions and find optimized v_half for each cost function and region,
    Compute statistics per region. Save calibration results to CSV.
    Requires CALIB 2 and 3. Required for CALIB 5."""
    print('\n...........................................................\n')
    print('Calibration step CALIB 4')
    print('.............................................................')

    # get required input data, i.e. th results of CALIB 3:
    if not (keep_variables and 'results_CALIB3' in locals()):
        try:
            results_CALIB3 = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB3_v_half', 'all', 'csv')),\
                                          encoding="ISO-8859-1", header=0)
        except FileNotFoundError:
            print('ERROR: File not found: ' + os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB3_v_half', 'all', 'csv')))

    print('Post-processing CALIB3 results...')
    # Clean up and potential corrections in calibration regions:
    print('events in df: %i' %(int(results_CALIB3.shape[0]/v_half_range.size)))
    if not 'None' in list(results_CALIB3.cal_region2):
        for cntry_ex in exclude_cntries: # change cal_region2 of excluded countries to 'None'
            results_CALIB3.loc[results_CALIB3['country']==cntry_ex, 'cal_region2'] = 'None'
        results_CALIB3.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB3_v_half', 'all', 'csv')), index=False)
    if new_reg_names and not 'OC' in list(results_CALIB3.cal_region2):
        results_CALIB3 = tc_cal.update_regions(results_CALIB3, region_ids_cal)[0]
        # results_CALIB3.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB3_v_half', 'all', 'csv')), index=False)
        results_CALIB3 = tc_cal.fix_data(results_CALIB3, mismatches=mismatched_emdat_events, year_range=YEAR_RANGE, \
                                          save_path=os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB3_v_half', 'all', 'csv')), rm_mismatches=False)
    print('events in df: %i' %(int(results_CALIB3.shape[0]/v_half_range.size)))
    results_CALIB3 = tc_cal.fix_data(results_CALIB3, mismatches=mismatched_emdat_events, year_range=YEAR_RANGE, \
                                      save_path=None, rm_mismatches=True)
    print('events in df after removing mismatches: %i' %(int(results_CALIB3.shape[0]/v_half_range.size)))
    if drop_excluded_cntries:
        results_CALIB3 = results_CALIB3[results_CALIB3.cal_region2!='None']
    print('events in df after dropping excluded countries: %i' %(int(results_CALIB3.shape[0]/v_half_range.size)))
    print('Find best fit to log(ratio)=0 for each event/country combo...')    
    # drop events with log of ratio (EDR) with default IF above or below threshold
    if def_log_ratio_threshold:
        rm_list_dict = dict()
        results_CALIB3_default_IF = results_CALIB3.loc[np.round(results_CALIB3.v_half, decimals=2)==74.70]
        rm_list_dict['EM_ID'] = list(results_CALIB3_default_IF.loc[(results_CALIB3_default_IF.log_ratio>=def_log_ratio_threshold) | (results_CALIB3_default_IF.log_ratio<=-def_log_ratio_threshold), 'EM_ID'].values)
        rm_list_dict['country'] = list(results_CALIB3_default_IF.loc[(results_CALIB3_default_IF.log_ratio>=def_log_ratio_threshold) | (results_CALIB3_default_IF.log_ratio<=-def_log_ratio_threshold), 'country'].values)
        for idx, _ in enumerate(rm_list_dict['EM_ID']):
            results_CALIB3 = results_CALIB3[(results_CALIB3.EM_ID!=rm_list_dict['EM_ID'][idx]) | (results_CALIB3.country!=rm_list_dict['country'][idx])]
    print('events in df after applying ratio threshold: %i' %(int(results_CALIB3.shape[0]/v_half_range.size)))

    ### Start of OPTIMIZATION:
    # call function to get v_half with EDR closest to one (SED=NRD) for each event:
    results_closest_ratio = tc_cal.closest_ratio(results_CALIB3, v0, v_half_range, scale, target=0)
    results_closest_ratio.reset_index(drop=True)
    for cntry_ex in exclude_cntries: # change cal_region2 of excluded countries to 'None'
        results_closest_ratio.loc[results_closest_ratio['country']==cntry_ex, 'cal_region2'] = 'None'
    results_closest_ratio = tc_cal.get_associated_disasters(results_closest_ratio, EMDAT_CSV)
    results_closest_ratio = tc_cal.fix_data(results_closest_ratio, mismatches=mismatched_emdat_events, \
                                        year_range=YEAR_RANGE, rm_mismatches=True, \
                                        save_path=os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                        REF_YEAR, HAZ, 'event_CALIB4_closest_log_ratio', 'all', 'csv')))

    # RMSF Optimization (RMSF=minimum):
    rmsf_results, min_rmsf = tc_cal.compute_metric_min(results_CALIB3, metric='RMSF', save=True)

    min_rmsf.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                  'CALIB4_%s_best_vhalf' % ('RMSF'), 'all', 'csv')), index=False)
    rmsf_results.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                  'CALIB4_%s_full_results_vhalf' % ('RMSF'), 'all', 'csv')), index=False)

    # TDR Optimization (TDR=1):
    tot_results, best_tot_vhalf = tc_cal.compute_vhalf_total_impact(results_CALIB3)
    tot_results_1_5, best_tot_vhalf_1_5 = tc_cal.compute_vhalf_total_impact(results_CALIB3, scaling_emdat=1.5)

    tot_results.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                    'CALIB4_tot_impact_res_scaling_%1.2f' % (1), 'all', 'csv')), index=False)
    best_tot_vhalf.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                        'CALIB4_tot_impact_best_vhalf_scaling_%1.2f' % (1), 'all', 'csv')), index=False)
    tot_results_1_5.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                    'CALIB4_tot_impact_res_scaling_%1.2f' % (1.5), 'all', 'csv')), index=False)
    best_tot_vhalf_1_5.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                        'CALIB4_tot_impact_best_vhalf_scaling_%1.2f' % (1.5), 'all', 'csv')), index=False)

    # prepare and save table with result overview:
    results_CALIB3_default_IF = results_CALIB3.loc[np.round(results_CALIB3.v_half, decimals=2)==74.70]
    glob_rmsf_rmsf, rmsf_control = tc_cal.compute_global_rmsf(results_CALIB3, min_rmsf, \
                                        list(np.round(min_rmsf.v_half.loc[min_rmsf.cal_region2=='GLB'], \
                                        decimals=3)) + [74.7])
    glob_tdr_rmsf, _ = tc_cal.compute_global_rmsf(results_CALIB3, best_tot_vhalf, \
                                        list(np.round(min_rmsf.v_half.loc[min_rmsf.cal_region2=='GLB'], \
                                        decimals=3)) + [74.7])
    result_table, countries_list = tc_cal.fill_result_table(regions_short, rmsf_results, min_rmsf, \
                    tot_results, best_tot_vhalf, results_CALIB3_default_IF, glob_rmsf_rmsf, glob_tdr_rmsf)

    countries_list.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, '_results_countries_list', 'all', 'csv')), index=False)
    result_table.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, '_results_table', 'all', 'csv')), index=False)

    # Export data for Philippines:
    # results_CALIB3_PHL_RMSF_IF = results_CALIB3.loc[np.round(results_CALIB3.v_half, decimals=2)\
    #             == np.round(min_rmsf.loc[min_rmsf.cal_region2=='WP2', 'v_half'].values[0], decimals=2)]
    # results_CALIB3_PHL_RMSF_IF.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_CALIB4_PHL_RMSF_IF', 'all', 'csv')), index=False)

    results_CALIB3_regional_RMSF_IF = pd.DataFrame(columns=list(results_CALIB3.columns))
    for idx, reg in enumerate(regions_short):
        df_tmp = results_CALIB3.loc[(results_CALIB3.cal_region2==reg) &\
              (np.round(results_CALIB3.v_half, decimals=2) \
              == np.round(min_rmsf.loc[min_rmsf.cal_region2==reg, 'v_half'].values[0], decimals=2))]
        results_CALIB3_regional_RMSF_IF = pd.concat([results_CALIB3_regional_RMSF_IF, df_tmp], ignore_index=True, sort =False)
    results_CALIB3_regional_RMSF_IF.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                           'results_CALIB4_regional_RMSF_IF', 'all', 'csv')), index=False)
    results_CALIB3_regional_TDR_IF = pd.DataFrame(columns=list(results_CALIB3.columns))
    for idx, reg in enumerate(regions_short):
        df_tmp = results_CALIB3.loc[(results_CALIB3.cal_region2==reg) &\
              (np.round(results_CALIB3.v_half, decimals=2) \
              == np.round(best_tot_vhalf.loc[best_tot_vhalf.cal_region2==reg, 'v_half'].values[0], decimals=2))]
        results_CALIB3_regional_TDR_IF = pd.concat([results_CALIB3_regional_TDR_IF, df_tmp], \
                                                   ignore_index=True, sort =False)
    results_CALIB3_regional_TDR_IF.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                          'results_CALIB4_regional_TDR_IF', 'all', 'csv')), index=False)

    # compute and save basic sctatistics (linear regression):
    linregress_df = pd.DataFrame(index=['slope', 'intercept', 'r_value', 'p_value', 'std_err'], \
                                  columns=['default', 'RMSF_opt', 'TDR_opt'])
    r2_df = pd.DataFrame(index=regions_short + ['GLB'], \
                                  columns=['default', 'RMSF_opt', 'TDR_opt'])
    std_err_df = pd.DataFrame(index=regions_short + ['GLB'], \
                                  columns=['default', 'RMSF_opt', 'TDR_opt'])
    col = 'default'
    linregress_df.loc['slope', col], linregress_df.loc['intercept', col], \
        linregress_df.loc['r_value', col], linregress_df.loc['p_value', col], \
        linregress_df.loc['std_err', col] = stats.linregress(results_CALIB3_default_IF.emdat_impact_scaled, \
                                results_CALIB3_default_IF.climada_impact)
    col = 'RMSF_opt'
    linregress_df.loc['slope', col], linregress_df.loc['intercept', col], \
        linregress_df.loc['r_value', col], linregress_df.loc['p_value', col], \
        linregress_df.loc['std_err', col] = stats.linregress(results_CALIB3_regional_RMSF_IF.emdat_impact_scaled, \
                                results_CALIB3_regional_RMSF_IF.climada_impact)
    col = 'TDR_opt'
    linregress_df.loc['slope', col], linregress_df.loc['intercept', col], \
        linregress_df.loc['r_value', col], linregress_df.loc['p_value', col], \
        linregress_df.loc['std_err', col] = stats.linregress(results_CALIB3_regional_TDR_IF.emdat_impact_scaled, \
                                results_CALIB3_regional_TDR_IF.climada_impact)
    linregress_df.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_CALIB4_linregress_cal', 'all', 'csv')), index=True)

    for reg in regions_short:
        col = 'default'
        _, _, r2_df.loc[reg, col], _, std_err_df.loc[reg, col] = stats.linregress(results_CALIB3_default_IF.loc[results_CALIB3_default_IF.cal_region2==reg].emdat_impact_scaled, \
                                    results_CALIB3_default_IF.loc[results_CALIB3_default_IF.cal_region2==reg].climada_impact)
        col = 'RMSF_opt'
        _, _, r2_df.loc[reg, col], _, std_err_df.loc[reg, col] = stats.linregress(results_CALIB3_regional_RMSF_IF.loc[results_CALIB3_regional_RMSF_IF.cal_region2==reg].emdat_impact_scaled, \
                                    results_CALIB3_regional_RMSF_IF.loc[results_CALIB3_regional_RMSF_IF.cal_region2==reg].climada_impact)
        col = 'TDR_opt'
        l_, _, r2_df.loc[reg, col], _, std_err_df.loc[reg, col] = stats.linregress(results_CALIB3_regional_TDR_IF.loc[results_CALIB3_regional_TDR_IF.cal_region2==reg].emdat_impact_scaled, \
                                    results_CALIB3_regional_TDR_IF.loc[results_CALIB3_regional_TDR_IF.cal_region2==reg].climada_impact)
    r2_df.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_CALIB4_linregress_r2_reg_cal', 'all', 'csv')), index=True)
    std_err_df.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_CALIB4_linregress_std_err_df_reg_cal', 'all', 'csv')), index=True)
    linregress_df.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_CALIB4_linregress_cal', 'all', 'csv')), index=True)

    reg_stats = pd.DataFrame(index=np.arange(0, len(region_ids_cal)), \
                              columns=['cal_region2', 'N_country_events', 'N_events', 'quantile .25', 'median', 'quantile .75', 'IQR', 'min', 'max', 'lim_min', 'lim_max'])
    for idx, cal_reg in enumerate(region_ids_cal): # ['CAR']: # region_ids_cal:
        results_closest_ratio_tmp = results_closest_ratio.loc[results_closest_ratio.cal_region2==cal_reg]
        reg_stats.loc[idx, 'cal_region2'] = cal_reg
        reg_stats.loc[idx, 'N_country_events'] = results_closest_ratio_tmp.shape[0]
        reg_stats.loc[idx, 'N_events'] = len(results_closest_ratio_tmp.EM_ID.unique())
        reg_stats.loc[idx, 'min'] = results_closest_ratio_tmp['v_half'].min()
        reg_stats.loc[idx, 'max'] = results_closest_ratio_tmp['v_half'].max()
        reg_stats.loc[idx, 'median'] = results_closest_ratio_tmp['v_half'].median()
        reg_stats.loc[idx, 'quantile .25'] = results_closest_ratio_tmp['v_half'].quantile(q=0.25)
        reg_stats.loc[idx, 'quantile .75'] = results_closest_ratio_tmp['v_half'].quantile(q=0.75)
        reg_stats.loc[idx, 'lim_min'] = v_half_range.min()
        reg_stats.loc[idx, 'lim_max'] = v_half_range.max()
    reg_stats['IQR'] = reg_stats['quantile .75'] - reg_stats['quantile .25']
    reg_stats.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB4_reg_stats', 'all', 'csv')), index=False)

if 5 in CALIB or 6 in CALIB:
    """prepare data and calibrated impact functions for CALIB 5 and 6"""

    # get input data (from CALIB 4):
    if (not 'results_closest_ratio' in locals()) or not keep_variables:
        try:
            results_closest_ratio = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB4_closest_log_ratio', 'all', 'csv')),\
                                      encoding="ISO-8859-1", header=0)
        except:
            print('ERROR: File not found: ' + os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB4_closest_log_ratio', 'all', 'csv')))
    if (not 'min_rmsf' in locals()) or not keep_variables:
        try:
            min_rmsf = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                  'CALIB4_%s_best_vhalf' % ('RMSF'), 'all', 'csv')),\
                                      encoding="ISO-8859-1", header=0)
        except:
            print('ERROR: File not found: ' + os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                      'CALIB4_%s_full_results_vhalf' % ('RMSF'), 'all', 'csv')))
    if (not 'best_tot_vhalf' in locals()) or not keep_variables:
        try:
            best_tot_vhalf = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                        'CALIB4_tot_impact_best_vhalf_scaling_%1.2f' % (1), 'all', 'csv')),\
                                      encoding="ISO-8859-1", header=0)
        except:
            print('ERROR: File not found: ' + os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, \
                                      'CALIB4_tot_impact_res_scaling_%1.2f' % (1), 'all', 'csv')))

    if_TC_all = np.ones(exp_coast.shape[0], int)*(1+len(regions_short))
    results_CALIB4 = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'event_CALIB4_reg_stats', 'all', 'csv')))
    
    if new_reg_names and (not 'OC' in list(results_closest_ratio.cal_region2) \
                          or not 'OC' in list(results_CALIB4.cal_region2)):
        results_CALIB4 = tc_cal.update_regions(results_CALIB4, region_ids_cal)[0]
        results_closest_ratio = tc_cal.update_regions(results_closest_ratio, region_ids_cal)[0]

    regions_ids_matched = dict()
    # Set impact function ID (if_) in exposure for each country to link 
    # to assigned region:
    for idx, region in enumerate(regions_short):
        regions_ids_matched[region] = list(np.unique(results_closest_ratio[results_closest_ratio.cal_region2==region].country))
        for cntry in region_ids_cal[region]:
            try:
                cntry_num = int(iso_cntry.get(cntry).numeric)
            except:
                cntry_num = 0
            if_TC_all[exp_coast.region_id==cntry_num] = idx+1
    exp_coast['if_'] = if_TC_all
    exp_coast['if_TC'] = if_TC_all

    # reduce hazard to year range:
    hazard_year_range = hazard.select(date=(str(YEAR_RANGE[0])+'-01-01', \
                                      str(YEAR_RANGE[-1])+'-12-31'), \
                                      reset_frequency=True)

    # Init Impact Function sets (IFS) (calibrated and uncalibrated):
    quantiles = [.05, .25, .5, .75, .95]
    IFSs = dict()
    v_half = dict()
    IFS_default = ImpactFuncSet()
    if_tc = IFTropCyclone()
    if_tc.haz_type = HAZ
    if_tc.set_emanuel_usa(v_thresh=25.7, v_half=74.7, scale=1)
    if_tc.id = 1
    IFS_default.append(if_tc)

    event_count = pd.DataFrame(index=['GLB'] + regions_short, columns=['CLIMADA', 'EM-DAT', 'matched'])
    event_country_count = pd.DataFrame(index=['GLB'] + regions_short, columns=['CLIMADA', 'EM-DAT', 'matched'])
    impact_count = pd.DataFrame(index=['GLB'] + regions_short, columns=['CLIMADA', 'CLIMADA matched', 'EM-DAT scaled', 'EM-DAT', 'EM-DAT matched', 'EM-DAT matched scaled'])


    for q in quantiles:
        IFSs['%.2f' % (q)] = IFSTropCyclone()
        IFSs['%.2f' % (q)].set_calibrated_regional_IFs(calibration_approach='EDR', q=q, \
                                       input_file_path=results_closest_ratio, version=1)
    IFSs['RMSF'] = IFSTropCyclone()
    IFSs['TDR'] = IFSTropCyclone()
    IFSs['RMSF'].set_calibrated_regional_IFs(calibration_approach='RMSF', q=q, \
                                             input_file_path=min_rmsf, version=1)
    IFSs['TDR'].set_calibrated_regional_IFs(calibration_approach='TDR', q=q, \
                                            input_file_path=best_tot_vhalf, version=1)

if 5 in CALIB:
    print('\n...........................................................\n')
    print('Calibration step CALIB 5')
    print('.............................................................')
    """Compute AAD per country for comparison with GAR 2013 (calibrated and uncalibrated)"""
    countries_list = pd.read_csv(os.path.join(RES_DIR, RES_STR % \
                                (RES_ARCSEC, REF_YEAR, HAZ, '_results_countries_list', 'all', 'csv')),\
                                  encoding="ISO-8859-1", header=0)
    gar13_df = pd.read_csv(os.path.join(DATA_DIR, 'GAR2013_COUNTRY_RISK.csv'),\
                                  encoding="ISO-8859-1", header=0)
    countries_TC_emdat = emdat_countries_by_hazard(HAZ, EMDAT_CSV, \
                                  year_range=YEAR_RANGE, target_version=2018)[0]
    countries_exposure = pd.read_csv(os.path.join(DATA_DIR, '202004_metadata_countries_v1_2.csv'),\
                                  encoding="ISO-8859-1", header=0)

    # Init DataFrame for AAD results per country:
    res_cntrs = pd.DataFrame(index=np.arange(len(list(countries_exposure.iso3.unique()))), \
                             columns=['country_id', 'alpha3', 'name', 'cal', 'region', 'SIDS', \
                                        'AAD_emdat', 'PC_gar13', 'AAD_gar13', 'AAD/PC gar13 [permil]', \
                                        'PC_2014', 'AAD_def', 'AAD_rmsf', 'AAD_tdr',\
                                        'AAD/PC def [permil]', 'AAD/PC rmsf [permil]', 'AAD/PC tdr [permil]'])
    res_cntrs.alpha3 = list(countries_exposure.iso3.unique())

    # Loop through countries, compute AAD and fill res_cntrs:
    for idx, iso in enumerate(res_cntrs.alpha3):
        res_cntrs.loc[idx, 'country_id'] = int(iso_cntry.get(iso).numeric)
        res_cntrs.loc[idx, 'name'] = iso_cntry.get(iso).name

        if iso in countries_list.country.values:
            res_cntrs.loc[idx, 'cal'] = 'yes'
            res_cntrs.loc[idx, 'region'] = countries_list[countries_list.country.values==iso].region.values[0]
        else:
            res_cntrs.loc[idx, 'cal'] = 'no'
            res_cntrs.loc[idx, 'region'] = 'ROW'
            for reg in region_ids_cal:
                if iso in region_ids_cal[reg]:
                    res_cntrs.loc[idx, 'region'] = reg
                    #print(idx + ' ' + iso + ': ' + exp_coast)
                    break
        if iso in gar13_df['ISO code'].values:
            res_cntrs.loc[idx, 'SIDS'] = gar13_df.loc[gar13_df['ISO code']==iso].SIDS.values[0]
            res_cntrs.loc[idx, 'PC_gar13'] = gar13_df.loc[gar13_df['ISO code']==iso, 'Produced Capital [PC] TOTAL [USD million]'].values[0]
            res_cntrs.loc[idx, 'AAD_gar13'] = gar13_df.loc[gar13_df['ISO code']==iso, 'AAL [USD million]'].values[0]
            res_cntrs.loc[idx, 'AAD/PC gar13 [permil]'] = gar13_df.loc[gar13_df['ISO code']==iso, 'AAL/PC [permil]'].values[0]
        if iso in countries_TC_emdat:
            imp_emdat, _ = emdat_to_impact(EMDAT_CSV, \
                                        year_range=YEAR_RANGE, 
                                        countries=[iso],\
                                        hazard_type_climada=HAZ, \
                                        reference_year=REF_YEAR, \
                                        target_version=2018)
            res_cntrs.loc[idx, 'AAD_emdat'] = imp_emdat.aai_agg * 1e-6 # convert to million USD
        else:
            res_cntrs.loc[idx, 'AAD_emdat'] = 0
        res_cntrs.loc[idx, 'PC_2014'] = 1e-6 * countries_exposure.loc[countries_exposure.iso3==iso, 'total_value [USD]'].values[0]
        
        exp_coast['if_TC'] = if_TC_all
        exp_coast['if_'] = if_TC_all
        exp_country = exp_coast.loc[exp_coast['region_id']==res_cntrs.loc[idx, 'country_id']]
        if exp_country.size>0:
            imp_c = Impact()
            imp_c.calc(exp_country, IFSs['RMSF'], hazard_year_range) 
            res_cntrs.loc[idx, 'AAD_rmsf'] = 1e-6 * imp_c.aai_agg
            imp_c = Impact()
            imp_c.calc(exp_country, IFSs['TDR'], hazard_year_range) 
            res_cntrs.loc[idx, 'AAD_tdr'] = 1e-6 * imp_c.aai_agg
    
            if_TC_tmp = exp_country.if_TC.unique()[0]
            imp_c = Impact()
            exp_country.if_TC = 1
            exp_country.if_ = 1
            imp_c.calc(exp_country, IFS_default, hazard_year_range) 
            res_cntrs.loc[idx, 'AAD_def'] = 1e-6 * imp_c.aai_agg
            exp_country.if_TC = if_TC_tmp
            exp_country.if_ = if_TC_tmp
        else:
            res_cntrs.loc[idx, 'AAD_rmsf'] = 0
            res_cntrs.loc[idx, 'AAD_tdr'] = 0
            res_cntrs.loc[idx, 'AAD_def'] = 0
        
    res_cntrs['AAD/PC def [permil]'] = 1e3 * res_cntrs.AAD_def / res_cntrs.PC_2014
    res_cntrs['AAD/PC rmsf [permil]'] = 1e3 * res_cntrs.AAD_rmsf / res_cntrs.PC_2014
    res_cntrs['AAD/PC tdr [permil]'] = 1e3 * res_cntrs.AAD_tdr / res_cntrs.PC_2014
    # save AAD results per country in one CSV:
    res_cntrs.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_per_country', 'GAR13_4', 'csv')), index=False)

    # aggregate AAD results per calibration region:
    res_regs = pd.DataFrame(index=np.arange(len(regions_short)+1), columns=['region', \
                                        'AAD_emdat', 'PC_gar13', 'AAD_gar13', 'AAD/PC gar13 [permil]', \
                                        'PC_2014', 'AAD_def', 'AAD_rmsf', 'AAD_tdr', \
                                        'AAD/PC def [permil]', 'AAD/PC rmsf [permil]', 'AAD/PC tdr [permil]'])
    res_regs.region = regions_short + ['GLB']
    for idx, reg in enumerate(res_regs.region):
        for var in list(res_regs.columns[1:-3]):
            if reg=='GLB':
                res_regs.loc[idx, var] =  res_cntrs[var].sum()
            else:
                res_regs.loc[idx, var] =  res_cntrs.loc[(res_cntrs.region==reg) & (res_cntrs.cal=='yes'), var].sum()
    res_regs['AAD/PC gar13 [permil]'] = 1e3 * res_regs.AAD_gar13 / res_regs.PC_gar13
    res_regs['AAD/PC def [permil]'] = 1e3 * res_regs.AAD_def / res_regs.PC_2014
    res_regs['AAD/PC rmsf [permil]'] = 1e3 * res_regs.AAD_rmsf / res_regs.PC_2014
    res_regs['AAD/PC tdr [permil]'] = 1e3 * res_regs.AAD_tdr / res_regs.PC_2014
    # save AAD results per region in one CSV:
    res_regs.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, REF_YEAR, HAZ, 'results_per_region', 'GAR13_4', 'csv')), index=False)

if 6 in CALIB:
    """ Compute damage time series, trend and significance, as well as
        standard deviation of annual damage for EM-DAT and CLIMADA
        per country and region."""
    print('\n...........................................................\n')
    print('Calibration step CALIB 6')
    print('.............................................................')
    countries_list = pd.read_csv(os.path.join(RES_DIR, RES_STR % \
                                (RES_ARCSEC, REF_YEAR, HAZ, '_results_countries_list', 'all', 'csv')),\
                                 encoding="ISO-8859-1", header=0)

    # init dictionaries with yearly impact (yi) and statistics (stats):
    yi_climada = dict()
    yi_climada_Q3 = dict()
    yi_climada_Q1 = dict()
    yi_climada_RMSF = dict()
    yi_climada_tot = dict()
    yi_climada_def = dict()
    stats_climada = dict()
    stats_climada_Q3 = dict()
    stats_climada_Q1 = dict()
    stats_climada_RMSF = dict()
    stats_climada_tot = dict()
    stats_climada_def = dict()
    stats_log_climada = dict()
    stats_log_climada_Q3 = dict()
    stats_log_climada_Q1 = dict()
    stats_log_climada_RMSF = dict()
    stats_log_climada_tot = dict()
    stats_log_climada_def = dict()
    yi_emdat = dict()
    stats_emdat = dict()
    stats_log_emdat = dict()
    plt_axes = dict()
    plt_figs = dict()
    
    try:
        # try to load yearly impacts from files:
        for region in regions_short + ['GLB', 'GLB_all']:
            yi_climada[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_climada_Q3[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada_Q3', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_climada_Q1[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada_Q1', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_climada_RMSF[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada_RMSF', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_climada_tot[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada_tot', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_climada_def[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_climada_def', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
            yi_emdat[region] = pd.read_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                  REF_YEAR, HAZ, 'yi_emdat', region, 'csv')),\
                                              encoding="ISO-8859-1", header=0)
        print('Loaded yearly impacts from CSV')
    except FileNotFoundError:
        # simulate yearly impacts if not found in files:
        yi_climada['GLB'] = pd.DataFrame()
        yi_climada_Q3['GLB'] = pd.DataFrame()
        yi_climada_Q1['GLB'] = pd.DataFrame()
        yi_climada_RMSF['GLB'] = pd.DataFrame()
        yi_climada_tot['GLB'] = pd.DataFrame()
        yi_climada_def['GLB'] = pd.DataFrame()
        yi_emdat['GLB'] = pd.DataFrame()
        # loop over regions, compute yearly impacts and liner regression statstics:
        aai_agg = pd.DataFrame(index=regions_short+['GLB'], columns=['def', 'rmsf', 'tot'])
        
        
        for reg_num, region in enumerate(regions_short + ['WP_all']):
            countries = list()
            if region == 'WP_all':
                regions_list = region_ids_cal['WP1']+region_ids_cal['WP2']+region_ids_cal['WP3']+region_ids_cal['WP4']
            else:
                regions_list = region_ids_cal[region]
            for cntry in regions_list:
                if cntry in countries_list.country.values:
                    countries.append(cntry)
                    
            exp_coast['if_TC'] = if_TC_all
            exp_coast['if_'] = if_TC_all
            yi_climada[region], stats_climada[region], stats_log_climada[region], _ = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.5)], \
                                         hazard_year_range)
            yi_climada_Q3[region], stats_climada_Q3[region], stats_log_climada_Q3[region], _ = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.75)], \
                                         hazard_year_range)
            yi_climada_Q1[region], stats_climada_Q1[region], stats_log_climada_Q1[region], _ = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.25)], \
                                         hazard_year_range)
            yi_climada_RMSF[region], stats_climada_RMSF[region], stats_log_climada_RMSF[region], aai_agg.loc[region, 'rmsf'] = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['RMSF'], \
                                         hazard_year_range)
            yi_climada_tot[region], stats_climada_tot[region], stats_log_climada_tot[region], aai_agg.loc[region, 'tot'] = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['TDR'], \
                                         hazard_year_range)
            yi_climada_def[region], stats_climada_def[region], stats_log_climada_def[region], aai_agg.loc[region, 'def'] = \
                tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFS_default, \
                                         hazard_year_range)
    
            yi_emdat[region], stats_emdat[region], stats_log_emdat[region] = \
                tc_cal.trend_by_country_EMDAT(countries, year_range=YEAR_RANGE, reference_year=REF_YEAR, \
                               hazard_type_climada=HAZ)
    
            exp_coast['if_TC'] = if_TC_all
            exp_coast['if_'] = if_TC_all
            yi_climada[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada', region, 'csv')))
            yi_climada_Q3[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_Q3', region, 'csv')))
            yi_climada_Q1[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_Q1', region, 'csv')))
            yi_climada_RMSF[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_RMSF', region, 'csv')))
            yi_climada_tot[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_tot', region, 'csv')))
            yi_climada_def[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_def', region, 'csv')))
            yi_emdat[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_emdat', region, 'csv')))
            exp_coast['if_TC'] = if_TC_all
            exp_coast['if_'] = if_TC_all
    
            ddff=pd.DataFrame.from_dict(stats_climada[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_climada_Q3[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada_Q3[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada_Q3', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_climada_Q1[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada_Q1[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada_Q1', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_climada_RMSF[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada_RMSF[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada_RMSF', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_climada_tot[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada_tot[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada_tot', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_climada_def[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_climada_def[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_climada_def', region, 'csv')))
            ddff=pd.DataFrame.from_dict(stats_emdat[region], orient='index', columns=[region])
            ddff['log']=pd.DataFrame.from_dict(stats_log_emdat[region], orient='index', columns=[region])
            ddff.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'stats_log_emdat', region, 'csv')))

            if reg_num == 0:
                yi_climada['GLB']['year'] = yi_climada[region]['year']
                yi_climada_Q3['GLB']['year'] = yi_climada_Q3[region]['year']
                yi_climada_Q1['GLB']['year'] = yi_climada_Q1[region]['year']
                yi_climada_RMSF['GLB']['year'] = yi_climada_RMSF[region]['year']
                yi_climada_tot['GLB']['year'] = yi_climada_tot[region]['year']
                yi_climada_def['GLB']['year'] = yi_climada_def[region]['year']
                yi_emdat['GLB']['year'] = yi_emdat[region]['year']
                yi_climada['GLB']['all'] = yi_climada[region]['all']
                yi_climada_Q3['GLB']['all'] = yi_climada_Q3[region]['all']
                yi_climada_Q1['GLB']['all'] = yi_climada_Q1[region]['all']
                yi_climada_RMSF['GLB']['all'] = yi_climada_RMSF[region]['all']
                yi_climada_tot['GLB']['all'] = yi_climada_tot[region]['all']
                yi_climada_def['GLB']['all'] = yi_climada_def[region]['all']
                yi_emdat['GLB']['all'] = yi_emdat[region]['all']
                yi_emdat['GLB']['all scaled'] = yi_emdat[region]['all scaled']
            else:          
                yi_climada['GLB']['all'] = yi_climada['GLB']['all']+yi_climada[region]['all']
                yi_climada_Q3['GLB']['all'] = yi_climada_Q3['GLB']['all']+yi_climada_Q3[region]['all']
                yi_climada_Q1['GLB']['all'] = yi_climada_Q1['GLB']['all']+yi_climada_Q1[region]['all']
                yi_climada_RMSF['GLB']['all'] = yi_climada_RMSF['GLB']['all']+yi_climada_RMSF[region]['all']
                yi_climada_tot['GLB']['all'] = yi_climada_tot['GLB']['all']+yi_climada_tot[region]['all']
                yi_climada_def['GLB']['all'] = yi_climada_def['GLB']['all']+yi_climada_def[region]['all']
                yi_emdat['GLB']['all'] = yi_emdat['GLB']['all']+yi_emdat[region]['all'] 
                yi_emdat['GLB']['all scaled'] = yi_emdat['GLB']['all scaled']+yi_emdat[region]['all scaled']

        region = 'GLB'
        yi_climada[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada', region, 'csv')))
        yi_climada_Q3[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_Q3', region, 'csv')))
        yi_climada_Q1[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_Q1', region, 'csv')))
        yi_climada_RMSF[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_RMSF', region, 'csv')))
        yi_climada_tot[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_tot', region, 'csv')))
        yi_climada_def[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_climada_def', region, 'csv')))
        yi_emdat[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'yi_emdat', region, 'csv')))

        region = 'GLB_all'
        countries='GLB_all'
        exp_coast['if_TC'] = if_TC_all
        exp_coast['if_'] = if_TC_all
        yi_climada[region], stats_climada[region], stats_log_climada[region], _ = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.5)], \
                                      hazard_year_range)
        yi_climada_Q3[region], stats_climada_Q3[region], stats_log_climada_Q3[region], _ = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.75)], \
                                      hazard_year_range)
        yi_climada_Q1[region], stats_climada_Q1[region], stats_log_climada_Q1[region], _ = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['%.2f' % (.25)], \
                                      hazard_year_range)
        yi_climada_RMSF[region], stats_climada_RMSF[region], stats_log_climada_RMSF[region], aai_agg.loc[region, 'rmsf'] = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['RMSF'], \
                                      hazard_year_range)
        yi_climada_tot[region], stats_climada_tot[region], stats_log_climada_tot[region], aai_agg.loc[region, 'tot'] = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFSs['TDR'], \
                                      hazard_year_range)
        yi_climada_def[region], stats_climada_def[region], stats_log_climada_def[region], aai_agg.loc[region, 'def'] = \
            tc_cal.trend_by_country_CLIMADA(countries, exp_coast, IFS_default, \
                                      hazard_year_range)
        exp_coast['if_TC'] = if_TC_all
        exp_coast['if_'] = if_TC_all
        yi_emdat[region], stats_emdat[region], stats_log_emdat[region] = \
            tc_cal.trend_by_country_EMDAT(countries, year_range=YEAR_RANGE, reference_year=REF_YEAR, \
                            hazard_type_climada=HAZ)
        yi_climada[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada', region, 'csv')))
        yi_climada_Q3[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada_Q3', region, 'csv')))
        yi_climada_Q1[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada_Q1', region, 'csv')))
        yi_climada_RMSF[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada_RMSF', region, 'csv')))
        yi_climada_tot[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada_tot', region, 'csv')))
        yi_climada_def[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_climada_def', region, 'csv')))
        yi_emdat[region].to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                     REF_YEAR, HAZ, 'yi_emdat', region, 'csv')))

        aai_agg.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'aai_agg', '', 'csv')))

    # compute correlations per region:
    corr_factors = pd.DataFrame(index=yi_climada.keys(), columns=['Pearson', 'Pearson_ignore_0', 'Spearman', 'Spearman_ignore_0'])
    for region in yi_climada.keys():
        if 'GLB_all' in yi_climada[region].columns:
            yi_climada[region]['all'] = yi_climada[region]['GLB_all']
        if 'GLB_all' in yi_climada_Q3[region].columns:
            yi_climada_Q3[region]['all'] = yi_climada_Q3[region]['GLB_all']
        if 'GLB_all' in yi_climada_Q1[region].columns:
            yi_climada_Q1[region]['all'] = yi_climada_Q1[region]['GLB_all']
        if 'GLB_all' in yi_climada_RMSF[region].columns:
            yi_climada_RMSF[region]['all'] = yi_climada_RMSF[region]['GLB_all']
        if 'GLB_all' in yi_climada_tot[region].columns:
            yi_climada_tot[region]['all'] = yi_climada_tot[region]['GLB_all']
        if 'GLB_all' in yi_climada_def[region].columns:
            yi_climada_def[region]['all'] = yi_climada_def[region]['GLB_all']
        corr_factors.loc[region, 'Pearson'] = yi_climada[region]['all'].corr(yi_emdat[region]['all scaled'], method='pearson')
        corr_factors.loc[region, 'Pearson_ignore_0'] = yi_climada[region]['all'].replace(0, np.NaN).corr(yi_emdat[region]['all scaled'].replace(0, np.NaN), method='pearson')
        corr_factors.loc[region, 'Spearman'] = yi_climada[region]['all'].corr(yi_emdat[region]['all scaled'], method='spearman')
        corr_factors.loc[region, 'Spearman_ignore_0'] = yi_climada[region]['all'].replace(0, np.NaN).corr(yi_emdat[region]['all scaled'].replace(0, np.NaN), method='spearman')
    corr_factors.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'corr', 'factors', 'csv')))

    # compute annual average damage (AAD) with standard deviation (std):
    aad_std = tc_cal.aai_from_yi(yi_climada)
    aad_std['climada_median'] = tc_cal.aai_from_yi(yi_climada)['value']
    aad_std['climada_Q3'] = tc_cal.aai_from_yi(yi_climada_Q3)['value']
    aad_std['climada_Q1'] = tc_cal.aai_from_yi(yi_climada_Q1)['value']
    aad_std['climada_RMSF'] = tc_cal.aai_from_yi(yi_climada_RMSF)['value']
    aad_std['climada_tot'] = tc_cal.aai_from_yi(yi_climada_tot)['value']
    aad_std['climada_def'] = tc_cal.aai_from_yi(yi_climada_def)['value']
    aad_std['emdat'] = tc_cal.aai_from_yi(yi_emdat, column_name='all scaled')['value']
    aad_std['climada_median_std'] = tc_cal.aai_from_yi(yi_climada)['std']
    aad_std['climada_Q3_std'] = tc_cal.aai_from_yi(yi_climada_Q3)['std']
    aad_std['climada_Q1_std'] = tc_cal.aai_from_yi(yi_climada_Q1)['std']
    aad_std['climada_RMSF_std'] = tc_cal.aai_from_yi(yi_climada_RMSF)['std']
    aad_std['climada_tot_std'] = tc_cal.aai_from_yi(yi_climada_tot)['std']
    aad_std['climada_def_std'] = tc_cal.aai_from_yi(yi_climada_def)['std']
    aad_std['emdat_std'] = tc_cal.aai_from_yi(yi_emdat, column_name='all scaled')['std']
    aad_std = aad_std.drop(columns=['value', 'std'])

    aad_std.to_csv(os.path.join(RES_DIR, RES_STR % (RES_ARCSEC, \
                                 REF_YEAR, HAZ, 'AAD', 'std', 'csv')))

print('\n...........................................................\n')
print('All done. \nCheck out results directory (RES_DIR) for output: \n %s\n' %(RES_DIR))
print('.............................................................')
