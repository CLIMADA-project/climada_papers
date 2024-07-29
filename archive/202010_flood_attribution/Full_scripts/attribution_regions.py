#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 21:45:27 2019

@author: insauer
"""
import pandas as pd
import numpy as np
import statsmodels.api as sm
from pyts.decomposition import SingularSpectrumAnalysis
import pymannkendall as mk
from scipy import stats





def rel_time_attr_MK(dataFrame71):
    """
    Theil-Sen-Slope estimation and Mann-Kendall-Test to estimate the
    contribution of each driver!

    Parameters
    ----------
    dataFrame71 : time series
        Time series

    Returns
    -------
    regH : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 1980 fixed exposure)(TS_Haz) 1980-2010
    regHE : List MK-output
        Sen_slope and MK-test result with uncertainty range of TS_HazExp
        1980-2010
    regF : List MK-output
        Sen_slope and MK-test result with uncertainty range of TS_Full
        1980-2010.
    regH7 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 1980 fixed exposure)(TS_Haz) 1971-2010
    regH107 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 2010 fixed exposure)(TS_Haz) 1971-2010
    regH10 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 2010 fixed exposure)(TS_Haz) 1980-2010
    regE : List MK-output
        Sen_slope and MK-test result with uncertainty range of exposure
        difference function (TS_HazExp - TS_Haz) 1980-2010 (not used)
    regE7 : List MK-output
        Sen_slope and MK-test result with uncertainty range of exposure
        difference function (TS_HazExp - TS_Haz) 1971-2010 (not used)
    regV : List MK-output
        Sen_slope and MK-test result with uncertainty range of vulnerability
        difference function (TS_full - TS_Haz_Exp)(not used)
    regI : List MK-output
        Sen_slope and MK-test result with uncertainty range of modeled damges
        (including vulnerability)
    regN : List MK-output
        Sen_slope and MK-test result with uncertainty range of observed damages

    """

    dataFrame = dataFrame71[dataFrame71['Year'] > 1979]


    regLHazExp = mk.original_test(dataFrame['Norm_Impact_2y_trend'], alpha=0.1)

    slopeLHazExp = stats.theilslopes(dataFrame['Norm_Impact_2y_trend'],
                                     alpha=0.1)

    regHE = [regLHazExp.slope, regLHazExp.p, slopeLHazExp[2], slopeLHazExp[3]]

    regLFull = mk.original_test(dataFrame['Norm_Impact_Pred'], alpha=0.1)

    slopeLFull = stats.theilslopes(dataFrame['Norm_Impact_Pred'], alpha=0.1)

    regF = [regLFull.slope, regLFull.p, slopeLFull[2], slopeLFull[3]]

    regHaz = mk.original_test(dataFrame['Norm_ImpFix_2y_trend'], alpha=0.1)

    slopeHaz = stats.theilslopes(dataFrame['Norm_ImpFix_2y_trend'],
                                 alpha=0.1)

    regH = [regHaz.slope, regHaz.p, slopeHaz[2], slopeHaz[3]]

    regHaz7 = mk.original_test(dataFrame71['Norm_ImpFix_2y_trend'],
                               alpha=0.1)

    slopeHaz7 = stats.theilslopes(dataFrame71['Norm_ImpFix_2y_trend'],
                                  alpha=0.1)

    regH7 = [regHaz7.slope, regHaz7.p, slopeHaz7[2], slopeHaz7[3]]

    regHaz107 = mk.original_test(dataFrame71['Norm_Imp2010_2y_trend'],
                                 alpha=0.1)

    slopeHaz107 = stats.theilslopes(dataFrame71['Norm_Imp2010_2y_trend'],
                                    alpha=0.1)

    regH107 = [regHaz107.slope, regHaz107.p, slopeHaz107[2], slopeHaz107[3]]

    regHaz10 = mk.original_test(dataFrame['Norm_Imp2010_2y_trend'],
                                alpha=0.1)

    slopeHaz10 = stats.theilslopes(dataFrame['Norm_Imp2010_2y_trend'],
                                   alpha=0.1)

    regH10 = [regHaz10.slope, regHaz10.p, slopeHaz10[2], slopeHaz10[3]]

    regNat = mk.original_test(dataFrame['natcat_flood_damages_2005_CPI'], alpha=0.1)

    slopeNat = stats.theilslopes(dataFrame['natcat_flood_damages_2005_CPI'], alpha=0.1)

    regN = [regNat.slope, regNat.p, slopeNat[2], slopeNat[3]]

    return regH, regHE, regH7, regH107, regH10, regF, regN



def normalise(dataFrame):
    """
    Normalisation of time series. Normalisation for trends estimation
    normalises everything to total observed damage. Normalisation for plotting
    shifts time series to the same starting point in 1980.
    ----------
    region : string
        Abbrevation of region
    dataFrame : DataFrame
        Time series

    Returns
    -------
    dataFrame : DataFrame
        Time series + normalised time series

    """

    obs_adj = dataFrame.loc[dataFrame['Year'] > 1979,
                            'natcat_flood_damages_2005_CPI'].mean() / \
        dataFrame.loc[dataFrame['Year'] > 1979, 'Impact_Pred'].mean()

    dataFrame['Norm_Impact_Pred'] = dataFrame['Impact_Pred']*obs_adj

    offsetExp = dataFrame.loc[dataFrame['Year'] == 1980,
                              'Norm_Impact_Pred'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980,
                      'Impact_2y_Flopros'].sum()
    offsetHaz = dataFrame.loc[dataFrame['Year'] == 1980,
                              'Norm_Impact_Pred'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980, 'ImpFix_2y_Flopros'].sum()

    # normalisation for plotting

    dataFrame['Norm_Impact_2y_offset'] = dataFrame['Impact_2y_Flopros'] * offsetExp
    dataFrame['Norm_ImpFix_2y_offset'] = dataFrame['ImpFix_2y_Flopros'] * offsetHaz

    dataFrame['Norm_Imp2010_2y_offset'] = dataFrame['Imp2010_2y_Flopros'] * offsetHaz

    # modelspread for plotting

    dataFrame['Norm_Impact_Pred_1thrd_offset'] = \
        dataFrame['Impact_Pred_1thrd']*obs_adj
    dataFrame['Norm_Impact_Pred_2thrd_offset'] = \
        dataFrame['Impact_Pred_2thrd']*obs_adj

    dataFrame['Norm_Impact_2y_onethird_offset'] = \
        dataFrame['model_flood_damages_onethird_quantile']*offsetExp
    dataFrame['Norm_Impact_2y_twothird_offset'] = \
        dataFrame['model_flood_damages_twothird_quantile']*offsetExp

    dataFrame['Norm_ImpFix_2y_onethird_offset'] = \
        dataFrame['model_flood_damages_onethird_quantile_1980flospros']*offsetHaz
    dataFrame['Norm_ImpFix_2y_twothird_offset'] = \
        dataFrame['model_flood_damages_twothird_quantile_1980flospros']*offsetHaz

    dataFrame['Norm_Imp2010_2y_onethird_offset'] =\
        dataFrame['model_flood_damages_onethird_quantile_2010flospros'] * offsetHaz
    dataFrame['Norm_Imp2010_2y_twothird_offset'] =\
        dataFrame['model_flood_damages_twothird_quantile_2010flospros'] * offsetHaz

    # normalisation for trend estimation

    trendNormExp = dataFrame.loc[dataFrame['Year'] >= 1980,
                                 'Norm_Impact_Pred'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Impact_2y_Flopros'].mean()
    trendNormHaz = dataFrame.loc[dataFrame['Year'] >= 1980,
                                 'Norm_Impact_Pred'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'ImpFix_2y_Flopros'].mean()
    trendNormHaz2010 = dataFrame.loc[dataFrame['Year'] >= 1980,
                                     'Norm_Impact_Pred'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Imp2010_2y_Flopros'].mean()

    dataFrame['Norm_Impact_2y_trend'] = dataFrame['Impact_2y_Flopros'] * trendNormExp
    dataFrame['Norm_ImpFix_2y_trend'] = dataFrame['ImpFix_2y_Flopros'] * trendNormHaz

    dataFrame['Norm_Imp2010_2y_trend'] = dataFrame['Imp2010_2y_Flopros'] * trendNormHaz2010

    return dataFrame


def prep_table_timeMK(region, dat, regLHaz, regLHazExp, regH7, regH107,
                      regH10, regI, regN):
    """
    Prepare output table for attribution done with Mann-Kendall and Theil-Sen slope

    Parameters
    ----------
    region : string
        Abbrevation of region
    dat : DataFrame
        Time series
    regMeta : DataFrame
        Fitting info (not used)
    regLHaz : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 1980 fixed exposure)(TS_Haz) 1980-2010
    regLHazExp :  List MK-output
        Sen_slope and MK-test result with uncertainty range of TS_HazExp
        1980-2010
    regLFull : List MK-output
        Sen_slope and MK-test result with uncertainty range of TS_Full
        1980-2010..
    regH7 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 1980 fixed exposure)(TS_Haz) 1971-2010
    regH107 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 2010 fixed exposure)(TS_Haz) 1971-2010
    regH10 : List MK-output
        Sen_slope and MK-test result with uncertainty range of hazard
        (with 2010 fixed exposure)(TS_Haz) 1980-2010
    regE : List MK-output
        Sen_slope and MK-test result with uncertainty range of exposure
        difference function (TS_HazExp - TS_Haz) 1980-2010 (not used)
    regE7 : List MK-output
        Sen_slope and MK-test result with uncertainty range of exposure
        difference function (TS_HazExp - TS_Haz) 1971-2010 (not used)
    regV : List MK-output
        Sen_slope and MK-test result with uncertainty range of vulnerability
        difference function (TS_full - TS_Haz_Exp)(not used)
    regI : List MK-output
        Sen_slope and MK-test result with uncertainty range of modeled damges
        (including vulnerability)
    regN : List MK-output
        Sen_slope and MK-test result with uncertainty range of observed damages

    Returns
    -------
    table1 : DataFrame
        final output for one region

    """

    nat_param = np.nanmean(dat.loc[(dat['Year'] > 1979) &
                        (dat['Year'] < 1996), 'natcat_flood_damages_2005_CPI'])
    # nat_2010 =  dat.loc[(dat['Year']==2010),'natcat_flood_damages_2005_CPI'].mean()

    # norm_param = regMeta.loc[regMeta['Region'] == region,'Predicted_damages'].sum()/31.

    cH7_norm = regH7[0]*100/nat_param  # dam_param

    cH_norm = regLHaz[0]*100/nat_param  # dam_param

    cH7_normCL = regH7[0]*100/nat_param  # dam_param

    cH7up_normCL = regH7[3]*100/nat_param - cH7_normCL  # dam_param

    cH7bot_normCL = cH7_normCL - regH7[2]*100/nat_param

    cH_normCL = regLHaz[0]*100/nat_param  # dam_param

    cHup_normCL = regLHaz[3]*100/nat_param - cH_normCL

    cHbot_normCL = cH_normCL - regLHaz[2]*100/nat_param

    cH107_normCL = regH107[0]*100/nat_param  # dam_param

    cH107up_normCL = regH107[3]*100/nat_param - cH107_normCL

    cH107bot_normCL = (cH107_normCL - regH107[2]*100/nat_param)

    cH10_normCL = regH10[0]*100/nat_param  # dam_param

    cH10up_normCL = regH10[3]*100/nat_param - cH10_normCL

    cH10bot_normCL = (cH10_normCL - regH10[2]*100/nat_param)

    cE_norm = (regLHazExp[0]-regLHaz[0])*100/nat_param  # dam_param
     
    cE10_norm = (regLHazExp[0]-regH10[0]) * 100 / nat_param

    cV_norm = (regI[0]-regLHazExp[0])*100/nat_param

    cI_norm = regI[0]*100/nat_param


    cN_norm = regN[0]*100/nat_param

    table1 = pd.DataFrame({'Region': region,
                           'Change H': regLHaz[0],  # develop_taylor(regH, 1995)[0],
                           'Change Hn': cH_norm,
                           'Change HnCL': cH_normCL,
                           'Change HnCLup': cHup_normCL,
                           'Change HnCLbot': cHbot_normCL,
                           'Change H10': regH10[0],  # develop_taylor(regH10, 1995)[0],
                           'Change H7': regH7[0],  # develop_taylor(regH7, 1991)[0],
                           'Change H10nCL': cH10_normCL,
                           'Change H10nCLup': cH10up_normCL,
                           'Change H10nCLbot': cH10bot_normCL,
                           'Change H107nCL': cH107_normCL,
                           'Change H107nCLup': cH107up_normCL,
                           'Change H107nCLbot': cH107bot_normCL,
                           'Change H7n': cH7_norm,
                           'Change H7nCL': cH7_normCL,
                           'Change H7nCLup': cH7up_normCL,
                           'Change H7nCLbot': cH7bot_normCL,
                           'Change H107': regH107[0],  
                           'Change En': cE_norm,
                           'Change En10': cE10_norm,
                           'Change Vn': cV_norm,  # develop_taylor(regV, 1995)[0],
                           'Change I': regI[0],  # develop_taylor(regI, 1995)[0],
                           'Change In': cI_norm,  # develop_taylor(regI, 1995)[0],
                           'Change N': regN[0],  # develop_taylor(regN, 1995)[0],
                           'Change Nn': cN_norm,
                           'Sign H': regLHaz[1],
                           'Sign H7': regH7[1],
                           'Sign H10': regH10[1],
                           'Sign H107': regH107[1],
                           'Sign I': regI[1],
                           'Sign N': regN[1],
                           '2010_haz_loss80': regH10[0]*31,
                           '2010_haz_loss71': regH107[0]*40,
                           '2010_haz_loss80_rel': regH10[0]*31/nat_param,
                           '2010_haz_loss71_rel': regH107[0]*40/nat_param,
                           '2010_haz_loss80_rel10': regH10[0]/regN[0],
                           '2010_haz_loss71_rel10': regH107[0]/regN[0],
                   }, index =[0])
    return table1



def attr_regr(dataFrame):
    """
    Wrapper function that coordinates attribution and prepares data output
    for each region.

    Parameters
    ----------
    dataFrame : DataFrame
        Time series
    metaData : DataFrame
        Vulnerability fit info
    normalisation : string
        parameter used for normalistion (not really used)
    regr : string
        regression method used for trend estimation

    Returns
    -------
    normData : DataFrame
        normalised time series
    attrTable : DataFrame
        attribution output

    """

    attrTable = pd.DataFrame()
    normData = pd.DataFrame()
    for i, test_region in enumerate(test_regions):

        DATA_region = dataFrame[(dataFrame['Region'] == test_region) &
                                (dataFrame['Year'] < 2011) &
                                (dataFrame['Year'] > 1970)]

        DATA_region = DATA_region.reset_index()

        DATA_region = normalise(DATA_region)

        regLHaz, regLHazExp, regH7, regH107, regH10, regI, regN = rel_time_attr_MK(DATA_region)
    
        attrReg = prep_table_timeMK(test_region, DATA_region,
                                    regLHaz, regLHazExp, regH7, regH107,
                                    regH10, regI, regN)

        attrTable = attrTable.append(attrReg, ignore_index=True)
        normData = normData.append(DATA_region, ignore_index=True)
    return normData, attrTable


DATA = pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/VulnerabilityAdjustmentTimeSeriesRegions.csv')

region_names = {
               'NAM': 'North America',
               'AUS': 'Oceania',
               'LAM': 'Central America',
               'EUR': 'Western Europe',
               'NAFARA': 'North Africa + Middle East',
               'SSAF': 'SSA + Southern Africa',
               'CAS': 'Central Asia + Eastern Europe',
               'SWEA': 'Southern Asia + South-East Asia',
               'CHN': 'Eastern Asia',
               'GLB': 'Global',
                }
test_regions = list(region_names)

regr = 'MK'


norm_names = {'Predicted_damages': 'RelPred',
              'Observed_damages': 'RelObs',
              'Absolute_damages': 'Abs'}

normData, attrTable = attr_regr(DATA)

attrTable.to_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AttributionMetaDataRegions.csv', index=False)

normData.to_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AttributionTimeSeriesRegions.csv', index=False)