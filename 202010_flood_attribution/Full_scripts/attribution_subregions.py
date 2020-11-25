#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 12:32:07 2020

@author: insauer
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from pyts.decomposition import SingularSpectrumAnalysis
import pymannkendall as mk
from scipy import stats


def rel_time_attr_MK(dataFrame71, disch):
    """
    Theil-Sen-Slope estimation and Mann-Kendall-Test to estimate the
    contribution of each driver!

    Parameters
    ----------
    dataFrame71 : time series
        Time series.
    disch : string
        discharge group

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

    regLHazExp = mk.original_test(dataFrame['NormExp_Impact_2y{}_trend'.format(disch)],
                                  alpha=0.1)
    slopeLHazExp = stats.theilslopes(dataFrame['NormExp_Impact_2y{}_trend'.format(disch)],
                                     alpha=1/3)

    regHE = [regLHazExp.slope, regLHazExp.p, slopeLHazExp[2], slopeLHazExp[3]]

    regLFull = mk.original_test(dataFrame['Norm_Impact_Pred_{}'.format(disch)],
                                alpha=0.1)

    slopeLFull = stats.theilslopes(dataFrame['Norm_Impact_Pred_{}'.format(disch)],
                                   alpha=1/3)

    regF = [regLFull.slope, regLFull.p, slopeLFull[2], slopeLFull[3]]

    regHaz = mk.original_test(dataFrame['NormHaz_ImpFix_2y{}_trend'.format(disch)],
                              alpha=0.1)

    slopeHaz = stats.theilslopes(dataFrame['NormHaz_ImpFix_2y{}_trend'.format(disch)],
                                 alpha=1/3)

    regH = [regHaz.slope, regHaz.p, slopeHaz[2], slopeHaz[3]]

    regHaz7 = mk.original_test(dataFrame71['NormHaz_ImpFix_2y{}_trend'.format(disch)], alpha=0.1)

    slopeHaz7 = stats.theilslopes(dataFrame71['NormHaz_ImpFix_2y{}_trend'.format(disch)],
                                  alpha=1/3)

    regH7 = [regHaz7.slope, regHaz7.p, slopeHaz7[2], slopeHaz7[3]]

    regHaz107 = mk.original_test(dataFrame71['NormHaz_Imp2010_2y{}_trend'.format(disch)],
                                 alpha=0.1)

    slopeHaz107 = stats.theilslopes(dataFrame71['NormHaz_Imp2010_2y{}_trend'.format(disch)],
                                    alpha=1/3)

    regH107 = [regHaz107.slope, regHaz107.p, slopeHaz107[2], slopeHaz107[3]]

    regHaz10 = mk.original_test(dataFrame['NormHaz_Imp2010_2y{}_trend'.format(disch)],
                                alpha=0.1)

    slopeHaz10 = stats.theilslopes(dataFrame['NormHaz_Imp2010_2y{}_trend'.format(disch)],
                                   alpha=1/3)

    regH10 = [regHaz10.slope, regHaz10.p, slopeHaz10[2], slopeHaz10[3]]

    regE = mk.original_test(dataFrame['NormExp_Impact_2y{}_trend'.format(disch)], alpha=0.1)
    regE7 = mk.original_test(dataFrame['NormExp_Impact_2y{}_trend'.format(disch)], alpha=0.1)

    regV = mk.original_test(dataFrame['Norm_Impact_Pred_{}'.format(disch)], alpha=0.1)

    regI = mk.original_test(dataFrame['Norm_Impact_Pred_{}'.format(disch)],
                            alpha=0.1)
    regI = regF

    regNat = mk.original_test(dataFrame['natcat_damages_2005_CPI_{}'.format(disch)],
                              alpha=0.1)

    slopeNat = stats.theilslopes(dataFrame['natcat_damages_2005_CPI_{}'.format(disch)],
                                 alpha=1/3)

    regN = [regNat.slope, regNat.p, slopeNat[2], slopeNat[3]]

    return regH, regHE, regF, regH7, regH107, regH10, regE, regE7, regV, regI, regN


def normalise(region, dataFrame):
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

    obs_adj_neg = dataFrame.loc[dataFrame['Year'] > 1979, 'natcat_damages_2005_CPI_Neg'].mean() / \
        dataFrame.loc[dataFrame['Year'] > 1979, 'Impact_Pred_Neg'].mean()
    obs_adj_pos = dataFrame.loc[dataFrame['Year'] > 1979, 'natcat_damages_2005_CPI_Pos'].mean() / \
        dataFrame.loc[dataFrame['Year'] > 1979, 'Impact_Pred_Pos'].mean()

    dataFrame['Norm_Impact_Pred_Neg'] = dataFrame['Impact_Pred_Neg'] * obs_adj_neg
    dataFrame['Norm_Impact_Pred_Pos'] = dataFrame['Impact_Pred_Pos'] * obs_adj_pos

    # normalisation for plotting
    offsetExpNeg = dataFrame.loc[dataFrame['Year'] == 1980, 'Norm_Impact_Pred_Neg'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980, 'Impact_2yNeg'].sum()
    dataFrame['NormExp_Impact_2yNeg_offset'] = dataFrame['Impact_2yNeg'] * offsetExpNeg

    offsetExpPos = dataFrame.loc[dataFrame['Year'] == 1980, 'Norm_Impact_Pred_Pos'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980, 'Impact_2yPos'].sum()
    dataFrame['NormExp_Impact_2yPos_offset'] = dataFrame['Impact_2yPos'] * offsetExpPos

    offsetHazNeg = dataFrame.loc[dataFrame['Year'] == 1980,
                                 'Norm_Impact_Pred_Neg'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980, 'ImpFix_2yNeg'].sum()
    dataFrame['NormHaz_ImpFix_2yNeg_offset'] = dataFrame['ImpFix_2yNeg'] * offsetHazNeg

    offsetHazPos = dataFrame.loc[dataFrame['Year'] == 1980, 'Norm_Impact_Pred_Pos'].sum() / \
        dataFrame.loc[dataFrame['Year'] == 1980, 'ImpFix_2yPos'].sum()
    dataFrame['NormHaz_ImpFix_2yPos_offset'] = dataFrame['ImpFix_2yPos'] * offsetHazPos

    dataFrame['NormHaz_Imp2010_2yPos_offset'] = dataFrame['Imp2010_2yPos']*offsetExpPos
    dataFrame['NormHaz_Imp2010_2yNeg_offset'] = dataFrame['Imp2010_2yNeg']*offsetExpNeg

    # modelspread for plotting

    dataFrame['Norm_Impact_Pred_1thrd_Neg'] = dataFrame['Impact_Pred_1thrd_Neg'] * obs_adj_neg

    dataFrame['Norm_Impact_Pred_1thrd_Pos'] = dataFrame['Impact_Pred_1thrd_Pos'] * obs_adj_pos

    dataFrame['Norm_Impact_Pred_2thrd_Neg'] = dataFrame['Impact_Pred_2thrd_Neg'] * obs_adj_neg

    dataFrame['Norm_Impact_Pred_2thrd_Pos'] = dataFrame['Impact_Pred_2thrd_Pos'] * obs_adj_pos

    dataFrame['Norm_Impact_2yNeg_onethird_quantile'] = \
        dataFrame['Impact_2yNeg_onethird_quantile'] * offsetExpNeg

    dataFrame['Norm_Impact_2yPos_onethird_quantile'] = \
        dataFrame['Impact_2yPos_onethird_quantile'] * offsetExpPos

    dataFrame['Norm_Impact_2yNeg_twothird_quantile'] = \
        dataFrame['Impact_2yNeg_twothird_quantile'] * offsetExpNeg

    dataFrame['Norm_Impact_2yPos_twothird_quantile'] = \
        dataFrame['Impact_2yPos_twothird_quantile'] * offsetExpPos

    dataFrame['Norm_ImpFix_2yNeg_onethird_quantile'] = \
        dataFrame['ImpFix_2yNeg_onethird_quantile'] * offsetExpNeg

    dataFrame['Norm_ImpFix_2yPos_onethird_quantile'] = \
        dataFrame['ImpFix_2yPos_onethird_quantile'] * offsetExpPos

    dataFrame['Norm_ImpFix_2yNeg_twothird_quantile'] = \
        dataFrame['ImpFix_2yNeg_twothird_quantile'] * offsetExpNeg

    dataFrame['Norm_ImpFix_2yPos_twothird_quantile'] = \
        dataFrame['ImpFix_2yPos_twothird_quantile'] * offsetExpPos

    dataFrame['Norm_Imp2010_2yNeg_onethird_quantile'] = \
        dataFrame['Imp2010_2yNeg_onethird_quantile'] * offsetExpNeg

    dataFrame['Norm_Imp2010_2yPos_onethird_quantile'] = \
        dataFrame['Imp2010_2yPos_onethird_quantile'] * offsetExpPos

    dataFrame['Norm_Imp2010_2yNeg_twothird_quantile'] = \
        dataFrame['Imp2010_2yNeg_twothird_quantile'] * offsetExpNeg

    dataFrame['Norm_Imp2010_2yPos_twothird_quantile'] = \
        dataFrame['Imp2010_2yPos_twothird_quantile'] * offsetExpPos

    # normalisation for trend estimation

    trendNormExpNeg = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Neg'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Impact_2yNeg'].mean()
    dataFrame['NormExp_Impact_2yNeg_trend'] = dataFrame['Impact_2yNeg'] * trendNormExpNeg

    trendNormExpPos = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Pos'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Impact_2yPos'].mean()
    dataFrame['NormExp_Impact_2yPos_trend'] = dataFrame['Impact_2yPos'] * trendNormExpPos

    trendNormHazNeg = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Neg'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'ImpFix_2yNeg'].mean()
    dataFrame['NormHaz_ImpFix_2yNeg_trend'] = dataFrame['ImpFix_2yNeg'] * trendNormHazNeg

    trendNormHazPos = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Pos'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'ImpFix_2yPos'].mean()

    dataFrame['NormHaz_ImpFix_2yPos_trend'] = dataFrame['ImpFix_2yPos'] * trendNormHazPos

    trendNormHaz10Neg = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Neg'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Imp2010_2yNeg'].mean()
    dataFrame['NormHaz_Imp2010_2yNeg_trend'] = dataFrame['Imp2010_2yNeg'] * trendNormHaz10Neg

    trendNormHaz10Pos = dataFrame.loc[dataFrame['Year'] >= 1980, 'Norm_Impact_Pred_Pos'].mean() / \
        dataFrame.loc[dataFrame['Year'] >= 1980, 'Imp2010_2yPos'].mean()
    dataFrame['NormHaz_Imp2010_2yPos_trend'] = dataFrame['Imp2010_2yPos'] * trendNormHaz10Pos

    return dataFrame


def prep_table_timeMK(region, data, regLHaz, regLHazExp, regLFull, regH7,
                      regH107, regH10, regE, regE7, regV, regI, regN, dis):
    """
    Prepare output table for attribution done with Mann-Kendall and Theil-Sen slope

    Parameters
    ----------
    region : string
        Abbrevation of region
    dat : DataFrame
        Time series
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
    dis : string
        discharge group

    Returns
    -------
    table1 : DataFrame
        final output for one region

    """

    dam_norm = np.nanmean(data.loc[(data['Year'] > 1979) & (data['Year'] < 1991),
                        'Impact_Pred_{}'.format(dis)])
    nat_norm = np.nanmean(data.loc[(data['Year'] > 1979) & (data['Year'] < 1996),
                        'natcat_damages_2005_CPI_{}'.format(dis)])

    nat_2010 = data.loc[(data['Year'] == 2010), 'natcat_damages_2005_CPI_{}'.format(dis)].mean()

    cH_norm = regLHaz[0]*100/nat_norm

    cH_bot = (cH_norm - regLHaz[2]*100/nat_norm)

    cH_up = (regLHaz[3]*100/nat_norm) - cH_norm

    cH7_normCL = regH7[0]*100/nat_norm

    cH7up_normCL = (regH7[3]*100/nat_norm) - cH7_normCL

    cH7bot_normCL = cH7_normCL - regH7[2]*100/nat_norm

    cH_normCL = regLHaz[0]*100/nat_norm

    cH10_normCL = regH10[0]*100/nat_norm

    cH10up_normCL = (regH10[3]*100/nat_norm)-cH10_normCL

    cH10bot_normCL = (cH10_normCL - regH10[2]*100/nat_norm)

    cH107_normCL = regH107[0]*100/nat_norm

    cH107up_normCL = (regH107[3]*100/nat_norm) - cH107_normCL

    cH107bot_normCL = (cH107_normCL - regH107[2]*100/nat_norm)

    # cH_norm = (regH10.params[0]*100)/dam_norm
    cH7_norm = regH7[0]  # *100/nat_norm
    cE_norm = (regLHazExp[0]-regLHaz[0])*100/nat_norm
    
    cE10_norm = (regLHazExp[0]-regH10[0])*100/nat_norm
    cE7_norm = regE7.slope*100/dam_norm
    cV_norm = (regLFull[0]-regLHazExp[0])*100/nat_norm

    cI_norm = regLFull[0]*100/nat_norm
    cN_norm = regN[0]*100/nat_norm


    table1 = pd.DataFrame({'Region': region+'_'+dis,
                           'Change H': regLHaz[0],  # develop_taylor(regH, 1995)[0],
                           'Change Hn': cH_norm,
                           'Change HnCLup': cH_up,
                           'Change HnCLbot': cH_bot,
                           'Change HnCL': cH_normCL,
                           'Change H10': regH10[0],  # develop_taylor(regH10, 1995)[0],
                           'Change H10nCL': cH10_normCL,
                           'Change H10nCLup': cH10up_normCL,
                           'Change H10nCLbot': cH10bot_normCL,
                           'Change H7': regH7[0],  # develop_taylor(regH7, 1991)[0]
                           'Change H7n': cH7_norm,
                           'Change H7nCLup': cH7up_normCL,
                           'Change H7nCLbot': cH7bot_normCL,
                           'Change H7nCL': cH7_normCL,
                           'Change H107': regH107[0],  # develop_taylor(regH107, 1991)[0],
                           'Change H107nCLup': cH107up_normCL,
                           'Change H107nCLbot': cH107bot_normCL,
                           'Change H107nCL': cH107_normCL,
                           'Change E7': regE7.slope,  # develop_taylor(regE7, 1991)[0],
                           'Change E': regE.slope,  # develop_taylor(regE, 1995)[0],
                           'Change En': cE_norm,
                           'Change En10': cE10_norm,
                           'Change V': regV.slope,  # develop_taylor(regV, 1995)[0],
                           'Change Vn': cV_norm,
                           'Change I': regI[0],  # develop_taylor(regI, 1995)[0],
                           'Change In': cI_norm,
                           'Change N': regN[0],  # develop_taylor(regN, 1995)[0],
                           'Change Nn': cN_norm,
                           'Sign H': regLHaz[1],
                           'Sign H10': regH10[1],
                           'Sign H7': regH7[1],
                           'Sign H107': regH107[1],
                           'Sign E': regE.p,
                           'Sign E7': regE7.p,
                           'Sign V': regV.p,
                           'Sign I': regI[1],
                           'Sign N': regN[1],
                           '2010_haz_loss80': regH10[0]*31,
                           '2010_haz_loss71': regH107[0]*40,
                           '2010_haz_loss80_rel': regH10[0]*31/nat_norm,
                           '2010_haz_loss71_rel': regH107[0]*40/nat_norm,
                           '2010_haz_loss80_rel10': regH10[0]*31/nat_2010,
                           '2010_haz_loss71_rel10': regH107[0]*40/nat_2010
                           },
                          index=[0])
    return table1


def attr_regr(dataFrame):

    attrTable = pd.DataFrame()
    normData = pd.DataFrame()

    for i, test_region in enumerate(test_regions):

        DATA_region = dataFrame[(dataFrame['Region'] == test_region) &
                                (dataFrame['Year'] < 2011) & (dataFrame['Year'] > 1970)]
        DATA_region = DATA_region.reset_index()

        DATA_region = normalise(test_region, DATA_region)

        regLHazPos, regLHazExpPos, regLFullPos, regH7Pos, regH107Pos,\
            regH10Pos, regEPos, regE7Pos, regVPos, regIPos, regNPos =\
            rel_time_attr_MK(DATA_region, 'Pos')
        regLHazNeg, regLHazExpNeg, regLFullNeg, regH7Neg, regH107Neg,\
            regH10Neg, regENeg, regE7Neg, regVNeg, regINeg, regNNeg =\
            rel_time_attr_MK(DATA_region, 'Neg')
        attrRegPos = prep_table_timeMK(test_region, DATA_region, regLHazPos,
                                       regLHazExpPos, regLFullPos, regH7Pos,
                                       regH107Pos, regH10Pos, regEPos, regE7Pos,
                                       regVPos, regIPos, regNPos, 'Pos')
        attrRegNeg = prep_table_timeMK(test_region, DATA_region, regLHazNeg,
                                       regLHazExpNeg, regLFullNeg, regH7Neg,
                                       regH107Neg, regH10Neg, regENeg, regE7Neg,
                                       regVNeg, regINeg, regNNeg, 'Neg')

        attrTable = attrTable.append(attrRegPos, ignore_index=True)

        attrTable = attrTable.append(attrRegNeg, ignore_index=True)
        normData = normData.append(DATA_region, ignore_index=True)

    return normData, attrTable


DATA = pd.read_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/VulnerabilityAdjustmentTimeSeriesSubregions.csv')


region_names = {'NAM': 'North America',
                'LAM': 'Latin America',
                'EUR': 'Western Europe',
                'NAFARA': 'North Africa + Middle East',
                'SSAF': 'SSA+ Southern Africa',
                'CAS': 'Central Asia + Eastern Europe',
                'SWEA': 'Southern Asia + South-East Asia',
                'CHN': 'Eastern Asia',
                'AUS': 'Oceania',
                'GLB': 'Global'}

test_regions = list(region_names)


normalisation = 'Predicted_damages'
norm_names = {'Predicted_damages': 'RelPred',
              'Observed_damages': 'RelObs',
              'Absolute_damages': 'Abs'}

regr = 'MK'


normData, attrTable = attr_regr(DATA)

attrTable.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/AttributionMetaDataSubregions.csv', index=False)

normData.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/AttributionTimeSeriesSubregions.csv', index=False)