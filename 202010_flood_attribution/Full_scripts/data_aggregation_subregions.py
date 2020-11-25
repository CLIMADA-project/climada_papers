#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This file aggregates multi-model flood damage output on country level
to regional model medians, taking into account subregional discharge trends.
Additional variables such as GDP, Pop, Capital Stock, GDP_pc and recorded
damages are added.
"""

import numpy as np
import pandas as pd
import os
from astropy.convolution import convolve


def runmean(data, halfwin):
    """
    Simple running mean.
    CAUTION: Data is *extended* at the edges by repeating the
    edge values; thereby any trend present in the data will
    become attenuated at the edges!
    """
    window = 2*halfwin + 1
    if window > len(data):
        print('Error: window too large!')
        import sys
        sys.exit(0)
    weights = np.repeat(1.0, window) / window
    # Treat edges: Extend data
    extended_data = np.hstack([[data[0]] * (halfwin), data, [data[len(data)-1]]
                               * (halfwin)])
    # rm = np.convolve(extended_data, weights, 'valid')
    rm = convolve(extended_data, weights, boundary=None, nan_treatment='fill',
                  preserve_nan=False)
    return rm[halfwin:-halfwin]


def aggregation_regions(x):
    """
    This function aggregates country-level damages and variables to
    regional level.
    Parameters
    ----------
    x : DataFrame
        country-level damages and other indicators for all model combinations

    Returns
    -------
    DataFrame
        regionally aggregated damages and other indicators
    """
    aggregated_model_damages_pos = x['Impact_2yPosFlopros'].sum()
    aggregated_model_damages_neg = x['Impact_2yNegFlopros'].sum()
    aggregated_model_damages_1980_pos = x['ImpFix_2yPosFlopros'].sum()
    aggregated_model_damages_1980_neg = x['ImpFix_2yNegFlopros'].sum()
    aggregated_model_damages_2010_pos = x['Imp2010_2yPosFlopros'].sum()
    aggregated_model_damages_2010_neg = x['Imp2010_2yNegFlopros'].sum()
    aggregated_observed_damages_pos = (x['natcat_damages_2005_CPI_pos']).sum()
    aggregated_observed_damages_neg = (x['natcat_damages_2005_CPI_neg']).sum()
    # Use the population-weighted GDP per cap
    aggregated_GDP_per_cap_pos = (x['gdp_Pos_RM'].sum())/(x['pop_Pos_RM'].sum())
    aggregated_GDP_per_cap_neg = (x['gdp_Neg_RM'].sum())/(x['pop_Neg_RM'].sum())
    aggregated_population_pos = x['pop_Pos_RM'].sum()
    aggregated_population_neg = x['pop_Neg_RM'].sum()
    aggregated_capstock_pos = x['CapStock_Pos'].sum()
    aggregated_capstock_neg = x['CapStock_Neg'].sum()
    aggregated_flooded_area_pos = x['FloodedAreaPosFlopros'].sum()
    aggregated_flooded_area_neg = x['FloodedAreaNegFlopros'].sum()
    aggregated_flood_vol_pos = x['FloodVolumePosFlopros'].sum()
    aggregated_flood_vol_neg = x['FloodVolumePosFlopros'].sum()
    # aggregated_flooded_vol = x['FloodVol_Flopros'].sum()
    return pd.Series([aggregated_model_damages_pos,
                      aggregated_model_damages_neg,
                      aggregated_model_damages_1980_pos,
                      aggregated_model_damages_1980_neg,
                      aggregated_model_damages_2010_pos,
                      aggregated_model_damages_2010_neg,
                      aggregated_observed_damages_pos,
                      aggregated_observed_damages_neg,
                      aggregated_GDP_per_cap_pos,
                      aggregated_GDP_per_cap_neg,
                      aggregated_population_pos, aggregated_population_neg,
                      aggregated_flooded_area_pos, aggregated_flooded_area_neg,
                      aggregated_flood_vol_pos, aggregated_flood_vol_neg,
                      aggregated_capstock_pos, aggregated_capstock_neg],
                     index=['Impact_2yPosFlopros', 'Impact_2yNegFlopros',
                            'ImpFix_2yPosFlopros', 'ImpFix_2yNegFlopros',
                            'Imp2010_2yPosFlopros', 'Imp2010_2yNegFlopros',
                            'natcat_damages_2005_CPI_pos',
                            'natcat_damages_2005_CPI_neg',
                            'GDP_pcPos', 'GDP_pcNeg', 'PopPos', 'PopNeg',
                            'FloodedAreaPosFlopros', 'FloodedAreaNegFlopros',
                            'FloodVolumePosFlopros', 'FloodVolumeNegFlopros',
                            'CapStock_Pos', 'CapStock_Neg'])


def func_median(x):
    """
    This function aggregates the damages and other indicators from the
    different model runs to the model median and adds basic statistics such as
    the one-third and two-third quantiles.

    Parameters
    ----------
    x : DataFrame
        regionally aggregated damages and other indicators for all model
        combinations

    Returns
    -------
    DataFrame
         model medians of regionally aggregated damages and other indicators
    """
    # identify the median of the model data:
    median_model_damages_pos = x['Impact_2yPosFlopros'].median()  # =quantile(0.5)
    median_model_damages_neg = x['Impact_2yNegFlopros'].median()
    median_model_damages_1980_pos = x['ImpFix_2yPosFlopros'].median()  # =quantile(0.5)
    median_model_damages_1980_neg = x['ImpFix_2yNegFlopros'].median()
    median_model_damages_2010_pos = x['Imp2010_2yPosFlopros'].median()  # =quantile(0.5)
    median_model_damages_2010_neg = x['Imp2010_2yNegFlopros'].median()
    median_observed_damages_pos = (x['natcat_damages_2005_CPI_pos']).mean()  # all the same value
    median_observed_damages_neg = (x['natcat_damages_2005_CPI_neg']).mean()  # all the same value
    flood_area_pos = x['FloodedAreaPosFlopros'].median()
    flood_area_neg = x['FloodedAreaNegFlopros'].median()
    flood_vol_pos = x['FloodVolumePosFlopros'].median()
    flood_vol_neg = x['FloodVolumeNegFlopros'].median()

    median_capstock_pos = x['CapStock_Pos'].mean()
    median_capstock_neg = x['CapStock_Neg'].mean()

    median_GDP_per_cap_pos = x['GDP_pcPos'].mean()  # all the same value
    median_population_pos = x['PopPos'].mean()  # all the same value
    median_GDP_per_cap_neg = x['GDP_pcNeg'].mean()  # all the same value
    median_population_neg = x['PopNeg'].mean()

    one_third_quantile_flood_area_pos = x['FloodedAreaPosFlopros'].quantile(0.3)  # 30
    two_third_quantile_flood_area_pos = x['FloodedAreaPosFlopros'].quantile(0.7)  # 70
    one_third_quantile_flood_area_neg = x['FloodedAreaNegFlopros'].quantile(0.3)  # 30
    two_third_quantile_flood_area_neg = x['FloodedAreaNegFlopros'].quantile(0.7)  # 70

    one_third_quantile_flood_vol_pos = x['FloodVolumePosFlopros'].quantile(0.3)  # 30
    two_third_quantile_flood_vol_pos = x['FloodVolumePosFlopros'].quantile(0.7)  # 70
    one_third_quantile_flood_vol_neg = x['FloodVolumeNegFlopros'].quantile(0.3)  # 30
    two_third_quantile_flood_vol_neg = x['FloodVolumeNegFlopros'].quantile(0.7)  # 70

    one_third_quantile_model_damages_pos = x['Impact_2yPosFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_pos = x['Impact_2yPosFlopros'].quantile(0.7)
    one_third_quantile_model_damages_neg = x['Impact_2yNegFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_neg = x['Impact_2yNegFlopros'].quantile(0.7)
    one_third_quantile_model_damages_1980_pos = x['ImpFix_2yPosFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_1980_pos = x['ImpFix_2yPosFlopros'].quantile(0.7)  # 70
    one_third_quantile_model_damages_1980_neg = x['ImpFix_2yNegFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_1980_neg = x['ImpFix_2yNegFlopros'].quantile(0.7)  # 70

    one_third_quantile_model_damages_2010_pos = x['Imp2010_2yPosFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_2010_pos = x['Imp2010_2yPosFlopros'].quantile(0.7)  # 70
    one_third_quantile_model_damages_2010_neg = x['Imp2010_2yNegFlopros'].quantile(0.3)  # 30
    two_third_quantile_model_damages_2010_neg = x['Imp2010_2yNegFlopros'].quantile(0.7)  # 70

    # flood_vol = x['FloodVol_Flopros'].median()
    return pd.Series([median_model_damages_pos,
                      median_model_damages_neg,
                      median_model_damages_1980_pos,
                      median_model_damages_1980_neg,
                      median_model_damages_2010_pos,
                      median_model_damages_2010_neg,
                      median_observed_damages_pos,
                      median_observed_damages_neg,
                      flood_area_pos,
                      flood_area_neg,
                      flood_vol_pos,
                      flood_vol_neg,
                      median_GDP_per_cap_pos,
                      median_population_pos,
                      median_GDP_per_cap_neg,
                      median_population_neg,
                      median_capstock_pos,
                      median_capstock_neg,
                      one_third_quantile_flood_area_pos,
                      two_third_quantile_flood_area_pos,
                      one_third_quantile_flood_area_neg,
                      two_third_quantile_flood_area_neg,
                      one_third_quantile_flood_vol_pos,
                      two_third_quantile_flood_vol_pos,
                      one_third_quantile_flood_vol_neg,
                      two_third_quantile_flood_vol_neg,
                      one_third_quantile_model_damages_pos,
                      two_third_quantile_model_damages_pos,
                      one_third_quantile_model_damages_neg,
                      two_third_quantile_model_damages_neg,
                      one_third_quantile_model_damages_1980_pos,
                      two_third_quantile_model_damages_1980_pos,
                      one_third_quantile_model_damages_1980_neg,
                      two_third_quantile_model_damages_1980_neg,
                      one_third_quantile_model_damages_2010_pos,
                      two_third_quantile_model_damages_2010_pos,
                      one_third_quantile_model_damages_2010_neg,
                      two_third_quantile_model_damages_2010_neg],
                     index=['Impact_2yPos',
                            'Impact_2yNeg',
                            'ImpFix_2yPos',
                            'ImpFix_2yNeg',
                            'Imp2010_2yPos',
                            'Imp2010_2yNeg',
                            'natcat_damages_2005_CPI_Pos',
                            'natcat_damages_2005_CPI_Neg',
                            'FloodedAreaPos',
                            'FloodedAreaNeg',
                            'FloodVolumePos',
                            'FloodVolumeNeg',
                            'GDP_pcPos',
                            'PopPos',
                            'GDP_pcNeg',
                            'PopNeg',
                            'CapStockPos',
                            'CapStockNeg',
                            'flood_area_onethird_quantile_Pos',
                            'flood_area_twothird_quantile_Pos',
                            'flood_area_onethird_quantile_Neg',
                            'flood_area_twothird_quantile_Neg',
                            'flood_vol_onethird_quantile_Pos',
                            'flood_vol_twothird_quantile_Pos',
                            'flood_vol_onethird_quantile_Neg',
                            'flood_vol_twothird_quantile_Neg',
                            'Impact_2yPos_onethird_quantile',
                            'Impact_2yPos_twothird_quantile',
                            'Impact_2yNeg_onethird_quantile',
                            'Impact_2yNeg_twothird_quantile',
                            'ImpFix_2yPos_onethird_quantile',
                            'ImpFix_2yPos_twothird_quantile',
                            'ImpFix_2yNeg_onethird_quantile',
                            'ImpFix_2yNeg_twothird_quantile',
                            'Imp2010_2yPos_onethird_quantile',
                            'Imp2010_2yPos_twothird_quantile',
                            'Imp2010_2yNeg_onethird_quantile',
                            'Imp2010_2yNeg_twothird_quantile'])


def add_GDP_NatCat(megaDataFrame, years, gdp_resc):
    """
    This function inserts national annual variables in the MegaDataFrame
    containing all the data. Damages are all converted to capital stock by
    applying a corresponding annual national conversion factor.
    Inserted are:
        GDP (not relevant for final paper)
        GDP 10 yr running mean (not relevant for final paper)
        GDP per capita (not relevant for final paper)
        Population (not relevant for final paper)
        Capital Stock (not relevant for final paper)
        GMT (not relevant for final paper)
        recorded damages (NatCat Munich Re)

    Parameters
    ----------
    megaDataFrame : DataFrame
        big data set containing all the data from all model runs
    years : int array
        years to be considered
    gdp_resc : bool
        gdp to capital stock conversion

    Returns
    -------
    DataFrame
         model medians of regionally aggregated damages and other indicators
    """
    # gdp_pc  and population information is available at 
    # Geiger, Tobias; Frieler, Katja (2017): Continuous national
    # Gross Domestic Product (GDP) time series for 195 countries:
    # past observations (1850-2005) harmonized with future projections according
    # the Shared Socio-economic Pathways (2006-2100). GFZ Data Services.
    # https://doi.org/10.5880/pik.2017.003
    # adapt path
    gdp_pc = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/Income-PPP2005_ISIMIP_merged_Maddison-pwt81_1850-2015_extended_WDI-1.csv')
    pop = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/Population_ISIMIP-pwt81_1850-2009_extended_WDI_filled-with-Hyde.csv')
    gdpop_posneg = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/GDP_pc_PosNeg.csv')
    # provide dataset with observational data assigned to each subregion...this needs to be processed previously
    natcat = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Review_Data/NatCat_PPP_conv/NatCatAssigned_PPP_bas_geo_corr_test3.csv')
    # natcat=pd.read_excel('/home/insauer/projects/Attribution/Floods/Data/Input_Data/1980-2016_Masterdatei_NatCat_worldwide_no-pw_2005conversions.xlsx', index_col=0)
    # asset rescaling to correct for the ssp transition and to convert GDP to capital stock
    # datasets can be generated with a separate code discribed in the readme.txt
    if gdp_resc:
        resc_factors = pd.read_csv('/home/insauer/projects/Attribution/Data/RescalingFactors_GDPobs_GDPjpnClean.csv')
        # cap_factors = pd.read_csv('/home/insauer/projects/Attribution/Data/PostProcess_GDP2CapStock_cgdpoClean.csv')
        cap_factors = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/PostProcess_GDP2CapStock_rgdpnaCleanRM.csv')
        conv_cap_stock = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/GDP2CapStock_rgdpnaClean.csv')
    countries = list(set(megaDataFrame['Country']).intersection(gdp_pc.iloc[:, 0]))

    megaDataFrame['GDP_pc'] = np.nan
    megaDataFrame['GDP'] = np.nan
    megaDataFrame['Pop'] = np.nan
    megaDataFrame['GDP_10y_runmean'] = np.nan
    megaDataFrame['gdp_Pos_RM'] = np.nan

    megaDataFrame['gdp_Neg_RM'] = np.nan
    megaDataFrame['pop_Pos_RM'] = np.nan
    megaDataFrame['pop_Neg_RM'] = np.nan

    megaDataFrame['natcat_damages_2005_CPI_pos'] = np.nan
    megaDataFrame['natcat_damages_2005_CPI_neg'] = np.nan
    megaDataFrame['GMT'] = np.nan

    gdp_transpose = gdp_pc.T
    gdp_transpose.columns = gdp_transpose.iloc[0]
    gdp_transpose.drop(gdp_transpose.index[0], inplace=True)
    gdp_runmean = gdp_transpose.rolling(window=10, min_periods=5).mean()

    # skip_cntr = ['KGZ','TJK']
    for country in countries:
        for year in years:
            print(str(year) + ' ' + country)
            
            gdp_pc_cnt_yr = gdp_pc.loc[gdp_pc.iloc[:, 0] == country, str(year)]
            resc_fac_cnt_yr = resc_factors.loc[resc_factors['ISO'] == country, str(year)].sum()
            capst_fac_cnt_yr = cap_factors.loc[cap_factors['ISO'] == country, str(year)].sum()

            capst_cnt_yr = conv_cap_stock.loc[conv_cap_stock['ISO'] == country, str(year)].sum()

            resc_fac_cnt_yr_1980 = resc_factors.loc[resc_factors['ISO'] == country, str(1980)].sum()
            capst_fac_cnt_yr_1980 = cap_factors.loc[cap_factors['ISO'] == country, str(1980)].sum()

            resc_fac_cnt_yr_2010 = resc_factors.loc[resc_factors['ISO'] == country, str(2010)].sum()
            capst_fac_cnt_yr_2010 = cap_factors.loc[cap_factors['ISO'] == country, str(2010)].sum()

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Impact_2yPosFlopros'] *= resc_fac_cnt_yr * capst_fac_cnt_yr
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Impact_2yNegFlopros'] *= resc_fac_cnt_yr * capst_fac_cnt_yr
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'ImpFix_2yPosFlopros'] *= resc_fac_cnt_yr_1980 * capst_fac_cnt_yr_1980
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'ImpFix_2yNegFlopros'] *= resc_fac_cnt_yr_1980 * capst_fac_cnt_yr_1980

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Imp2010_2yPosFlopros'] *= resc_fac_cnt_yr_2010 * capst_fac_cnt_yr_2010
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Imp2010_2yNegFlopros'] *= resc_fac_cnt_yr_2010 * capst_fac_cnt_yr_2010

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'GDP_pc'] = np.float(gdp_pc_cnt_yr)
            pop_cnt_yr = pop.loc[gdp_pc.iloc[:, 0] == country, str(year)]

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Pop'] = np.float(pop_cnt_yr)
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'GDP'] = np.float(pop_cnt_yr) * np.float(gdp_pc_cnt_yr)

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'GDP_10y_runmean'] = gdp_runmean.loc[str(year), country]

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'gdp_Pos_RM'] = gdpop_posneg.loc[(gdpop_posneg['Country'] == country) &
                                                               (gdpop_posneg['Year']==year),
                                                               'GDP_RM_Pos'].mean()

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'CapStock_Pos'] = gdpop_posneg.loc[(gdpop_posneg['Country'] == country) &
                                                                 (gdpop_posneg['Year'] == year),
                                                                 'GDP_RM_Pos'].mean() * np.float(capst_cnt_yr)

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'gdp_Neg_RM'] = gdpop_posneg.loc[(gdpop_posneg['Country'] == country) &
                                                               (gdpop_posneg['Year'] == year),
                                                               'GDP_RM_Neg'].mean()

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'CapStock_Neg'] = gdpop_posneg.loc[(gdpop_posneg['Country'] == country) &
                                                                 (gdpop_posneg['Year'] == year),
                                                                 'GDP_RM_Neg'].mean() * np.float(capst_cnt_yr)

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'pop_Pos_RM'] = gdpop_posneg.loc[(gdpop_posneg['Country']==country) &
                                                               (gdpop_posneg['Year']==year),
                                                               'POP_RM_Pos'].mean()
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'pop_Neg_RM'] = gdpop_posneg.loc[(gdpop_posneg['Country'] == country) &
                                                               (gdpop_posneg['Year'] == year),
                                                               'POP_RM_Neg'].mean()

            if year > 1979:

                natcat_dam_pos = natcat.loc[(natcat['Country'] == country) &
                                            (natcat['Year'] == year), 'pos_risk'].sum()

                natcat_dam_neg = natcat.loc[(natcat['Country'] == country) &
                                            (natcat['Year'] == year), 'neg_risk'].sum()

                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_damages_2005_CPI_pos'] = natcat_dam_pos*1.0e6

                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_damages_2005_CPI_neg'] = natcat_dam_neg*1.0e6

            else:
                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_damages_2005_CPI_pos'] = np.nan
                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_damages_2005_GDP_neg'] = np.nan

    return megaDataFrame


def aggregate_new_region(dataFrame):
    """
    This function combines world regions to bigger regions

    Parameters
    ----------
    dataFrame : DataFrame
        DataFrame containig smaller regions

    Returns
    -------
    DataFrame
         DataFrame with manipulated regions
    """
    
    dataFrame.loc[dataFrame['Country'] == 'RUS',
                  'Region'] = 'CAS'

    dataFrame.loc[(dataFrame['Region'] == 'CAR') |
                  (dataFrame['Region'] == 'LAS') |
                  (dataFrame['Region'] == 'LAN'),
                  'Region'] = 'LAM'

    dataFrame.loc[(dataFrame['Region'] == 'NAF') |
                  (dataFrame['Region'] == 'ARA'), 'Region'] = 'NAFARA'

    dataFrame.loc[(dataFrame['Region'] == 'SSA') |
                  (dataFrame['Region'] == 'SAF'), 'Region'] = 'SSAF'

    dataFrame.loc[(dataFrame['Region'] == 'EUR') |
                  (dataFrame['Region'] == 'EUA'), 'Region'] = 'EUR'

    dataFrame.loc[(dataFrame['Region'] == 'SWA') |
                  (dataFrame['Region'] == 'SEA'),
                  'Region'] = 'SWEA'

    dataFrame.loc[(dataFrame['Region'] == 'PIS1') |
                  (dataFrame['Region'] == 'PIS2') |
                  (dataFrame['Region'] == 'AUS'),
                  'Region'] = 'AUS'

    return dataFrame


def region_aggregation(out_cols, dataFrame):
    """
    This function is a wrapper for the aggregation and selects the columns to
    be aggregated to regional level.

    Parameters
    ----------
    out_cols : string list
        Columns to be aggregated

    Returns
    -------
    DataFrame
         regionally aggregated damages and other indicators regions
    """
    data_region = dataFrame.groupby(['Year', 'rf_model', 'cl_model', 'Region'])\
                                     [out_cols].apply(aggregation_regions)\
                                     .reset_index()  # groupby year and model

    return data_region


def model_aggregation(out_cols, dataFrame, years, select_model):
    """
    This function is a wrapper for the multi-model aggregation and provides
    the model median for each region of all variables.

    Parameters
    ----------
    out_cols : string list
        Columns to be aggregated

    Returns
    -------
    DataFrame
         regionally aggregated model medians
    """

    if select_model:

        dataFrame = dataFrame[dataFrame['rf_model'] == select_model]

    data_models = dataFrame[(dataFrame['Year'] <= np.max(years)) &
                            (dataFrame['Year'] >= np.min(years))]
    # Get the median for model and datasets
    data_models = data_models.groupby(['Year', 'Region'])\
                              [out_cols].apply(func_median).reset_index()

    return data_models

def assemble_data_frame(path, sort, years):
    """
    This function gathers all the data from all model_runs and
    provides one big data set containing all model runs

    Parameters
    ----------
    path : string
        Path to directory on cluster where all runs are stored
    sort : string list
        columns used for sorting
    years : int array

    Returns
    -------
    DataFrame
         full data
    """
    megaDataFrame = pd.DataFrame()
    list_of_model_output = os.listdir(path)

    for i, model_output in enumerate(list_of_model_output):  # loop over all model output 

        [n1, n2, n3, n4, rf_model, cl_model, n7, n8, n9] = model_output.split('_')

        print('Loading ' + model_output)
        temp = pd.read_csv(path+model_output)
        temp['rf_model'] = rf_model
        temp['cl_model'] = cl_model
        temp = temp[temp['Year'] >= 1971]
        megaDataFrame = megaDataFrame.append(temp, ignore_index=True)

    megaDataFrame = megaDataFrame.sort_values(by=sort)

    return megaDataFrame

def get_rm(col_names, dataFrame, rm_window):
    """
    This function calculates running means for selected columns.
    Not applied in the paper!

    Parameters
    ----------
    col_names : string list
        selected columns
    dataFrame : DataFrame
        full dataset
    rm_window : int
         window size

    Returns
    -------
    DataFrame
         DataFrame containing the running means for selected columns
    """
    for col in col_names:
        dataFrame[col] = dataFrame[col].replace(
                                            [-np.inf, np.inf],
                                            [np.nan, np.nan])
        dataFrame[col] = runmean(dataFrame[col], rm_window)

    return dataFrame

path = '/home/insauer/mnt/ebm/results/floodDam/Paper_Inga_NC_20_06_20_review/dam_basin_geo/'
sort = ['Year', 'Country']
years = np.arange(1971, 2012)

in_cols = ['Impact_2yPosFlopros',
           'Impact_2yNegFlopros',
           'ImpFix_2yPosFlopros',
           'ImpFix_2yNegFlopros',
           'Imp2010_2yPosFlopros',
           'Imp2010_2yNegFlopros',
           'FloodedAreaPosFlopros',
           'FloodedAreaNegFlopros',
           'FloodVolumePosFlopros',
           'FloodVolumeNegFlopros',
           'natcat_damages_2005_CPI_pos',
           'natcat_damages_2005_CPI_neg',
           'gdp_Pos_RM',
           'gdp_Neg_RM',
           'pop_Pos_RM',
           'pop_Neg_RM',
           'CapStock_Pos',
           'CapStock_Neg']

# groupby year and model...

out_cols = ['Impact_2yPosFlopros', 'Impact_2yNegFlopros',
            'ImpFix_2yPosFlopros', 'ImpFix_2yNegFlopros',
            'Imp2010_2yPosFlopros', 'Imp2010_2yNegFlopros',
            'natcat_damages_2005_CPI_pos', 'natcat_damages_2005_CPI_neg',
            'GDP_pcPos', 'GDP_pcNeg', 'PopPos',
            'PopNeg', 'FloodedAreaPosFlopros',
            'FloodedAreaNegFlopros', 'FloodVolumePosFlopros',
            'FloodVolumeNegFlopros', 'CapStock_Pos', 'CapStock_Neg']


# Put all different model outputs in one dataframe
# and add GDP and NatCat Damages for each year and country
assDataFrame = assemble_data_frame(path, sort, years)
assDataFrame = add_GDP_NatCat(assDataFrame, years, gdp_resc=True)
assDataFrame.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/AssembledDataSubregions.csv', index=False)
assDataFrame = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/AssembledDataSubregions.csv')

#  rename regions
assDataFrame = aggregate_new_region(assDataFrame)

# aggregate all the country based data to regions
regDataFrame = region_aggregation(in_cols, assDataFrame)
#regDataFrame.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/RegionalAggregationDataSubregions.csv', index=False)

#  building model median
modDataFrame = model_aggregation(out_cols, regDataFrame, years, None)
modDataFrame.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/ModelMedianSubregions.csv', index=False)
