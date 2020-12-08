# -*- coding: utf-8 -*-
"""
Spyder Editor

This file aggregates multi-model flood damage output on country level 
to regional model medians. Additional variables such as GDP, Pop, Capital Stock,
GDP_pc and recorded damages are added.
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
    # Aggreagted damages is given as the sum over country-wise damages
    aggregated_model_damages = x['Impact_2y_Flopros'].sum()
    aggregated_model_damages0 = x['Impact_2y_0'].sum()
    aggregated_model_damages_1980flospros = x['ImpFix_2y_Flopros'].sum()
    aggregated_model_damages_1980_0 = x['ImpFix_2y_0'].sum()
    aggregated_model_damages_2010flospros = x['Imp2010_2y_Flopros'].sum()
    aggregated_model_damages_2010_0 = x['Imp2010_2y_0'].sum()
    aggregated_observed_damages = (x['natcat_flood_damages_2005_CPI']).sum()
    # Use the population-weighted GDP per cap
    aggregated_flooded_area = x['FloodedAreaFlopros'].sum()
    aggregated_flooded_area0 = x['FloodedArea0'].sum()
    # aggregated_flooded_vol = x['FloodVol_Flopros'].sum()
    return pd.Series([aggregated_model_damages, aggregated_model_damages0,
                      aggregated_model_damages_1980flospros,
                      aggregated_model_damages_1980_0,
                      aggregated_model_damages_2010flospros,
                      aggregated_model_damages_2010_0,
                      aggregated_flooded_area,
                      aggregated_flooded_area0, aggregated_observed_damages],
                     index=['Impact_2y_Flopros', 'Impact_2y_0',
                            'ImpFix_2y_Flopros', 'ImpFix_2y_0',
                            'Imp2010_2y_Flopros', 'Imp2010_2y_0',
                            'FloodedAreaFlopros', 'FloodedArea0',
                            'natcat_flood_damages_2005_CPI'])

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
    median_model_damages = x['Impact_2y_Flopros'].median()  # =quantile(0.5)
    median_model_damages0 = x['Impact_2y_0'].median()
    median_model_damages_1980flospros = x['ImpFix_2y_Flopros'].median()
    # =quantile(0.5)
    median_model_damages_1980_0 = x['ImpFix_2y_0'].median()
    median_model_damages_2010flospros = x['Imp2010_2y_Flopros'].median()
    # =quantile(0.5)
    median_model_damages_2010_0 = x['Imp2010_2y_0'].median()
    median_observed_damages = (x['natcat_flood_damages_2005_CPI']).mean()

    one_third_quantile_flood_area = x['FloodedAreaFlopros'].quantile(0.3)  # 30
    two_third_quantile_flood_area = x['FloodedAreaFlopros'].quantile(0.7)  # 70
    one_third_quantile_model_damages = x['Impact_2y_Flopros'].quantile(0.3)
    # 30
    one_third_quantile_model_damages_1980flospros = x['ImpFix_2y_Flopros'].quantile(0.3)  # 30
    one_third_quantile_model_damages_2010flospros = x['Imp2010_2y_Flopros'].quantile(0.3)
    two_third_quantile_model_damages = x['Impact_2y_Flopros'].quantile(0.7)  # 70
    two_third_quantile_model_damages_1980flospros = x['ImpFix_2y_Flopros'].quantile(0.7)  # 70
    two_third_quantile_model_damages_2010flospros = x['Imp2010_2y_Flopros'].quantile(0.7)  # 70
    flood_area_flospros = x['FloodedAreaFlopros'].median()
    flood_area_0 = x['FloodedArea0'].median()
    # flood_vol = x['FloodVol_Flopros'].median()
    return pd.Series([median_model_damages,
                      median_model_damages0,
                      median_model_damages_1980flospros,
                      median_model_damages_1980_0,
                      median_model_damages_2010flospros,
                      median_model_damages_2010_0,
                      median_observed_damages,
                      one_third_quantile_flood_area,
                      two_third_quantile_flood_area,
                      one_third_quantile_model_damages,
                      one_third_quantile_model_damages_1980flospros,
                      one_third_quantile_model_damages_2010flospros,
                      two_third_quantile_model_damages,
                      two_third_quantile_model_damages_1980flospros,
                      two_third_quantile_model_damages_2010flospros,
                      flood_area_flospros,
                      flood_area_0],
                     index=['Impact_2y_Flopros',
                            'Impact_2y_0',
                            'ImpFix_2y_Flopros',
                            'ImpFix_2y_0',
                            'Imp2010_2y_Flopros',
                            'Imp2010_2y_0',
                            'natcat_flood_damages_2005_CPI',
                            'flood_area_onethird_quantile',
                            'flood_area_twothird_quantile',
                            'model_flood_damages_onethird_quantile', 
                            'model_flood_damages_onethird_quantile_1980flopros',
                            'model_flood_damages_onethird_quantile_2010flopros',
                            'model_flood_damages_twothird_quantile',
                            'model_flood_damages_twothird_quantile_1980flopros',
                            'model_flood_damages_twothird_quantile_2010flopros',
                            'FloodedAreaFlopros',
                            'FloodedArea0'
                            #  'FloodedArea0'
                            ])


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

    # provide dataset with observational data
    natcat = pd.read_excel('/home/insauer/projects/Attribution/Floods/Paper_NC_Review_Data/Input_PPP_conversion/1980-2016_Masterdatei_NatCat_worldwide_no-pw_2005conversions_PPP.xlsx', index_col=0)
    countries = pd.read_csv('/home/insauer/projects/NC_Submission/Data/final_country_list.csv')
    
    # asset rescaling to correct for the ssp transition and to convert GDP to capital stock
    # datasets can be generated with a separate code discribed in the readme.txt
    if gdp_resc:
        resc_factors = pd.read_csv('/home/insauer/projects/Attribution/Data/RescalingFactors_GDPobs_GDPjpnClean.csv')
        # cap_factors = pd.read_csv('/home/insauer/projects/Attribution/Data/PostProcess_GDP2CapStock_cgdpoClean.csv')
        cap_factors = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/Input_Data/PostProcess_GDP2CapStock_rgdpnaCleanRM.csv')

    countries = list(set(megaDataFrame['Country']).intersection(countries.iloc[:, 0]))
    natcat_short = natcat.iloc[:,[4,7,8,9,13,12,11,21,20,25,29,33,-17,-3,-2,-1]] # selection of NatCat columns, containing all rows

    cols = ['year', 'event', 'type', 'subtype', 'ISO', 'country', 'Continent',
            'tot_loss', 'ins_loss', 'tot_loss_GDP', 'tot_loss_GCP',
            'tot_loss_CPI', 'Fatalities', 'CPI_con', 'GDP_con', 'Reg_name']
    natcat_short.columns = cols
    natcat_hydro = natcat_short.loc[natcat_short['subtype'] == 'gf:General flood']
    megaDataFrame['natcat_flood_damages_2005_CPI'] = np.nan


    for country in countries:
        
        resc_fac_cnt_yr_1980 = resc_factors.loc[resc_factors['ISO'] == country, str(1980)].sum()
        capst_fac_cnt_yr_1980 = cap_factors.loc[cap_factors['ISO'] == country, str(1980)].sum()
        resc_fac_cnt_yr_2010 = resc_factors.loc[resc_factors['ISO'] == country, str(2010)].sum()
        capst_fac_cnt_yr_2010 = cap_factors.loc[cap_factors['ISO'] == country, str(2010)].sum()


        for year in years:
            print(str(year) + ' ' + country)
            resc_fac_cnt_yr = resc_factors.loc[resc_factors['ISO'] == country, str(year)].sum()
            capst_fac_cnt_yr = cap_factors.loc[cap_factors['ISO'] == country, str(year)].sum()

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Impact_2y_Flopros'] *= resc_fac_cnt_yr * capst_fac_cnt_yr
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Impact_2y_0'] *= resc_fac_cnt_yr * capst_fac_cnt_yr
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'ImpFix_2y_Flopros'] *= resc_fac_cnt_yr_1980 * capst_fac_cnt_yr_1980
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'ImpFix_2y_0'] *= resc_fac_cnt_yr_1980 * capst_fac_cnt_yr_1980

            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Imp2010_2y_Flopros'] *= resc_fac_cnt_yr_2010 * capst_fac_cnt_yr_2010
            megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                              (megaDataFrame['Year'] == year),
                              'Imp2010_2y_0'] *= resc_fac_cnt_yr_2010 * capst_fac_cnt_yr_2010

            if year > 1979:
                tmp_natcat_05_cpi = natcat_hydro.loc[(natcat_hydro['year']==year) & (natcat_hydro['ISO']==country), 'CPI_con'].sum()

                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_flood_damages_2005_CPI'] = tmp_natcat_05_cpi*1.0e6

            else:
                megaDataFrame.loc[(megaDataFrame['Country'] == country) &
                                  (megaDataFrame['Year'] == year),
                                  'natcat_flood_damages_2005_CPI'] = np.nan

    return megaDataFrame


def region_aggregation(cols, dataFrame):
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
    data_region = dataFrame.groupby(['Year', 'GHM', 'clim_forc', 'Region'])\
                                    [cols].apply(aggregation_regions)\
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

        dataFrame = dataFrame[dataFrame['GHM'] == select_model]

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

        [n1, n2, n3, n4, ghm, clim_forc, n7, n8, n9] = model_output.split('_') 
#        if cl_model == 'watch':
#            continue

        print('Loading ' + model_output)
        temp = pd.read_csv(path+model_output)
        temp['GHM'] = ghm
        temp['clim_forc'] = clim_forc
        temp = temp[temp['Year'] >= 1971]
        megaDataFrame = megaDataFrame.append(temp, ignore_index=True)

    megaDataFrame = megaDataFrame.sort_values(by=sort)

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
#  Cluster path where data from damage assessment calculated with CLIMADA
#  is stored
path = '/home/insauer/mnt/ebm/inga/paper_results/paper_resubmission_1_12/'
sort = ['Year', 'Country']
years = np.arange(1971, 2012)

# Columns used for region aggregation
in_cols = ['Impact_2y_Flopros',
           'Impact_2y_0',
           'ImpFix_2y_Flopros',
           'ImpFix_2y_0',
           'Imp2010_2y_Flopros',
           'Imp2010_2y_0',
           'FloodedAreaFlopros',
           'FloodedArea0',
           'natcat_flood_damages_2005_CPI']  

# groupby year and model...




#Building one big data set with all the data from all model runs

#assDataFrame = assemble_data_frame(path, sort, years)
# assDataFrame.rename(columns={"Impact_2yFlopros": "Impact_2y_Flopros",
#                               "Impact_2y0": "Impact_2y_0",
#                               'ImpFix_2yFlopros': 'ImpFix_2y_Flopros',
#                               'Imp2010_2yFlopros': 'Imp2010_2y_Flopros',
#                               'ImpFix_2y0': 'ImpFix_2y_0',
#                               'Imp2010_2y0': 'Imp2010_2y_0'}, inplace=True)
# # adding additional data
#assDataFrame = add_GDP_NatCat(assDataFrame, years, gdp_resc = True)
# # storing giant dataframe containing all the data for all model runs
#assDataFrame.to_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AssembledDataRegions.csv', index=False)

assDataFrame = pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AssembledDataRegions.csv')

#  rename regions
assDataFrame = aggregate_new_region(assDataFrame)

# aggregate all the country based data to regions
regDataFrame = region_aggregation(in_cols, assDataFrame)
#regDataFrame.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/RegionalAggregationDataRegions.csv', index=False)

#  building model median
modDataFrame = model_aggregation(in_cols, regDataFrame, years, None)

modDataFrame.to_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/ModelMedianRegions.csv', index=False)
