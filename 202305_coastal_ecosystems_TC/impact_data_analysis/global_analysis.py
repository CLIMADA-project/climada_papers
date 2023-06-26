"""
Script to analyse impact data on a global level
Author: Sarah Hülsen
"""

from pathlib import Path
import pandas as pd
import copy
import numpy as np

# dataset characteristics
scenario = 'hist_STORM'
pop_years = ['2000',
             '2020']
pop_data = 'wp'
rhab_years = ['1992',
              '2020']

path = Path('../results/final/')
rhab_path = Path('../results/intermediate/')

# Merging population impact data (per exposure points) for all basins
for pop_yr in pop_years:
    # load basin impact files
    AP_path = f'{path}AP_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    IO_path = f'{path}IO_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    WP_path = f'{path}WP_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    SH_path = f'{path}SH_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'

    AP_df = pd.read_csv(AP_path)
    IO_df = pd.read_csv(IO_path)
    WP_df = pd.read_csv(WP_path)
    SH_df = pd.read_csv(SH_path)

    # merge basin data
    dfs = [AP_df, IO_df, WP_df, SH_df]
    df_all = pd.concat(dfs)

    # calculate impact across all windspeed categories
    df_all['eai_exp_total'] = df_all['eai_exp_cat1'] \
                              + df_all['eai_exp_cat2']\
                              + df_all['eai_exp_cat3'] \
                              + df_all['eai_exp_cat4'] \
                              + df_all['eai_exp_cat5']

    # save merged dataframe
    df_all.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')


# Merging population impact per exposure point with file containing protective habitat information
for rhab_yr in rhab_years:
    for pop_yr in pop_years:
        # load population data
        pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')
        pop = pop.reset_index(drop=True)
        pop['longitude'] = pop['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
        pop['latitude'] = pop['latitude'].round(decimals=11)

        # load habitat data
        rhab = pd.read_csv(f'{rhab_path}{pop_data}{pop_yr}_rhab{rhab_yr}.csv')
        rhab = rhab.reset_index(drop=True)
        rhab['longitude'] = rhab['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
        rhab['latitude'] = rhab['latitude'].round(decimals=11)

        # joining population impact data and habitat data
        join = pd.merge(pop, rhab, how='left', on=['longitude', 'latitude'])
        cols = ['longitude',
                'latitude',
                'eai_exp_cat1',
                'eai_exp_cat2',
                'eai_exp_cat3',
                'eai_exp_cat4',
                'eai_exp_cat5',
                'eai_exp_total',
                'rhab']
        join = join[cols]
        join['rhab'] = join['rhab'].replace(np.nan, 5)  # assigning rhab value 5 to all unmatched population exposures

        # save file containing impacts for all exposure points
        join.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}_incl_5_rhab.csv')
        # save file containing only impacts for habitat values below 5
        join = join[join.rhab != 5]
        join.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv')

        # Calculate aggregated impact data:
        # total & relative number of people impacted per wind speed category
        # total & relative number of people protected per wind speed category

        # read population impact data of all basins per exposure point (calculated above)
        pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')

        # summing population impacts per category
        pop_cat1 = pop['eai_exp_cat1'].sum()
        pop_cat2 = pop['eai_exp_cat2'].sum()
        pop_cat3 = pop['eai_exp_cat3'].sum()
        pop_cat4 = pop['eai_exp_cat4'].sum()
        pop_cat5 = pop['eai_exp_cat5'].sum()
        pop_total = pop['eai_exp_total'].sum()

        # creating dataframe of sums
        pop_sums = pd.DataFrame([[pop_cat1, pop_cat2, pop_cat3, pop_cat4, pop_cat5, pop_total]],
                                columns=['eai_exp_cat1',
                                         'eai_exp_cat2',
                                         'eai_exp_cat3',
                                         'eai_exp_cat4',
                                         'eai_exp_cat5',
                                         'eai_exp_total'
                                         ])

        # saving file with population impacts per wind speed category
        pop_sums_rd = pop_sums.round()
        pop_sums_rd.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat.csv')

        # read file containing exposure point-wise information on habitats
        file = f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
        df = pd.read_csv(file)

        # select only necessary columns
        col_select = ['rhab',
                      'eai_exp_cat1',
                      'eai_exp_cat2',
                      'eai_exp_cat3',
                      'eai_exp_cat4',
                      'eai_exp_cat5',
                      'eai_exp_total'
                      ]
        df = df[col_select]

        # group  and sum impacts by habitat value
        df_sum = df.groupby(['rhab'], as_index=False)[['eai_exp_cat1',
                                                       'eai_exp_cat2',
                                                       'eai_exp_cat3',
                                                       'eai_exp_cat4',
                                                       'eai_exp_cat5',
                                                       'eai_exp_total']].sum()
        # save file containing numbers of people protected per windspeed category
        df_sum_rd = df_sum.round()
        df_sum_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat.csv')

        # relative impact (share of people protected)
        df_perc = copy.copy(df_sum)
        df_perc['eai_exp_cat1'] = df_perc['eai_exp_cat1'] / pop_cat1*100
        df_perc['eai_exp_cat2'] = df_perc['eai_exp_cat2'] / pop_cat2 * 100
        df_perc['eai_exp_cat3'] = df_perc['eai_exp_cat3'] / pop_cat3 * 100
        df_perc['eai_exp_cat4'] = df_perc['eai_exp_cat4'] / pop_cat4 * 100
        df_perc['eai_exp_cat5'] = df_perc['eai_exp_cat5'] / pop_cat5 * 100
        df_perc['eai_exp_total'] = df_perc['eai_exp_total'] / pop_total * 100

        # save file containing share of people protected
        df_perc_rd = df_perc.round(2)
        df_perc_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat.csv')

        ## Calculate total number of people protected by habitats and difference to total number of people impacted

        # read file with impacts at all points protected by habitat
        file = f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
        df = pd.read_csv(file)
        types = 'Total'
        rhab_cat1 = df['eai_exp_cat1'].sum()
        rhab_cat2 = df['eai_exp_cat2'].sum()
        rhab_cat3 = df['eai_exp_cat3'].sum()
        rhab_cat4 = df['eai_exp_cat4'].sum()
        rhab_cat5 = df['eai_exp_cat5'].sum()
        rhab_total = df['eai_exp_total'].sum()

        rhab_sums = pd.DataFrame([[types,
                                   rhab_cat1,
                                   rhab_cat2,
                                   rhab_cat3,
                                   rhab_cat4,
                                   rhab_cat5,
                                   rhab_total]], columns=['rhab',
                                                          'eai_exp_cat1',
                                                          'eai_exp_cat2',
                                                          'eai_exp_cat3',
                                                          'eai_exp_cat4',
                                                          'eai_exp_cat5',
                                                          'eai_exp_total'])
        types = 'None'
        diff_cat1 = pop_cat1 - rhab_cat1
        diff_cat2 = pop_cat2 - rhab_cat2
        diff_cat3 = pop_cat3 - rhab_cat3
        diff_cat4 = pop_cat4 - rhab_cat4
        diff_cat5 = pop_cat5 - rhab_cat5
        diff_total = pop_total - rhab_total

        rhab_diff = pd.DataFrame([[types,
                                   diff_cat1,
                                   diff_cat2,
                                   diff_cat3,
                                   diff_cat4,
                                   diff_cat5,
                                   diff_total]], columns=['rhab',
                                                          'eai_exp_cat1',
                                                          'eai_exp_cat2',
                                                          'eai_exp_cat3',
                                                          'eai_exp_cat4',
                                                          'eai_exp_cat5',
                                                          'eai_exp_total'])

        # save file
        df_comp = pd.concat([rhab_sums, rhab_diff])
        df_comp_rd = df_comp.round()
        df_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp.csv')
        df_comp_total = pd.concat([df_sum_rd, df_comp_rd])
        df_comp_total.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp_total.csv')

        # difference in relative impact
        df_perc_comp = copy.copy(df_comp)
        df_perc_comp['eai_exp_cat1'] = df_perc_comp['eai_exp_cat1'] / pop_cat1*100
        df_perc_comp['eai_exp_cat2'] = df_perc_comp['eai_exp_cat2'] / pop_cat2 * 100
        df_perc_comp['eai_exp_cat3'] = df_perc_comp['eai_exp_cat3'] / pop_cat3 * 100
        df_perc_comp['eai_exp_cat4'] = df_perc_comp['eai_exp_cat4'] / pop_cat4 * 100
        df_perc_comp['eai_exp_cat5'] = df_perc_comp['eai_exp_cat5'] / pop_cat5 * 100
        df_perc_comp['eai_exp_total'] = df_perc_comp['eai_exp_total'] / pop_total * 100

        # save file
        df_perc_comp_rd = df_perc_comp.round(2)
        df_perc_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_comp.csv')
        df_perc_comp_total = pd.concat([df_perc_rd, df_perc_comp_rd])
        df_perc_comp_total.to_csv(f'{path}'
                                  f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_comp_total.csv')

## Run for climate models
clim_mods = ['EC-Earth3P-HR',
             'CMCC-CM2-VHR4',
             'CNRM-CM6-1-HR',
             'HadGEM3-GC31-HM']

for model in clim_mods:
    scenario = f'clim_STORM_{model}'
    rhab_yr = '2020'
    pop_yr = '2020'

    # load basin impact files
    AP_path = f'{path}AP_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    IO_path = f'{path}IO_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    WP_path = f'{path}WP_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'
    SH_path = f'{path}SH_{scenario}_{pop_yr}_{pop_data}_all_imp.csv'

    AP_df = pd.read_csv(AP_path)
    IO_df = pd.read_csv(IO_path)
    WP_df = pd.read_csv(WP_path)
    SH_df = pd.read_csv(SH_path)

    # merge basin data
    dfs = [AP_df, IO_df, WP_df, SH_df]
    df_all = pd.concat(dfs)

    # calculate impact across all windspeed categories
    df_all['eai_exp_total'] = df_all['eai_exp_cat1'] \
                              + df_all['eai_exp_cat2']\
                              + df_all['eai_exp_cat3'] \
                              + df_all['eai_exp_cat4'] \
                              + df_all['eai_exp_cat5']

    # save merged dataframe
    df_all.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')

    # load population data
    pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')
    pop = pop.reset_index(drop=True)
    pop['longitude'] = pop['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
    pop['latitude'] = pop['latitude'].round(decimals=11)

    # load habitat data
    rhab = pd.read_csv(f'{rhab_path}{pop_data}{pop_yr}_rhab{rhab_yr}.csv')
    rhab = rhab.reset_index(drop=True)
    rhab['longitude'] = rhab['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
    rhab['latitude'] = rhab['latitude'].round(decimals=11)

    # joining population impact data and habitat data
    join = pd.merge(pop, rhab, how='left', on=['longitude', 'latitude'])
    cols = ['longitude',
            'latitude',
            'eai_exp_cat1',
            'eai_exp_cat2',
            'eai_exp_cat3',
            'eai_exp_cat4',
            'eai_exp_cat5',
            'eai_exp_total',
            'rhab']
    join = join[cols]
    join['rhab'] = join['rhab'].replace(np.nan, 5)  # assigning habitat value 5 to all unmatched population exposures

    # save file containing impacts for all exposure points
    join.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}_incl_5_rhab.csv')
    # save file containing only impacts for habitat values below 5
    join = join[join.rhab != 5]
    join.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv')

    # Calculate aggregated impact data:
    # total & relative number of people impacted per wind speed category
    # total & relative number of people protected per wind speed category

    # read population impact data of all basins per exposure point (calculated above)
    pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')

    # summing population impacts per category
    pop_cat1 = pop['eai_exp_cat1'].sum()
    pop_cat2 = pop['eai_exp_cat2'].sum()
    pop_cat3 = pop['eai_exp_cat3'].sum()
    pop_cat4 = pop['eai_exp_cat4'].sum()
    pop_cat5 = pop['eai_exp_cat5'].sum()
    pop_total = pop['eai_exp_total'].sum()

    # creating dataframe of sums
    pop_sums = pd.DataFrame([[pop_cat1, pop_cat2, pop_cat3, pop_cat4, pop_cat5, pop_total]],
                            columns=['eai_exp_cat1',
                                     'eai_exp_cat2',
                                     'eai_exp_cat3',
                                     'eai_exp_cat4',
                                     'eai_exp_cat5',
                                     'eai_exp_total'
                                     ])

    # saving file with population impacts per wind speed category
    pop_sums_rd = pop_sums.round()
    pop_sums_rd.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat.csv')

    # read file containing exposure point-wise information on habitats
    file = f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
    df = pd.read_csv(file)

    # select only necessary columns
    col_select = ['rhab',
                  'eai_exp_cat1',
                  'eai_exp_cat2',
                  'eai_exp_cat3',
                  'eai_exp_cat4',
                  'eai_exp_cat5',
                  'eai_exp_total'
                  ]
    df = df[col_select]

    # group  and sum impacts by habitat value
    df_sum = df.groupby(['rhab'], as_index=False)[['eai_exp_cat1',
                                                   'eai_exp_cat2',
                                                   'eai_exp_cat3',
                                                   'eai_exp_cat4',
                                                   'eai_exp_cat5',
                                                   'eai_exp_total']].sum()
    # save file containing numbers of people protected per windspeed category
    df_sum_rd = df_sum.round()
    df_sum_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat.csv')

    # relative impact (share of people protected)
    df_perc = copy.copy(df_sum)
    df_perc['eai_exp_cat1'] = df_perc['eai_exp_cat1'] / pop_cat1 * 100
    df_perc['eai_exp_cat2'] = df_perc['eai_exp_cat2'] / pop_cat2 * 100
    df_perc['eai_exp_cat3'] = df_perc['eai_exp_cat3'] / pop_cat3 * 100
    df_perc['eai_exp_cat4'] = df_perc['eai_exp_cat4'] / pop_cat4 * 100
    df_perc['eai_exp_cat5'] = df_perc['eai_exp_cat5'] / pop_cat5 * 100
    df_perc['eai_exp_total'] = df_perc['eai_exp_total'] / pop_total * 100

    # save file containing share of people protected
    df_perc_rd = df_perc.round(2)
    df_perc_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat.csv')

    ## Calculate total number of people protected by habitats and difference to total number of people impacted

    # read file with impacts at all points protected by habitat
    file = f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
    df = pd.read_csv(file)
    types = 'Total'
    rhab_cat1 = df['eai_exp_cat1'].sum()
    rhab_cat2 = df['eai_exp_cat2'].sum()
    rhab_cat3 = df['eai_exp_cat3'].sum()
    rhab_cat4 = df['eai_exp_cat4'].sum()
    rhab_cat5 = df['eai_exp_cat5'].sum()
    rhab_total = df['eai_exp_total'].sum()

    rhab_sums = pd.DataFrame([[types,
                               rhab_cat1,
                               rhab_cat2,
                               rhab_cat3,
                               rhab_cat4,
                               rhab_cat5,
                               rhab_total]], columns=['rhab',
                                                      'eai_exp_cat1',
                                                      'eai_exp_cat2',
                                                      'eai_exp_cat3',
                                                      'eai_exp_cat4',
                                                      'eai_exp_cat5',
                                                      'eai_exp_total'])
    types = 'None'
    diff_cat1 = pop_cat1 - rhab_cat1
    diff_cat2 = pop_cat2 - rhab_cat2
    diff_cat3 = pop_cat3 - rhab_cat3
    diff_cat4 = pop_cat4 - rhab_cat4
    diff_cat5 = pop_cat5 - rhab_cat5
    diff_total = pop_total - rhab_total

    rhab_diff = pd.DataFrame([[types,
                               diff_cat1,
                               diff_cat2,
                               diff_cat3,
                               diff_cat4,
                               diff_cat5,
                               diff_total]], columns=['rhab',
                                                      'eai_exp_cat1',
                                                      'eai_exp_cat2',
                                                      'eai_exp_cat3',
                                                      'eai_exp_cat4',
                                                      'eai_exp_cat5',
                                                      'eai_exp_total'])

    # save file
    df_comp = pd.concat([rhab_sums, rhab_diff])
    df_comp_rd = df_comp.round()
    df_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp.csv')
    df_comp_total = pd.concat([df_sum_rd, df_comp_rd])
    df_comp_total.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp_total.csv')

    # difference in relative impact
    df_perc_comp = copy.copy(df_comp)
    df_perc_comp['eai_exp_cat1'] = df_perc_comp['eai_exp_cat1'] / pop_cat1 * 100
    df_perc_comp['eai_exp_cat2'] = df_perc_comp['eai_exp_cat2'] / pop_cat2 * 100
    df_perc_comp['eai_exp_cat3'] = df_perc_comp['eai_exp_cat3'] / pop_cat3 * 100
    df_perc_comp['eai_exp_cat4'] = df_perc_comp['eai_exp_cat4'] / pop_cat4 * 100
    df_perc_comp['eai_exp_cat5'] = df_perc_comp['eai_exp_cat5'] / pop_cat5 * 100
    df_perc_comp['eai_exp_total'] = df_perc_comp['eai_exp_total'] / pop_total * 100

    # save file
    df_perc_comp_rd = df_perc_comp.round(2)
    df_perc_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_comp.csv')
    df_perc_comp_total = pd.concat([df_perc_rd, df_perc_comp_rd])
    df_perc_comp_total.to_csv(
        f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_comp_total.csv')

## calculate median of all climate models
# variables
scenario = 'clim_STORM'
pop_data = 'wp'
rhab_yr = '2020'
pop_yr = '2020'

## calculate median per exposure point (population only)
CMCC_df = pd.read_csv(f'{path}CMCC-CM2-VHR4/all_basins_{scenario}_CMCC-CM2-VHR4_{pop_yr}_{pop_data}.csv')
CNRM_df = pd.read_csv(f'{path}CNRM-CM6-1-HR/all_basins_{scenario}_CNRM-CM6-1-HR_{pop_yr}_{pop_data}.csv')
HadGEM_df = pd.read_csv(f'{path}HadGEM3-GC31-HM/all_basins_{scenario}_HadGEM3-GC31-HM_{pop_yr}_{pop_data}.csv')
EC_Earth_df = pd.read_csv(f'{path}EC-Earth3P-HR/all_basins_{scenario}_EC-Earth3P-HR_{pop_yr}_{pop_data}.csv')

# merge basin data
dfs = [CMCC_df, CNRM_df, HadGEM_df, EC_Earth_df]
df_all = pd.concat(dfs)

# sum eai_exp of duplicate exposures
df_clean = df_all.groupby(['longitude', 'latitude'],
                          as_index=False)['eai_exp_cat1',
                                          'eai_exp_cat2',
                                          'eai_exp_cat3',
                                          'eai_exp_cat4',
                                          'eai_exp_cat5',
                                          'eai_exp_total'].median()

# save merged dataframe
df_clean.to_csv(f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}.csv')

## calculate median per exposure point including rhab
# joining population impact data and habitat data
pop = pd.read_csv(f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}.csv')
pop = pop.reset_index(drop=True)
pop['longitude'] = pop['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
pop['latitude'] = pop['latitude'].round(decimals=11)

rhab_path = '/Users/sarah/Jobs/Work/ETH_research_assistant/manuscript_calculations/new_data/exposures_rhab/'
rhab = pd.read_csv(f'{rhab_path}{pop_data}{pop_yr}_rhab{rhab_yr}.csv')
rhab = rhab.reset_index(drop=True)
rhab['longitude'] = rhab['longitude'].round(decimals=11)  # rounding decimals to ensure 1:1 matching
rhab['latitude'] = rhab['latitude'].round(decimals=11)

join = pd.merge(pop, rhab, how='left', on=['longitude', 'latitude'])
cols = ['longitude',
        'latitude',
        'eai_exp_cat1',
        'eai_exp_cat2',
        'eai_exp_cat3',
        'eai_exp_cat4',
        'eai_exp_cat5',
        'eai_exp_total',
        'rhab']
join = join[cols]
join['rhab'] = join['rhab'].replace(np.nan, 5)  # assigning habitat value 5 to all unmatched population exposures

# save file containing impacts for all exposure points
join.to_csv(f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}_rhab{rhab_yr}_incl_5_rhab.csv')
# save merged dataframe
join = join[join.rhab != 5]
join.to_csv(f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv')

# caculate impact per category
# read population impact data of all basins per exposure point (calculated above)
pop = pd.read_csv(f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}.csv')

# summing population impacts per category
pop_cat1 = pop['eai_exp_cat1'].sum()
pop_cat2 = pop['eai_exp_cat2'].sum()
pop_cat3 = pop['eai_exp_cat3'].sum()
pop_cat4 = pop['eai_exp_cat4'].sum()
pop_cat5 = pop['eai_exp_cat5'].sum()
pop_total = pop['eai_exp_total'].sum()

# creating dataframe of sums
pop_sums = pd.DataFrame([[pop_cat1, pop_cat2, pop_cat3, pop_cat4, pop_cat5, pop_total]],
                        columns=['eai_exp_cat1',
                                 'eai_exp_cat2',
                                 'eai_exp_cat3',
                                 'eai_exp_cat4',
                                 'eai_exp_cat5',
                                 'eai_exp_total'
                                 ])

# saving file with population impacts per wind speed category
pop_sums_rd = pop_sums.round()
pop_sums_rd.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_median_imp_per_cat.csv')

# read file containing exposure point-wise information on habitats
file = f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
df = pd.read_csv(file)

# select only necessary columns
col_select = ['rhab',
              'eai_exp_cat1',
              'eai_exp_cat2',
              'eai_exp_cat3',
              'eai_exp_cat4',
              'eai_exp_cat5',
              'eai_exp_total'
              ]
df = df[col_select]

# group  and sum impacts by habitat value
df_sum = df.groupby(['rhab'], as_index=False)[['eai_exp_cat1',
                                               'eai_exp_cat2',
                                               'eai_exp_cat3',
                                               'eai_exp_cat4',
                                               'eai_exp_cat5',
                                               'eai_exp_total']].sum()
# save file containing numbers of people protected per windspeed category
df_sum_rd = df_sum.round()
df_sum_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_imp_per_cat.csv')

# relative impact (share of people protected)
df_perc = copy.copy(df_sum)
df_perc['eai_exp_cat1'] = df_perc['eai_exp_cat1'] / pop_cat1 * 100
df_perc['eai_exp_cat2'] = df_perc['eai_exp_cat2'] / pop_cat2 * 100
df_perc['eai_exp_cat3'] = df_perc['eai_exp_cat3'] / pop_cat3 * 100
df_perc['eai_exp_cat4'] = df_perc['eai_exp_cat4'] / pop_cat4 * 100
df_perc['eai_exp_cat5'] = df_perc['eai_exp_cat5'] / pop_cat5 * 100
df_perc['eai_exp_total'] = df_perc['eai_exp_total'] / pop_total * 100

# save file containing share of people protected
df_perc_rd = df_perc.round(2)
df_perc_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_relative_imp_per_cat.csv')

## Calculate total number of people protected by habitats and difference to total number of people impacted

# read file with impacts at all points protected by habitat
file = f'{path}all_basins_{scenario}_median_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv'
df = pd.read_csv(file)
types = 'Total'
rhab_cat1 = df['eai_exp_cat1'].sum()
rhab_cat2 = df['eai_exp_cat2'].sum()
rhab_cat3 = df['eai_exp_cat3'].sum()
rhab_cat4 = df['eai_exp_cat4'].sum()
rhab_cat5 = df['eai_exp_cat5'].sum()
rhab_total = df['eai_exp_total'].sum()

rhab_sums = pd.DataFrame([[types,
                           rhab_cat1,
                           rhab_cat2,
                           rhab_cat3,
                           rhab_cat4,
                           rhab_cat5,
                           rhab_total]], columns=['rhab',
                                                  'eai_exp_cat1',
                                                  'eai_exp_cat2',
                                                  'eai_exp_cat3',
                                                  'eai_exp_cat4',
                                                  'eai_exp_cat5',
                                                  'eai_exp_total'])
types = 'None'
diff_cat1 = pop_cat1 - rhab_cat1
diff_cat2 = pop_cat2 - rhab_cat2
diff_cat3 = pop_cat3 - rhab_cat3
diff_cat4 = pop_cat4 - rhab_cat4
diff_cat5 = pop_cat5 - rhab_cat5
diff_total = pop_total - rhab_total

rhab_diff = pd.DataFrame([[types,
                           diff_cat1,
                           diff_cat2,
                           diff_cat3,
                           diff_cat4,
                           diff_cat5,
                           diff_total]], columns=['rhab',
                                                  'eai_exp_cat1',
                                                  'eai_exp_cat2',
                                                  'eai_exp_cat3',
                                                  'eai_exp_cat4',
                                                  'eai_exp_cat5',
                                                  'eai_exp_total'])

# save file
df_comp = pd.concat([rhab_sums, rhab_diff])
df_comp_rd = df_comp.round()
df_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_imp_per_cat_comp.csv')
df_comp_total = pd.concat([df_sum_rd, df_comp_rd])
df_comp_total.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_imp_per_cat_comp_total.csv')

# difference in relative impact
df_perc_comp = copy.copy(df_comp)
df_perc_comp['eai_exp_cat1'] = df_perc_comp['eai_exp_cat1'] / pop_cat1 * 100
df_perc_comp['eai_exp_cat2'] = df_perc_comp['eai_exp_cat2'] / pop_cat2 * 100
df_perc_comp['eai_exp_cat3'] = df_perc_comp['eai_exp_cat3'] / pop_cat3 * 100
df_perc_comp['eai_exp_cat4'] = df_perc_comp['eai_exp_cat4'] / pop_cat4 * 100
df_perc_comp['eai_exp_cat5'] = df_perc_comp['eai_exp_cat5'] / pop_cat5 * 100
df_perc_comp['eai_exp_total'] = df_perc_comp['eai_exp_total'] / pop_total * 100

# save file
df_perc_comp_rd = df_perc_comp.round(2)
df_perc_comp_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_relative_imp_per_cat_comp.csv')
df_perc_comp_total = pd.concat([df_perc_rd, df_perc_comp_rd])
df_perc_comp_total.to_csv(
    f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_median_relative_imp_per_cat_comp_total.csv')
