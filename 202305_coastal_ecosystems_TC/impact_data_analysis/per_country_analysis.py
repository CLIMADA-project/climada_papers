"""
Script to perform country-level analysis of impact data
Author: Sarah Hülsen
"""

import pandas as pd
import copy

# data variables
model = 'median'
scenario = f'hist_STORM'  # 'clim_STORM_{model}'
pop_data = 'wp'
pop_yr = '2020'
rhab_yr = '2020'

# data paths
path = f'../results/final/'

## for population impacts (without rhab)
# load impact data
pop_imp = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_regions.csv')

# keep only relevant columns
cols = ['eai_exp_cat1',     # number of people impacted per wind speed category & in total
        'eai_exp_cat2',
        'eai_exp_cat3',
        'eai_exp_cat4',
        'eai_exp_cat5',
        'eai_exp_total',
        'ISO3',             # ISO3 country code (e.g. GBR - Great Britain)
        'country_name',
        'custom_region']
pop_imp = pop_imp[cols]

# group by countries and sum impacts
pop_imp_sum = pop_imp.groupby(['country_name',
                               'ISO3',
                               'custom_region'], as_index=False)[['eai_exp_cat1',
                                                                  'eai_exp_cat2',
                                                                  'eai_exp_cat3',
                                                                  'eai_exp_cat4',
                                                                  'eai_exp_cat5',
                                                                  'eai_exp_total']].sum()
# save file
pop_imp_sum_rd = pop_imp_sum.round()
pop_imp_sum_rd.to_csv(f'{path}per_country_{pop_data}_{pop_yr}_{scenario}_imp_per_cat.csv')

## for rhab dataframes
rhab_imp = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab_{rhab_yr}_regions.csv')

cols = ['eai_exp_cat1',     # number of people impacted per wind speed category & in total
        'eai_exp_cat2',
        'eai_exp_cat3',
        'eai_exp_cat4',
        'eai_exp_cat5',
        'eai_exp_total',
        'rhab',             # habitat ranking
        'ISO3',             # ISO3 country code (e.g. GBR - Great Britain)
        'country_name',
        'custom_region']

rhab_imp = rhab_imp[cols]
rhab_imp_sum = rhab_imp.groupby(['country_name',
                                 'ISO3',
                                 'custom_region',
                                 'rhab'], as_index=False)[['eai_exp_cat1',
                                                           'eai_exp_cat2',
                                                           'eai_exp_cat3',
                                                           'eai_exp_cat4',
                                                           'eai_exp_cat5',
                                                           'eai_exp_total']].sum()
# save file
rhab_imp_sum_rd = rhab_imp_sum.round()
rhab_imp_sum_rd.to_csv(f'{path}per_country_{pop_data}_{pop_yr}_rhab{rhab_yr}_{scenario}_imp_per_cat.csv')

## relative numbers (share of rhab as percentage of total numbers)
# rename columns before merging pop impact dataframe with rhab impact datafram
cols = ['country_name',
        'eai_exp_cat1',
        'eai_exp_cat2',
        'eai_exp_cat3',
        'eai_exp_cat4',
        'eai_exp_cat5',
        'eai_exp_total']
pop_imp_sum_rd = pop_imp_sum_rd[cols]
pop_imp_sum_rd = pop_imp_sum_rd.rename(columns={'eai_exp_cat1': 'pop_eai_exp_cat1',
                                                'eai_exp_cat2': 'pop_eai_exp_cat2',
                                                'eai_exp_cat3': 'pop_eai_exp_cat3',
                                                'eai_exp_cat4': 'pop_eai_exp_cat4',
                                                'eai_exp_cat5': 'pop_eai_exp_cat5',
                                                'eai_exp_total': 'pop_eai_exp_total'})
# merge dataframes on country level
join = pd.merge(rhab_imp_sum_rd, pop_imp_sum_rd, how='left', on=['country_name'])

# calculate relative numbers
join['relative_eai_exp_cat1'] = join['eai_exp_cat1'] / join['pop_eai_exp_cat1'] * 100
join['relative_eai_exp_cat2'] = join['eai_exp_cat2'] / join['pop_eai_exp_cat2'] * 100
join['relative_eai_exp_cat3'] = join['eai_exp_cat3'] / join['pop_eai_exp_cat3'] * 100
join['relative_eai_exp_cat4'] = join['eai_exp_cat4'] / join['pop_eai_exp_cat4'] * 100
join['relative_eai_exp_cat5'] = join['eai_exp_cat5'] / join['pop_eai_exp_cat5'] * 100
join['relative_eai_exp_total'] = join['eai_exp_total'] / join['pop_eai_exp_total'] * 100

join = join.reset_index(drop=True)
# drop rhab impact columns OR save only one file containing both absolute and relative numbers
# add row per country containing total rhab percentage
# copy df
join_total = copy.copy(join)
# in copied df drop rhab column
join_total = join_total.drop('rhab', axis=1)
# group by country and sum absolute and relative impacts
join_total = join_total.groupby(['country_name',
                                 'ISO3',
                                 'custom_region',
                                 'pop_eai_exp_cat1',
                                 'pop_eai_exp_cat2',
                                 'pop_eai_exp_cat3',
                                 'pop_eai_exp_cat4',
                                 'pop_eai_exp_cat5',
                                 'pop_eai_exp_total'], as_index=False)[['eai_exp_cat1',
                                                                        'eai_exp_cat2',
                                                                        'eai_exp_cat3',
                                                                        'eai_exp_cat4',
                                                                        'eai_exp_cat5',
                                                                        'eai_exp_total',
                                                                        'relative_eai_exp_cat1',
                                                                        'relative_eai_exp_cat2',
                                                                        'relative_eai_exp_cat3',
                                                                        'relative_eai_exp_cat4',
                                                                        'relative_eai_exp_cat5',
                                                                        'relative_eai_exp_total']].sum()
# add new column and assign rhab='Total'
join_total = join_total.assign(rhab='Total')
# concatenate both dataframes
join_concat = pd.concat([join, join_total])
join_concat = join_concat.sort_values(by=['country_name', 'rhab'], ascending=[True, True])
join_concat = join_concat.reset_index(drop=True)

# save dataframe
join = join.round(2)
join_concat = join_concat.round(2)
join.to_csv(f'{path}per_country_{pop_data}_{pop_yr}_rhab{rhab_yr}_{scenario}_imp_per_cat_relative.csv')
join_concat.to_csv(f'{path}'
                   f'per_country_{pop_data}_{pop_yr}_rhab{rhab_yr}_{scenario}_imp_per_cat_relative_total.csv')
