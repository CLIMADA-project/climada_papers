"""
Script to perform regional analysis of impact data
Author: Sarah Hülsen
"""

import pandas as pd
import copy

path = f'../results/final/'
region_path = '../data/'

## Split dataframes containing global impact and region information into individual df for each region
# dictionary with regions
reg_index = {
    'Australia and New Zealand': 1,
    'Caribbean': 2,
    'Central America': 3,
    'Central Asia': 4,
    'Eastern Africa': 5,
    'Eastern Asia': 6,
    'Eastern Europe': 7,
    'Europe': 8,
    'Middle Africa': 9,
    'Northern Africa': 10,
    'Northern America': 11,
    'Pacific Islands': 12,
    'South America': 13,
    'South-Eastern Asia': 14,
    'Southern Africa': 15,
    'Southern Asia': 16,
    'Western Africa': 17,
    'Western Asia': 18
}


# define function for splitting dataframes into regions
def region_split(impact_df, region):
    """ Split impact dataframes into regions for analysis """
    # get region subset
    subset = impact_df.loc[impact_df['custom_region'] == region]
    region_index = reg_index[region]
    subset = subset.assign(index=region_index)  # values from reg_index
    return subset


## hist and baseline datasets
# variables to determine which dataset is analysed
scenario = 'hist_STORM'
pop_data = 'wp'
pop_years = ['2000', '2020']
rhab_years = ['1992', '2020']

for pop_yr in pop_years:
    for rhab_yr in rhab_years:
        # read pop and rhab impact files
        imp_pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')
        imp_rhab = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv')

        imp_pop['longitude'] = imp_pop['longitude'].round(decimals=11)
        imp_pop['latitude'] = imp_pop['latitude'].round(decimals=11)

        imp_rhab['longitude'] = imp_rhab['longitude'].round(decimals=11)
        imp_rhab['latitude'] = imp_rhab['latitude'].round(decimals=11)

        # load file containing region and country information for each exposure point
        # prepared by matching each exposure point to the closest polygon of region layer
        reg_info = pd.read_csv(f'{region_path}custom_regions_analysis_file_{pop_yr}.csv')
        reg_info['longitude'] = reg_info['longitude'].round(decimals=11)  # round to ensure 1:1 matching
        reg_info['latitude'] = reg_info['latitude'].round(decimals=11)
        reg_info = reg_info.drop_duplicates(['longitude', 'latitude'], keep='first')  # remove duplicate matches

        # join impact files with files containing region and country information
        pop_regions = pd.merge(imp_pop, reg_info, how='left', on=['longitude', 'latitude'])
        rhab_regions = pd.merge(imp_rhab, reg_info, how='left', on=['longitude', 'latitude'])

        # save joined file
        # these files contain region, country, impact, and protection information per exposure point across all basins
        pop_regions.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_regions.csv')
        rhab_regions.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab_{rhab_yr}_regions.csv')

        ## Calculate total impact, protection, and relative protection per region
        # initialize empty lists to append regional dataframes
        df_list_pop = []        # for population exposure impacts
        df_list_abs = []        # for absolute rhab numbers
        df_list_rel = []        # for relative rhab numbers
        df_list_abs_comp = []   # for absolute rhab numbers including total and none rhab values
        df_list_rel_comp = []   # for relative rhab numbers including total and none rhab values

        for rgn in reg_index:
            # split the different dataframes (all exposure points per region)
            pop_split = region_split(pop_regions, rgn)
            rhab_split = region_split(rhab_regions, rgn)
            pop_split.to_csv(f"{path}all_basins_{scenario}_{pop_yr}_{pop_data}_{rgn}.csv")
            rhab_split.to_csv(f"{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}_{rgn}.csv")

            # sum impacts per category
            pop_cat1 = pop_split['eai_exp_cat1'].sum()
            pop_cat2 = pop_split['eai_exp_cat2'].sum()
            pop_cat3 = pop_split['eai_exp_cat3'].sum()
            pop_cat4 = pop_split['eai_exp_cat4'].sum()
            pop_cat5 = pop_split['eai_exp_cat5'].sum()
            pop_total = pop_split['eai_exp_total'].sum()
            region_name = rgn
            region_index = reg_index[rgn]

            # dataframe containing summed impacts per category
            pop_sums = pd.DataFrame([[region_name,
                                      region_index,
                                      pop_cat1,
                                      pop_cat2,
                                      pop_cat3,
                                      pop_cat4,
                                      pop_cat5,
                                      pop_total]],
                                    columns=['region',
                                             'region_index',
                                             'eai_exp_cat1',
                                             'eai_exp_cat2',
                                             'eai_exp_cat3',
                                             'eai_exp_cat4',
                                             'eai_exp_cat5',
                                             'eai_exp_total'])
            # round numbers
            pop_sums_rd = pop_sums.round()
            # save file
            pop_sums_rd.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_{rgn}.csv')
            df_list_pop.append(pop_sums_rd)

            # read file containing number of people protected per region
            col_select = ['rhab',
                          'eai_exp_cat1',
                          'eai_exp_cat2',
                          'eai_exp_cat3',
                          'eai_exp_cat4',
                          'eai_exp_cat5',
                          'eai_exp_total']
            rhab_split = rhab_split[col_select]
            rhab_sum = rhab_split.groupby(['rhab'], as_index=False)[['eai_exp_cat1',
                                                                     'eai_exp_cat2',
                                                                     'eai_exp_cat3',
                                                                     'eai_exp_cat4',
                                                                     'eai_exp_cat5',
                                                                     'eai_exp_total']].sum()
            # save file
            rhab_sum = rhab_sum.assign(region=region_name)
            rhab_sum = rhab_sum.assign(region_index=region_index)
            # figure out issue with region_index name shadowing
            rhab_sum_rd = rhab_sum.round()
            rhab_sum_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_{rgn}.csv')
            df_list_abs.append(rhab_sum_rd)

            ## calculate protection across rhab categories 1-4 (total) and 5 (none)
            # absolute numbers
            types = 'Total'
            rhab_cat1 = rhab_split['eai_exp_cat1'].sum()
            rhab_cat2 = rhab_split['eai_exp_cat2'].sum()
            rhab_cat3 = rhab_split['eai_exp_cat3'].sum()
            rhab_cat4 = rhab_split['eai_exp_cat4'].sum()
            rhab_cat5 = rhab_split['eai_exp_cat5'].sum()
            rhab_total = rhab_split['eai_exp_total'].sum()

            rhab_total_sums = pd.DataFrame([[types,
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
            df_comp = pd.concat([rhab_total_sums, rhab_diff])
            df_comp = df_comp.assign(region=region_name)
            df_comp = df_comp.assign(region_index=region_index)
            df_comp_rd = df_comp.round()
            df_comp_total = pd.concat([rhab_sum_rd, df_comp_rd])
            df_comp_total.to_csv(f'{path}'
                                 f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp_total_{rgn}.csv')
            df_list_abs_comp.append(df_comp_total)

            # calculate share of people protected (relative numbers)
            df_perc = copy.copy(rhab_sum)
            df_perc['eai_exp_cat1'] = df_perc['eai_exp_cat1'] / pop_cat1*100
            df_perc['eai_exp_cat2'] = df_perc['eai_exp_cat2'] / pop_cat2 * 100
            df_perc['eai_exp_cat3'] = df_perc['eai_exp_cat3'] / pop_cat3 * 100
            df_perc['eai_exp_cat4'] = df_perc['eai_exp_cat4'] / pop_cat4 * 100
            df_perc['eai_exp_cat5'] = df_perc['eai_exp_cat5'] / pop_cat5 * 100
            df_perc['eai_exp_total'] = df_perc['eai_exp_total'] / pop_total * 100

            # round numbers
            df_perc_rd = df_perc.round(1)
            # save file
            df_perc_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_{rgn}.csv')
            df_list_rel.append(df_perc_rd)

            ## calculate protection across rhab categories 1-4 (total) and 5 (none)
            # relative numbers
            df_perc_comp = copy.copy(df_comp)
            df_perc_comp['eai_exp_cat1'] = df_perc_comp['eai_exp_cat1'] / pop_cat1 * 100
            df_perc_comp['eai_exp_cat2'] = df_perc_comp['eai_exp_cat2'] / pop_cat2 * 100
            df_perc_comp['eai_exp_cat3'] = df_perc_comp['eai_exp_cat3'] / pop_cat3 * 100
            df_perc_comp['eai_exp_cat4'] = df_perc_comp['eai_exp_cat4'] / pop_cat4 * 100
            df_perc_comp['eai_exp_cat5'] = df_perc_comp['eai_exp_cat5'] / pop_cat5 * 100
            df_perc_comp['eai_exp_total'] = df_perc_comp['eai_exp_total'] / pop_total * 100

            # save file
            df_perc_comp = df_perc_comp.assign(region=region_name)
            df_perc_comp = df_perc_comp.assign(region_index=region_index)
            df_perc_comp_rd = df_perc_comp.round(2)
            df_perc_comp_total = pd.concat([df_perc_rd, df_perc_comp_rd])
            df_perc_comp_total.to_csv(f'{path}'
                                      f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}'
                                      f'_relative_imp_per_cat_comp_total_{rgn}.csv')
            df_list_rel_comp.append(df_perc_comp_total)

        ## Concatenate individual regions dataframes collected in lists above
        final_df_pop = pd.concat(df_list_pop)
        final_df_abs = pd.concat(df_list_abs)
        final_df_rel = pd.concat(df_list_rel)
        final_df_abs_comp = pd.concat(df_list_abs_comp)
        final_df_rel_comp = pd.concat(df_list_rel_comp)

        # save final dataframes to csv
        final_df_pop.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_all_regions.csv')
        final_df_abs.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_all_regions.csv')
        final_df_rel.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_all_regions.csv')
        final_df_abs_comp.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}'
                                 f'_imp_per_cat_comp_total_all_regions.csv')
        final_df_rel_comp.to_csv(f'{path}'
                                 f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_all_regions.csv')


## clim change dataset
model = 'median'
scenario = f'clim_STORM_{model}'
pop_yr = '2020'
rhab_yr = '2020'

# read pop and rhab impact files
imp_pop = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}.csv')
imp_rhab = pd.read_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}.csv')

imp_pop['longitude'] = imp_pop['longitude'].round(decimals=11)
imp_pop['latitude'] = imp_pop['latitude'].round(decimals=11)

imp_rhab['longitude'] = imp_rhab['longitude'].round(decimals=11)
imp_rhab['latitude'] = imp_rhab['latitude'].round(decimals=11)

# load file containing region and country information for each exposure point
# prepared by matching each exposure point to the closest polygon of region layer
reg_info = pd.read_csv(f'{region_path}custom_regions_analysis_file_{pop_yr}.csv')
reg_info['longitude'] = reg_info['longitude'].round(decimals=11)  # round to ensure 1:1 matching
reg_info['latitude'] = reg_info['latitude'].round(decimals=11)
reg_info = reg_info.drop_duplicates(['longitude', 'latitude'], keep='first')  # remove duplicate matches

# join impact files with files containing region and country information
pop_regions = pd.merge(imp_pop, reg_info, how='left', on=['longitude', 'latitude'])
rhab_regions = pd.merge(imp_rhab, reg_info, how='left', on=['longitude', 'latitude'])

# save joined file
# these files now contain region, country, impact, and protection information per exposure point across all basins
pop_regions.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_regions.csv')
rhab_regions.to_csv(f'{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab_{rhab_yr}_regions.csv')

## Calculate total impact, protection, and relative protection per region
# initialize empty lists to append regional dataframes
df_list_pop = []        # for population exposure impacts
df_list_abs = []        # for absolute rhab numbers
df_list_rel = []        # for relative rhab numbers
df_list_abs_comp = []   # for absolute rhab numbers including total and none rhab values
df_list_rel_comp = []   # for relative rhab numbers including total and none rhab values

for rgn in reg_index:
    # split the different dataframes (all exposure points per region)
    pop_split = region_split(pop_regions, rgn)
    rhab_split = region_split(rhab_regions, rgn)
    pop_split.to_csv(f"{path}all_basins_{scenario}_{pop_yr}_{pop_data}_{rgn}.csv")
    rhab_split.to_csv(f"{path}all_basins_{scenario}_{pop_yr}_{pop_data}_rhab{rhab_yr}_{rgn}.csv")

    # sum impacts per category
    pop_cat1 = pop_split['eai_exp_cat1'].sum()
    pop_cat2 = pop_split['eai_exp_cat2'].sum()
    pop_cat3 = pop_split['eai_exp_cat3'].sum()
    pop_cat4 = pop_split['eai_exp_cat4'].sum()
    pop_cat5 = pop_split['eai_exp_cat5'].sum()
    pop_total = pop_split['eai_exp_total'].sum()
    region_name = rgn
    region_index = reg_index[rgn]

    # dataframe containing summed impacts per category
    pop_sums = pd.DataFrame([[region_name,
                              region_index,
                              pop_cat1,
                              pop_cat2,
                              pop_cat3,
                              pop_cat4,
                              pop_cat5,
                              pop_total]],
                            columns=['region',
                                     'region_index',
                                     'eai_exp_cat1',
                                     'eai_exp_cat2',
                                     'eai_exp_cat3',
                                     'eai_exp_cat4',
                                     'eai_exp_cat5',
                                     'eai_exp_total'])
    # round numbers
    pop_sums_rd = pop_sums.round()
    # save file
    pop_sums_rd.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_{rgn}.csv')
    df_list_pop.append(pop_sums_rd)

    # read file containing number of people protected per region
    col_select = ['rhab',
                  'eai_exp_cat1',
                  'eai_exp_cat2',
                  'eai_exp_cat3',
                  'eai_exp_cat4',
                  'eai_exp_cat5',
                  'eai_exp_total']
    rhab_split = rhab_split[col_select]
    rhab_sum = rhab_split.groupby(['rhab'], as_index=False)[['eai_exp_cat1',
                                                             'eai_exp_cat2',
                                                             'eai_exp_cat3',
                                                             'eai_exp_cat4',
                                                             'eai_exp_cat5',
                                                             'eai_exp_total']].sum()
    # save file
    rhab_sum = rhab_sum.assign(region=region_name)
    rhab_sum = rhab_sum.assign(region_index=region_index)
    # figure out issue with region_index name shadowing
    rhab_sum_rd = rhab_sum.round()
    rhab_sum_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_{rgn}.csv')
    df_list_abs.append(rhab_sum_rd)

    ## calculate protection across rhab categories 1-4 (total) and 5 (none)
    # absolute numbers
    types = 'Total'
    rhab_cat1 = rhab_split['eai_exp_cat1'].sum()
    rhab_cat2 = rhab_split['eai_exp_cat2'].sum()
    rhab_cat3 = rhab_split['eai_exp_cat3'].sum()
    rhab_cat4 = rhab_split['eai_exp_cat4'].sum()
    rhab_cat5 = rhab_split['eai_exp_cat5'].sum()
    rhab_total = rhab_split['eai_exp_total'].sum()

    rhab_total_sums = pd.DataFrame([[types,
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
    df_comp = pd.concat([rhab_total_sums, rhab_diff])
    df_comp = df_comp.assign(region=region_name)
    df_comp = df_comp.assign(region_index=region_index)
    df_comp_rd = df_comp.round()
    df_comp_total = pd.concat([rhab_sum_rd, df_comp_rd])
    df_comp_total.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_comp_total_{rgn}.csv')
    df_list_abs_comp.append(df_comp_total)

    # calculate share of people protected (relative numbers)
    df_perc = copy.copy(rhab_sum)
    df_perc['eai_exp_cat1'] = df_perc['eai_exp_cat1'] / pop_cat1*100
    df_perc['eai_exp_cat2'] = df_perc['eai_exp_cat2'] / pop_cat2 * 100
    df_perc['eai_exp_cat3'] = df_perc['eai_exp_cat3'] / pop_cat3 * 100
    df_perc['eai_exp_cat4'] = df_perc['eai_exp_cat4'] / pop_cat4 * 100
    df_perc['eai_exp_cat5'] = df_perc['eai_exp_cat5'] / pop_cat5 * 100
    df_perc['eai_exp_total'] = df_perc['eai_exp_total'] / pop_total * 100

    # round numbers
    df_perc_rd = df_perc.round(1)
    # save file
    df_perc_rd.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_{rgn}.csv')
    df_list_rel.append(df_perc_rd)

    ## calculate protection across rhab categories 1-4 (total) and 5 (none)
    # relative numbers
    df_perc_comp = copy.copy(df_comp)
    df_perc_comp['eai_exp_cat1'] = df_perc_comp['eai_exp_cat1'] / pop_cat1 * 100
    df_perc_comp['eai_exp_cat2'] = df_perc_comp['eai_exp_cat2'] / pop_cat2 * 100
    df_perc_comp['eai_exp_cat3'] = df_perc_comp['eai_exp_cat3'] / pop_cat3 * 100
    df_perc_comp['eai_exp_cat4'] = df_perc_comp['eai_exp_cat4'] / pop_cat4 * 100
    df_perc_comp['eai_exp_cat5'] = df_perc_comp['eai_exp_cat5'] / pop_cat5 * 100
    df_perc_comp['eai_exp_total'] = df_perc_comp['eai_exp_total'] / pop_total * 100

    # save file
    df_perc_comp = df_perc_comp.assign(region=region_name)
    df_perc_comp = df_perc_comp.assign(region_index=region_index)
    df_perc_comp_rd = df_perc_comp.round(2)
    df_perc_comp_total = pd.concat([df_perc_rd, df_perc_comp_rd])
    df_perc_comp_total.to_csv(f'{path}rrhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}'
                              f'_relative_imp_per_cat_comp_total_{rgn}.csv')
    df_list_rel_comp.append(df_perc_comp_total)

## Concatenate individual regions dataframes collected in lists above
final_df_pop = pd.concat(df_list_pop)
final_df_abs = pd.concat(df_list_abs)
final_df_rel = pd.concat(df_list_rel)
final_df_abs_comp = pd.concat(df_list_abs_comp)
final_df_rel_comp = pd.concat(df_list_rel_comp)

# save final dataframes to csv
final_df_pop.to_csv(f'{path}pop_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_all_regions.csv')
final_df_abs.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_imp_per_cat_all_regions.csv')
final_df_rel.to_csv(f'{path}'
                    f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_all_regions.csv')
final_df_abs_comp.to_csv(f'{path}rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}'
                         f'_imp_per_cat_comp_total_all_regions.csv')
final_df_rel_comp.to_csv(f'{path}'
                         f'rhab_{rhab_yr}_{pop_data}_{pop_yr}_{scenario}_relative_imp_per_cat_all_regions.csv')
