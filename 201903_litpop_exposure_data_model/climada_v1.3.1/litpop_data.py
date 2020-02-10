"""
This file is part of CLIMADA-papers.

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019. 

Reproduces LitPop based gridded asset exposure data as in the ETH data repository:
LitPop: Global Exposure Data for Disaster Risk Assessment
DOI: 10.3929/ethz-b-000331316
https://www.research-collection.ethz.ch/handle/20.500.11850/331316

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0
or later

The required gridded population data GPWv4.10 is available from SEDAC's Beta site, please see
https://beta.sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse

For more guidance on the LitPop module please refer to the CLIMADA tutorial:
https://climada-python.readthedocs.io/en/latest/tutorial/climada_entity_LitPop.html

@author: Samuel Eberenz
"""
import os
import pandas as pd
import numpy as np
from iso3166 import countries as iso_cntry

from climada.entity.exposures import litpop as lp
from climada.util.constants import DATA_DIR
from climada.util.finance import world_bank_wealth_account

# set output directory:
output_path = os.path.join(DATA_DIR, 'results')

# Data export choices:
export_litpop_norm = True # export gridded normalised Lit and Pop for (14) selected countries
export_litpop_all = True # export gridded LitPop for all (227) countries
save_csv = True # save entities as CSV-file
save_hdf5 = False # save entities as HDF5-file

version = 'v1'
yyyy = 2014 # reference year for input data
fin_mode = 'pc' # produced capital as prefered asset base
resolution = 30 # resolution of exposure map. set to 30 for best results (slow)

countries = ['AUS', 'BRA', 'CAN', 'CHE', 'CHN', 'DEU', 'FRA', 'IDN', 'IND', \
             'JPN', 'MEX', 'TUR', 'USA', 'ZAF'] # selected countries for extra data
countries = sorted(countries)
countries_sel = np.arange(0, len(countries)) # indices of countries in countries list

if not os.path.isdir(output_path): # create output directory if it does not exist
    os.mkdir(output_path)

print('CREATE AND EXPORT GRIDDED DATA:')
for i in countries_sel:
    print('*** ' + countries[i] + ' *** ')
    ent = lp.LitPop()
    # Gridded LitPop for selected countries:
    ent.set_country(countries[i], res_arcsec=resolution, fin_mode=fin_mode, \
                    exponents=[1, 1], reference_year=yyyy)
    if save_csv:
        ent.to_csv(os.path.join(output_path, 'LitPop_pc_' + str(resolution) + \
                                'arcsec_' + countries[i] +  '_' + version + '.csv'))
    if save_hdf5:
        ent.write_hdf5(os.path.join(output_path, 'LitPop_pc_' + str(resolution) + \
                                'arcsec_' + countries[i] +  '_' + version + '.hdf5'))
    if export_litpop_norm:
        # Normalized gridded Lit and Pop for selected countries:
        ent.set_country(countries[i], res_arcsec=resolution, fin_mode='norm', \
                        exponents=[1, 0], reference_year=yyyy)
        if save_csv:
            ent.to_csv(os.path.join(output_path, 'Lit_norm_' + str(resolution) + \
                                    'arcsec_' + countries[i] +  '_' + version + '.csv'))
        if save_hdf5:
            ent.write_hdf5(os.path.join(output_path, 'Lit_norm_' + str(resolution) + \
                                    'arcsec_' + countries[i] +  '_' + version + '.hdf5'))
        ent.set_country(countries[i], res_arcsec=resolution, fin_mode='norm', \
                        exponents=[0, 1], reference_year=yyyy)
        if save_csv:
            ent.to_csv(os.path.join(output_path, 'Pop_norm_' + str(resolution) + \
                                    'arcsec_' + countries[i] + '_' + version + '.csv'))
        if save_hdf5:
            ent.write_hdf5(os.path.join(output_path, 'Pop_norm_' + str(resolution) + \
                                    'arcsec_' + countries[i] +  '_' + version + '.hdf5'))

if export_litpop_all:
    success = list()
    failure = list()
    yyyy = 2014
    for c in iso_cntry:
        print('*** ' + c.alpha3 + ' ***')
        if not (c.alpha3 in countries):
            # Try computing and exporting gridded LitPop for rest of countries:
            try:
                ent = lp.LitPop()
                ent.set_country(c.alpha3, res_arcsec=resolution, fin_mode=fin_mode, \
                            exponents=[1, 1], reference_year=yyyy)
                if save_csv:
                    ent.to_csv(os.path.join(output_path, 'LitPop_pc_' + \
                                            str(resolution) + 'arcsec_' + \
                                            c.alpha3 + '_' + version + '.csv'))
                if save_hdf5:
                    ent.write_hdf5(os.path.join(output_path, 'LitPop_pc_' + \
                                            str(resolution) + 'arcsec_' + \
                                            c.alpha3 + '_' + version + '.hdf5'))
                success.append(c.alpha3)
                print('Success: ' + c.name)
            except:
                failure.append(c.alpha3)
                print('Fail: ' + c.name)
        else:
            print('Skipped: ' + c.name)
    print(failure)
    print(success)
"""
# Lists of countries for which LitPop computation worked (2019/03):
if not 'countries' in locals():
        countries = ['AUS', 'BRA', 'CAN', 'CHE', 'CHN', 'DEU', 'FRA', 'IDN', \
                     'IND', 'JPN', 'MEX', 'TUR', 'USA', 'ZAF']
if not 'success' in locals():
        success = ['AFG', 'ALB', 'DZA', 'ASM', 'AND', 'AGO', 'AIA', 'ATG', \
                   'ARG', 'ARM', 'ABW', 'AUT', 'AZE', 'BHS', 'BHR', 'BGD', \
                   'BRB', 'BLR', 'BEL', 'BLZ', 'BEN', 'BMU', 'BTN', 'BOL', \
                   'BIH', 'BWA', 'IOT', 'BRN', 'BGR', 'BFA', 'BDI', 'KHM', \
                   'CMR', 'CPV', 'CYM', 'CAF', 'TCD', 'CHL', 'COL', 'COM', \
                   'COG', 'COD', 'COK', 'CRI', 'CIV', 'HRV', 'CUB', 'CUW', \
                   'CYP', 'CZE', 'DNK', 'DJI', 'DMA', 'DOM', 'ECU', 'EGY', \
                   'SLV', 'GNQ', 'ERI', 'EST', 'ETH', 'FLK', 'FRO', 'FJI', \
                   'FIN', 'PYF', 'ATF', 'GAB', 'GMB', 'GEO', 'GHA', 'GIB', \
                   'GRC', 'GRL', 'GRD', 'GUM', 'GTM', 'GGY', 'GIN', 'GNB', \
                   'GUY', 'HTI', 'HND', 'HKG', 'HUN', 'ISL', 'IRN', 'IRQ', \
                   'IRL', 'IMN', 'ISR', 'ITA', 'JAM', 'JEY', 'JOR', 'KAZ', \
                   'KEN', 'KIR', 'PRK', 'KOR', 'KWT', 'KGZ', 'LAO', 'LVA', \
                   'LBN', 'LSO', 'LBR', 'LBY', 'LIE', 'LTU', 'LUX', 'MAC', \
                   'MKD', 'MDG', 'MWI', 'MYS', 'MDV', 'MLI', 'MLT', 'MHL', \
                   'MRT', 'MUS', 'FSM', 'MDA', 'MCO', 'MNG', 'MNE', 'MSR', \
                   'MAR', 'MOZ', 'MMR', 'NAM', 'NRU', 'NPL', 'NLD', 'NCL', \
                   'NZL', 'NIC', 'NER', 'NGA', 'NIU', 'NFK', 'MNP', 'NOR', \
                   'OMN', 'PAK', 'PLW', 'PAN', 'PNG', 'PRY', 'PER', 'PHL', \
                   'POL', 'PRT', 'PRI', 'QAT', 'ROU', 'RUS', 'RWA', 'BLM', \
                   'SHN', 'KNA', 'LCA', 'MAF', 'SPM', 'VCT', 'WSM', 'SMR', \
                   'STP', 'SAU', 'SEN', 'SRB', 'SYC', 'SLE', 'SGP', 'SXM', \
                   'SVK', 'SVN', 'SLB', 'SOM', 'SGS', 'ESP', 'LKA', 'SDN', \
                   'SUR', 'SWZ', 'SWE', 'SYR', 'TWN', 'TJK', 'TZA', 'THA', \
                   'TLS', 'TGO', 'TON', 'TTO', 'TUN', 'TKM', 'TCA', 'TUV', \
                   'UGA', 'UKR', 'ARE', 'GBR', 'URY', 'UZB', 'VUT', 'VEN', \
                   'VNM', 'VGB', 'VIR', 'WLF', 'YEM', 'ZMB', 'ZWE']

"""
# write countries' metadata to CSV:
countries_metadata = pd.DataFrame(columns=['country_name','iso3',\
                                             'region_id', 'total_value', 'pc'])
if not 'success' in locals():
    success=list()
for c in iso_cntry:
    print(c.name)
    if (c.alpha3 in success) or (c.alpha3 in countries):
        countries_metadata = countries_metadata.append({'country_name' : c.name, \
                               'iso3' : c.alpha3, \
                               'region_id' : c.numeric, \
                               'total_value' : world_bank_wealth_account(c.alpha3, 2014)[1] , \
                               'pc' : world_bank_wealth_account(c.alpha3, 2014)[2]} , \
                               ignore_index=True)
    else:
        countries_metadata = countries_metadata.append({'country_name' : c.name, \
                               'iso3' : c.alpha3, \
                               'region_id' : c.numeric, \
                               'total_value' : float('NaN')} , \
                               ignore_index=True)

countries_metadata.to_csv(os.path.join(output_path, \
                                       'LitPop_countries_META_pc_' + str(yyyy) + \
                                       '_' + version + '.csv'))
