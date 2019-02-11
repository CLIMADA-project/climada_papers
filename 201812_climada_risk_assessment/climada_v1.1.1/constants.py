"""
"""

RESOL = 0.5
""" Approx. resolution in km """

YEAR = 2016
""" Year of exposure data """

CNTRIES = ['Saint Barthelemy', 'Saint Martin', 'Sint Maarten', 'Anguilla',
           'British Virgin Islands', 'United States Virgin Islands',
           'Turks And Caicos Islands', 'Saint Kitts And Nevis',
           'Antigua And Barbuda', 'Netherlands']
""" Country (island groups) names """

CNTRIES_ISO = ['BLM', 'MAF', 'SXM', 'AIA', 'VGB', 'VIR', 'TCA', 'KNA', 'ATG', 'NLD']
""" Country (island groups) ISO3 codes """

GDP = {'BLM': 414710000, 'MAF': 614258169, 'SXM': 1081577185, \
       'AIA': 337201995, 'VGB': 971237110, \
       'VIR': 3765000000, 'TCA': 917550492, \
       'KNA': 909854630, 'ATG': 1460144703, 'NLD': ''}
""" GDP at YEAR per island group """

GDP_NLD_ISL = 48.0e6 + 100.0e6
""" GDP Saba and St. Eustatius """

INC_GRP_DEF = 4
INC_GRP = {'BLM': INC_GRP_DEF, 'MAF': INC_GRP_DEF, 'SXM': INC_GRP_DEF,
           'AIA': INC_GRP_DEF, 'VGB': INC_GRP_DEF, 'VIR': INC_GRP_DEF,
           'TCA': INC_GRP_DEF, 'KNA': INC_GRP_DEF, 'ATG': INC_GRP_DEF,
           'NLD': INC_GRP_DEF}
""" income group level at YEAR per island group """

POLY_VAL = [0, 0, 1]
""" Polygonal transformation in night lights """