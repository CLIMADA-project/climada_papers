"""
This file is part of CLIMADA-papers.

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019. 

Plot scatter plots for evaluation countries.
Sections 3.4
Figure 5

Requires litpop_evaluation.py to be run first.

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.3.1
or later

The required gridded population data GPWv4.10 is available from SEDAC's Beta site, please see
https://beta.sedac.ciesin.columbia.edu/data/collection/gpw-v4/sets/browse

For more guidance on the LitPop module please refer to the CLIMADA tutorial:
https://climada-python.readthedocs.io/en/latest/tutorial/climada_entity_LitPop.html

@author: Samuel Eberenz
"""
import os
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import stats

from climada.entity.exposures import litpop as lp
from climada.util.constants import DATA_DIR
from climada.util.finance import income_group

# set output directory:
output_path = os.path.join(DATA_DIR, 'results')
input_path = output_path
if not os.path.isdir(output_path):
    os.mkdir(output_path)

experiment_name = 'v1'

nGRP_mod_fn = '%s_%s30_adm0.csv'
nGRP_ref_fn = '%s_%s30_adm1_ref.csv'


countries = ['AUS', 'BRA', 'CAN', 'CHE', 'CHN', 'DEU', 'FRA', 'IDN', 'IND', \
             'JPN', 'MEX', 'TUR', 'USA', 'ZAF']
methods = ['Lit', 'Lit2', 'Lit3', 'Lit4', 'Lit5', 'Pop', 'Pop2', 'Lit3Pop', 'Lit2Pop', 'LitPop']
meth_labels = [r'$Lit^1$', r'$Lit^2$', r'$Lit^3$', r'$Lit^4$', r'$Lit^5$', \
                   r'$Pop^1$', r'$Pop^2$', r'$Lit^3Pop^1$', r'$Lit^2Pop^1$', r'$Lit^1Pop^1$']
markers_list = ['o', 's', '^', 'd', 'p', 'o', 's', '^', 's', 'o']
colors3 = ['#1b9e77', '#7570b3', '#d95f02']
c3_10 = [0, 0, 0, 0, 0, 1, 1, 2, 2, 2]
methods_plot = [3, 5, 9]



nGRP = pd.DataFrame()
for idx, country in enumerate(countries):
    
    nGRP_country_ref = pd.read_csv(os.path.join(input_path, nGRP_ref_fn %(experiment_name, country)), \
                                   encoding="ISO-8859-1", header=0)
    nGRP_country = pd.read_csv(os.path.join(input_path, nGRP_mod_fn %(experiment_name, country)), \
                                   encoding="ISO-8859-1", header=0)
    nGRP_country['reference'] = nGRP_country_ref['LitPop']
    nGRP_country['country'] = country
    
    nGRP = pd.concat([nGRP, nGRP_country], ignore_index=True)
    nGRP = nGRP.dropna(subset=['reference'])
    
    print('%s: %i' %(country, nGRP.loc[nGRP['country']==country].shape[0]))
nGRP = nGRP.reset_index()
nGRP = nGRP.drop(columns=['index', 'Unnamed: 0'])
nGRP.to_csv(os.path.join(output_path, '%s_all30_nGRP.csv' %(experiment_name)))

meth_labels_plot = []
plt.figure(figsize=[8,8]) # Scatter plot per country
axes_sc = list()
for i in methods_plot:
    meth_labels_plot.append(meth_labels[i])
    axes_sc.append(plt.scatter(nGRP['reference'], nGRP[methods[i]], \
                   c=colors3[c3_10[i]], marker='.', \
                   alpha=1, s=50))

plt.plot([0, np.max([plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]])],
          [0, np.max([plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]])],\
          ls="--", c=".3")

plt.legend(axes_sc, meth_labels_plot)
plt.xlabel('Reference nGRP')
plt.ylabel('Modelled nGRP')

plt.savefig(os.path.join(output_path, '%s_all30_nGRP.pdf' %(experiment_name)), \
            dpi=300, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format='pdf',
            transparent=False, bbox_inches=None, pad_inches=0.1,
            frameon=None, metadata=None)

statistics = pd.DataFrame(index=['rho', 'beta', 'p_val', 'RMSF'], columns=methods)
for idx, meth in enumerate(methods):
    statistics.loc['rho', meth] = stats.pearsonr(nGRP['reference'], nGRP[meth])[0]
    statistics.loc['beta', meth] = stats.linregress(nGRP['reference'], nGRP[meth])[0]
    statistics.loc['RMSF', meth] = np.exp(np.sqrt(np.sum((np.log(nGRP[meth]/ \
                                        nGRP['reference']))**2)/nGRP.shape[0]))

    statistics.loc['p_val', meth] = stats.linregress(nGRP['reference'], nGRP[meth])[3]

statistics.to_csv(os.path.join(output_path, '%s_all30_nGRP_statistics.csv' %(experiment_name)))
