"""
This file is part of CLIMADA-papers.

Eberenz, S., Stocker, D., Röösli, T., and Bresch, D. N.:
Exposure data for global physical risk assessment,
Earth Syst. Sci. Data Discuss., https://doi.org/10.5194/essd-2019-189, in review, 2019. 

LitPop exposure data model evaluation for 14 countries and plotting of scatter and box plots
Sections 2.6, 3.2, 3.3
Figures 3, 5
Tables (A1), A2, A3

Requires https://github.com/CLIMADA-project/climada_python/releases/tag/v1.2.0
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
if not os.path.isdir(output_path):
    os.mkdir(output_path)

experiment_name = 'v1'

# SWITCH FEATURES ON/OFF:
"""
- compute_validation:
Compute full validation statistics for all selected countries:
Pearson correlation coefficient (rho), linear regression slope beta,
and root mean squared fraction RMSF  for variations of Lit^n * Pop^m.
# warning: computational intensive. This can take several hours.
Plots scatter plots per country (Figure 4).

- validation_plots:
Make and save box plots (Figure 3).
This requires compute_validation to be run for all selected countries first.
"""
compute_validation = True   # default: True
validation_plots = True     # default: True

# quick_test: make quick test to check whether engine works
quick_test = False

# countries: list of countries taken into account.
# each country requires regional gross regional product (GRP) data file as XLS in data > system > GSDP
countries = ['AUS', 'BRA', 'CAN', 'CHE', 'CHN', 'DEU', 'FRA', 'IDN', 'IND', \
             'JPN', 'MEX', 'TUR', 'USA', 'ZAF']

countries = sorted(countries)
if quick_test:
    countries_sel = [3] # Switzerland only
    resolution = 120 # reduced resolution
    experiment_name = 'test'
else:
    countries_sel = np.arange(0, len(countries)) # all countries in list
    # set resolution of exposure in arcsec. set to 30 for best results (slow):
    resolution = 30

# name per method evaluated:
methods = ['Lit', 'Lit2', 'Lit3', 'Lit4', 'Lit5', 'Pop', 'Pop2', 'Lit3Pop', 'Lit2Pop', 'LitPop']
# exponents per method:
exponents_list = [[1, 0], [2, 0], [3, 0], [4, 0], [5, 0], [0, 1], [0, 2], [3, 1], [2, 1], [1, 1]]
# choose which method to include in plots + marker style per method:
methods_show = [True, False, True, False, True, True, False, False, False, True]
markers_list = ['o', 's', '^', 'd', 'p', 'o', 's', '^', 's', 'o']

# coefficients to be calculated per country and method:
coeff_types = ['rp', 'rs', 'rmse', 'rmsf']

# initiating variables...
cc = 0
all_coeffs = list()
rmsf_coeffs = list()
slope_coeffs = list()
pval_coeffs = list()
for meth in methods:
    rmsf_coeffs.append('rmsf_' + meth)
    slope_coeffs.append('slope_' + meth)
    pval_coeffs.append('pval_' + meth)
    for coeff in coeff_types:
        all_coeffs.append(coeff + '_' + meth)
        cc = cc + 1

# indicies of each skill metric:
rp_i = np.arange(0, cc, len(coeff_types))
rs_i = np.arange(1, cc, len(coeff_types))
rmse_i = np.arange(2, cc, len(coeff_types))
rmsf_i = np.arange(3, cc, len(coeff_types))

colors3 = ['#1b9e77', '#7570b3', '#d95f02']
c3_10 = [0, 0, 0, 0, 0, 1, 1, 2, 2, 2]

income_groups = list()
for cntry in countries:
    income_groups.append(income_group(cntry, 2016)[1])


if compute_validation:
    rho = dict()
    adm0 = dict()
    adm1 = dict()

    for i in countries_sel:
        print('*** ' + countries[i] + ' *** ')
        start_time_c = time.time()
        rho[countries[i]], adm0[countries[i]], adm1[countries[i]] =\
            lp.admin1_validation(countries[i], methods, exponents_list, \
                                 res_arcsec=resolution, check_plot=False)

        plt.figure() # Scatter plot per country
        lit3_scatter = plt.scatter(adm1[countries[i]]['Lit3'], \
                                   adm0[countries[i]]['Lit3'], c=colors3[0], marker='^')
        pop_scatter = plt.scatter(adm1[countries[i]]['Pop'], \
                                  adm0[countries[i]]['Pop'], c=colors3[1])
        litpop_scatter = plt.scatter(adm1[countries[i]]['LitPop'], \
                                     adm0[countries[i]]['LitPop'], c=colors3[2])
        plt.plot([0, np.max([plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]])],
                 [0, np.max([plt.gca().get_xlim()[1], plt.gca().get_ylim()[1]])],\
                 ls="--", c=".3")

        plt.legend((litpop_scatter, lit3_scatter, pop_scatter),\
                   (r'$LitPop$', r'$Lit^3$', r'$Pop$',))
        plt.xlabel('Reference nGRP')
        plt.ylabel('Modelled nGRP')

        plt.savefig(os.path.join(output_path, experiment_name + '_' + countries[i] + str(resolution) + '_.pdf'), \
                    dpi=600, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format='pdf',
                    transparent=False, bbox_inches=None, pad_inches=0.1,
                    frameon=None, metadata=None)
        plt.show()

        df = pd.DataFrame(adm0[countries[i]])
        df.to_csv(os.path.join(output_path, experiment_name + '_' + countries[i] + \
                               str(resolution) + '_adm0.csv'))
        df = pd.DataFrame(adm1[countries[i]])
        df.to_csv(os.path.join(output_path, experiment_name + '_' +  countries[i] + \
                               str(resolution) + '_adm1_ref.csv'))

        df_r = pd.DataFrame(rho)
        df_r['COEFF'] = all_coeffs
        cols = df_r.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_r = df_r[cols]
        # df_r.to_csv(os.path.join(output_path, experiment_name + '_' + countries[i] + str(resolution) + \
        #                          '_corr_coeffs_zwschnspchr.csv'))
        print('Computing admin1-validation for ' + countries[i] + ' took ' + \
              str(round(time.time()-start_time_c, 2)) + 's')

    df_r = pd.DataFrame(rho)
    df_r['COEFF'] = all_coeffs
    cols = df_r.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df_r = df_r[cols]
    df_r.to_csv(os.path.join(output_path, experiment_name + '_' + str(resolution) + '_corr_coeffs.csv'))

def plot_skillpercountry(data_df, **args):
    """
    Make plot of skill scores with countries on x-axis, methods in legend
    """

    name = args.get('name', 'VARIABLE NAME')
    idx = args.get('idx', data_df.index.values)
    dd = args.get('dd', 5.8) # 3.3
    wdth = args.get('wdth', 8) # 7
    hght = args.get('hght', 4)
    markersize = 60
    target_y = args.get('target_y', 1)
    label_y = args.get('label_y', r'$\rho$')
    meth_labels = [r'$Lit$', r'$Lit^2$', r'$Lit^3$', r'$Lit^4$', r'$Lit^5$', \
                   r'$Pop$', r'$Pop^2$', r'$Lit^3Pop$', r'$Lit^2Pop$', r'$LitPop$']

    plt.figure(facecolor='w', figsize=(wdth, hght))

    for i in np.arange(0, len(methods_show)):
        if not methods_show[i]:
            markers_list[i] = ''
        else:
            plt.scatter([], [], marker=markers_list[i], lw=1, c=colors3[c3_10[i]], \
                        s=markersize, edgecolor='black', linewidth='.4', label=meth_labels[i])
    plt.legend()
    # legendspace:
    plt.scatter([0, len(idx)+dd], [0.7, 0.7], marker='.', lw=1, c='white')

    # actual plotting:
    for i in countries_sel: # country
        for j in np.arange(0, len(idx)):
            # rp - pearson correlation:
            plt.scatter([i], data_df[countries[i]][idx[j]], marker=markers_list[j], \
                        c=colors3[c3_10[j]],\
                        s=markersize, edgecolor='black', linewidth='.5',\
                        alpha=1., zorder=j+10)
    if not target_y == 'none':
        plt.plot([0, i], [target_y, target_y], c='#d3d3d3', lw=5, ls='-', zorder=1)

    plt.xticks(countries_sel, [countries[i] for i in countries_sel], color='black')
    plt.grid(axis='y')
    plt.xlabel('Country')
    plt.ylabel(label_y)
    plt.title(name)

    plt.savefig(os.path.join(output_path, experiment_name + '_' + 'allcountries_v4_' + name + '.pdf'),\
                dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format='pdf',
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    plt.show()

def plot_countryperskill(data_df, **args):
    """
    Make plot of skill scores with method on x-axis, countries in legend
    """
    name = args.get('name', 'VARIABLE NAME')
    idx = args.get('idx', data_df.index.values)
    order = args.get('order', np.array([9, 0, 1, 2, 3, 4, 5, 6, 8, 7], int))
    dd = args.get('dd', .7) # 3.3
    wdth = args.get('wdth', 8) # 7
    hght = args.get('hght', 4)
    markersize = 60
    target_y = args.get('target_y', 1)
    label_y = args.get('label_y', r'$\rho$')
    colors14 = args.get('colors14', ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', \
                                     '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', \
                                     '#cab2d6', '#6a3d9a', '#ffff99', '#b15928', \
                                     '#dd1c77', '#8dd3c7'])
    plt.figure(facecolor='w', figsize=(wdth, hght))
    meth_labels = [r'$Lit$', r'$Lit^2$', r'$Lit^3$', r'$Lit^4$', r'$Lit^5$', \
                   r'$Pop$', r'$Pop^2$', r'$Lit^3Pop$', r'$Lit^2Pop$', r'$LitPop$']
    idx = idx[order]
    meth_labels = [meth_labels[i] for i in order]
    # empty plots for legend handlers:
    for i in np.arange(0, len(countries_sel)): # country
        plt.scatter([], [], marker='o', s=markersize, edgecolor='black', linewidth='.4',\
                    c=colors14[i], label=countries[countries_sel[i]])
    plt.legend()

    plt.scatter([0, len(idx)+dd], [0.7, 0.7], marker='.', lw=1, c='white') # legendspace

    # actual plotting:
    for i in np.arange(0, len(countries_sel)): # country
        for j in np.arange(0, len(idx)):
            # rp - pearson correlation:
            plt.scatter([j], data_df[countries[countries_sel[i]]][idx[j]], marker='o', \
                        s=markersize, edgecolor='black', linewidth='.4',\
                        alpha=1., c=colors14[i], zorder=j+10)
    if not target_y == 'none':
        plt.plot([0, j], [target_y, target_y], c='#d3d3d3', lw=5, ls='-', zorder=1)

    plt.xticks(np.arange(0, len(idx)), meth_labels, color='black', rotation=30)
    plt.grid(axis='y')
    # plt.xlabel('Method')
    plt.ylabel(label_y)
    plt.title(name)

    plt.savefig(os.path.join(output_path, experiment_name + '_' + 'allcountries_perScore_v4_' + name + '.pdf'),\
                dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format='pdf',
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    plt.show()

def boxplot_skillpermethod(data_df, **args):
    """
    Make boxplot of skill scores with method on x-axis
    """
    name = args.get('name', 'VARIABLE NAME')
    idx = args.get('idx', data_df.index.values)
    order = args.get('order', np.array([9, 0, 1, 2, 3, 4, 5, 6, 8, 7], int))
    # dd = args.get('dd', .7) # 3.3
    wdth = args.get('wdth', 6) # 7
    hght = args.get('hght', 3.5)
    target_y = args.get('target_y', 1)
    label_y = args.get('label_y', r'$\rho$')
    ticks_y = args.get('ticks_y', None)

    meth_labels = [r'$Lit^1$', r'$Lit^2$', r'$Lit^3$', r'$Lit^4$', r'$Lit^5$', \
                   r'$Pop^1$', r'$Pop^2$', r'$Lit^3Pop^1$', r'$Lit^2Pop^1$', r'$Lit^1Pop^1$']
    idx = idx[order]
    meth_labels = [meth_labels[i] for i in order]
    data_df = data_df.set_index('COEFF')

    f_h = plt.figure(facecolor='w', figsize=(wdth, hght))
    ax_h = f_h.add_subplot(1,1,1)

    if not target_y == 'none':
        ax_h.plot([0, len(idx)+1], [target_y, target_y], c='black', alpha=.25, lw=3, ls='-', zorder=1)

    data_df.iloc[np.array(idx)].T.boxplot(ax=ax_h)
    plt.xticks(np.arange(1, len(idx)+1), meth_labels, color='black', rotation=30)
    if not ticks_y is None:
        ax_h.set_yscale('log')
        ax_h.yaxis.set_ticks(ticks_y)
        ax_h.yaxis.set_ticklabels( ['%1.1f' % i for i in ticks_y] )

    plt.grid(axis='x')
    plt.ylabel(label_y)
    # plt.title(name)
    f_h.tight_layout()

    f_h.savefig(os.path.join(output_path, experiment_name + '_' + 'allcountries_BOX_v4_' + name + '.pdf'),\
                dpi=600, facecolor='w', edgecolor='w',
                orientation='portrait', papertype=None, format='pdf',
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None, metadata=None)
    f_h.show()



if validation_plots:

    adm1_gdp_share_all = pd.DataFrame()
    # load GDP-share for each country and combine into 1 dataframe:
    rmsf = dict()
    slopes = dict()
    p_vals = dict()
    for i in countries_sel:
        adm1_gdp_share = pd.read_csv(os.path.join(output_path, experiment_name + '_' + countries[i] + \
                                    str(resolution) + '_adm0.csv'), index_col=0)
        adm1_reference = pd.read_csv(os.path.join(output_path, experiment_name + '_' + countries[i] + \
                                    str(resolution) + '_adm1_ref.csv'), index_col=0)

        adm1_gdp_share['Reference'] = adm1_reference['LitPop']
        adm1_gdp_share = adm1_gdp_share[adm1_gdp_share > 1e-12]
        adm1_gdp_share['country_num'] = i
        adm1_gdp_share['country'] = countries[i]
        adm1_gdp_share_all = pd.concat([adm1_gdp_share_all, adm1_gdp_share])


    # compute RMSF (Root mean squared fraction) for each method:
        rmsf[countries[i]] = list()
        slopes[countries[i]] = list()
        p_vals[countries[i]] = list()
        for i_meth in np.arange(0, len(methods)):
            rmsf[countries[i]].append(np.exp(np.sqrt(np.sum( \
                (np.log(adm1_gdp_share[methods[i_meth]]/ \
                adm1_gdp_share['Reference']))**2)/adm1_gdp_share.shape[0])))

            val1 = ~np.isnan(adm1_gdp_share['Reference'])
            val2 = ~np.isnan(adm1_gdp_share[methods[i_meth]])
            slopes[countries[i]].append(stats.linregress( \
                  adm1_gdp_share['Reference'][val1 & val2],\
                  adm1_gdp_share[methods[i_meth]][val1 & val2])[0])
            p_vals[countries[i]].append(stats.linregress( \
                  adm1_gdp_share['Reference'][val1 & val2],\
                  adm1_gdp_share[methods[i_meth]][val1 & val2])[3])
            #rmsf = np.exp(np.sqrt(np.sum((np.log(adm0_data/adm1_data))**2)/ \
            #                            adm0_data.shape[0]))

    df_rmsf = pd.DataFrame(rmsf)
    df_rmsf['COEFF'] = rmsf_coeffs
    df_rmsf.to_csv(os.path.join(output_path, experiment_name + '_' + str(resolution) + '_rmsf.csv'))
    df_slope = pd.DataFrame(slopes)
    df_slope['COEFF'] = slope_coeffs
    df_slope.to_csv(os.path.join(output_path, experiment_name + '_' + str(resolution) + '_slope.csv'))
    df_p = pd.DataFrame(p_vals)
    df_p['COEFF'] = pval_coeffs
    df_p.to_csv(os.path.join(output_path, experiment_name + '_' + str(resolution) + '_p_val.csv'))

    r = pd.read_csv(os.path.join(output_path, experiment_name + '_' + str(resolution) + '_corr_coeffs.csv'))
    r.__delitem__('Unnamed: 0')
#    r['all countries']=rho_all_df['ALL']

    # PLOTTING:
    order=np.array([9, 0, 1, 2, 3, 4, 5, 6, 8, 7], int)
    plot_skillpercountry(r, idx=rp_i, name='Pearson Correlation', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\rho$', \
                         methods_show=methods_show)
    plot_skillpercountry(df_slope, name='Linear Regression Slope', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\beta$', methods_show=methods_show)
    plot_skillpercountry(r, idx=rmsf_i, name='Root Mean Squared Fraction', countries=countries, \
                         countries_sel=countries_sel, label_y='RMSF', methods_show=methods_show)

    plot_countryperskill(r, idx=rp_i, order=order, name='Pearson Correlation', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\rho$', \
                         methods_show=methods_show)
    plot_countryperskill(df_slope, order=order, name='Linear Regression Slope', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\beta$', methods_show=methods_show)
    plot_countryperskill(r, idx=rmsf_i, order=order, name='Root Mean Squared Fraction', countries=countries, \
                         countries_sel=countries_sel, label_y='RMSF', methods_show=methods_show)
    boxplot_skillpermethod(r, idx=rp_i, order=order, name='Pearson Correlation', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\rho$', \
                         methods_show=methods_show)
    boxplot_skillpermethod(df_slope, order=order, name='Linear Regression Slope', countries=countries, \
                         countries_sel=countries_sel, label_y=r'$\beta$', methods_show=methods_show)
    boxplot_skillpermethod(r, idx=rmsf_i, order=order, name='Root Mean Squared Fraction', countries=countries, \
                         countries_sel=countries_sel, label_y='RMSF', methods_show=methods_show, \
                         ticks_y = np.array([1, 2, 5, 10, 20]))

    # Rearrange Data frames for nice CSV output
    methods = [methods[i] for i in order]
    countries_ = [countries[i] for i in countries_sel]
    r_rp = r[countries_].iloc[rp_i[order]]
    r_rp['COEFF'] = r['COEFF'].iloc[rp_i[order]]

    r_slope = df_slope[countries_].iloc[order]
    r_slope['COEFF'] = df_slope['COEFF'].iloc[order]

    r_rmsf = r[countries_].iloc[rmsf_i[order]]
    r_rmsf['COEFF'] = r['COEFF'].iloc[rmsf_i[order]]

    cols = r_rp.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    r_rp = r_rp[cols]
    r_rp.set_index('COEFF')
    r_slope = r_slope[cols]
    r_slope.set_index('COEFF')
    r_rmsf = r_rmsf[cols]
    r_rmsf.set_index('COEFF')
    r_rp = r_rp.reset_index(drop=True)
    r_slope = r_slope.reset_index(drop=True)
    r_rmsf = r_rmsf.reset_index(drop=True)

    #Create a DataFrame
    statistics_ = {'Method':methods}
    r_stat = pd.DataFrame(statistics_)
    r_stat['rp_median'] = r_rp.quantile(q=0.5, axis=1, numeric_only=True, interpolation='linear')
    r_stat['rp_IQR'] = r_rp.quantile(q=0.75, axis=1, numeric_only=True, interpolation='linear')\
                      -r_rp.quantile(q=0.25, axis=1, numeric_only=True, interpolation='linear')
    r_stat['slope_median'] = r_slope.quantile(q=0.5, axis=1, numeric_only=True, interpolation='linear')
    r_stat['slope_IQR'] = r_slope.quantile(q=0.75, axis=1, numeric_only=True, interpolation='linear')\
                         -r_slope.quantile(q=0.25, axis=1, numeric_only=True, interpolation='linear')
    r_stat['rmsf_median'] = r_rmsf.quantile(q=0.5, axis=1, numeric_only=True, interpolation='linear')
    r_stat['rmsf_IQR'] = r_rmsf.quantile(q=0.75, axis=1, numeric_only=True, interpolation='linear')\
                        -r_rmsf.quantile(q=0.25, axis=1, numeric_only=True, interpolation='linear')
    r_stat = r_stat.round(2)
    r_rp = r_rp.round(2)
    r_slope = r_slope.round(2)
    r_rmsf = r_rmsf.round(2)

    # save to CSV:
    r_stat.to_csv(os.path.join(output_path, experiment_name + '_STAT_' + str(resolution) + '.csv'))
    r_rp.to_csv(os.path.join(output_path, experiment_name + '_RP_' + str(resolution) + '.csv'))
    r_slope.to_csv(os.path.join(output_path, experiment_name + '_SLOPE_' + str(resolution) + '.csv'))
    r_rmsf.to_csv(os.path.join(output_path, experiment_name + '_RMSF_' + str(resolution) + '.csv'))
