#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 19:34:03 2020

This script test damage time series on the influence of teleconnections and
GMT. Tested are the teleconnections with ENSO, NAO, PDO and the effect of GMT or AMO.
This script is for disaggregated regions.

@author: insauer
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import itertools
import pymannkendall as mk
from scipy.stats import shapiro
import numpy.ma as ma
import matplotlib.pyplot as plt
DATA_TS = pd.read_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/AttributionTimeSeriesSubregions.csv')
DATA_ATTR = pd.read_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/AttributionMetaDataSubregions.csv')

teleData = pd.read_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/InputData/teleconnections_lag.csv')

output_path = '/home/insauer/projects/NC_Submission/Climada_papers/Test/Lag_ENSO_GMT_PDO_NAO_Loo_Subregions.csv'

normClim = True

teleData = teleData.loc[teleData['Year'] >= 1971]
teleData80 = teleData.loc[teleData['Year'] >= 1980]
telecon80 = teleData80[['ENSO', 'ENSO_lag', 'GMT', 'PDO', 'PDO_lag', 'NAO', 'NAO_lag']]

telecon = teleData[['ENSO', 'ENSO_lag', 'GMT', 'PDO', 'PDO_lag', 'NAO', 'NAO_lag']]

region_names = {'SSAF': 'SSA + Southern Africa',
                'EUR': 'Western Europe',
                'GLB': 'Global',
                'LAM': 'Central America',
                'NAFARA': 'North Africa + Middle East',
                'CAS': 'Central Asia + Eastern Europe',
                'SWEA': 'Southern Asia + South-East Asia',
                'CHN': 'Eastern Asia',
                'AUS': 'Oceania',
                'NAM': 'North America'
                }

predictors = ['ENSO', 'ENSO_lag', 'GMT', 'PDO', 'PDO_lag', 'NAO', 'NAO_lag']

link_fnc_list = [sm.families.links.log(), sm.families.links.identity(),
                 sm.families.links.inverse_power()]

test_regions = list(region_names)


def get_pearson(pred, climdat):
    """
    pearson correlation of model predicted data and damage time series

    Parameters
    ----------
    pred : GLM
        model
    climdat : np.array
        damage time series

    Returns
    -------
    float
        Pearson correlation coefficient

    """

    a = ma.masked_invalid(climdat)
    b = ma.masked_invalid(pred.predict())
    msk = (~a.mask & ~b.mask)
    corrcoef = ma.corrcoef(a[msk], b[msk])

    return corrcoef[0, 1]


def looCV(clim, predic, fnc):

    err = 0
    for lo_index in range(len(clim)):

        clim_mask = np.ma.array(clim, mask=False)
        clim_mask.mask[lo_index] = True
        clim_lo = clim_mask.compressed()

        predic_lo = predic.reset_index().drop(lo_index, axis=0).drop('index', 1)

        model_res = sm.GLM(clim_lo, predic_lo,
                           family=sm.families.Gamma(fnc)).fit(maxiter=5000, scale=1.)

        value_pred = model_res.predict(predic).iloc[lo_index]

        err = err + (clim[lo_index] - value_pred)**2

    return err/len(clim)


def pred_double(comb):

    pred_names = ['GMT', 'ENSO', 'NAO', 'PDO']

    for p in pred_names:
        if p in comb:
            if p + '_lag' in comb:
                return True

    return False


def find_best_model(climateDat, telecon):
    """
    Wrapper function to select the best model. Function provides all possible
    combination of predictors and a constant and evaluates the model applying
    the LooCV. It selects the model with the smallest out-of sample error.

    Parameters
    ----------
    climateDat : np.array
        damage time series
    telecon : DataFrame
        teleconnections and GMT

    Returns
    -------
    best_model: GLMObject
        full GLM object of the best model
    best_model_indices[0]: int
        x-index of the best model (indicates combination of predictors)
    best_model_indices[1]: int
        y-index of the best model (indicates link function)
    iter_max: int
        number iterations needed for convergence of the best model
    pearson_corr: float
        pearson correlation of model and data
    best_loo: float
        out-of-sample-error of the best model
    looICs_lf: np.array
        all out-off-sample errors
    """

    max_i = 5000

    if test_region == 'AUS':

        climateDat = np.nan_to_num(climateDat)

    models_lf = []
    deviances_lf = []
    chi2_lf = []
    iter_max = 0
    looICs_lf = []
    comb_list_lf = []
    for link_fnc in link_fnc_list:
        models = []
        deviances = []
        chi2 = []
        looICs = []
        comb_list = []
        for n_preds in range(0, 5):
            for comb in list(itertools.combinations(predictors, n_preds)):
                print(list(comb))
                if pred_double(comb):
                    print('skip combination')
                    continue
                data_exog = sm.add_constant(telecon[list(comb)])
                try:

                    model_result = sm.GLM(climateDat, data_exog,
                                          family=sm.families.Gamma(link_fnc)).fit(maxiter=max_i,
                                                                                  scale=1.)
                    looIC = looCV(climateDat, data_exog, link_fnc)

                    models.append(model_result)
                    looICs.append(looIC)
                    deviances.append(model_result.aic)
                    chi2.append(model_result.pearson_chi2)
                    comb_list.append(comb)
                except ValueError:

                    models.append(sm.GLM(data_exog, np.ones(len(climateDat))))
                    deviances.append(1e10)
                    looICs.append(1e10)
                    chi2.append(1e10)
                    comb_list.append(comb)
                if model_result.fit_history['iteration'] == max_i:
                    iter_max += 1

                if n_preds == 4:
                    print('stop')

        looICs_lf.append(looICs)
        models_lf.append(models)
        deviances_lf.append(deviances)
        chi2_lf.append(chi2)
        comb_list_lf.append(comb_list)

    best_model_indices = np.array(np.unravel_index(np.argmin(np.array(looICs_lf), axis=None),
                                                   np.array(looICs_lf).shape))

    best_model = models_lf[best_model_indices[0]][best_model_indices[1]]

    best_loo = looICs_lf[best_model_indices[0]][best_model_indices[1]]

    pearson_corr = get_pearson(best_model, climateDat)

    return best_model, best_model_indices[0], best_model_indices[1],\
        iter_max, pearson_corr, best_loo, looICs_lf, comb_list_lf


def test_residuals(model, timeperiod,reg,dis):
    """
    Test for a residual trend, applying a Mann-Kendall-test

    Parameters
    ----------
    model : GLMObject
        Best model
    timeperiod : np.array
        considered years (not used here)

    Returns
    -------
    float
        slope in residuals
    float
        p-value

    """
    res_trend = mk.original_test(model.resid_response, alpha=0.1)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    sm.graphics.tsa.plot_acf(model.resid_response, lags=39, ax = ax)
    ax.set_xlabel('lag')
    ax.set_title('Autocorrelation {}_{}'.format(reg,dis))
    #fig.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Response_letter/Autocorr_Residuals/AutocorrResidualsGMT_{}_{}.png'.format(reg,dis),bbox_inches = 'tight',dpi =600)

    alt_trend_test = mk.hamed_rao_modification_test(model.resid_response)

    return res_trend.slope, res_trend.p, alt_trend_test.trend, alt_trend_test.p


def test_autocorrelation(time_series):
    """
    Test for a residual trend, applying a Mann-Kendall-test

    Parameters
    ----------
    model : GLMObject
        Best model
    timeperiod : np.array
        considered years (not used here)

    Returns
    -------
    float
        slope in residuals
    float
        p-value

    """
    auto = mk.original_test(time_series, alpha=0.1)

    return auto.Tau


def extract_model_coefs(region, model, link, model10, link10,
                        model80, link80, model8010, link8010, dis):
    """
    Reads the coefficients and p-values for each predictor of the best model
    and saves data in a csv file. To achieve comparability between coefficients of
    models with different link functions the partial devaritive in a centric
    development point is calculated.

    Parameters
    ----------
    region : string
        region abbreviation
    model : GLMObject
        best model (1971-2010, 1980 fixed exposure)
    link : int
        index of link function (1971-2010, 1980 fixed exposure)
    model10 : GLMObject
        best model (1971-2010, 2010 fixed exposure)
    link10 : int
        index of link function (1971-2010, 2010 fixed exposure)
    model80 : GLMObject
        best model (1980-2010, 1980 fixed exposure)
    link80 : int
        index of link function (1980-2010, 1980 fixed exposure)
    model8010 : GLMObject
        best model (1980-2010, 2010 fixed exposure)
    link8010 : int
        index of link function (1980-2010, 1980 fixed exposure)
    dis : string
        discharge group
    """
    shortage = ['', '10', '80', '8010']
    dev_index = [20, 20, 15, 15]
    mods = [model,  model10, model80, model8010]
    coefs = ['ENSO', 'AMO', 'PDO', 'NAO', 'GMT']
    lnks = [link,  link10, link80, link8010]

    for m, mod in enumerate(mods):

        coeff_sum = 0

        for c, coef in enumerate(coefs):
            try:
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'_'+shortage[m]] = mod.params[coef]
                coef_deriv = link_fnc_list[lnks[m]].inverse_deriv(teleData[coef])\
                    * mod.params[coef]
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'dv_'+shortage[m]] = np.array(coef_deriv)[dev_index[m]]
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'pval_'+shortage[m]] = mod.pvalues[coef]

                coeff_sum += mod.params[coef]
            except KeyError:
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'_'+shortage[m]] = np.nan
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'pval_'+shortage[m]] = np.nan
                DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                              coef+'dv_'+shortage[m]] = np.nan
        DATA_ATTR.loc[DATA_ATTR['Region'] == region+'_'+dis,
                      'CoefSum_'+shortage[m]] = coeff_sum

    DATA_ATTR.to_csv(output_path)


def unexplainable_Trends(pval, slope, test_region, change, sign):
    """
    Check for the presence of an unexplainable trend after adjustment for
    teleconnections

    Parameters
    ----------
    pval : float
        p-value of residual trend
    slope : float
        slope of residual trend
    test_region : string
        region
    change : slope
        slope of trend in damages
    sign : float
        significance of slope in damages

    Returns
    -------
    bool
    """

    if pval > 0.1:
        return False

    haz_slope = DATA_ATTR.loc[DATA_ATTR['Region'] == test_region, change].sum()

    if (haz_slope < 0) and (slope > 0):
        return False

    if (haz_slope > 0) and (slope < 0):
        return False

    return True


for i, test_region in enumerate(test_regions):
    
    print(test_region)
    DATA_region = DATA_TS[(DATA_TS['Region'] == test_region) &
                          (DATA_TS['Year'] < 2011) & (DATA_TS['Year'] > 1970)]
    # DATA_region = DATA_region.reset_index()

    DATA_region80 = DATA_TS[(DATA_TS['Region'] == test_region) &
                            (DATA_TS['Year'] < 2011) & (DATA_TS['Year'] > 1979)]
    # DATA_region80 = DATA_region80.reset_index()

    if test_region != 'NAM':
        climateDataPos = np.array(DATA_region['NormHaz_ImpFix_2yPos_offset'])
        auto_corrPos = test_autocorrelation(climateDataPos)
        
        if normClim is True:
            climateDataPos = climateDataPos/np.nanmax(climateDataPos)
            
        t, shap_logPos = shapiro(np.log(climateDataPos))
        t, shap_normPos = shapiro(climateDataPos)
        
        best_modelPos, yPos, xPos, maxiPos, pearson_corrPos,\
            best_looPos, loosPos, combPos = find_best_model(climateDataPos, telecon)
        
        comb_df_pos = pd.DataFrame(combPos)
        comb_df_pos = comb_df_pos.T
        
        loo_dfPos = pd.DataFrame(loosPos)
        loo_df_pos = loo_dfPos.T
        loo_df_pos.columns = ['log', 'identity', 'inverse-power']
        loo_df_pos['combination'] = comb_df_pos.iloc[:, 0]
        
        extract_model_coefs(test_region,  best_modelPos, yPos,  'Pos')
        coefPos, pvalPos, alt_trendPos, alt_pvalPos = test_residuals(best_modelPos, np.arange(1971, 2011),test_region,'Pos')
        
        unexTPos = unexplainable_Trends(pvalPos, coefPos, test_region+'_Pos', 'Change H7', 'Sign H7')
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Res_Sig'] = pvalPos
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Res_Slope'] = coefPos
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Alt_Sig'] = alt_pvalPos
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Alt_Slope'] = alt_trendPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Unexplained Haz'] = unexTPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'CorrCoef'] = pearson_corrPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'BestLoo'] = best_looPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Best_link'] =\
        best_modelPos.fit_history['iteration']
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Link_func'] = yPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Combi'] = xPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Norm_dist'] = shap_normPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'LogNorm_dist'] = shap_logPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Maxi'] = maxiPos
        
        DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Pos', 'Autocorr'] = auto_corrPos
    
    
    climateDataNeg = np.array(DATA_region['NormHaz_ImpFix_2yNeg_offset'])
    auto_corrNeg = test_autocorrelation(climateDataNeg)

    if normClim is True:

        climateDataNeg = climateDataNeg/np.nanmax(climateDataNeg)

    t, shap_logNeg = shapiro(np.log(climateDataNeg))
    t, shap_normNeg = shapiro(climateDataNeg)

    best_modelNeg, yNeg, xNeg, maxiNeg, pearson_corrNeg,\
        best_looNeg, loosNeg, combNeg = find_best_model(climateDataNeg, telecon)

    comb_df_neg = pd.DataFrame(combNeg)
    comb_df_neg = comb_df_neg.T
    # store all out-of-sample-errors
    loo_df_pos.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/LooIC_ENSO_GMT_PDO_NAO_Subregions_{}Pos.csv'.format(test_region))

    loo_dfNeg = pd.DataFrame(loosNeg)
    loo_df_neg = loo_dfNeg.T
    loo_df_neg.columns = ['log', 'identity', 'inverse-power']
    loo_df_neg['combination'] = comb_df_neg.iloc[:, 0]
    # store all out-of-sample-errors
    loo_df_neg.to_csv('/home/insauer/projects/NC_Submission/Climada_papers/Test/LooIC_ENSO_GMT_PDO_NAO_Subregions_{}Neg.csv'.format(test_region))

    extract_model_coefs(test_region,  best_modelNeg,  yNeg, 'Neg')

    coefNeg, pvalNeg, alt_trendNeg, alt_pvalNeg = test_residuals(best_modelNeg, np.arange(1971, 2011),test_region,'Neg')

    unexTNeg = unexplainable_Trends(pvalNeg, coefNeg, test_region+'_Neg', 'Change H7', 'Sign H7')

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Res_Sig'] = pvalNeg
    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Alt_Sig'] = alt_pvalNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Res_Slope'] = coefNeg
    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Alt_Slope'] =  alt_trendNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Unexplained Haz'] = unexTNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'CorrCoef'] = pearson_corrNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'BestLoo'] = best_looNeg
    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Autocorr'] = auto_corrNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Best_link'] =\
        best_modelNeg.fit_history['iteration']

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Link_func'] = yNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Combi'] = xNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Norm_dist'] = shap_normNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'LogNorm_dist'] = shap_logNeg

    DATA_ATTR.loc[DATA_ATTR['Region'] == test_region+'_Neg', 'Maxi'] = maxiNeg

    DATA_ATTR.to_csv(output_path)


DATA_ATTR.to_csv(output_path)
