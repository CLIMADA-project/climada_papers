"""
@author: Samuel Luethi WCR ETHZ
@date: 05.07.2021

BiasCorrection

"""

import numpy as np
import pandas as pd


"""Contains functions for bias correction.

Bias correction is performed using a quantile mapping approach.
This code is a slightly simplified version of the method published by
Rajczak et al. (2016). doi:10.1002/joc.4417 and is available in R under
https://github.com/SvenKotlarski/qmCH2018

"""

def bias_correct_time_series(dat_mod, dat_obs, dat_mod_all,
                             minq=0.001, maxq=1.000, incq=0.001):
    """ Wrapper function for bias correction. Calculates quantiles for mapping,
        estimates cumulative distribution function (CDF) of observational
        (dat_obs) and modelled (dat_mod) data to estimate quantile specific
        correction function. CDF of these two timeseries need to correspond to
        the same time period. The estimated correction function is then applied
        to dat_mod_all which is model output data the same model than dat_mod
        but can cover any time range.

        Parameters
        ----------
        dat_mod : pd.Series
            Model data series as reference for bias adjustment
        dat_obs : pd.Series
            Observational data series as ground truth
        dat_mod_all : pd.Series
            Data series to be bias adjusted
        minq : float
            Minimum quantile for correction function
        maxq: float
            Maximum quantile for correction function
        incq : float
            Quantile increment for correction function (bin size)

        Returns
        -------
        dat_mod_all_corrected : pd.Series
            bias corrected dat_mod_all

        """
    # define quantiles used for mapping
    q = np.arange(minq, maxq, incq)
    # (1) calc cdf of observational and modeled data
    cdf_obs = _calc_cdf(dat_obs, q)
    cdf_mod = _calc_cdf(dat_mod, q)

    # (2) estimate correction function
    cdf_dif = cdf_mod - cdf_obs

    # (3) perform quantile mapping to data
    dat_mod_all_corrected = _map_quantile(dat_mod_all, cdf_mod, cdf_dif, q)

    return dat_mod_all_corrected

def bias_correct_ensemble(dat_mod, dat_obs, dat_mod_all,
                          minq=0.001, maxq=1.000, incq=0.001):
    """ Wrapper function for bias correction of a large ensemble.
         Calculates quantiles for mapping, estimates cumulative distribution
         function (CDF) of observational (dat_obs) and modelled (dat_mod) data
         to estimate one single quantile specific correction function for the
         whole ensemble. CDF of the observational data and the ensemble
         DataFrame need to correspond to the same time period. The estimated
         correction function is then applied to dat_mod_all which is ensemble
         model output data of the same model than dat_mod but can cover
         any time range.

        Parameters
        ----------
        dat_mod : pd.DataFrame
            DataFrame with climate data from large ensemble. Needs to cover
            same range as dat_obs
        dat_obs : pd.Series
            Observational data series as ground truth. Needs to cover same
            range as dat_mod
        dat_mod_all : pd.DataFrame
            DataFrame to be bias adjusted
        minq : float
            Minimum quantile for correction function
        maxq: float
            Maximum quantile for correction function
        incq : float
            Quantile increment for correction function (bin size)

        Returns
        -------
        dat_mod_all_corrected : pd.DataFrame
            bias corrected dat_mod_all

        """
    # define quantiles used for mapping
    q = np.arange(minq, maxq, incq)
    # (1) calc cdf of observational and modeled data
    cdf_obs = _calc_cdf(dat_obs, q)
    cdf_mod = _calc_cdf_ens(dat_mod, q)

    # (2) estimate correction function
    cdf_dif = cdf_mod - cdf_obs

    # (3) perform quantile mapping to data
    dat_mod_corrected = _map_quantile_ens(dat_mod_all, cdf_mod, cdf_dif, q)

    return dat_mod_corrected

def _calc_cdf(data_series, q):
    """ Calculates cumulative distribution function (CDF) of any time series.
        Takes no assumption on distribution.

        Parameters
        ----------
        data_series : pd.Series
            Data series
        q : np.array
            quantiles of cdf to be calculated

        Returns
        -------
        cdf : np.array
            cdf of data_series on quantiles q

        """

    # sort data
    dat_sorted = np.sort(data_series.values)
    # calculate the proportional values of samples
    p = 1. * np.arange(len(data_series)) / (len(data_series) - 1)
    # map to percentiles
    cdf = np.interp(q, p, dat_sorted)

    return cdf

def _map_quantile(dat_mod_all, cdf_mod_orig, cdf_dif, q):
    """ Performs bias correction using quantile mapping

        Parameters
        ----------
        dat_mod_all : pd.Series
            Data series to be bias adjusted
        cdf_mod_orig : np.array
            original cdf of model used for bias correction
        cdf_dif : np.array
            cdf correction function
        q : np.array
            quantiles of cdf to be calculated

        Returns
        -------
        dat_mod_adj : pd.Series
            bias corrected data series

        """
    # calc percentile value of each temperature value in modelled time series
    perc_mod = np.interp(dat_mod_all, cdf_mod_orig, q)
    # correction term for each temperature value in modelled time series
    cor_term = np.interp(perc_mod, q, cdf_dif)
    # adjust for bias
    dat_mod_adj = dat_mod_all-cor_term

    return pd.Series(dat_mod_adj)

def _calc_cdf_ens(dat_mod, q):
    """ Calculates cumulative distribution function (CDF) an ensemble.
        Ensemble CDF is calculated as the mean over all CDF of the ensemble
        members.
        Takes no assumption on distribution.

        Parameters
        ----------
        dat_mod : pd.DataFrame
            DataFrame with climate data from large ensemble
        q : np.array
            quantiles of cdf to be calculated

        Returns
        -------
        cdf : np.array
            mean cdf over all ensemble members on quantiles q

        """

    # array to store cdfs
    cdf_array = np.zeros((dat_mod.shape[1], len(q)))
    for i in range(dat_mod.shape[1]):
        # calc cdf per member
        cdf_array[i,:] = _calc_cdf(dat_mod.iloc[:,i], q)

    # average cdf
    cdf_ens = np.mean(cdf_array, axis=0)

    return cdf_ens

def _map_quantile_ens(dat_mod_all, cdf_mod, cdf_dif, q):
    """ Performs bias correction for each ensemble member.

        Parameters
        ----------
        dat_mod_all : pd.DataFrame
            DataFrame to be bias adjusted
        cdf_mod_orig : np.array
            original cdf of model used for bias correction
        cdf_dif : np.array
            cdf correction function
        q : np.array
            quantiles of cdf to be calculated

        Returns
        -------
        dat_mod_adj : pd.DataFrame
            bias corrected data series

        """
    ens_array = np.zeros(dat_mod_all.shape)
    for i in range(dat_mod_all.shape[1]):
        ens_array[:,i] = _map_quantile(dat_mod_all.iloc[:,i], cdf_mod, cdf_dif, q)
    dat_mod_corrected = pd.DataFrame(ens_array)
    dat_mod_corrected.columns = dat_mod_all.columns

    return dat_mod_corrected
