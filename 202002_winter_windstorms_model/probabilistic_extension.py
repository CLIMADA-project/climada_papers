# -*- coding: utf-8 -*-
"""
Created on Fri Feb  28  2020
    
description: creation of all "WISC probabilistic extension" hazard dataset

@author: ThomasRoosli, thomas.roeoesli@usys.ethz.ch
"""




import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pandas as pd
import datetime
import statsmodels.api as sm

from climada.hazard import StormEurope


# All these variables need to be defined correctly for the code to work.
# The strings in these variables should point to existing folders on your computer.

project_folder = 'C:\\Users\\ThomasRoosli\\Documents\\PhD\\WISC_phd\\paper GVZ\\jupyter'
file_identifier = '_v01' # this string is added to all files written by this code


wisc_hist_filename = os.path.join(project_folder, 'WISC_hist' + file_identifier + '.hdf5')
wisc_prob_filename = os.path.join(project_folder, 'WISC_prob_CH' + file_identifier + '.hdf5')


### tryout GEV with 5 year maxima
wisc_hist = StormEurope()
wisc_hist.read_hdf5(wisc_hist_filename)
#get 5 year maxima
block_size_years = 5
quinquenniums = np.arange(0,80,block_size_years) + 1940
ssi_5y_maxima = np.zeros_like(quinquenniums[:-1], dtype=float)
for ind_i, start_year_i in enumerate(quinquenniums[:-1]):
    start_i = datetime.datetime(start_year_i,1,1).toordinal()
    end_i = datetime.datetime(quinquenniums[ind_i+1],1,1).toordinal()
    selection_i = (wisc_hist.date >= start_i) & (wisc_hist.date < end_i)
    ssi_5y_maxima[ind_i] = wisc_hist.ssi[selection_i].max()
    
    
#params_first = using R script 
params_final_gev = scipy.stats.genextreme.fit(ssi_5y_maxima,-0.1154771, loc=21360556255,scale= 6491392290)

return_periods = 3200/np.power(np.sqrt(2),np.arange(21))
frequencies = 1-(1/return_periods) # in [per year]
#quantiles = 1-(1/(return_periods*wisc_hist.frequency.sum())) # in quantiles of the distribution
quantiles_gev = 1-(1/(return_periods*wisc_hist.frequency.sum()/block_size_years)) # in quantiles of the distribution
ssi_quantiles_distribution = scipy.stats.genextreme.ppf(quantiles_gev, params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2])


return_period_hist = 1/np.cumsum(wisc_hist.frequency[np.argsort(-wisc_hist.ssi)])
ssi_quantiles_hist = wisc_hist.ssi[np.argsort(-wisc_hist.ssi)]
plt.figure()
plt.plot(return_periods,ssi_quantiles_distribution)
plt.plot(return_period_hist,ssi_quantiles_hist,'r.')
plt.xlim([1,200])
plt.ylim([1/10,10**11])
plt.xscale('log')


### resampling
number_of_resamplings = 1000
number_of_events = len(ssi_5y_maxima)
resampling_sample = np.empty([number_of_resamplings,number_of_events])
resampling_params = np.empty([number_of_resamplings,3])
resampling_ppf = np.empty([number_of_resamplings,len(return_periods)])
for ind_i in np.arange(number_of_resamplings):
    #random sample
    sample_i = scipy.stats.genextreme.rvs(params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2],size=number_of_events)
    #fit back and save the new parameters and the pdf
    params_i = scipy.stats.genextreme.fit(sample_i,-0.1154771, loc=21360556255,scale= 6491392290)
    resampling_ppf[ind_i,:] = scipy.stats.genextreme.ppf(quantiles_gev, params_i[0],loc=params_i[1],scale=params_i[2])
    resampling_sample[ind_i,:] = sample_i
    resampling_params[ind_i,:] = params_i
    


plt.figure()
plt.fill_between(return_periods,
                 np.quantile(resampling_ppf,[0.05],axis=0)[0],
                 np.quantile(resampling_ppf,[0.95],axis=0)[0],
                 facecolor = 'b',
                 alpha= 0.3)
plt.plot(return_periods,ssi_quantiles_distribution)
plt.plot(return_period_hist,ssi_quantiles_hist,'r.')
#plt.plot(return_period_prob,ssi_quantiles_prob,'b.')
plt.xlim([1,3000])
plt.ylim([0,10**11])
plt.xlabel('Returnperiod [years]')
plt.ylabel('Storm Severity Index []')
plt.title('GEV fit on 5y-maxima')
plt.xscale('log')


### dependence of GEV parameters on blocksize... assymptotic zone... stability
def compare_blocksize(wisc_hist):
    #get 5 year maxima
    blocksize_test = np.arange(4,20)
    blocksize_test_params = np.empty([len(blocksize_test),3])
    for block_ind_i in np.arange(len(blocksize_test)):
        block_size_years = blocksize_test[block_ind_i]
        quinquenniums = np.arange(0,80,block_size_years) + 1940
        ssi_5y_maxima = np.zeros_like(quinquenniums[:-1], dtype=float)
        for ind_i, start_year_i in enumerate(quinquenniums[:-1]):
            start_i = datetime.datetime(start_year_i,1,1).toordinal()
            end_i = datetime.datetime(quinquenniums[ind_i+1],1,1).toordinal()
            selection_i = (wisc_hist.date >= start_i) & (wisc_hist.date < end_i)
            ssi_5y_maxima[ind_i] = wisc_hist.ssi[selection_i].max()
            
            
        #params_first = using R script 
        params_final_gev = scipy.stats.genextreme.fit(ssi_5y_maxima,-0.1154771, loc=21360556255,scale= 6491392290)
        blocksize_test_params[block_ind_i,:] = params_final_gev
    
    plt.figure()
    plt.plot(blocksize_test,blocksize_test_params[:,0])
    plt.xlabel('blocksize')
    plt.ylabel('shape parameter of GEV fit')
    plt.figure()
    plt.plot(blocksize_test,blocksize_test_params[:,1]-blocksize_test_params[:,2]/blocksize_test_params[:,0])
    plt.xlabel('blocksize')
    plt.ylabel('modified scale parameter of GEV fit')
    return None


compare_blocksize(wisc_hist) #-> block size of 5 is ok

### create prob events
power_list = []
scale_list = []
ks_extrapolate = []


for power_i in [1.05,1.1,1.15,1.2,1.5]: #[1.05,1.1,1.25,1.5,2]:
    for scale_i in np.arange(0.0025,0.033,0.0025): #[0.05,0.1,0.15,0.2]:
        #scale_i=np.arange(0.01,0.08,0.005)[8] #0.05
        #power_i=1.2
#        scale_i=np.arange(0.0025,0.026,0.0025)[-1] #0.05
#        power_i=1.1
        hazard_name = ('WISC_prob_CH_p' +
                       str(power_i).replace('.','_') +
                       '_s' +
                       str(scale_i).replace('.','_') +
                       file_identifier +
                       '.hdf5')
        try:
            wisc_prob_CH_i = StormEurope()
            wisc_prob_CH_i.read_hdf5(os.path.join(project_folder,hazard_name))
        except:
            ssi_args = {
                        'threshold': 25,
                        }
            wisc_prob_CH_i = wisc_hist.generate_prob_storms(reg_id = 756, 
                                                          ssi_args=ssi_args,
                                                          scale=scale_i,
                                                          power=power_i)
            wisc_prob_CH_i.write_hdf5(os.path.join(project_folder,'create_final_hazards',hazard_name))
        
        return_period_prob = 1.0/np.cumsum(wisc_prob_CH_i.frequency[np.argsort(-wisc_prob_CH_i.ssi_full_area)])
        ssi_quantiles_prob = wisc_prob_CH_i.ssi_full_area[np.argsort(-wisc_prob_CH_i.ssi_full_area)]
        
        
        plt.figure()
        plt.fill_between(return_periods,
                         np.quantile(resampling_ppf,[0.05],axis=0)[0],
                         np.quantile(resampling_ppf,[0.95],axis=0)[0],
                         facecolor = 'b',
                         alpha= 0.3)
        plt.plot(return_periods,ssi_quantiles_distribution)
        plt.plot(return_period_hist,ssi_quantiles_hist,'r.')
        plt.plot(return_period_prob,ssi_quantiles_prob,'b.')
        plt.xlim([1,3000])
        plt.ylim([0,10**11])
        plt.xscale('log')
        plt.title('power {} and scale {}'.format(power_i,scale_i))
        plt.xlabel('Returnperiod [years]')
        plt.ylabel('Storm Severity Index []')
        
        

        #plot cdfs to get ks test statistic
        cdf_gev_i = scipy.stats.genextreme.cdf(ssi_quantiles_prob, params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2])
        plt.figure()
        plt.plot(ssi_quantiles_prob,cdf_gev_i,'.k')
#        plt.plot(ssi_quantiles_distribution,(1-1/return_periods)**5,'.b')
        cdf_wisc_prob_i = (1-1/return_period_prob)**block_size_years
        plt.plot(ssi_quantiles_prob,cdf_wisc_prob_i,'.b')
        plt.title('power {} and scale {}'.format(power_i,scale_i))
        plt.ylabel('cummulative distribution function')
        plt.xlabel('Storm Severity Index []')

        ks_extrapolate_i = np.max(np.absolute((cdf_gev_i-cdf_wisc_prob_i)[cdf_wisc_prob_i>(1-1/(1/(wisc_hist.frequency[0])))**block_size_years]))
        
        plt.draw()
        plt.pause(0.05)
#        D, p_value = scipy.stats.kstest(wisc_prob_CH_i.ssi_full_area,'genextreme',args= params_final_gev)
#        D, p_value = scipy.stats.kstest(wisc_prob_CH_i.ssi_full_area, lambda x: 1-(1/scipy.stats.genextreme.cdf(x, params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2]))**block_size_years)
#        D, p_value = scipy.stats.kstest(wisc_prob_CH_i.ssi_full_area, lambda x: 1.0-1.0/(1.0/(1.0-scipy.stats.genextreme.cdf(x, params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2]))*block_size_years))
        power_list.append(power_i)
        scale_list.append(scale_i)
        ks_extrapolate.append(ks_extrapolate_i)
        print('for power {} and scale {} kstest reveals an adapted ks teststatistic {} .'.format(power_i,scale_i,ks_extrapolate_i))

results_df = pd.DataFrame(list(zip(power_list, scale_list, ks_extrapolate)),columns=['power','scale','ks_ts_extr'])


plt.figure()
plt.scatter(x=results_df['power'],y=results_df['scale'],c=results_df['ks_ts_extr'],s=64,cmap='plasma_r')
plt.xlabel('power')
plt.ylabel('scale')
cbar = plt.colorbar()
cbar.set_label('ks teststat. only in extrapolation range')

## selected parameters
power_final = results_df['power'][results_df['ks_ts_extr']==results_df['ks_ts_extr'].min()].values[0]
scale_final = results_df['scale'][results_df['ks_ts_extr']==results_df['ks_ts_extr'].min()].values[0]


# how much do the maximum intensities per event change with the chosen parameters?
plt.boxplot((scale_final*np.power(wisc_hist.intensity.max(axis=1).todense()[:],power_final))[:])
plt.ylabel('max delta wind gust per storm [m/s]')

# load the probabilistic hazard
wisc_prob_CH = StormEurope()
wisc_prob_CH.read_hdf5(wisc_prob_filename)

# get the ssi values and their return periods from the probabilistic hazard
return_period_prob = 1.0/np.cumsum(wisc_prob_CH.frequency[np.argsort(-wisc_prob_CH.ssi_full_area)])
ssi_quantiles_prob = wisc_prob_CH.ssi_full_area[np.argsort(-wisc_prob_CH.ssi_full_area)]

# plot the ssi with the fitted GEV distribution
plt.figure()
plt.fill_between(return_periods,
                 np.quantile(resampling_ppf,[0.05],axis=0)[0],
                 np.quantile(resampling_ppf,[0.95],axis=0)[0],
                 facecolor = 'b',
                 alpha= 0.3)
plt.plot(return_periods,ssi_quantiles_distribution)
plt.plot(return_period_hist,ssi_quantiles_hist,'r.')
plt.plot(return_period_prob,ssi_quantiles_prob,'b.')
plt.xlim([1,3000])
plt.ylim([0,10**11])
plt.xscale('log')
#plt.title('power {} and scale {}'.format(power_final,scale_final))
plt.xlabel('Returnperiod [years]')
plt.ylabel('Storm Severity Index []')


#plot cdfs to get ks test statistic
cdf_gev_i = scipy.stats.genextreme.cdf(ssi_quantiles_prob, params_final_gev[0],loc=params_final_gev[1],scale=params_final_gev[2])
plt.figure()
plt.plot(ssi_quantiles_prob,cdf_gev_i,'.k')
#        plt.plot(ssi_quantiles_distribution,(1-1/return_periods)**5,'.b')
cdf_wisc_prob_i = (1-1/return_period_prob)**block_size_years
plt.plot(ssi_quantiles_prob,cdf_wisc_prob_i,'.b')
plt.xlim(xmin=10^10)
#plt.title('power {} and scale {}'.format(power_final,scale_final))
plt.ylabel('cummulative distribution function')
plt.xlabel('Storm Severity Index []')


# show qq plot of historic and probabilistic ssis
sm.qqplot_2samples(wisc_hist.ssi,
                   wisc_prob_CH.ssi_full_area,
                   line='45',
                   )
#                   xlabel='"WISC historic" pan-European SSI',
#                   ylabel='"WISC probabilistic extension" pan-European SSI'
plt.xlabel('"WISC probabilistic extension" pan-European SSI')
plt.ylabel('"WISC historic" pan-European SSI')

