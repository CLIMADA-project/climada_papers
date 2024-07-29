#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:28:17 2021

@author: sam
"""

import os
from climada.hazard.wildfire import WildFire

from climada.hazard import Centroids
from climada.entity.exposures.base import Exposures
from climada.entity.exposures.litpop import LitPop

from climada.entity.impact_funcs.wildfire import ImpfWildfire
from climada.entity.impact_funcs import ImpactFuncSet
from climada.engine import Impact

from matplotlib import colors
import numpy as np
import pandas as pd
import geopandas as gp
import seaborn as sns
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import gridspec

%matplotlib inline

''' ------------- LOAD DATA -----------------------------------------------'''
# Data folder
DP = "/path/to/data"

# Read events file
events_comp_res = pd.read_csv(os.path.join(DP, 'EMDAT_WF_MATCHED_COMPARE_RESOLUTION.csv'))
events_comp_res['MODIS_ID'] = events_comp_res['MODIS_ID'].apply(str)

# results of the calibration as well as estimated damages are made available
# within the data folder

# please load them as pandas dataframes (pd.read_csv(...))

''' ------------- Plotting specs ------------------------------------------'''

# set colour
cm = plt.get_cmap('viridis')
scheme = [cm(i / results_CV.shape[0]) for i in range(results_CV.shape[0])]
dcol = scheme[0]
ccol = scheme[3]
sns.set_style("whitegrid")
sns.set(font_scale = 1.2)

''' ------------- Scatter plots -------------------------------------------'''

# scatter plot
def plot_scatter_damage_continents(dat, events_clean, lim_low=10**6,
                                   lim_up=10**11, col='inferno', numticks=7,
                                   xlabel=True, ylabel=True, legend=True):
    
    data = pd.DataFrame(np.vstack((dat, events_clean.Continent.values)).transpose(),
                          columns= ['y_true', 'y_est', 'Continent'])
    sns.set_style("whitegrid")
    g = sns.scatterplot(data=data, x="y_true", y="y_est", hue="Continent",
                        palette=col, s=50, legend=legend)
    g.set(xscale="log", yscale="log")
    g.set_ylim(lim_low,lim_up)
    g.set_xlim(lim_low,lim_up)    
    g.yaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=numticks))
    g.xaxis.set_major_locator(plt.LogLocator(base=10.0, numticks=numticks))
    if xlabel:
        g.set_xlabel("Reported damage [USD]", fontsize=15)
    else:
        g.set_xlabel("")
    if ylabel:
        g.set_ylabel("Estimated damage [USD]", fontsize=15)
    else:
        g.set_ylabel("")
    if legend:
        g.legend(ncol=int(len(data.Continent.unique())),
             bbox_to_anchor=(0.4, -0.5), loc='lower center',
             fontsize=15, frameon=False)
    plt.plot([lim_low, lim_up], [lim_low, lim_up], 'k-',linewidth=1)
    plt.plot([lim_low*10, lim_up*10], [lim_low, lim_up], 'k--',linewidth=1)
    plt.plot([lim_low/10, lim_up/10], [lim_low, lim_up], 'k--',linewidth=1)

 
''' ------------- Line plots ----------------------------------------------'''

def calc_MDD(i_half=535.8, intensity=np.arange(295, 501, 5)):

    i_thresh = 295
    i_n = (intensity-i_thresh)/(i_half-i_thresh)
    mdd = i_n**3/(1+i_n**3)

    return(mdd)

def prep_IF_space(results_CV, intensity=np.arange(295, 501, 5)):
    
    if_space = np.zeros((len(intensity), results_CV.shape[0]))
    for i in range(results_CV.shape[0]):
        if_space[:,i] = calc_MDD(results_CV.i_half.values[i])
        
    return(if_space)

def plot_IF_space(results_CV, intensity=np.arange(295, 501, 5), col=dcol,
                  ax=False, legend=True, ylabel=True):
    
    if_space = prep_IF_space(results_CV)
    idx_max = results_CV.i_half.idxmax()
    idx_min = results_CV.i_half.idxmin()
    
    if not ax:
        fig, ax = plt.subplots()
    ax.plot(intensity, if_space[:,0], '-', color=dcol)
    ax.fill_between(intensity, if_space[:,idx_max], if_space[:,idx_min],
                    color=dcol, alpha=0.2)
    p1 = ax.plot(np.NaN, np.NaN, color=dcol)
    p2 = ax.fill(np.NaN, np.NaN, color=dcol, alpha=0.2)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    if legend:
        plt.legend([(p2[0], p1[0]), ], ['Calibrated Impact function'],
                   bbox_to_anchor=(0.4, -0.5), loc='lower center',
                   fontsize=15, frameon=False)
    plt.ylim([0,1])
    plt.xlim(275,500)    
    plt.xlabel("Intensity [K]", fontsize=15)
    if ylabel:
        plt.ylabel("Affected assets [%]", fontsize=15)
    else:
        plt.ylabel("")

def plot_IF_all(results_CV, intensity=np.arange(295, 501, 5)):
    
    n_CV = results_CV.shape[0]
    cm = plt.get_cmap('viridis')
    scheme = [cm(i / n_CV) for i in range(n_CV)]
    

    for i in range(n_CV):
        
        y =  calc_MDD(results_CV.i_half.values[i])
        
        if i == 0:
            g = sns.lineplot(intensity, y, color=scheme[i], linewidth = 4)
        else: g = sns.lineplot(intensity, y, color=scheme[i], linewidth = 2, dashes=True)
        
    g.set_ylim(0,1)
    g.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    g.set_xlim(275,500)
    g.set_xlabel("Fire intensity [K]")
    g.set_ylabel("Precentage of affected assets")
    
''' ------------- Box Plots -----------------------------------------------'''

def plot_box_RMSF(results_CV, train_col, test_col, ylabel=True):
    
    sns.set_style("whitegrid")
    g = sns.boxplot(data=[results_CV.train_RMSF.values, results_CV.test_RMSF.values])
    g = sns.swarmplot(data=[results_CV.train_RMSF.values, results_CV.test_RMSF.values], color="0.01")
    g.set_ylim(0,125)
    if ylabel:
        g.set_ylabel("RMSF of 10 fold CV")
    else:
        g.set_ylabel("")
    #g.set_xlabel("Error scores")
    box = g.artists[0]
    box.set_facecolor(train_col) 
    box = g.artists[1]
    box.set_facecolor(test_col)
    
    g.set_xticklabels(['Training data', 'Test data'])

def plot_box_i_half(results_CV, box_col):
    
    sns.set_style("whitegrid")
    g = sns.boxplot(data=results_CV.i_half.values)
    g = sns.swarmplot(data=results_CV.i_half.values, color="0.01")
    g.set_ylim(250,700)
    g.set_xlabel("Calibration results")
    g.set_ylabel("i-half of 10 fold CV")
    box = g.artists[0]
    box.set_facecolor(box_col) 
    g.set_xticklabels([' '])
    
def plot_boxplots(results_CV, box_col):

    fig = plt.figure(figsize=(16, 6))
    sns.set(font_scale = 1.5)
            
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2], wspace=0.2, hspace=0.4)
    
    plt.subplot(gs[:,0])
    g1 = plot_box_i_half(results_CV, box_col)
        
    plt.subplot(gs[0,1])
    g2 = plot_box_RMSF(results_CV, box_col)




''' ------- Figure 1 -------------------- '''

fig = plt.figure(figsize=(16, 9))
sns.set(font_scale = 1.5)
sns.set_style("whitegrid")

gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1], wspace=0.3, hspace=0.6)
plt.subplot(gs[0,0])
plot_scatter_damage_continents(dat_haz1_exp30, events_comp_haz1_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False)
plt.subplot(gs[0,1])
plot_scatter_damage_continents(dat_haz4_exp120, events_comp_haz1_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, ylabel=False)
plt.subplot(gs[0,2])
plot_scatter_damage_continents(dat_haz10_exp300, events_comp_haz1_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)

plot_IF_space(results_CV_haz1_exp30,ax=plt.subplot(gs[1,0]), legend=False)
plot_IF_space(results_CV_haz4_exp120,ax=plt.subplot(gs[1,1]), ylabel=False)
plot_IF_space(results_CV_haz10_exp300,ax=plt.subplot(gs[1,2]), legend=False, ylabel=False)

''' ------- Figure Appendix Scatter -------------------- '''

fig = plt.figure(figsize=(16, 12))
sns.set(font_scale = 1.5)
sns.set_style("whitegrid")

gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], wspace=0.3, hspace=0.5)

plt.subplot(gs[0,0])
plot_scatter_damage_continents(dat_haz1_exp30, events_comp_haz1_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False)
plt.subplot(gs[0,1])
plot_scatter_damage_continents(dat_haz1_exp120, events_comp_haz1_exp120, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)
plt.subplot(gs[0,2])
plot_scatter_damage_continents(dat_haz1_exp300, events_comp_haz1_exp300, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)
plt.subplot(gs[1,0])
plot_scatter_damage_continents(dat_haz4_exp30, events_comp_haz4_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False)
plt.subplot(gs[1,1])
plot_scatter_damage_continents(dat_haz4_exp120, events_comp_haz4_exp120, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)
plt.subplot(gs[1,2])
plot_scatter_damage_continents(dat_haz4_exp300, events_comp_haz4_exp300, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)
plt.subplot(gs[2,0])
plot_scatter_damage_continents(dat_haz10_exp30, events_comp_haz10_exp30, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False)
plt.subplot(gs[2,1])
plot_scatter_damage_continents(dat_haz10_exp120, events_comp_haz10_exp120, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, ylabel=False)
plt.subplot(gs[2,2])
plot_scatter_damage_continents(dat_haz10_exp300, events_comp_haz10_exp300, lim_low=10**4, lim_up=10**11,
                               col='inferno', numticks=11, legend=False, ylabel=False)

''' ------- Figure Appendix IF curves -------------------- '''

fig = plt.figure(figsize=(16, 12))
sns.set(font_scale = 1.5)
sns.set_style("whitegrid")

gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], wspace=0.3, hspace=0.4)

plot_IF_space(results_CV_haz1_exp30,ax=plt.subplot(gs[0,0]), legend=False)
plot_IF_space(results_CV_haz1_exp120,ax=plt.subplot(gs[0,1]), legend=False, ylabel=False)
plot_IF_space(results_CV_haz1_exp300,ax=plt.subplot(gs[0,2]), legend=False, ylabel=False)

plot_IF_space(results_CV_haz4_exp30,ax=plt.subplot(gs[1,0]), legend=False)
plot_IF_space(results_CV_haz4_exp120,ax=plt.subplot(gs[1,1]), legend = False, ylabel=False)
plot_IF_space(results_CV_haz4_exp300,ax=plt.subplot(gs[1,2]), legend=False, ylabel=False)

plot_IF_space(results_CV_haz10_exp30,ax=plt.subplot(gs[2,0]), legend=False)
plot_IF_space(results_CV_haz10_exp120,ax=plt.subplot(gs[2,1]), ylabel=False)
plot_IF_space(results_CV_haz10_exp300,ax=plt.subplot(gs[2,2]), legend=False, ylabel=False)

''' ------- Figure Appendix Boxplots curves -------------------- '''

fig = plt.figure(figsize=(16, 12))
sns.set(font_scale = 1.5)
sns.set_style("whitegrid")

gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], wspace=0.3, hspace=0.25)
plt.subplot(gs[0,0])
plot_box_RMSF(results_CV_haz1_exp30, dcol, ccol)
plt.subplot(gs[0,1])
plot_box_RMSF(results_CV_haz1_exp120, dcol, ccol, ylabel=False)
plt.subplot(gs[0,2])
plot_box_RMSF(results_CV_haz1_exp300, dcol, ccol, ylabel=False)

plt.subplot(gs[1,0])
plot_box_RMSF(results_CV_haz4_exp30, dcol, ccol)
plt.subplot(gs[1,1])
plot_box_RMSF(results_CV_haz4_exp120, dcol, ccol, ylabel=False)
plt.subplot(gs[1,2])
plot_box_RMSF(results_CV_haz4_exp300, dcol, ccol, ylabel=False)

plt.subplot(gs[2,0])
plot_box_RMSF(results_CV_haz10_exp30, dcol, ccol)
plt.subplot(gs[2,1])
plot_box_RMSF(results_CV_haz10_exp120, dcol, ccol, ylabel=False)
plt.subplot(gs[2,2])
plot_box_RMSF(results_CV_haz10_exp300, dcol, ccol, ylabel=False)

