# -*- coding: utf-8 -*-
""" 
Plot Figure 5 : Performance diagramm
Uses the 'Tkagg' backend of matplotlib to manually place labels
of the bias plot.

Created on Tue Feb  7 13:22:01 2023

@author: Raphael Portmann
"""
import sys, os
#change path to path of current skript, where also utility.py lies
sys.path.append('C:/Users/F80840370/projects/scClim/climada/scClim/subproj_D/papers/NHESS/code_and_data/')

import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib 
from utility import data_dir, read_at_centroid_data, compute_verification_stats, plot_performance_diagram, plot_diagram_points

## uncomment to manually place labels in performance diagramm
matplotlib.use('TkAgg')

#%% Directories
datadir = f"{data_dir}/data_at_centroid/"
figdir = 'C:/Users/F80840370/projects/scClim/climada/scClim/subproj_D/papers/NHESS/fig/'
#croptypes (use ['Weizen','Mais','Gerste','Raps'] for field crops)
croptypes_lists=[['wheat','maize','barley','rapeseed'],['grapevine']]#
variable='MESHS'
units={'MESHS':'mm','E_kin':r'Jm$^{-2}$','POH':'%'}
#resolutions to be plotted
opt_resolutions=['1km','4km','8km']

#thresholds to be plotted
opt_threshs=[20,30,40]

#exposure threshold(s)
exp_threshs=[1]
score='1-FAR' #1-FAR for Performance diagram POFD for ROC curve
#colorbar label (CSI for Performance diagram and HK for ROC curve)
clab='CSI'

#%% read data
at_centroid_data_crops={}
for croptypes in croptypes_lists:
    at_centroid_data, croptype = read_at_centroid_data(datadir,croptypes,variable=variable)
    at_centroid_data_crops[croptype] = at_centroid_data

#%% Plot Performance diagram

# setup plot 
# Colors Markers
colors=['k','lightgrey','wheat']
markers=['.','d','o']
#fontsizes
font_s=16
font_l=20
matplotlib.rcParams.update({'font.size': font_l, 'axes.labelsize':font_l})
labels=['a)','b)','c)','d)','e)','f)','g)']
#figsize
figsize=(14,10)
th=0 #exposure threshold to label
titles=['field crops', 'grapevine']

fig,axs=plt.subplots(figsize=(14,6),nrows=1,ncols=2,sharey=True)
if not isinstance(axs,np.ndarray):
    axs=[axs]
plt.subplots_adjust(wspace=0.1)
hs=[] #list of handles
#loop over crops
for j,croptype in enumerate(at_centroid_data_crops.keys()):

    ax=axs[j]
    ax.set_title(titles[j],fontweight='bold',fontsize=font_l,pad=40)   

    #plot raw performance diagram
    ax,h,hl=plot_performance_diagram(ax)
    hs.append(hl)

    #add colorbar
    if j==1:
        # Adding the colorbar
        cbaxes = fig.add_axes([0.25,-0.08, 0.5, 0.05]) 
        fig.colorbar(h,ax=ax,orientation='horizontal',label=clab,cax=cbaxes)
        ax.set_ylabel('')
    
    #plot points in diagram
    for exp_thresh in exp_threshs:
        if exp_thresh==th or len(exp_threshs)==1:
            if croptype=='Reben':
                ax=plot_diagram_points(at_centroid_data_crops[croptype],variable,opt_resolutions,exp_thresh,opt_threshs,units[variable],score,colors,markers,font_s,ax,
                                       labeled=True,
                                       label_resolutions=False,
                                       loc_thresh_labels='r')
            else:
                ax=plot_diagram_points(at_centroid_data_crops[croptype],variable,opt_resolutions,exp_thresh,opt_threshs,units[variable],score,colors,markers,font_s,ax,
                                       labeled=True,
                                       label_resolutions=True) 
        else:
            ax=plot_diagram_points(at_centroid_data_crops[croptype],variable,opt_resolutions,exp_thresh,opt_threshs,units[variable],score,colors,markers,font_s,ax,
                                       labeled=False,
                                       label_resolutions=False)
        #format axis
        ax.spines[['right', 'top']].set_visible(False)
        ax.text(0,1.05,labels[j],
                    transform=ax.transAxes,fontweight='bold')
        ax.set_xlabel(score)
        ax.set_ylim([0,1])
        ax.set_xlim([0,1])

        #add legend
        if j==1:
            leg=ax.legend(ncol=1,loc='center', title='resolution',bbox_to_anchor=(1.2,0.5),edgecolor='none')
            leg._legend_box.align = "left"
            
#uncomment to manually place bias labels
plt.clabel(hs[0],manual=True,fontsize=font_s)
plt.clabel(hs[1],manual=True,fontsize=font_s)

#save figure
#fig.savefig(f'{figdir}/Figure_Performance_diagram_field_crops_grapevine_{exp_threshs[0]}.png',dpi=300,bbox_inches='tight')
#fig.savefig(f'{figdir}/Figure_Performance_diagram_field_crops_grapevine_{exp_threshs[0]}.pdf',dpi=300,bbox_inches='tight')
