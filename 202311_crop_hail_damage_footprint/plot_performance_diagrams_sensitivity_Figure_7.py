# -*- coding: utf-8 -*-
""" 
Plot Figure 7: Performance diagramm (Sensitivity)
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
#croptype
croptypes=['wheat','maize','barley','rapeseed']
units={'MESHS':'mm','E_kin':r'Jm$^{-2}$','POH':'%'}

#%% read data
at_centroid_data={}
# load data for MESHS  and POH
for var in ['MESHS','POH']:
    data, croptype = read_at_centroid_data(datadir,croptypes,variable=var)
    at_centroid_data[var]=data
    
#%% Plot Performance diagram and ROC diagram Paper

# figure setup
colors=['k','lightgrey','wheat']
markers=['.','d','o']
opt_resolutions=['1km','4km','8km'] #resolutions to be plotted
opt_threshs={'POH':[70,85,100], 'MESHS': [20,30,40]} #thresholds to be plotted
exp_threshs=[0,20] #exposure thresholds to be plotted
th=20 #exposure threshold to label
score='1-FAR'  #score on horizontal axis
clab='CSI' #score shown in shading
labels=['a)','b)']
titles=[r'exposure (n$_{thresh}$=20)','hazard (POH)']
font_s=16 #fontsizes
font_l=20
figsize=(14,10)

#update fontsizes
matplotlib.rcParams.update({'font.size': font_l, 'axes.labelsize':font_l})

#%% Create Figure - Performance Diagram only (2 Panels, POH, exposure)

fig,axs=plt.subplots(figsize=(14,6),nrows=1,ncols=2,sharey=True)
if not isinstance(axs,np.ndarray):
    axs=[axs]
plt.subplots_adjust(wspace=0.1)
hs=[] #handles

#loop over axis (length of axs is 1 if only performance diagram plotted)
for j,ax in enumerate(axs):

    ax.set_title(titles[j],fontweight='bold',fontsize=font_l,pad=40) 

    #setup performance diagram
    ax,h,hl=plot_performance_diagram(ax)
    hs.append(hl)

    #add colorbar
    if j==1:
        # Adding the colorbar
        cbaxes = fig.add_axes([0.25,-0.08, 0.5, 0.05]) 
        fig.colorbar(h,ax=ax,orientation='horizontal',label=clab,cax=cbaxes)
        ax.set_ylabel('')
    
    #plot points in diagram (POH/Hazard sensitivity) PANEL A
    if j==1:
        for variable in ['MESHS','POH']:
            exp_thresh=0
            if variable=='POH':
                ax=plot_diagram_points(at_centroid_data[variable],variable,opt_resolutions,exp_thresh,opt_threshs[variable],units[variable],score,colors,markers,font_s,ax,
                                       labeled=True,
                                       label_resolutions=True)
            elif variable=='MESHS':
                ax=plot_diagram_points(at_centroid_data[variable],variable,opt_resolutions,exp_thresh,opt_threshs[variable],units[variable],score,colors,markers,font_s,ax,
                                       labeled=False,
                                       label_resolutions=False)
    #plot points in diagram (exposure threshold sensitivity) PANEL B
    if j==0:
        variable='MESHS'
        for exp_thresh in exp_threshs:
            if exp_thresh==th:
                ax=plot_diagram_points(at_centroid_data[variable],variable,opt_resolutions,exp_thresh,opt_threshs[variable],units[variable],score,colors,markers,font_s,ax,
                                       labeled=True,
                                       label_resolutions=True)            
            else:
                ax=plot_diagram_points(at_centroid_data[variable],variable,opt_resolutions,exp_thresh,opt_threshs[variable],units[variable],score,colors,markers,font_s,ax,
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
#fig.savefig(f'{figdir}/Figure_Performance_diagram_field_crops_POH_exp.png',dpi=300,bbox_inches='tight')
#fig.savefig(f'{figdir}/Figure_Performance_diagram_field_crops_POH_exp.pdf',dpi=300,bbox_inches='tight')
