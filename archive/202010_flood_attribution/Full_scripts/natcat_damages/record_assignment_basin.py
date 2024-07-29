#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 14:34:34 2020

@author: insauer
"""

import numpy as np
import pandas as pd
import sys
import os
import copy
import matplotlib.pyplot as plt

#from climada.entity.exposures.exp_people import ExpPop

from climada.hazard.river_flood import RiverFlood
from climada.hazard.centroids import Centroids
from climada.entity import ImpactFuncSet
from climada.util.constants import RIVER_FLOOD_REGIONS_CSV
from climada.entity.exposures.litpop import LitPop
from sklearn.neighbors import NearestNeighbors

def find_dis_group(lon, lat, dis_pos, dis_neg):
    
    centroids = np.zeros((dis_pos.centroids.lat.shape[0], 2))
    centroids[:, 0] = dis_pos.centroids.lat
    centroids[:, 1] = dis_pos.centroids.lon
    if dis_pos.centroids.lat.shape[0] < 100:
        n_b = dis_pos.centroids.lat.shape[0]
    else:
        n_b = 100
    nbrs = NearestNeighbors(n_neighbors=n_b, algorithm='ball_tree').fit(centroids)
    distances, indices = nbrs.kneighbors(np.array([[lat, lon]]))
    
    pos_group = dis_pos.intensity[0, indices[0][0]]
    neg_group = dis_neg.intensity[0, indices[0][0]]
    if pos_group == 1:
        return 'pos_risk'
    elif neg_group ==1:
        return 'neg_risk'
    else:
        pos_size = np.where(dis_pos.intensity[0, indices[0, :]].todense()==1)[0].shape[0]
        if pos_size == 0:
            pos_dist = 100000
        else:
            pos_dist_ind = np.where(dis_pos.intensity[0, indices[0, :]].todense()==1)[1][0]
            pos_dist = distances[0, pos_dist_ind]

        neg_size = np.where(dis_neg.intensity[0, indices[0, :]].todense()==1)[0].shape[0]
        if neg_size == 0:
            neg_dist = 100000
        else:
            neg_dist_ind = np.where(dis_neg.intensity[0, indices[0, :]].todense()==1)[1][0]
            neg_dist = distances[0, neg_dist_ind]
            
        if pos_dist < neg_dist:
            
            return 'pos_risk'
        
        elif neg_dist < pos_dist:
            
            return 'neg_risk'
        else:
            
            if dis_pos.centroids.lat.shape[0] < 10000:
                n_b = dis_pos.centroids.lat.shape[0]
            else:
                n_b = 10000
            nbrs = NearestNeighbors(n_neighbors=n_b, algorithm='ball_tree').fit(centroids)
            distances, indices = nbrs.kneighbors(np.array([[lat, lon]]))
            
            pos_size = np.where(dis_pos.intensity[0, indices[0, :]].todense()==1)[0].shape[0]
            if pos_size == 0:
                pos_dist = 100000
            else:
                pos_dist_ind = np.where(dis_pos.intensity[0, indices[0, :]].todense()==1)[1][0]
                pos_dist = distances[0, pos_dist_ind]

            neg_size = np.where(dis_neg.intensity[0, indices[0, :]].todense()==1)[0].shape[0]
            if neg_size == 0:
                neg_dist = 100000
            else:
                neg_dist_ind = np.where(dis_neg.intensity[0, indices[0, :]].todense()==1)[1][0]
                neg_dist = distances[0, neg_dist_ind]
                
            if pos_dist < neg_dist:
            
                return 'pos_risk'
        
            elif neg_dist < pos_dist:
            
                return 'neg_risk'
            else:
                
                return 'unclear'


country_info = pd.read_csv(RIVER_FLOOD_REGIONS_CSV)

dph_path = '/home/insauer/mnt/ebm/data/hazard/floods/isimip2a/gswp3/clm40/depth-150arcsec/flddph_annual_max_gev_0.1mmpd_protection-flopros.nc'
frc_path = '/home/insauer/mnt/ebm/data/hazard/floods/isimip2a/gswp3/clm40/area-150arcsec/fldfrc_annual_max_gev_0.1mmpd_protection-flopros.nc'


isos = country_info['ISO'].tolist()

#isos = ['BRB']

natcat=pd.read_excel('/home/insauer/projects/Attribution/Floods/Paper_NC_Review_Data/Input_PPP_conversion/1980-2016_Masterdatei_NatCat_worldwide_no-pw_2005conversions_PPP.xlsx', index_col=0)

trend_path = '/home/insauer/projects/RiverDischarge/Data/basin_trends_geo.nc'

#trend_path = '/home/insauer/projects/Attribution/Floods/Paper_NC_20_06_Data/Input_Damage_Assessment/TrendsMedianDischarge_MK.nc'

natcat_assigned = pd.DataFrame(columns = ['Country', 'Year', 'pos_risk', 'neg_risk', 'unclear'])

for iso in isos:
    print(iso)
    if iso =='GIB' or iso == 'MCO':
        continue
    
    dis_pos = RiverFlood()
    dis_pos.set_from_nc(dph_path=dph_path, frc_path=frc_path, countries=[iso], ISINatIDGrid=True)
    dis_neg = copy.copy(dis_pos)
    dis_pos.get_dis_mask(trend_path, 'pos')
    dis_neg.get_dis_mask(trend_path, 'neg')

    # dis_pos = FloodTrend()
    # dis_pos.set_from_nc(dph_path=trend_path, countries=[iso])
    # dis_neg = copy.copy(dis_pos)
    # dis_pos.get_dis_mask(dis = 'pos')
    # dis_neg.get_dis_mask(dis = 'neg')
    
    natcat_cnt_all = natcat[(natcat['Country ISO3']==iso) &
                        (natcat['Subevent']=='gf:General flood')]
    
    years = list(set(natcat_cnt_all['Year']))
    
    for year in years:
        print(year)
        natcat_row = pd.DataFrame(columns = ['Country', 'Year', 'pos_risk', 'neg_risk', 'unclear'])
        natcat_row.loc[0,'Country']= iso
        natcat_row.loc[0,'Year'] = year
        natcat_cnt = natcat_cnt_all[natcat_cnt_all['Year']==year]
        n_events = natcat_cnt.shape[0]
        natcat_row.loc[0,'pos_risk'] = 0
        natcat_row.loc[0,'neg_risk'] = 0
        natcat_row.loc[0,'pos_events'] = 0
        natcat_row.loc[0,'neg_events'] = 0
        natcat_row.loc[0,'unclear_events'] = 0
        natcat_row.loc[0,'unclear'] = 0
        
        for ev in range(n_events):
            yearDF = natcat_cnt.loc[:,['Longitude', 'Latitude', 'CPI_con']]
            
            event = yearDF.iloc[ev,:]
            lon = event['Longitude']
            lat = event['Latitude']
            cpi_con = event['CPI_con']
            risk_group = find_dis_group(lon, lat, dis_pos, dis_neg)
            natcat_row.loc[0, risk_group] += cpi_con
            if risk_group == 'pos_risk':
                natcat_row.loc[0,'pos_events'] +=1
            elif risk_group == 'neg_risk':
                natcat_row.loc[0,'neg_events'] +=1
            else:
                natcat_row.loc[0,'unclear_events'] += 1
        natcat_assigned= natcat_assigned.append(natcat_row,ignore_index=True)
    
    natcat_assigned.to_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Review_Data/NatCat_PPP_conv/natcat_subregions.csv', index=False)
        
natcat_assigned.to_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Review_Data/NatCat_PPP_conv/natcat_subregions.csv', index=False)
