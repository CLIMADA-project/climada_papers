#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 22:09:13 2020

@author: insauer
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import matplotlib.patches as mpatches

fig3 = plt.figure(constrained_layout=True, figsize=(8.3, 11.7))
gs = fig3.add_gridspec(19, 14)
plt.subplots_adjust(wspace=0., hspace=0)


DATA_ATTR_Full = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Teleconnections/Lag_ENSO_GMT_PDO_NAO_Loo.csv')

DATA_ATTR = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Teleconnections/Lag_PosNeg_ENSO_PDO_NAO_GMT_Loo.csv')

teleDataNorm = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTest/tele_norm.csv')

AMO_DATA_ATTR_Full = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Teleconnections/Lag_ENSO_AMO_PDO_NAO_Loo.csv')

AMO_DATA_ATTR = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Teleconnections/Lag_PosNeg_ENSO_AMO_PDO_NAO_Loo.csv')

region_names={ 'GLB': 'Global (GLB)',
              'NAM':'North America (NAM)',
              'CHN':'Eastern Asia (EAS)',
              'AUS':'Oceania (OCE)',
              'LAM':'Latin America (LAM)',
              'EUR':'Europe (EUR)',
              'SWEA':'South & South-East Asia (SEA)',
              'CAS':'Central Asia & Russia (CAS)',
              }

region_abs={'NAM':'NAM', 
          'LAM':'LAM', 
          'EUR':'EUR',
          'NAFARA':'NAR',
          'SSAF':'SSA',
          'CAS':'CAS',
          'SWEA':'SEA', 
          'CHN':'EAS', 
          'AUS':'OCE',
          'GLB': 'GLB'}

regions = list(region_names)

three = [2]
two = [0,1,3,4]
one = [5,6,7]

r =0

for i in range(1):
    #r =0
    
    z =0
    for j in range(8):
        
        if r > 8:
            continue
        if r in one:
            f3_ax1 = fig3.add_subplot(gs[5*(i+1):5*(i+1)+4,z:z+1])
            z = z+1
        elif r in two:
            f3_ax1 = fig3.add_subplot(gs[5*(i+1):5*(i+1)+4,z:z+2])
            z = z+2
        else:
            f3_ax1 = fig3.add_subplot(gs[5*(i+1):5*(i+1)+4,z:z+3])
            z = z+3
        #f3_ax2 = fig3.add_subplot(gs[5*i+2:5*i+4,j*5:(j*5)+4])



        data_attr_reg = AMO_DATA_ATTR_Full[AMO_DATA_ATTR_Full['Region'] == regions[r]]
        data_attr_reg_pos = AMO_DATA_ATTR[AMO_DATA_ATTR['Region'] == regions[r]+'_Pos']
        data_attr_reg_neg = AMO_DATA_ATTR[AMO_DATA_ATTR['Region'] == regions[r]+'_Neg']
        
        unexp = data_attr_reg[['Unexplained Haz']]
        
        unexp_pos = data_attr_reg_pos[['Unexplained Haz']]
        
        unexp_neg = data_attr_reg_neg[['Unexplained Haz']]


        tele = data_attr_reg[['ENSOdv_','PDOdv_','NAOdv_', 'AMOdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        tele_pos =  data_attr_reg_pos[['ENSOdv_','PDOdv_','NAOdv_', 'AMOdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        tele_neg =  data_attr_reg_neg[['ENSOdv_','PDOdv_','NAOdv_', 'AMOdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        
        
        tele_pv = data_attr_reg[['ENSOpval_','PDOpval_','NAOpval_', 'AMOpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        tele_pos_pv =  data_attr_reg_pos[['ENSOpval_','PDOpval_','NAOpval_', 'AMOpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        tele_neg_pv =  data_attr_reg_neg[['ENSOpval_','PDOpval_','NAOpval_', 'AMOpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        
        
        
        
        tele_pv = [tele_pos_pv, tele_pv, tele_neg_pv]
        
        tele_cons = [tele_pos, tele, tele_neg]
        
        df_teles=pd.DataFrame(columns= ['limit','tele_pos', 'tele', 'tele_neg'])
                              
        tele_names = ['tele_pos', 'tele', 'tele_neg']
        
        df_teles['limit'] = ['ENSO_up', 'ENSO_bot','ENSO_pv',\
                 'PDO_up', 'PDO_bot','PDO_pv','NAO_up', 'NAO_bot', 'NAO_pv',\
                 'AMO_up', 'AMO_bot','AMO_pv',\
                 'ENSO_lag_up', 'ENSO_lag_bot','ENSO_lag_pv',\
                 'PDO_lag_up', 'PDO_lag_bot','PDO_lag_pv','NAO_lag_up', 'NAO_lag_bot', 'NAO_lag_pv',\
                 'Unex_up', 'Unex_bot']
        
        preds = ['ENSO','PDO','NAO', 'AMO', 'ENSO_lag','PDO_lag','NAO_lag']
        
        for t,tel in enumerate(tele_cons):
            
            abs_tel = np.nan_to_num(np.abs(tel))
            sum_tel = np.sum(np.array(abs_tel))

            perc = np.nan_to_num(tel/sum_tel)
            
            bottom = np.array(perc)[np.less(perc,0)].sum(axis =0)
            
            
            for m in np.argsort(perc)[0]:
                
                df_teles.loc[df_teles['limit']==preds[m]+'_bot', tele_names[t]] = bottom
                
                df_teles.loc[df_teles['limit']==preds[m]+'_up', tele_names[t]] = np.abs(perc[0,m])
            
            
                bottom += np.abs(perc[0,m])
                
        unex_order =  ['Unexplained Haz']
        t_iter = 0
        
        shortage=['']
        
        shor = 0
        
        for t,tel_pv in enumerate(tele_pv):
            
            if t in [3,6,9]:
                shor+=1
            
            
            for p in preds:
                
                if tel_pv[p+'pval_'+shortage[shor]].values[0]<0.1:
                    
                    alpha = 'k'
                else:
                    alpha = 'w'
                    
                df_teles.loc[df_teles['limit']==p+'_pv', tele_names[t]] = alpha
                
        
        for unex in unex_order:
            
            
            if unexp_pos[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            t_iter+=1
            
            if unexp[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            
            t_iter+=1
            
            if unexp_neg[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            
            t_iter+=1
         

        inds= [[0,1],[1,2],[0,1,2],[1,2], [0,1], [1], [0],[0]]
        
        if r in three:
        
            x = [0,1,2]
        
        elif r in two:
        
            x = [0,1]
        
        else:
        
            x = [0]
        
        
        colour_code = ['#5ab4ac', '#4575b4', '#d8b365', '#5ab4ac', '#4575b4', '#d8b365']
        
        df_teles_a1 = df_teles.iloc[:,0:4]
        
        cmap = plt.cm.get_cmap('tab20b')
        
        rgba_oce_d = cmap(0.02)
        rgba_oce_l = cmap(0.17)
        
        col_enso = rgba_oce_d
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_bot'])[0,1:][inds[r]], color = 'cornflowerblue',
        edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='PDO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_bot'])[0,1:][inds[r]], color = 'gold', label = 'PDO',
        edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='PDO_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='NAO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_bot'])[0,1:][inds[r]], color = 'gray',
         edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='NAO_pv'])[0,1:][inds[r]], label = 'NAO', linewidth = 1.5)
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='AMO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_bot'])[0,1:][inds[r]], color = 'seagreen',
        label ='AMO', edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='AMO_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_bot'])[0,1:][inds[r]], color = 'steelblue',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_pv'])[0,1:][inds[r]], linewidth = 1.5 )
        
        # f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='AMO_up'])[0,1:][inds[r]], 
        # bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_bot'])[0,1:][inds[r]], color = 'mediumseagreen',
        # width = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_pv'])[0,1:][inds[r]] )
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_bot'])[0,1:][inds[r]], color = 'gold',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_pv'])[0,1:][inds[r]] )
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_bot'])[0,1:][inds[r]], color = 'gray',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_pv'])[0,1:][inds[r]] )
        
        f3_ax1.bar(np.array(x), np.nan_to_num(np.array(df_teles_a1[df_teles_a1['limit']=='Unex_up'])[0,1:].astype(float))[inds[r]], 
        bottom = np.nan_to_num(np.array(df_teles_a1[df_teles_a1['limit']=='Unex_bot'])[0,1:].astype(float))[inds[r]], color = 'white', hatch='///', alpha =0.2, label = 'unexplained Trend')
        
        handles, labels = f3_ax1.get_legend_handles_labels()
        
        ax1_lims_up = [1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4]
        
        ax1_lims_low = [-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3]
        
        if r in [0,4]:

            f3_ax1.set_xlim(-1.,2)
            
            f3_ax1.set_xticks([0,1])
            
            f3_ax1.set_xticklabels([ '$R_{+}$', '$R$'], fontsize = 8)
            
        elif r == 2:
            f3_ax1.set_xlim(-1.,3)
            
            f3_ax1.set_xticks([0,1,2])
            
            f3_ax1.set_xticklabels(['$R_{+}$', '$R$', '$R_{-}$'], fontsize = 8)
            
        elif r in [1,3]:
            f3_ax1.set_xlim(-1.,2)
            
            f3_ax1.set_xticks([0, 1])
            
            f3_ax1.set_xticklabels([ '$R$', '$R_{-}$'], fontsize = 8)
        
        elif r ==5:
            f3_ax1.set_xlim(-1,1)
            
            f3_ax1.set_xticks([0])
            
            f3_ax1.set_xticklabels([ '$R$'], fontsize = 8)
        else:
            f3_ax1.set_xlim(-1,1)
            
            f3_ax1.set_xticks([0])
            
            f3_ax1.set_xticklabels(['$R_{+}$'], fontsize = 8)
        
        ax1_ticks_up = [1,1,1,1,1,1,1,1]
        
        ax1_ticks_low = [-1,-1,-1,-1,-1,-1,-1,-1]
        
        ax1_labels_up = ['1','','','','','', '','']
        
        ax1_labels_low = ['-1','','','','','','','']
        
        # ax2_labels_up = ['1','1','0.5','1','1', '1','0.5']
        
        # ax2_labels_low = ['-0.5','-1','-1','-1','-0.5','-0.5','-1']
       
        # f3_ax2.set_xlim(-1.,7)
        
        f3_ax1.set_ylim(ax1_lims_low[r],ax1_lims_up[r])
        
        f3_ax1.set_yticks([ax1_ticks_low[r] , 0, ax1_ticks_up[r]])
        
        f3_ax1.set_yticklabels([ax1_labels_low[r], ax1_labels_up[r] ],fontsize =6)
        
        if r==0 :
            f3_ax1.set_yticklabels([ax1_labels_low[r],'0', ax1_labels_up[r] ],fontsize =6)
        
        # f3_ax2.set_ylim(ax2_lims_low[r],ax2_lims_up[r])
        
        # f3_ax2.set_yticks([ax2_ticks_low[r] , 0, ax2_ticks_up[r]])
        
        # f3_ax2.set_yticklabels([ax2_labels_low[r],'0', ax2_labels_up[r] ],fontsize =6)
        
        

        f3_ax1.axhline(y=0,linewidth=0.3, color='k', linestyle = '-', alpha = 0.5)
        # f3_ax2.axvline(x=3,linewidth=0.3, color='k', linestyle = '-', alpha = 0.5)
        
        
        
        
        
        if r ==0 :
            f3_ax1.set_ylabel('Coefs in $D_{Haz-1980}$',  fontsize = 10, labelpad=+1)
            

        
        
        f3_ax1.set_title(' '+ region_abs[regions[r]], position = (0.5,0.83), fontsize = 9)
            
        if r ==0:
        
            f3_ax1.text(-0.55, 0.98, 'b', transform=f3_ax1.transAxes, 
                size=15, weight='bold')
            
            f3_ax1.text(-0.55, 2.2, 'a', transform=f3_ax1.transAxes, 
                size=15, weight='bold')
       
        r+=1


r =0

for i in range(1):
    #r =0
    
    z =0
    for j in range(8):
        
        if r in one:
            f3_ax1 = fig3.add_subplot(gs[5*(i+2):5*(i+2)+4,z:z+1])
            z = z+1
        elif r in two:
            f3_ax1 = fig3.add_subplot(gs[5*(i+2):5*(i+2)+4,z:z+2])
            z = z+2
        else:
            f3_ax1 = fig3.add_subplot(gs[5*(i+2):5*(i+2)+4,z:z+3])
            z = z+3
        #f3_ax2 = fig3.add_subplot(gs[5*i+2:5*i+4,j*5:(j*5)+4])

        
        
        
            
            

        data_attr_reg = DATA_ATTR_Full[DATA_ATTR_Full['Region'] == regions[r]]
        data_attr_reg_pos = DATA_ATTR[DATA_ATTR['Region'] == regions[r]+'_Pos']
        data_attr_reg_neg = DATA_ATTR[DATA_ATTR['Region'] == regions[r]+'_Neg']
        
        unexp = data_attr_reg[['Unexplained Haz']]
        
        unexp_pos = data_attr_reg_pos[['Unexplained Haz']]
        
        unexp_neg = data_attr_reg_neg[['Unexplained Haz']]


        tele = data_attr_reg[['ENSOdv_','PDOdv_','NAOdv_', 'GMTdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        tele_pos =  data_attr_reg_pos[['ENSOdv_','PDOdv_','NAOdv_', 'GMTdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        tele_neg =  data_attr_reg_neg[['ENSOdv_','PDOdv_','NAOdv_', 'GMTdv_','ENSO_lagdv_','PDO_lagdv_','NAO_lagdv_']]
        
        
        tele_pv = data_attr_reg[['ENSOpval_','PDOpval_','NAOpval_', 'GMTpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        tele_pos_pv =  data_attr_reg_pos[['ENSOpval_','PDOpval_','NAOpval_', 'GMTpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        tele_neg_pv =  data_attr_reg_neg[['ENSOpval_','PDOpval_','NAOpval_', 'GMTpval_','ENSO_lagpval_','PDO_lagpval_','NAO_lagpval_']]
        

        tele_pv = [tele_pos_pv, tele_pv, tele_neg_pv]
        
        tele_cons = [tele_pos, tele, tele_neg]
        
        df_teles=pd.DataFrame(columns= ['limit','tele_pos', 'tele', 'tele_neg'])
                              
        tele_names = ['tele_pos', 'tele', 'tele_neg']
        
        df_teles['limit'] = ['ENSO_up', 'ENSO_bot','ENSO_pv',\
                 'PDO_up', 'PDO_bot','PDO_pv','NAO_up', 'NAO_bot', 'NAO_pv',\
                 'GMT_up', 'GMT_bot','GMT_pv',\
                 'ENSO_lag_up', 'ENSO_lag_bot','ENSO_lag_pv',\
                 'PDO_lag_up', 'PDO_lag_bot','PDO_lag_pv','NAO_lag_up', 'NAO_lag_bot', 'NAO_lag_pv',\
                 'Unex_up', 'Unex_bot']
        
        preds = ['ENSO','PDO','NAO', 'GMT', 'ENSO_lag','PDO_lag','NAO_lag']
        
        for t,tel in enumerate(tele_cons):
            
            abs_tel = np.nan_to_num(np.abs(tel))
            sum_tel = np.sum(np.array(abs_tel))

            perc = np.nan_to_num(tel/sum_tel)
            
            bottom = np.array(perc)[np.less(perc,0)].sum(axis =0)
            
            
            for m in np.argsort(perc)[0]:
                
                df_teles.loc[df_teles['limit']==preds[m]+'_bot', tele_names[t]] = bottom
                
                df_teles.loc[df_teles['limit']==preds[m]+'_up', tele_names[t]] = np.abs(perc[0,m])
            
            
                bottom += np.abs(perc[0,m])
                
        unex_order =  ['Unexplained Haz']
        t_iter = 0
        
        shortage=['']
        
        shor = 0
        
        for t,tel_pv in enumerate(tele_pv):
            
            if t in [3,6,9]:
                shor+=1
            
            
            for p in preds:
                
                if tel_pv[p+'pval_'+shortage[shor]].values[0]<0.1:
                    
                    alpha = 'k'
                else:
                    alpha = 'w'
                    
                df_teles.loc[df_teles['limit']==p+'_pv', tele_names[t]] = alpha
                
        
        for unex in unex_order:
            
            
            if unexp_pos[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            t_iter+=1
            
            if unexp[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            
            t_iter+=1
            
            if unexp_neg[unex].sum() >0:
                df_teles.loc[df_teles['limit']=='Unex_up', tele_names[t_iter]] = 1
                df_teles.loc[df_teles['limit']=='Unex_bot', tele_names[t_iter]] = 0
            
            t_iter+=1
         

        inds= [[0,1],[1,2],[0,1,2],[1,2], [0,1], [1], [0],[0]]
        
        if r in three:
        
            x = [0,1,2]
        
        elif r in two:
        
            x = [0,1]
        
        else:
        
            x = [0]
        
        colour_code = ['#5ab4ac', '#4575b4', '#d8b365', '#5ab4ac', '#4575b4', '#d8b365']
        
        df_teles_a1 = df_teles.iloc[:,0:4]
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_bot'])[0,1:][inds[r]], color = 'cornflowerblue',
        edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='PDO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_bot'])[0,1:][inds[r]], color = 'gold', label = 'PDO',
        edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='PDO_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='NAO_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_bot'])[0,1:][inds[r]], color = 'gray',
        edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='NAO_pv'])[0,1:][inds[r]], label = 'NAO', linewidth = 1.5)
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='GMT_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='GMT_bot'])[0,1:][inds[r]], color = 'sandybrown',
        label ='GMT', edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='GMT_pv'])[0,1:][inds[r]], linewidth = 1.5)
        
        # f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='AMO_up'])[0,1:][inds[r]], 
        # bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_bot'])[0,1:][inds[r]], color = 'mediumseagreen',
        # label ='AMO', edgecolor=np.array(df_teles_a1[df_teles_a1['limit']=='AMO_pv'])[0,1:][inds[r]])
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_bot'])[0,1:][inds[r]], color = 'steelblue',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_lag_pv'])[0,1:][inds[r]] )
        
        # f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='AMO_up'])[0,1:][inds[r]], 
        # bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_bot'])[0,1:][inds[r]], color = 'mediumseagreen',
        # width = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_pv'])[0,1:][inds[r]] )
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_bot'])[0,1:][inds[r]], color = 'gold',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_lag_pv'])[0,1:][inds[r]] )
        
        f3_ax1.bar(np.array(x), np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_up'])[0,1:][inds[r]], 
        bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_bot'])[0,1:][inds[r]], color = 'gray',
        edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_lag_pv'])[0,1:][inds[r]] )
        
        
        f3_ax1.bar(np.array(x), np.nan_to_num(np.array(df_teles_a1[df_teles_a1['limit']=='Unex_up'])[0,1:].astype(float))[inds[r]], 
        bottom = np.nan_to_num(np.array(df_teles_a1[df_teles_a1['limit']=='Unex_bot'])[0,1:].astype(float))[inds[r]], color = 'white', hatch='///', alpha =0.2, label = 'unexplained Trend')
        
        handles, labels = f3_ax1.get_legend_handles_labels()
        
        ax1_lims_up = [1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4]
        
        ax1_lims_low = [-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3,-1.3]
        
        if r in [0,4]:

            f3_ax1.set_xlim(-1.,2)
            
            f3_ax1.set_xticks([0,1])
            
            f3_ax1.set_xticklabels([ '$R_{+}$', '$R$'], fontsize = 8)
            
        elif r == 2:
            f3_ax1.set_xlim(-1.,3)
            
            f3_ax1.set_xticks([0,1,2])
            
            f3_ax1.set_xticklabels(['$R_{+}$', '$R$', '$R_{-}$'], fontsize = 8)
            
        elif r in [1,3]:
            f3_ax1.set_xlim(-1.,2)
            
            f3_ax1.set_xticks([0, 1])
            
            f3_ax1.set_xticklabels([ '$R$', '$R_{-}$'], fontsize = 8)
        
        elif r ==5:
            f3_ax1.set_xlim(-1,1)
            
            f3_ax1.set_xticks([0])
            
            f3_ax1.set_xticklabels([ '$R$'], fontsize = 8)
        else:
            f3_ax1.set_xlim(-1,1)
            
            f3_ax1.set_xticks([0])
            
            f3_ax1.set_xticklabels(['$R_{+}$'], fontsize = 8)
            
        
        ax1_ticks_up = [1,1,1,1,1,1,1,1]
        
        ax1_ticks_low = [-1,-1,-1,-1,-1,-1,-1,-1]
        
        ax1_labels_up = ['1','','','','','', '','']
        
        ax1_labels_low = ['-1','','','','','','','']
        
        # ax2_labels_up = ['1','1','0.5','1','1', '1','0.5']
        
        # ax2_labels_low = ['-0.5','-1','-1','-1','-0.5','-0.5','-1']
       
        # f3_ax2.set_xlim(-1.,7)
        
        f3_ax1.set_ylim(ax1_lims_low[r],ax1_lims_up[r])
        
        f3_ax1.set_yticks([ax1_ticks_low[r] , 0, ax1_ticks_up[r]])
        
        f3_ax1.set_yticklabels([ax1_labels_low[r], ax1_labels_up[r] ],fontsize =6)
        
        if r==0 :
            f3_ax1.set_yticklabels([ax1_labels_low[r],'0', ax1_labels_up[r] ],fontsize =6)
        
        # f3_ax2.set_ylim(ax2_lims_low[r],ax2_lims_up[r])
        
        # f3_ax2.set_yticks([ax2_ticks_low[r] , 0, ax2_ticks_up[r]])
        
        # f3_ax2.set_yticklabels([ax2_labels_low[r],'0', ax2_labels_up[r] ],fontsize =6)
        
        

        f3_ax1.axhline(y=0,linewidth=0.3, color='k', linestyle = '-', alpha = 0.5)
        # f3_ax2.axvline(x=3,linewidth=0.3, color='k', linestyle = '-', alpha = 0.5)
        
        
        
        
        
        if r ==0 :
            f3_ax1.set_ylabel('Coefs in $D_{Haz-1980}$',  fontsize = 10, labelpad=+1)
            

        
        
        f3_ax1.set_title(' '+ region_abs[regions[r]], position = (0.5,0.83), fontsize = 9)
            
        if r ==0:
        
            f3_ax1.text(-0.55, 0.98, 'c', transform=f3_ax1.transAxes, 
                size=15, weight='bold')
       
        r+=1
        
        
            


# enso =pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTotalRegion/22_05_ENSO_only_Loo.csv')
# ensoPosNeg= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionsPosNeg/22_05_PosNeg_ENSO_only_Loo.csv')
# amo=pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTotalRegion/22_05_AMO_only_Loo.csv')
# amoPosNeg= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionsPosNeg/22_05_PosNeg_AMO_only_Loo.csv')
# pdo = pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTotalRegion/22_05_PDO_only_Loo.csv')
# pdoPosNeg= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionsPosNeg/22_05_PosNeg_PDO_only_Loo.csv')
# nao= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTotalRegion/22_05_NAO_only_Loo.csv')
# naoPosNeg= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionsPosNeg/22_05_PosNeg_NAO_only_Loo.csv')
# gmt= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionTotalRegion/22_05_GMT_only_Loo.csv')
# gmtPosNeg= pd.read_csv('/home/insauer/projects/Attribution/Floods/Data/FinalPaper/TeleconnectionsPosNeg/22_05_PosNeg_GMT_only_Loo.csv')       

# h = 12
# r = 0
# z = 0
# for i in range(1):
#     #r =0
    

#     for j in range(7):
        
#         if r > 7:
#             continue
#         if r == 1 or r==2:
#             f3_ax1 = fig3.add_subplot(gs[h:h+3,z:z+6])
#             z = z+6
#         else:
#             f3_ax1 = fig3.add_subplot(gs[h:h+3,z:z+3])
#             z = z+3
#         #f3_ax2 = fig3.add_subplot(gs[5*i+2:5*i+4,j*5:(j*5)+4])

#         # if r== 7:
            
#         #     f3_ax1.legend(handles,labels, frameon=True, fontsize = 8, loc = 'upper left')
#         #     f3_ax1.axis('off')
#         #     r=r+1
#         #     continue
        
#         predictors = ['AMO','GMT','NAO','PDO','ENSO']
        
#         colors_preds = ['mediumseagreen','sandybrown','silver','khaki','cornflowerblue']
        
#         regs= ['R+', 'R', 'R-']
        
#         stack_data = pd.DataFrame()
#         stack_data['info'] = ['AMO','AMOpv','AMOunex','GMT','GMTpv','GMTunex',
#                               'PDO','PDOpv','PDOunex','NAO','NAOpv','NAOunex',
#                               'ENSO','ENSOpv','ENSOunex']
        
#         data_sets = [amo,gmt,nao,pdo,enso]
        
#         for da,pr in enumerate(predictors):
            
#             if np.isnan(data_sets[da].loc[data_sets[da]['Region'] == regions[r], pr+'dv_'].values[0]):
#                 stack_data.loc[stack_data['info']==pr,'R'] = 0
#             else:
#                 stack_data.loc[stack_data['info']==pr,'R'] = 1
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r], pr+'pval_'].values[0] < 0.1:
#                 stack_data.loc[stack_data['info']==pr+'pv','R'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'pv','R'] = 0
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r], 'Unexplained Haz'].sum()==1:
#                 stack_data.loc[stack_data['info']==pr+'unex','R'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'unex','R'] = 0
            
#         data_sets = [amoPosNeg,gmtPosNeg,naoPosNeg,pdoPosNeg,ensoPosNeg]

#         for da,pr in enumerate(predictors):
            
#             if np.isnan(data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Pos', pr+'dv_'].values[0]):
#                 stack_data.loc[stack_data['info']==pr,'R+'] = 0
#             else:
#                 stack_data.loc[stack_data['info']==pr,'R+'] = 1
                               
#             if np.isnan(data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Neg', pr+'dv_'].values[0]):
#                 stack_data.loc[stack_data['info']==pr,'R-'] = 0
#             else:
#                 stack_data.loc[stack_data['info']==pr,'R-'] = 1
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Pos', pr+'pval_'].sum() < 0.1:
#                 stack_data.loc[stack_data['info']==pr+'pv','R+'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'pv','R+'] = 0
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Neg', pr+'pval_'].sum() < 0.1:
#                 stack_data.loc[stack_data['info']==pr+'pv','R-'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'pv','R-'] = 0
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Pos', 'Unexplained Haz'].sum()==1:
#                 stack_data.loc[stack_data['info']==pr+'unex','R+'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'unex','R+'] = 0
                               
#             if data_sets[da].loc[data_sets[da]['Region'] == regions[r]+'_Neg', 'Unexplained Haz'].sum()==1:
#                 stack_data.loc[stack_data['info']==pr+'unex','R-'] = 1
#             else:
#                 stack_data.loc[stack_data['info']==pr+'unex','R-'] = 0
        
        
        
        
        
#         df_teles=pd.DataFrame(columns= ['limit','R+', 'R', 'R-'])

        
#         df_teles['limit'] = ['AMO_up', 'AMO_bot','AMO_sig', 'AMO_unex_bot', 'AMO_unex_up', 'AMO_col','GMT_up', 'GMT_bot', 'GMT_sig',
#                              'GMT_unex_bot', 'GMT_unex_up','GMT_col','NAO_up', 'NAO_bot', 'NAO_sig', 'NAO_unex_bot', 'NAO_unex_up','NAO_col',
#                              'PDO_up','PDO_bot', 'PDO_sig', 'PDO_unex_bot', 'PDO_unex_up', 'PDO_col',
#                              'ENSO_bot','ENSO_up', 'ENSO_sig', 'ENSO_unex_bot', 'ENSO_unex_up', 'ENSO_col' ]
        
#         bots = [0,0,0]
    
#         for t,tel in enumerate(predictors):
            
#             for s,sr in enumerate(regs):
                
#                 df_teles.loc[df_teles['limit']==tel+'_bot', sr] = bots[s]
                
#                 df_teles.loc[df_teles['limit']==tel+'_unex_bot', sr ]= bots[s]
                
                
                    
#                 if stack_data.loc[stack_data['info']==tel+'pv', sr].sum() == 1:
                
#                     df_teles.loc[df_teles['limit']==tel+'_sig', sr ]= 'k'

#                 else:
#                     df_teles.loc[df_teles['limit']==tel+'_sig', sr] = 'lightgray'
                    
#                 if stack_data.loc[stack_data['info']==tel, sr].sum() == 1:
                
#                     df_teles.loc[df_teles['limit']==tel+'_up', sr ]= 1
#                     df_teles.loc[df_teles['limit']==tel+'_col', sr ]= colors_preds[t]
                    
#                 else:
#                     df_teles.loc[df_teles['limit']==tel+'_up', sr] = 1
#                     df_teles.loc[df_teles['limit']==tel+'_col', sr ]= 'w'
#                     df_teles.loc[df_teles['limit']==tel+'_sig', sr] = 'lightgray'
                    
#                 if stack_data.loc[stack_data['info']==tel+'unex', sr].sum() == 1:
                    
#                     if stack_data.loc[stack_data['info']==tel, sr].sum() == 1:
                
#                         df_teles.loc[df_teles['limit']==tel+'_unex_up', sr ]= 1
#                     else:
                        
#                         df_teles.loc[df_teles['limit']==tel+'_unex_up', sr] = 0

#                 else:
#                     df_teles.loc[df_teles['limit']==tel+'_unex_up', sr] = 0
                    
#                 bots[s] += 1

                
#         print('test')
       
         

#         inds= [[1],[0,1,2],[0,1,2],[1], [0], [0], [2]]
        
#         x = [0,1,2]
        
#         colour_code = ['#5ab4ac', '#4575b4', '#d8b365', '#5ab4ac', '#4575b4', '#d8b365']
        
#         df_teles_a1 = df_teles.iloc[:,0:4]
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_bot'])[0,1:][inds[r]], color = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_col'])[0,1:][inds[r]],
#          label = 'ENSO', edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_sig'])[0,1:][inds[r]])
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_unex_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='ENSO_unex_bot'])[0,1:][inds[r]], color = 'white', hatch = '////', alpha =0.3)
        
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='PDO_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_bot'])[0,1:][inds[r]], color = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_col'])[0,1:][inds[r]], label = 'PDO',
#         edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_sig'])[0,1:][inds[r]])
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='PDO_unex_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='PDO_unex_bot'])[0,1:][inds[r]],  color = 'white', hatch = '////', alpha =0.3)
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='NAO_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_bot'])[0,1:][inds[r]], color = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_col'])[0,1:][inds[r]],
#          label = 'NAO', edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_sig'])[0,1:][inds[r]])
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='NAO_unex_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='NAO_unex_bot'])[0,1:][inds[r]],  color = 'white', hatch = '////', alpha =0.3)
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='GMT_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='GMT_bot'])[0,1:][inds[r]], color = np.array(df_teles_a1[df_teles_a1['limit']=='GMT_col'])[0,1:][inds[r]],
#         label ='GMT', edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='GMT_sig'])[0,1:][inds[r]])
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='GMT_unex_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='GMT_unex_bot'])[0,1:][inds[r]],  color = 'white', hatch = '////', alpha =0.3)
        
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='AMO_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_bot'])[0,1:][inds[r]], color = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_col'])[0,1:][inds[r]],
#         label ='AMO', edgecolor = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_sig'])[0,1:][inds[r]])
        
#         f3_ax1.bar(np.array(x)[inds[r]], np.array(df_teles_a1[df_teles_a1['limit']=='AMO_unex_up'])[0,1:][inds[r]], 
#         bottom = np.array(df_teles_a1[df_teles_a1['limit']=='AMO_unex_bot'])[0,1:][inds[r]],  color = 'white', hatch = '////', alpha =0.3)
        
        
#         f3_ax1.set_ylim(0.,6.)
        
#         if r==1 or r ==2:

#             f3_ax1.set_xlim(-1.,3)
            
#             f3_ax1.set_xticks([0,1,2])
            
#             f3_ax1.set_xticklabels([ '$R_{+}$', '$R$','$R_{-}$', '$R_{+}$',  '$R$', '$R_{-}$'], fontsize = 8)
            
#         elif r ==4 or r ==5:
#             f3_ax1.set_xlim(-1.,1)
            
#             f3_ax1.set_xticks([0])
            
#             f3_ax1.set_xticklabels([ '$R_{+}$'], fontsize = 8)
            
#         elif r ==6:
#             f3_ax1.set_xlim(1.,3)
            
#             f3_ax1.set_xticks([2])
            
#             f3_ax1.set_xticklabels([ '$R_{-}$'], fontsize = 8)
        
#         else:
#             f3_ax1.set_xlim(0,2)
            
#             f3_ax1.set_xticks([1])
            
#             f3_ax1.set_xticklabels([ '$R$'], fontsize = 8)
        
#         f3_ax1.set_yticks([0.5,1.5,2.5,3.5,4.5])
        
#         if r ==0:
#             f3_ax1.set_yticklabels(['AMO','GMT','NAO','PDO','ENSO'])
#         else:
#             f3_ax1.set_yticklabels(['','','','',''])
        
#         f3_ax1.set_title(' '+ region_abs[regions[r]], position = (0.5,0.83), fontsize = 9)
        
#         r+=1
        
enso_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='cornflowerblue')
pdo_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='gold')
nao_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='gray')
gmt_box = mpatches.Rectangle((0, 0), 1, 1, facecolor= 'sandybrown')
amo_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='seagreen')
sig_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='white', edgecolor ='k')
#unex_box = mpatches.Rectangle((0, 0), 1, 1, facecolor='white', hatch = '////')

labels=['ENSO','PDO','NAO','GMT','AMO','Significant predictor']
handles = [enso_box,pdo_box, nao_box,gmt_box, amo_box, sig_box]

f3_ax2 = fig3.add_subplot(gs[0:3,z-14:z])
          
f3_ax2.legend(handles,labels, frameon=True, fontsize = 8, loc = 'center right', edgecolor = 'k')
f3_ax2.axis('off')


       

f3_ax3 = fig3.add_subplot(gs[0:4,0:9])
years = np.arange(1951,2011)
f3_ax3.set_xlim(1950,2010)   
f3_ax3.plot(years,teleDataNorm['AMO'], label='AMO', color = 'mediumseagreen')
f3_ax3.plot(years,teleDataNorm['ENSO']+4,  color = 'cornflowerblue')
f3_ax3.plot(years,teleDataNorm['PDO']+3, color = 'gold')
f3_ax3.plot(years,teleDataNorm['NAO']+2, color = 'gray')
f3_ax3.plot(years,teleDataNorm['GMT']+1, color =  'sandybrown')

f3_ax3.set_yticks([0,1,2,3,4,5])
f3_ax3.set_yticklabels(['','','','','',''])
        

f3_ax3.axhline(y=1,linewidth=0.3, color='k', linestyle = '-', alpha = 0.8)
f3_ax3.axhline(y=2,linewidth=0.3, color='k', linestyle = '-', alpha = 0.8)
f3_ax3.axhline(y=3,linewidth=0.3, color='k', linestyle = '-', alpha = 0.8)
f3_ax3.axhline(y=4,linewidth=0.3, color='k', linestyle = '-', alpha = 0.8)
f3_ax3.axvline(x=1971,linewidth=0.3, color='k', linestyle = '-', alpha = 0.8)
f3_ax3.axvspan(1971,2010, facecolor='gainsboro')
f3_ax3.set_ylabel('Normalized indices',  fontsize = 8, labelpad=+1)
f3_ax3.set_xlabel('Time',  fontsize = 8, labelpad=+2)

f3_ax3.tick_params(axis='both', labelsize=7)
#plt.legend()

plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Figure4.png',bbox_inches = 'tight',dpi =600)
plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Figure4.svg',bbox_inches = 'tight', format = 'svg')
     
   
        
