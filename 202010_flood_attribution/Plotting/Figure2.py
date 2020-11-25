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
gs = fig3.add_gridspec(40, 15)
plt.subplots_adjust(wspace=0., hspace=0)

DATA_TSFull= pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/AttributionTimeSeriesRegions.csv')
DATA_TS= pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/AttributionTimeSeriesSubregions.csv')

DATA_FIT_Full= pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/VulnerabilityAdjustmentMetaDataRegions.csv')
DATA_FIT= pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/VulnerabilityAdjustmentMetaDataSubregions.csv')


DATA_ATTR_Full = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/AttributionMetaDataRegions.csv')

DATA_ATTR = pd.read_csv('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Aggregation_attribution/AttributionMetaDataSubregions.csv')

region_names={'GLB': 'Global (GLB)',
              'NAM':'North America (NAM)',
              'CHN':'Eastern Asia (EAS)',
              'AUS':'Oceania (OCE)',
              'LAM':'Latin America (LAM)',
              'EUR':'Europe (EUR)',
              'SSAF':'South & Sub-Sahara Africa (SSA)',
              'SWEA':'South & South-East Asia (SEA)',
              'CAS':'Central Asia & Russia (CAS)',
              'NAFARA':'North Africa & Middle East (NAF)',
              }

region_abs={'GLB': 'GLB',
            'NAM':'NAM', 
            'CHN':'EAS',
            'LAM':'LAM', 
            'EUR':'EUR',
            'AUS':'OCE',
            'CAS':'CAS',
            'SSAF':'SSA',
            'SWEA':'SEA', 
            'NAFARA': 'NAF'}

regions = list(region_names)
r =0

for i in range(10):
    for j in range(4):
        
        DATA_regionFull = DATA_TSFull[(DATA_TSFull['Region']==regions[r]) & 
                      (DATA_TSFull ['Year']<2011) & (DATA_TSFull ['Year']>1979)]
        
        DATA_region = DATA_TS[(DATA_TS['Region']==regions[r]) & 
                      (DATA_TS['Year']<2011) & (DATA_TS['Year']>1979)]
        
        if j<3:
        
            f3_ax1 = fig3.add_subplot(gs[4*i:4*i+4,j*5:(j*5)+5])
            
            if j ==0:
            
            
            
        
                # f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Impact_Pred_1thrd']), color='#8856a7', alpha = 0.5, linewidth = 1.)
                # f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Impact_Pred_2thrd']), color='#8856a7', alpha = 0.5, linewidth = 1.)
                f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['natcat_flood_damages_2005_CPI']), label='$D_{Obs}$', color='black', linewidth = 1.) 
                f3_ax1.scatter(DATA_regionFull['Year'], np.log10(DATA_regionFull['natcat_flood_damages_2005_CPI']), color='black', marker = '_', s = 3) 
                f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Norm_Impact_Pred']), label='$D_{Full}$', color='#8856a7', linewidth = 1.)
                
                f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Norm_Impact_2y_offset']), label='$D_{CliExp}$', color='#ff7f00', linewidth = 1.)
            
                f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Norm_ImpFix_2y_offset']), label='$D_{Cli}$', color='#4575b4', linewidth = 1.)
                
                #f3_ax1.plot(DATA_regionFull['Year'], np.log10(DATA_regionFull['Norm_Imp2010_2y_offset']), label='$Loss2010_{Haz}$', color='mediumseagreen', linewidth = 1., linestyle ='--', alpha = 0.5)
        
                
                
                #f3_ax1.fill_between(DATA_regionFull['Year'],np.log10(DATA_regionFull['Impact_Pred_1thrd']) , np.log10(DATA_regionFull['Impact_Pred_2thrd']), color='#8856a7', alpha=0.2, linewidth = 1.)
                if i ==5:
                
                    f3_ax1.set_title(' '+ region_names[regions[r]], position = (0.5,0.78), fontsize = 7)
                
                else:
                    
                    f3_ax1.set_title(' '+ region_names[regions[r]], position = (0.5,0.78), fontsize = 8)
                
                if i ==0 and j ==0:
                    handles, labels = f3_ax1.get_legend_handles_labels()
                    leg =f3_ax1.legend(handles[:2], labels[:2], loc ='lower left', labelspacing = 0.1, frameon=True, fontsize = 8, handlelength = 1.1 ) 
                    f3_ax1.legend(handles[2:], labels[2:], loc ='lower right', labelspacing = 0.1, frameon=True, fontsize = 8,  handlelength = 1.1)
                    f3_ax1.add_artist(leg)
                # f3_ax1.plot(DATA_regionFull['Year'], 
                #     np.log10(DATA_regionFull['model_flood_damages_onethird_quantile_2005flospros']),
                #     color='#4575b4', alpha = 0.5, linestyle = '--', linewidth = 0.5)
                
                
                # f3_ax1.plot(DATA_regionFull['Year'],
                #     np.log10(DATA_regionFull['model_flood_damages_twothird_quantile_2005flospros']),
                #     color='#4575b4', alpha = 0.5, linestyle = '--', linewidth = 0.5)
        
                # f3_ax1.fill_between(DATA_regionFull['Year'],np.log10(DATA_regionFull['model_flood_damages_onethird_quantile_2005flospros']) ,
                #              np.log10(DATA_regionFull['model_flood_damages_twothird_quantile_2005flospros']),
                #              color='#4575b4', alpha=0.2)
            
                # f3_ax1.plot(DATA_regionFull['Year'], 
                #     np.log10(DATA_regionFull['model_flood_damages_onethird_quantile']),
                #     color='#ff7f00', alpha = 0.5, linestyle = '--', linewidth = 0.5)
                # f3_ax1.plot(DATA_regionFull['Year'],
                #     np.log10(DATA_regionFull['model_flood_damages_twothird_quantile']),
                #     color='#ff7f00', alpha = 0.5, linestyle = '--', linewidth = 0.5)
        
                # f3_ax1.fill_between(DATA_regionFull['Year'],np.log10(DATA_regionFull['model_flood_damages_onethird_quantile']) ,
                #              np.log10(DATA_regionFull['model_flood_damages_twothird_quantile']),
                #              color='#ff7f00', alpha=0.2)
                r_rang = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'SP_corrcoef_pred_observed'].sum()

                
                r_log = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'P_LOG_ExpVar_pred_observed'].sum()
                r_lin = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'P_ExpVar_pred_observed'].sum()

                #r2 = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'New_explained_variance'].sum()
               
            else:
                
                if j ==1:
                    dis = 'Pos'
                    f3_ax1.set_title('{}'.format(region_abs[regions[r]])+'$_{+}$', position = (0.5,0.78), fontsize = 8)
                else:
                    dis = 'Neg'
                    f3_ax1.set_title('{}'.format(region_abs[regions[r]])+'$_{-}$', position = (0.5,0.78), fontsize = 8)
                
                
                #f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['Impact_Pred_1thrd_{}'.format(dis)]), color='#8856a7', alpha = 0.5, linewidth = 1.)
                #f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['Impact_Pred_2thrd_{}'.format(dis)]), color='#8856a7', alpha = 0.5, linewidth = 1.)
                
                f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['natcat_damages_2005_CPI_{}'.format(dis)]), label='Observed Flood Losses (NatCat)', color='black', linewidth = 1.) 
                f3_ax1.scatter(DATA_region['Year'], np.log10(DATA_region['natcat_damages_2005_CPI_{}'.format(dis)]), label='Observed Flood Losses (NatCat)', color='black', marker = '.', s = 3) 
            
                
                if not (i==1 and dis=='Pos'):  
                    
                 
                    f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['NormExp_Impact_2y{}_offset'.format(dis)]), label='$Loss_{HazExp}$', color='#ff7f00', linewidth = 1.)
                
                    f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['NormHaz_ImpFix_2y{}_offset'.format(dis)]), label='$Loss_{Haz}$', color='#4575b4', linewidth = 1.)
                    
                    #f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['NormHaz_Imp2010_2y{}_offset'.format(dis)]), label='$Loss2010_{Haz}$', color='mediumseagreen', linewidth = 1.,linestyle ='--', alpha = 0.5)
            
                    f3_ax1.plot(DATA_region['Year'], np.log10(DATA_region['Norm_Impact_Pred_{}'.format(dis)]), label='$Loss_{Full}$', color='#8856a7', linewidth = 1.)
                    
                # f3_ax1.fill_between(DATA_region['Year'],np.log10(DATA_region['Impact_Pred_1thrd_{}'.format(dis)]) , np.log10(DATA_region['Impact_Pred_2thrd_{}'.format(dis)]), color='#8856a7', alpha=0.2, linewidth = 1.)
                
        
                # f3_ax1.plot(DATA_region['Year'], 
                #     np.log10(DATA_region['ImpFix_2y{}_onethird_quantile'.format(dis)]),
                #     color='#4575b4', alpha = 0.5, linestyle = '--', linewidth = 0.5)
                
                
                # f3_ax1.plot(DATA_region['Year'],
                #     np.log10(DATA_region['ImpFix_2y{}_twothird_quantile'.format(dis)]),
                #     color='#4575b4', alpha = 0.5, linestyle = '--', linewidth = 0.5)
        
                # f3_ax1.fill_between(DATA_region['Year'],np.log10(DATA_region['ImpFix_2y{}_onethird_quantile'.format(dis)]) ,
                #              np.log10(DATA_region['ImpFix_2y{}_onethird_quantile'.format(dis)]),
                #              color='#4575b4', alpha=0.2)
            
                # f3_ax1.plot(DATA_region['Year'], 
                #     np.log10(DATA_region['Impact_2y{}_onethird_quantile'.format(dis)]),
                #     color='#ff7f00', alpha = 0.5, linestyle = '--', linewidth = 0.5)
                # f3_ax1.plot(DATA_region['Year'],
                #     np.log10(DATA_region['Impact_2y{}_twothird_quantile'.format(dis)]),
                #     color='#ff7f00', alpha = 0.5, linestyle = '--', linewidth = 0.5)
        
                # f3_ax1.fill_between(DATA_region['Year'],np.log10(DATA_region['Impact_2y{}_onethird_quantile'.format(dis)]) ,
                #              np.log10(DATA_region['Impact_2y{}_twothird_quantile'.format(dis)]),
                #              color='#ff7f00', alpha=0.2)
                
                
                
                r_lin = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_'+dis, 'ExpVar_model_pred_observed'].sum()
                r_log = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_'+dis, 'LOG_ExpVar_pred_observed'].sum()
                
                #r2 = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_'+dis, 'New_explained_variance'].sum()
                r_rang = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_'+dis, 'KT_corrcoef_pred_observed'].sum()

                
            #text_LOG = 'R²='+str(round(r_log*100,1))+ '% (LOG)'
            text_lin = 'R²='+str(round(r_lin*100,1))+ '%'
            text_log = 'R²='+str(round(r_log*100,1))+ '% (LOG)'
            #text_lin_rm = 'R²3yrm='+str(round(r_lin_rm*100,1))+ '% (LIN)'
            text_rang = 'r_rank='+str(round(r_rang,3))
            
            
            if r_lin> 0.2:
                
                f3_ax1.set_facecolor('gainsboro')
            
            f3_ax1.set_yticks([6, 8, 10])
            f3_ax1.set_yticklabels(['','',''])
            if i in [2]:
                f3_ax1.set_ylim((5.5, 12))
            
            elif i in [1]:
                f3_ax1.set_ylim((4.5, 11.5))
                
            elif i in [4]:
                f3_ax1.set_ylim((5, 11.5))
            
            elif i in [3]:
                f3_ax1.set_ylim((4.5, 11))
                
            elif i in [6]:
                f3_ax1.set_ylim((4, 10.5))
                
            elif i in [7]:
                f3_ax1.set_ylim((5.5, 11.5))
            
            elif i in [8]:
                f3_ax1.set_ylim((4.5, 11))
            
            elif i in [0]:
                f3_ax1.set_ylim((7, 12))
                f3_ax1.set_yticks([8, 10])
                
                if j == 0:
                    f3_ax1.set_yticklabels(['8','10'])
            else:
                f3_ax1.set_ylim((5., 11.5))
            
            f3_ax1.set_xlim((1978 ,2013))
            if not (i==1 and j==1):
                f3_ax1.annotate( xy=(1990, f3_ax1.get_ylim()[0]+0.08*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s=text_lin, fontsize=7 )
            #f3_ax1.annotate( xy=(1998, f3_ax1.get_ylim()[0]+0.15*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s=text_log, fontsize=7 )
            
            
            #f3_ax1.annotate( xy=(1980, f3_ax1.get_ylim()[0]+0.24*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s=text_rang, fontsize=7 )
            #f3_ax1.annotate( xy=(1980, f3_ax1.get_ylim()[0]+0.15*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s=text_lin, fontsize=7 )
            #f3_ax1.annotate( xy=(1980, f3_ax1.get_ylim()[0]+0.06*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s=text_r290, fontsize=6 )
            # if i==0:
            #     f3_ax1.annotate( xy=(1980, f3_ax1.get_ylim()[0]+0.88*(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) ) ,s='a', fontsize=8, fontweight = 'bold')
            
            
            
            f3_ax1.set_xticks([1980,1990,2000,2010])
            f3_ax1.set_xticklabels(['','','',''])
            
            
            if i == 9:
                f3_ax1.set_xticklabels(['1980','1990','2000','2010'], fontsize =7)
            
            if i ==4 and j ==0:
                f3_ax1.set_ylabel('LOG10 Damages (2005 USD)', fontsize=8, labelpad=-1)
            
            if j == 0 and i !=0:
                f3_ax1.set_yticklabels(['6','8','10'], fontsize=8)
            
            if i == 4 and j ==0:
                f3_ax1.set_xticklabels(['1980','1990','2000', '2010'],fontsize=8)
                f3_ax1.set_xlabel('Year', fontsize=9, labelpad=-2)
            
            if i == 4 and j ==0:
                f3_ax1.set_xticklabels(['1980','1990','2000', '2010'],fontsize=8)
                f3_ax1.set_xlabel('Year', fontsize=9, labelpad=-2)
            
            if i == 3 and j ==1:
                f3_ax1.set_xticklabels(['1980','1990','2000', '2010'],fontsize=8)
                f3_ax1.set_xlabel('Year', fontsize=9, labelpad=-2)
            
            f3_ax1.tick_params(axis="x", direction = 'in',length = 4)
            
            handles, labels = f3_ax1.get_legend_handles_labels() 
        
    else:
        

        
        r_linPos = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Pos', 'ExpVar_model_pred_observed'].sum()
        r_linNeg = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Neg', 'ExpVar_model_pred_observed'].sum()
        r_lin = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'P_ExpVar_pred_observed'].sum()
        
        
        
        r_rankPos = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Pos', 'KT_corrcoef_pred_observed'].sum()
        r_rankNeg = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Neg', 'KT_corrcoef_pred_observed'].sum()
        r_rank = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'SP_corrcoef_pred_observed'].sum()
        

        
        
    r+=1
        
        
        # indc = np.arange(4)
    
        
        # if data_attr_reg.loc[:,'Sign H'].sum() < 0.1:
        #     sigH = '*'
        #     if data_attr_reg.loc[:,'Sign H'].sum() < 0.05:
        #         sigH=' ** '
        #         if data_attr_reg.loc[:,'Sign H'].sum() < 0.01:
        #             sigH='***'
        # else:
        #     sigH=''
        
        # if data_attr_reg.loc[:,'Sign E'].sum() < 0.1:
        #     sigE = '*'
        #     if data_attr_reg.loc[:,'Sign E'].sum() < 0.05:
        #         sigE=' ** '
        #         if data_attr_reg.loc[:,'Sign E'].sum() < 0.01:
        #             sigE='***'
        # else:
        #     sigE=''
            
        
        # if data_attr_reg.loc[:,'Sign V'].sum() < 0.1:
        #     sigV = '*'
        #     if data_attr_reg.loc[:,'Sign V'].sum() < 0.05:
        #         sigV=' ** '
        #         if data_attr_reg.loc[:,'Sign V'].sum() < 0.01:
        #             sigV='***'
        # else:
        #     sigV=''
    
        
        # if data_attr_reg.loc[:,'Sign I'].sum() < 0.1:
        #     sigI = '*'
        #     if data_attr_reg.loc[:,'Sign I'].sum() < 0.05:
        #         sigI=' ** '
        #         if data_attr_reg.loc[:,'Sign I'].sum() < 0.01:
        #             sigI='***'
        # else:
        #     sigI=''
        # width = 0.5
        # f3_ax2.bar(indc,pl, color= ['#4575b4', '#fee090', '#d73027', '#8856a7'],tick_label= ['Haz', 'Exp', 'Vul', 'Full'],width =width, align = 'center') 
        # f3_ax2.tick_params(axis="both", direction = 'in',length = 3)
        # f3_ax2.set_yticks([-5,5])
        # f3_ax2.tick_params(axis = 'y', direction = 'in', length = 3)
        
        

        # #f3_ax2.set_ylabel('Change in damage',  fontsize = 5, labelpad=+2)
        # f3_ax2.yaxis.tick_right()
        # f3_ax2.yaxis.set_label_position("right")
        
        # if i == 2 and j == 1:
        #     f3_ax2.set_ylabel('Trends % / year', fontsize=8)
        
        
        #f3_ax2.axhline(y=f3_ax2.get_ylim()[0]+((f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0])*0.18),linewidth=0.5, color='k')
        
        # if i == 4 and j ==0:
        #     f3_ax2.set_xticklabels(['$Haz_{t}$', '$Exp_{t}$', '$Vul_{t}$', '$Full_{t}$'],rotation = 90., fontsize = 8)
        #     f3_ax2.tick_params(axis="x", direction = 'in',length = 3, pad = 5)
        
        # if i == 4 and j ==1:
        #     f3_ax2.set_xticklabels(['$Haz_{t}$', '$Exp_{t}$', '$Vul_{t}$', '$Full_{t}$'],rotation = 90., fontsize = 8)
        #     f3_ax2.tick_params(axis="x", direction = 'in',length = 3, pad = 5)
        
        # if pl[0]<0:
        #     hh = pl[0] - 1.2
        #     hh2= hh -0.5
        # else:
        #     hh = pl[0] - 0.5
        #     hh2= hh+ 0.5
        
        # if pl[1]<0:
        #     he = pl[1] - 1.2
        #     he2 = he - 0.5
        # else:
        #     he = pl[1] - 0.5
        #     he2 = he +0.5
            
        # if pl[2]<0:
        #     hv = pl[2] - 1.2
        #     hv2 = hv - 0.5
        # else:
        #     hv = pl[2] - 0.5
        #     hv2 = hv + 0.5
            
        # if pl[3]<0:
        #     hi = pl[3] - 1.2
        #     hi2 = hi - 0.5
        # else:
        #     hi = pl[3] - 0.5
        #     hi2 = hi + 0.5
        
        # #f3_ax2.annotate( xy=(0-width, f3_ax2.get_ylim()[0]+((f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0])*0.02)) ,s='Haz', fontsize=3.5, color= 'k',rotation = 70.)
        
        # if sigH == '***':   
        #     f3_ax2.annotate( xy=(0, hh),s=' ** ', fontsize=7, color= '#4575b4', ha='center')
        #     f3_ax2.annotate( xy=(0, hh2),s='*', fontsize=7, color= '#4575b4', ha='center')
        # else:   
        #     f3_ax2.annotate( xy=(0, hh),s=sigH, fontsize=7, color= '#4575b4', ha='center')
        
        # if sigE == '***':
        #     f3_ax2.annotate( xy=(1, he),s=' ** ', fontsize=7, color= '#fee090', ha = 'center')
        #     f3_ax2.annotate( xy=(1, he2),s=' * ', fontsize=7, color= '#fee090', ha = 'center')
        # else:
        #     f3_ax2.annotate( xy=(1, he),s=sigE, fontsize=7, color= '#fee090', ha='center')
        
        # if sigV == '***':
        #     f3_ax2.annotate( xy=(2, hv),s=' ** ', fontsize=7, color= '#d73027', ha = 'center')
        #     f3_ax2.annotate( xy=(2, hv2),s=' * ', fontsize=7, color= '#d73027', ha = 'center')
        # else:
        #     f3_ax2.annotate( xy=(2, hv),s=sigV, fontsize=7, color= '#d73027', ha='center')
        
        # if sigI == '***':
        #     f3_ax2.annotate( xy=(3, hi),s=' ** ', fontsize=7, color= '#8856a7', ha = 'center')
        #     f3_ax2.annotate( xy=(3, hi2),s=' * ', fontsize=7, color= '#8856a7', ha = 'center')
        # else:
        #     f3_ax2.annotate( xy=(3, hi),s=sigI, fontsize=7, color= '#8856a7', ha='center')
        
        #f3_ax2.annotate( xy=(1-width, f3_ax2.get_ylim()[0]+((f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0])*0.02)) ,s='Exp', fontsize=3.5, color= 'k',rotation = 70.)
        
        #f3_ax2.annotate( xy=(2-width/2, f3_ax2.get_ylim()[0]+((f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0])*0.02)) ,s='Vul', fontsize=3.5, color= 'k',rotation = 70.)

        #f3_ax2.annotate( xy=(3-width/2, f3_ax2.get_ylim()[0]+((f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0])*0.02)) ,s='Imp', fontsize=3.5, color= 'k',rotation = 70.)
#         if i==0:
#             f3_ax2.annotate( xy=(-0.7, f3_ax2.get_ylim()[0]+0.88*(f3_ax2.get_ylim()[1]-f3_ax2.get_ylim()[0]) ) ,s='b', fontsize=8, fontweight = 'bold')
    

plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Alt_Fig2a.png',bbox_inches = 'tight',dpi =600)
plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Alt_Fig2a.svg',bbox_inches = 'tight', format = 'svg')

#f3_ax1.set_title('gs[0, :]')
#f3_ax2 = fig3.add_subplot(gs[1, :-1])
#f3_ax2.set_title('gs[1, :-1]')
#f3_ax3 = fig3.add_subplot(gs[1:, -1])
#f3_ax3.set_title('gs[1:, -1]')
#f3_ax4 = fig3.add_subplot(gs[-1, 0])
#f3_ax4.set_title('gs[-1, 0]')
#f3_ax5 = fig3.add_subplot(gs[-1, -2])
#f3_ax5.set_title('gs[-1, -2]')