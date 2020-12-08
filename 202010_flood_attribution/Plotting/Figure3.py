#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:42:33 2020

@author: insauer
"""



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 22:09:13 2020

@author: insauer
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

fig3 = plt.figure(constrained_layout=True, figsize=(8.3, 11.7))
gs = fig3.add_gridspec(20, 14)
plt.subplots_adjust(wspace=-0.2, hspace=0)


DATA_ATTR_Full = pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AttributionMetaDataRegions.csv')

DATA_ATTR = pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/AttributionMetaDataSubregions.csv')

DATA_FIT_Full= pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/VulnerabilityAdjustmentMetaDataRegions.csv')
DATA_FIT= pd.read_csv('/home/insauer/projects/NC_Submission/Data/postprocessing/VulnerabilityAdjustmentMetaDataSubregions.csv')


region_names={'GLB': 'Global (GLB)',
              'NAM':'North America (NAM)',
              'CHN':'Eastern Asia (EAS)',
              'AUS':'Oceania (OCE)',
              'LAM':'Latin America (LAM)',
              'EUR':'Europe (EUR)',
              'CAS':'Central Asia & Russia (CAS)',
              'SSAF':'South & Sub-Sahara Africa (SSA)',
              'SWEA':'South & South-East Asia (SEA)', 
              'NAFARA':'North Africa & Middle East (NAF)'
              }

region_abs={'NAM':'NAM', 
          'LAM':'LAM', 
          'EUR':'WEU',
          'NAFARA':'NAR',
          'SSAF':'SSA',
          'CAS':'CAS',
          'SWEA':'SEA', 
          'CHN':'EAS', 
          'AUS':'OCE',
          'GLB': 'GLB'}

regions = list(region_names)

#order = [9,0,1,2,]
r =0

for i in range(4):
    for j in range(3):
        
        if i == 0 and (j ==0 or j ==2):
            continue
        
        
        f3_ax1 = fig3.add_subplot(gs[5*i:5*i+4,j*5:(j*5)+4])
        

 
        
        r_lin = DATA_FIT_Full.loc[DATA_FIT_Full['Region']==regions[r], 'P_ExpVar_pred_observed'].sum()     

        r_lin_pos = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Pos', 'ExpVar_model_pred_observed'].sum()
        r_lin_neg = DATA_FIT.loc[DATA_FIT['Region']==regions[r]+'_Neg', 'ExpVar_model_pred_observed'].sum()
        
        
        data_attr_reg = DATA_ATTR_Full[DATA_ATTR_Full['Region'] == regions[r]]
        data_attr_reg_pos = DATA_ATTR[DATA_ATTR['Region'] == regions[r]+'_Pos']
        data_attr_reg_neg = DATA_ATTR[DATA_ATTR['Region'] == regions[r]+'_Neg']
        
        h8 = data_attr_reg.loc[:,'Change HnCL'].sum()
        e8 = data_attr_reg.loc[:,'Change En'].sum()
        v8 = data_attr_reg.loc[:,'Change Vn'].sum()
        i8 = data_attr_reg.loc[:,'Change In'].sum()
        n8 = data_attr_reg.loc[:,'Change Nn'].sum()
        
        h8_up = data_attr_reg.loc[:,'Change HnCLup'].sum()

        
        h8_bot = data_attr_reg.loc[:,'Change HnCLbot'].sum()

        
        
        h_pos8 =  data_attr_reg_pos.loc[:,'Change HnCL'].sum()
        e_pos8 = data_attr_reg_pos.loc[:,'Change En'].sum()
        v_pos8 = data_attr_reg_pos.loc[:,'Change Vn'].sum()
        i_pos8 = data_attr_reg_pos.loc[:,'Change In'].sum()
        n_pos8 = data_attr_reg_pos.loc[:,'Change Nn'].sum()
        
        h_pos8_up =  data_attr_reg_pos.loc[:,'Change HnCLup'].sum()

        h_pos8_bot =  data_attr_reg_pos.loc[:,'Change HnCLbot'].sum()
 
        
        
        h_neg8 =  data_attr_reg_neg.loc[:,'Change Hn'].sum()
        e_neg8 =  data_attr_reg_neg.loc[:,'Change En'].sum()
        v_neg8 =  data_attr_reg_neg.loc[:,'Change Vn'].sum()
        i_neg8 =  data_attr_reg_neg.loc[:,'Change In'].sum()
        n_neg8 =  data_attr_reg_neg.loc[:,'Change Nn'].sum()
        

        
        h_sign = data_attr_reg.loc[:,'Sign H'].sum()
        h_sign_pos = data_attr_reg_pos.loc[:,'Sign H'].sum()
        h_sign_neg = data_attr_reg_neg.loc[:,'Sign H'].sum()

        x=[0,1,2,3,5,6,7,8,10,11,12,13]
        colour_code = ['#4575b4','orange','brown','#8856a7',
                       '#4575b4','orange','brown','#8856a7',
                       '#4575b4','orange','brown','#8856a7']
        y1= [h8,e8,v8,i8,h_pos8,e_pos8,v_pos8,i_pos8,h_neg8,e_neg8,v_neg8,i_neg8]


        
        # ax1_lims_up = [10,10,10,3,3,8,10,35,20,10]
        
        # ax1_lims_low = [-6,-10,-5,-2.2,-3,-8,-10,-5,-10,-10]
        
        # ax2_lims_up = [10,10,10,2.2,3,8,10,35,15,10]
        
        # ax2_lims_low = [-6,-10,-5,-2.2,-3,-8,-10,-5,-10,-10]
        
        # ax1_ticks_up = [5,5,5,1,2,5,5,20,10,5]
        
        # ax1_ticks_low = [-5,-5,-5,-1,-2,-5,-5,-5,-5,-5]
        
        # ax2_ticks_up = [5,5,5,1,2,5,5,20,10,5]
        
        # ax2_ticks_low = [-5,-5,-5,-1,-2,-5,-5,-5,-5, -5]
        
        # ax1_labels_up = ['5%', '5%','5%','1%','2%','5%', '5%', '20%','10%','5%']
        
        # ax1_labels_low = ['-5%','-5%','-5%','-1%','-2%','-5','-5%','','-5%','-5%']
        
        # ax2_labels_up = ['5%','5%','5%','5%','2%', '5%','5%', '20%','10%','5%']
        
        # ax2_labels_low = ['-5%','-5%','-5%','-1%','-2%','-5%','-5%','','-5%','-5%']
        
        if r_lin_pos > 0.2:
             f3_ax1.axvspan(4,9, facecolor='gainsboro')
        #     f3_ax2.axvspan(-1,0.5, facecolor='gainsboro')
        #     f3_ax1.axvspan(3,4.5, facecolor='gainsboro')
        #     f3_ax2.axvspan(3,4.5, facecolor='gainsboro')
        #     #f3_ax2.axvspan(4, 9, facecolor='gainsboro')
            
        if r_lin > 0.2:
             f3_ax1.axvspan(-1,4, facecolor='gainsboro')
        #     f3_ax2.axvspan(1.5,3., facecolor='gainsboro')
        #     f3_ax1.axvspan(5.5,7., facecolor='gainsboro')
        #     f3_ax2.axvspan(5.5,7., facecolor='gainsboro')
            
        if r_lin_neg > 0.2:
            
            f3_ax1.axvspan(9,14, facecolor='gainsboro')
        #     f3_ax2.axvspan(0.5,1.5, facecolor='gainsboro')
        #     f3_ax1.axvspan(4.5,5.5, facecolor='gainsboro')
        #     f3_ax2.axvspan(4.5,5.5, facecolor='gainsboro')
        
        
        for a in range(12):
            
            if not (r==1 and a>3 and a<8) :
                
                if ((a == 0)and (h_sign > 0.1)) or ((a == 4) and (h_sign_pos > 0.1)) or ((a == 8)and (h_sign_neg > 0.1)):
                    
                    f3_ax1.bar(x[a],y1[a], color = "#4575b4", alpha = 0.5 )
                
                else:
                    
                    f3_ax1.bar(x[a],y1[a], color = colour_code[a], #fmt='o',
                                #elinewidth = 0.5,
                            #yerr=np.reshape(np.array([y1err_bot[a],y1err_up[a]]),(2,1)),
                            ecolor='black')#,mfc =colour_code[a])


                
            

   
            # f3_ax1.axhline(y=0.85, xmin =0, xmax = 1/4,linewidth=0.5, color='white')
            # f3_ax1.axhline(y=0.75, xmin =0, xmax = 1/4,linewidth=0.5, color='white')
            # f3_ax1.axhline(y=0.85, xmin =1/2, xmax = 3/4,linewidth=0.5, color='white')
            # f3_ax1.axhline(y=0.75, xmin =1/2, xmax = 3/4,linewidth=0.5, color='white')
                    
            
        # y_max = f3_ax1.get_ylim()[1]*1.6
        # y_min = f3_ax1.get_ylim()[0]*1.4
        
        #f3_ax1.set_ylim(ax1_lims_low[r],ax1_lims_up[r])
        
        #f3_ax1.set_yticks([ax1_ticks_low[r] , 0, ax1_ticks_up[r]])
        
        #f3_ax1.set_yticklabels([ax1_labels_low[r],'0', ax1_labels_up[r] ],fontsize =6.5)
        
        f3_ax1.set_xticklabels(['','', '', '' ],fontsize =6.5)
        
        if r in [0,1,4,7]:
            f3_ax1.set_ylabel('Relative trend in %/year', fontsize=8, labelpad = 2.5)
            
        

        f3_ax1.set_xlim(-1.,14)
        
        if r <7:
            f3_ax1.set_ylim(-5, 12)
            f3_ax1.set_yticks([-5.,0, 5, 10 ])
            f3_ax1.set_yticklabels(['-5%','0%', '5%','10%'],fontsize =7)
        else:
            f3_ax1.set_ylim(-5, 33)
            f3_ax1.set_yticks([0, 10, 20, 30])
            f3_ax1.set_yticklabels(['0%','10%', '20%', '30%' ],fontsize =7)
        
        nes = [n8,n_pos8,n_neg8]
        
        xn = [4,9,14]
        
        xx= [3,8,13]
        
        for ns in range(3):
            if nes[ns]>0:
                lims = np.abs(f3_ax1.get_ylim()[0])/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) + nes[ns]/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0])
                zerlim = np.abs(f3_ax1.get_ylim()[0])/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0])
                f3_ax1.axhline(y=nes[ns], xmin=(xn[ns]-0.4)/15, xmax=(xn[ns]+0.4)/15,linewidth=1, color='k')
                f3_ax1.axvline(x=xx[ns]-0.4, ymin=zerlim, ymax=lims,linewidth=1., color='k')
                f3_ax1.axvline(x=xx[ns]+0.4, ymin=zerlim, ymax=lims,linewidth=1., color='k')
            else:
                lims = np.abs(f3_ax1.get_ylim()[0])/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0]) + nes[ns]/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0])
                zerlim = np.abs(f3_ax1.get_ylim()[0])/(f3_ax1.get_ylim()[1]-f3_ax1.get_ylim()[0])
                f3_ax1.axhline(y=nes[ns], xmin=(xn[ns]-0.4)/15, xmax=(xn[ns]+0.4)/15,linewidth=1, color='k')
                f3_ax1.axvline(x=(xx[ns]-0.4), ymin=lims, ymax=zerlim, linewidth=1, color='k')
                f3_ax1.axvline(x=(xx[ns]+0.4), ymin=lims, ymax=zerlim, linewidth=1, color='k')

        
        

        f3_ax1.axvline(x=4,linewidth=0.3, color='k', linestyle = '-', alpha = 1)
        f3_ax1.axvline(x=9,linewidth=0.3, color='k', linestyle = '-', alpha = 1)
        
        
        f3_ax1.set_xlabel('R            $R_{+}$          $R_{-}$',  fontsize = 10.5, labelpad=-10)
        f3_ax1.xaxis.set_label_position('top')
        #f3_ax2.axhline(y=-5,linewidth=0.2, color='k', linestyle = '-.')
        f3_ax1.axhline(y=0,linewidth=0.3, color='k', linestyle = '-')
        
        if (i == 1) and (j==1):
            f3_ax1.set_title(' '+ region_names[regions[r]], position = (0.5,1.), fontsize = 9)
            
        elif (i == 2) and ((j ==0) or (j==1)):
            
            f3_ax1.set_title(' '+ region_names[regions[r]], position = (0.5,1.), fontsize = 9)
        else:
            f3_ax1.set_title(' '+ region_names[regions[r]], position = (0.5,1.), fontsize = 9)
        
        
        
        r+=1

f3_ax4 = fig3.add_subplot(gs[0:5,10:14])
handles, labels = f3_ax4.get_legend_handles_labels()
f3_ax4.axis('off')





cli_box = mpatches.Circle((0.5,0.5), 1, facecolor="#4575b4", label ='Significant Climate contribution')
cli2_box = mpatches.Circle((0.5,0.5), 1, facecolor="#4575b4", alpha = 0.5, label ='Insignificant Climate contribution')
exp_box = mpatches.Circle((0.5,0.5), 1, color='orange', label ='Exposure contribution')
vul_box = mpatches.Circle((0.5,0.5), 1, facecolor='brown', label ='Vulnerability contribution')

imp_box = mpatches.Circle((0.5,0.5), 1, facecolor='#8856a7', label ='Total Damage Trend')





# circ = mpatches.Circle((0.5,0.5), 1,facecolor='#5ab4ac')
# sig1 = mpatches.Circle((0.5,0.5), 1,facecolor='#5ab4ac')
# sig2 = mpatches.Circle((0.5,0.5), 1,facecolor='#5ab4ac')
# sig3 = mpatches.Circle((0.5,0.5), 1,facecolor='#5ab4ac')

# labels=['Region with positive \n discharge trend $R_{+}$','Full world region $R$',
#         'Region with negative \n discharge trend $R_{-}$']
        #'circ','sig1','sig2','sig3']

f3_ax4.legend(handles = [cli_box, cli2_box, exp_box,vul_box, imp_box], frameon=True, fontsize = 8, loc = 'center', edgecolor = 'k')  

# plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Alt_Fig_2b.png',bbox_inches = 'tight',dpi =600)
# plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Alt_Fig_2b.svg',bbox_inches = 'tight', format = 'svg')

#f3_ax1.set_title('gs[0, :]')
#f3_ax2 = fig3.add_subplot(gs[1, :-1])
#f3_ax2.set_title('gs[1, :-1]')
#f3_ax3 = fig3.add_subplot(gs[1:, -1])
#f3_ax3.set_title('gs[1:, -1]')
#f3_ax4 = fig3.add_subplot(gs[-1, 0])
#f3_ax4.set_title('gs[-1, 0]')
#f3_ax5 = fig3.add_subplot(gs[-1, -2])
#f3_ax5.set_title('gs[-1, -2]')