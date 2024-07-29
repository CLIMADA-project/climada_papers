#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:49:30 2020

@author: insauer
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 22:09:32 2019

@author: insauer
"""
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy
import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import pandas as pd
import netCDF4 as nc4
import matplotlib.patches as mpatches

def smooth():
    return
    

sig = '/home/insauer/projects/RiverDischarge/Data/TrendsMedianDischarge_MK_pval.nc'
disch= '/home/insauer/projects/RiverDischarge/Data/TrendsMedianDischarge_MK.nc'

dis_bin = '/home/insauer/projects/RiverDischarge/Data/basin_trends_geo.nc'

natcat_id = '/home/insauer/projects/RiverDischarge/Data/NatID_0_25.nc'



reg_frame = pd.read_csv('/home/insauer/Climada/climada_python/data/system/NatRegIDs.csv')

reg_frame.loc[reg_frame['Reg_name']=='NAM', 'Reg_plot']= 509
reg_frame.loc[reg_frame['Reg_name']=='CHN', 'Reg_plot']= 507
reg_frame.loc[(reg_frame['Reg_name']=='EUR') |
              (reg_frame['Reg_name']=='EUA')
              , 'Reg_plot']= 505
reg_frame.loc[(reg_frame['Reg_name']=='CAR') |
              (reg_frame['Reg_name']=='LAN') |
              (reg_frame['Reg_name']=='LAS')
              , 'Reg_plot']= 503
reg_frame.loc[(reg_frame['Reg_name']=='CAS') |
              (reg_frame['ISO']=='RUS')
              , 'Reg_plot']= 501


reg_frame.loc[(reg_frame['Reg_name']=='SSA') |
              (reg_frame['Reg_name']=='SAF') 
              , 'Reg_plot']= 507

reg_frame.loc[(reg_frame['Reg_name']=='SWA') |
              (reg_frame['Reg_name']=='SEA') 
              , 'Reg_plot']= 505

reg_frame.loc[(reg_frame['Reg_name']=='ARA') |
              (reg_frame['Reg_name']=='NAF') 
              , 'Reg_plot']= 503

reg_frame.loc[reg_frame['Reg_name']=='AUS', 'Reg_plot']= 501

file = xr.open_dataset(disch)

lat  = file.lat.data
lon  = file.lon.data

#data = file.variables['p-value'][0,:,:]
#p-value =1
data_dis =file.Discharge.data[0,:,:]
#file.p-value[0,:,:].data
file.close()


fig = plt.figure(dpi=600)
fig.subplots_adjust(left=0.5, hspace = -0.2, wspace = 0.05)



ax = fig.add_subplot(211)
ax.set_anchor('W')

m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), \
  urcrnrlon=lon.max(),llcrnrlat=-60.,urcrnrlat=80., \
  resolution='c')

x, y = m(*np.meshgrid(lon,lat))

m.pcolormesh(x,y,data_dis,shading='flat',cmap=plt.cm.BrBG, vmin = -0.5, vmax = 0.5)
#m.pcolormesh(x,y,data_dis,shading='flat',cmap=plt.cm.tab10)
col = m.colorbar(location='right', size = 0.08)
col.set_label('Trend in discharge $m^{3}$/s', fontsize = 5)
col.ax.tick_params(axis="y", length = 4, labelsize =4)
# Add a coastline and axis values.

m.drawcoastlines(linewidth=0.2)
#m.drawcountries()
m.drawmapboundary(fill_color= 'white')
m.drawlsmask(ocean_color='aqua',lakes=True)

m.drawparallels(np.arange(-90.,90.,30.),labels=[1,0,0,0],fontsize=4, linewidth=0.2)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=4, linewidth=0.2)


nat_file = xr.open_dataset(natcat_id)

nat_lat  = nat_file.lat.data
nat_lon  = nat_file.lon.data

nat_grid =nat_file.NatIdGrid.data[:,:]

nat_file.close()

ax.text(-0.1, 0.95, 'a', transform=ax.transAxes, 
            size=11, weight='bold')


file = xr.open_dataset(dis_bin)

lat  = file.lat.data
lon  = file.lon.data

#data = file.variables['p-value'][0,:,:]
#p-value =1
data_dis =file.basin_trend.data[0,:,:]
#file.p-value[0,:,:].data
file.close()

data_dis[np.where(data_dis>0)]=0.2
data_dis[np.where(data_dis<0)]=1.8


for maske in range(501,511):
    
    ids = reg_frame.loc[reg_frame['Reg_plot']==maske, 'ID']
    
    sig_index = np.where(np.isin(nat_grid,ids.tolist()))
    
    data_dis[sig_index]+=maske

sig_index = np.where(data_dis<500)


data_dis[sig_index]+=501

#file_sig = xr.open_dataset(sig)




# file_sig = xr.open_dataset(sig)



# data_sig =file_sig.pvalues.data[0,:,:]
# #file.p-value[0,:,:].data
# file_sig.close()



#sig_index = np.where(data_sig<0.1)

ax = fig.add_subplot(212)
ax.set_anchor('W')
#fig, axes = plt.subplots(2, 1)
#ax = axes.flatten()[0]
#plt.title("Trends (Binary) in discharge 1971-2010")
#plt.subplots_adjust(wspace=0.2, hspace=0.08)
# Miller projection:
m=Basemap(projection='mill',lat_ts=10,llcrnrlon=lon.min(), \
  urcrnrlon=lon.max(),llcrnrlat=-60.,urcrnrlat=80., \
  resolution='c')
    


#m.fillcontinents()

# convert the lat/lon values to x/y projections.

x, y = m(*np.meshgrid(lon,lat))



# data_dis[np.where(data_dis>0)]=0.9
# data_dis[np.where(data_dis<0)]=-1.15


# data_dis[sig_index]*=1.5

#f = nc4.Dataset('/home/insauer/projects/RiverDischarge/Data/Binary_Trends.nc','w', format='NETCDF4')
#dis = xr.open_dataset(path)
#lat = dis.lat.data
#lon = dis.lon.data
##lat = lat[lat_inds]
##lon = lon[lon_inds]
#
#f.createDimension('lon', len(lon))
#f.createDimension('lat', len(lat))
#f.createDimension('time', None) 
#time = f.createVariable('time', 'f8', ('time',))
#longitude = f.createVariable('lon', 'f4', ('lon'))
#latitude = f.createVariable('lat', 'f4', ('lat')) 
#binary_trends = f.createVariable('binary_trends', 'f4', ('time','lat', 'lon'))
##pval = f.createVariable('p-values', 'f4', ('time','lat', 'lon'))
#longitude[:] = lon #The "[:]" at the end of the variable instance is necessary
#latitude[:] = lat
#binary_trends[0,:,:] = data_dis
##pval[0,:,:] = p_values
#f.close()
# plot the field using the fast pcolormesh routine 
# set the colormap to jet.
vmax = np.nanmax(data_dis)
vmin = np.nanmin(data_dis)



abs_max = np.max([np.abs(vmin),np.abs(vmax)])


m.pcolormesh(x,y,data_dis,shading='flat',cmap=plt.cm.tab20b, vmin = 501, vmax = 511)







#m.pcolormesh(x,y,data_dis,shading='flat',cmap=plt.cm.tab10)
#col = m.colorbar(location='right', size = 0.08)

# Add a coastline and axis values.

m.drawcoastlines(linewidth=0.2)
m.drawcountries(linewidth=0.2)
m.drawmapboundary(fill_color='white')
m.drawlsmask(ocean_color='aqua',lakes=True)


cmap = plt.cm.get_cmap('tab20b')
rgba_dark = cmap(0.82)
rgba_light = cmap(0.97)

rgba_lam_d = cmap(0.22)
rgba_lam_l = cmap(0.37)

rgba_weu_d = cmap(0.42)
rgba_weu_l = cmap(0.57)

rgba_ssa_d = cmap(0.62)
rgba_ssa_l = cmap(0.77)

rgba_naf_d = cmap(0.22)
rgba_naf_l = cmap(0.37)

rgba_oce_d = cmap(0.02)
rgba_oce_l = cmap(0.17)

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_dark, label ='$NAM_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_light, label ='$NAM_{-}$')
leg1 = plt.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.15, 0.62))  

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_lam_d, label ='$LAM_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_lam_l, label ='$LAM_{-}$')
leg2 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.28, 0.32))  
 
pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_weu_d, label ='$EUR_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_weu_l, label ='$EUR_{-}$')
leg3 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.47, 0.77))  

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_ssa_d, label ='$SSA_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_ssa_l, label ='$SSA_{-}$')
leg4 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.53, 0.35)) 

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_naf_d, label ='$NAF_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_naf_l, label ='$NAF_{-}$')
leg5 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.735, 0.455)) 

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_weu_d, label ='$SEA_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_weu_l, label ='$SEA_{-}$')
leg6 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.8, 0.35)) 

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_oce_d, label ='$OCE_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_oce_l, label ='$OCE_{-}$')
leg7 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.9, 0.15))

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_ssa_d, label ='$EAS_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_ssa_l, label ='$EAS_{-}$')
leg8 = ax.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.95, 0.55))  

pos_box = mpatches.Rectangle((0, 0), 1,1, color=rgba_oce_d, label ='$CAS_{+}$')
neg_box = mpatches.Rectangle((0, 0), 1, 1, color=rgba_oce_l, label ='$CAS_{-}$')
plt.legend(handles = [pos_box, neg_box], frameon=False, fontsize = 2.5,bbox_to_anchor=(0.9925, 0.67))

ax.add_artist(leg1)
ax.add_artist(leg2)
ax.add_artist(leg3)
ax.add_artist(leg4)
ax.add_artist(leg5)
ax.add_artist(leg6)
ax.add_artist(leg7)
ax.add_artist(leg8)

ax.text(-0.1, 0.95, 'b', transform=ax.transAxes, 
            size=11, weight='bold')
#plt.savefig('/home/insauer/projects/Attribution/Floods/Plots/FinalPaper/PreliminaryPlots/Region_map.svg', bbox_inches = 'tight',format = 'svg')
plt.savefig('/home/insauer/projects/Attribution/Floods/Paper_NC_Resubmission_data/Main_Figures/Figure1_geo.png',  bbox_inches = 'tight', resolution=600)