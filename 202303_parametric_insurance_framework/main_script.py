"""
This file is part of CLIMADA.

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.

---

Main script to design parametric insurance products and create plots
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import geopandas
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib.cm as cm_mp
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs

from climada.hazard.base import Hazard
from climada.hazard import TCTracks
from climada.hazard import TropCyclone
from climada.entity import ImpactFuncSet
from climada.entity.impact_funcs.trop_cyclone import ImpfSetTropCyclone
from climada.entity.impact_funcs.storm_europe import ImpfStormEurope
from climada.engine import Impact
import climada.util.coordinates as u_coord

import functions as fct


project_path = '/Users/carmenst/Documents/CLIMADA/own_projects/parametric_casestudy/data/' 
figures_path = '/Users/carmenst/Library/CloudStorage/Dropbox/Aplicaciones/Overleaf/WP1_Environment_Systems_and_Decisions/sn-article-template/art/'

#Plotting settings
plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = '11'
cm = 1/2.54  # centimeters in inches
fig_width = 17.4 #cm
max_fig_height = 23.4 #cm



"""Case study 1: Mozambique """
country_iso = 'MOZ'
country_name = 'mozambique'
basin = "SI"
impf_id = 5
bounds = [30, -26, 41, -10] # [lon_min , lat_min, lon_max, lat_max]
hazard_name = 'tc'
radius = np.array([10, 20, 30, 40, 50])*1000 #meter
steps = 25
nr_categories = 4
categories = np.array([43, 50, 59, 70])

hazard, hazard_present_folder, hospitals, insured = fct.wrapper_haz_exp(hazard_name, country_iso, country_name, impf_id, 
                                             project_path, basin, bounds)

fct.wrapper_generate_haz_in_a_circle(hazard, hospitals, hazard_name, hazard_present_folder, radius)

hospitals_plot = hospitals.gdf.set_crs(epsg=4326).to_crs(epsg=3857) #cartesian coordinate system
admin1_info, admin1_geo = u_coord.get_admin1_info(country_iso)
b = admin1_info[country_iso]

reg_names = []
for liste in b:
    reg_names.append((liste['name']))

#Comput impacts
if hazard_name == 'tc':
    #impact function
    impf = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet()
elif hazard_name == 'storm_europe':
    impf = ImpactFuncSet()
    impf.append(ImpfStormEurope.from_welker())

#calibration with synthetic impacts 
impact_syn = Impact()
impact_syn.calc(hospitals, impf, hazard.select(orig=False), save_mat=True)
#validation with historic tracks
impact_his = Impact()
impact_his.calc(hospitals, impf, hazard.select(orig=True), save_mat=True)

#Parametric options
payout_options = fct.generate_payout_options(steps, categories)

payouts_syn, list_radius, list_payout_structure = fct.compute_payout(hazard_name, payout_options, categories, 
                                                                    radius, hazard_present_folder, insured, impf_id, orig=False)

#Optimization
RMSE_all = fct.RMSE_payouts(payouts_syn, impact_syn)
idx_min_errors = np.where(RMSE_all == np.min(RMSE_all)) 
min_radius = np.asarray(list_radius)[idx_min_errors]
min_structure = np.asarray(list_payout_structure)[idx_min_errors]
one_best = idx_min_errors[0][0]

#Payouts for historic tracks
payouts_his, _, _ = fct.compute_payout(hazard_name, payout_options, categories, 
                                      radius, hazard_present_folder, insured, impf_id, orig=True)

"""Impacts and payouts for climate change scenario: Berechnungen Fig.5"""
hazard_future_folder = project_path + country_iso + '/Hazard/'+ hazard_name+'/future'

# Read or write (muss noch ausformuliert werden), auch, dass es für beide scenarien separat gemacht wird
# tc_knu = hazard.select(orig=False).apply_climate_scenario_knu(ref_year=Year, rcp_scenario=RCP)
# Path(hazard_future_folder_sc).mkdir(parents=True, exist_ok=True)
# tc_knu.write_hdf5(os.path.join(hazard_future_folder_sc, hazard_file))


# """Generate cats in circles"""
# for rad in radius:
#     circle_filename = 'TC_circle_'+str(int(rad/1000))+'km.h5'
#     circle_max_filename = 'TC_circle_max_'+str(int(rad/1000))+'km.h5'
    
#     if not os.path.exists((os.path.join(hazard_future_folder_sc, circle_filename))) or not os.path.exists(
#             (os.path.join(hazard_future_folder_sc, circle_max_filename))):
#         tc_knu.centroids.set_geometry_points()
#         list_idx = pf.centroids_in_a_circle(tc_knu.centroids.geometry, hospitals.gdf, rad)
#         tc_circle = pf.cat_in_a_circle(tc_knu, list_idx)
#         tc_circle.write_hdf5(os.path.join(hazard_future_folder_sc, circle_filename)) 
#         tc_circle_max = pf.max_in_circle(tc_circle, hospitals.gdf, list_idx)
#         tc_circle_max.write_hdf5(os.path.join(hazard_future_folder_sc, circle_max_filename))

RCP = 26
Year = 2100
hazard_future_folder_sc = hazard_future_folder+'_'+str(2100)+'_'+str(RCP)
hazard_file = hazard_name + '_'+country_iso+'.h5'
tc_knu = Hazard.from_hdf5(os.path.join(hazard_future_folder_sc, hazard_file))
impact_knu15 = Impact()
impact_knu15.calc(hospitals, impf, tc_knu, save_mat=True)
payouts_future15, _, _ = fct.compute_payout(hazard_name, payout_options, categories, radius, hazard_future_folder_sc, insured, impf_id)
RMSE_15 = fct.RMSE_payouts(payouts_future15, impact_knu15)
idx_min_errors15 = np.where( RMSE_15== np.min(RMSE_15)) 
min_structure15 = np.asarray(list_payout_structure)[idx_min_errors15]
option15 = min_structure15[-1]

RCP = 45
Year = 2100
hazard_future_folder_sc = hazard_future_folder+'_'+str(2100)+'_'+str(RCP)
hazard_file = hazard_name + '_'+country_iso+'.h5'
tc_knu = Hazard.from_hdf5(os.path.join(hazard_future_folder_sc, hazard_file))
impact_knu25 = Impact()
impact_knu25.calc(hospitals, impf, tc_knu, save_mat=True)
payouts_future25, _, _ = fct.compute_payout(hazard_name, payout_options, categories, radius, hazard_future_folder_sc, insured, impf_id)
RMSE_25 = fct.RMSE_payouts(payouts_future25, impact_knu25)
idx_min_errors25 = np.where(RMSE_25 == np.min(RMSE_25)) 
min_structure25 = np.asarray(list_payout_structure)[idx_min_errors25]
option25 = min_structure25[-1]

# option = min_structure[-1]
# impf_payout = pf.payout_from_TCcategory(option, impf_id)
# payout_knu = Impact()
# payout_knu.calc(hospitals, impf_payout, tc_knu, save_mat = True)

# payouts_future25, _, _ = pf.compute_payout(hazard_name, payout_options, categories, radius, hazard_future_folder_sc, insured, impf_id)


"""Calculations for Fig. 6"""
#per event
big_events = np.argsort(impact_his.at_event)[-3:]
event_ids = np.asarray(impact_his.event_name)[big_events]
idx_event = impact_his.event_name.index(event_ids[-2])

nr_events = 1
list_idx = []
imp_per_region = np.zeros((nr_events,len(reg_names)))
for idx, region in enumerate(admin1_geo['MOZ']):
    clipped = geopandas.clip(hospitals.gdf.geometry, region)
    list_idx.append(clipped.index)
    impact_region = np.squeeze(np.asarray((np.sum(impact_his.imp_mat[idx_event, list_idx[idx]], axis=1))).flatten())
    imp_per_region[:, idx] = impact_region
    impact_region[ impact_region==0 ] = np.nan

dense_his_payout = np.squeeze(np.asarray(payouts_syn[48].imp_mat.todense()).flatten())
dense_his_impact = np.squeeze(np.asarray(impact_syn.imp_mat.todense()).flatten())

nonzero_impacts = np.nonzero(dense_his_impact)

steps_impact = 0.1 #10% steps
bins_impact = np.arange(0,1+steps_impact, steps_impact)
bins_payout = np.arange(0, 1.25, 0.25)
nr, bins = np.histogram(dense_his_impact[nonzero_impacts],bins=bins_impact)
centroid_bin_his = np.digitize(dense_his_impact[nonzero_impacts], bins=bins_impact)

imp_bins_with_values = np.unique(centroid_bin_his)
mean_impact = np.zeros(len(imp_bins_with_values))

for idx, binn in enumerate(imp_bins_with_values):
    mean_impact[idx] = np.mean(dense_his_impact[nonzero_impacts][np.where(centroid_bin_his==binn)])


histogram = np.zeros((len(bins), len(bins_payout)-1))
for idx_bin in range(1,len(bins)+1):
    centroid_payout = np.where(centroid_bin_his == idx_bin)
    nr_payout, _ = np.histogram(dense_his_payout[nonzero_impacts][centroid_payout], bins=bins_payout)
    histogram[idx_bin-1, :] = nr_payout
    
plt.hist(histogram[0, ], bins_payout[:-1])

histogram_sum = np.sum(histogram, axis=1)
test = histogram.T/histogram_sum
mean_per_imp_bin = np.round(np.sum(test.T*bins_payout[:-1],axis=1)*100).astype(int)


"""Fig. 1: schematic plot on payout functions"""
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
fig1.set_figwidth(fig_width*cm)

ax1 = fct.plot_schematic_figure(ax1)
plt.savefig(figures_path+'Fig1.pdf', format='pdf', bbox_inches='tight')    

"""Fig 2. PLOT HOSPITALS AND IDAI TRACKS"""
event = '2019063S18038' # Idai
west, south, east, north = 30, -26, 41, -10
projection= ccrs.PlateCarree()

#Historic track
tracks = TCTracks.from_ibtracs_netcdf(basin=basin, storm_id=event)
#synthetic tracks
tracks.equal_timestep()
tracks.calc_perturbed_trajectories(nb_synth_tracks=3)

#PLOTTING
fig1, ax1 = plt.subplots(subplot_kw={'projection':projection},figsize=(max_fig_height*cm,fig_width*cm))
ax1.set_extent((west, east+2.4, south-1, north))
fig1.set_figwidth(fig_width*cm)

ax1.add_geometries(admin1_geo['MOZ'], crs=ccrs.epsg('3857'))
handles, labels = ax1.get_legend_handles_labels()
ax1 = tracks.plot(axis=ax1, legend=False)
ax1.scatter(hospitals_plot['longitude'], hospitals_plot['latitude'], color='k',  s=15, zorder=3)
#create legend
CAT_NAMES = {
    -1: 'Tropical depression',
    0: 'Tropical storm',
    1: 'Tropical cyclone Cat. 1',
    2: 'Tropical cyclone Cat. 2',
    3: 'Tropical cyclone Cat. 3',
    4: 'Tropical cyclone Cat. 4',
    5: 'Tropical cyclone Cat. 5',
}
SAFFIR_SIM_CAT = [34, 64, 83, 96, 113, 137, 1000]
CAT_COLORS = cm_mp.rainbow(np.linspace(0, 1, len(SAFFIR_SIM_CAT)))
leg_lines = [Line2D([0], [0], color=CAT_COLORS[i_col], lw=2)
             for i_col in range(len(SAFFIR_SIM_CAT))]
leg_lines.append(Line2D([0], [0], color='grey', lw=2, ls='solid'))
leg_lines.append(Line2D([0], [0], color='grey', lw=2, ls=':'))
new_leg_lines = [Line2D([], [], marker='o', markersize=5, color='k', markeredgecolor='k', linewidth=0)]+leg_lines
leg_names = [CAT_NAMES[i_col] for i_col in sorted(CAT_NAMES.keys())]
leg_names.append('IBTrACS')
leg_names.append('IBTrACS_p')
new_leg_names = ['Hospital']+leg_names
ax1.legend(new_leg_lines, new_leg_names, loc="lower right")
ax1.set_extent((west, east+2.4, south-1, north))
ax1.set_yticks(ax1.get_yticks()[1:])
ax1.set_xticks(ax1.get_xticks())
ax1.add_geometries(admin1_geo['MOZ'], crs=projection, edgecolor='k', alpha=0.2, facecolor='grey')
reg_names1 = ['Maputo', 'Maputo City',
 'Cabo Delgado',
 'Niassa',
 'Tete',
 'Manica',
 'Gaza',
 'Zambezia',
 'Inhambane',
 'Sofala',
 'Nampula']
reg_location_x = [32, 32.8, 38.1, 36, 32, 32.6, 32, 37.5, 33.4, 34.3, 38]
reg_location_y = [-25.4, -25.9, -12.7, -13, -15, -21, -23, -16.2, -22.6, -19, -14.7]
for idx, region in enumerate(admin1_geo['MOZ']):
    ax1.annotate(reg_names1[idx], (reg_location_x[idx], reg_location_y[idx]))
xlabel= ax1.get_xticks()
new_xlabel = []
for label in xlabel:
   new_xlabel.append(str(int(label))+'°E') 
ax1.set_xticklabels(new_xlabel)
ylabel= ax1.get_yticks()
new_ylabel = []
for label in ylabel:
   new_ylabel.append(str(int(-label))+'°S') 
ax1.set_yticklabels(new_ylabel)
plt.savefig(figures_path+'Fig2.pdf', format='pdf', bbox_inches='tight')


"""FIG 3. Idai over Central Mozambique, with hospital locations and radius of 30km"""
rad= 30*1000 #plot this radius
event = '2019063S18038' # Idai

cmap0 = LinearSegmentedColormap.from_list('', ['white', *plt.cm.GnBu(np.arange(255))])
projection= ccrs.PlateCarree()

fig3, ax3 = plt.subplots(subplot_kw={'projection':projection})
fig3.set_figwidth(fig_width*cm)
ax3.set_extent((30, 41, -22, -16))
# Idai track
idai = TCTracks.from_ibtracs_netcdf(basin=basin, storm_id=event)
idai.equal_timestep(time_step_h=0.08) #high temporal resolution for the tropical cyclone wind field
idai_wind = TropCyclone.from_tracks(idai, centroids=hazard.centroids)
ax3 = idai_wind.plot_intensity(event, cmap=cmap0, axis=ax3, smooth=False, vmin=0.0, vmax=60)

ax3.set_extent((30, 41, -22, -16))
ax3.set_title('')
ax3.set_yticks(ax3.get_yticks()[1:])
ax3.set_xticks(ax3.get_xticks())
ax3.scatter(hospitals_plot['longitude'], hospitals_plot['latitude'], color='k', s=25, label='Hospital')
ax3.set_extent((30, 41, -22, -16))
crs_orig = hospitals.gdf.geometry.crs
exp_gdf_metric = hospitals.gdf.geometry.to_crs(epsg=6933)
circles_metric = exp_gdf_metric.buffer(rad, cap_style=1)
circles  = circles_metric.to_crs(crs_orig)
circles.index = np.arange(0, circles.index.shape[0])

ax3.add_geometries(circles, crs=projection, edgecolor='k', facecolor='none', label='Circle (30 km radius)')

ax3.add_geometries(admin1_geo['MOZ'], crs=projection, edgecolor='k', alpha=0.2, facecolor='grey')

#annotation of regions
reg_names1 = ['Tete',
 'Manica',
 'Gaza',
 'Zambezia',
 'Inhambane',
 'Sofala',
 'Nampula']
reg_location_x = [32.4, 32.7, 32.2, 36, 33.5, 34.3, 38.5]
reg_location_y = [-16.3, -21, -21.8, -17, -21.8, -19, -16.3]
for idx, region in enumerate(reg_names1):
    # ax.annotate(reg_names[idx], (region.centroid.x, region.centroid.y))
    ax3.annotate(reg_names1[idx], (reg_location_x[idx], reg_location_y[idx]))

xlabel= ax3.get_xticks()
new_xlabel = []
for label in xlabel:
   new_xlabel.append(str(int(label))+'°E') 
ax3.set_xticklabels(new_xlabel)

ylabel= ax3.get_yticks()
new_ylabel = []
for label in ylabel:
   new_ylabel.append(str(int(-label))+'°S') 
ax3.set_yticklabels(new_ylabel)

ax3.legend(loc='lower right')
proxy_artist = mpatches.CirclePolygon((0, 0), radius=2, edgecolor='black', facecolor='none')
line4 = Line2D([], [], marker='o', markersize=15, markerfacecolor="none", markeredgecolor='k', linewidth=0, markeredgewidth=2)
handles, labels = ax3.get_legend_handles_labels()
handles.append(line4)
labels.append('Circle (30 km radius)')
ax3.legend(handles=handles, labels=labels, loc='lower right')

plt.savefig(figures_path+'Fig3.pdf', format='pdf', bbox_inches='tight')


"""Fig4: Plotting parametric products and basis risk"""
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
fig4.set_figwidth(fig_width*cm)
categories_plot = np.array([18, 33, 43, 50, 59, 70, 120]) #TC
ax4 = fct.plot_vul_his(hazard, categories_plot, impf, impf_id, min_structure[1])
plt.savefig(figures_path+'Fig4.pdf', format='pdf', bbox_inches='tight')


"""Fig.5: Parametric product both for present and future climate"""
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
fig5.set_figwidth(fig_width*cm)
ax5 = fct.plot_products(payouts_syn, impact_syn, payouts_his, impact_his, idx_min_errors)
ax5.scatter(impact_knu15.at_event, payouts_future15[48].at_event, c='#7570b3', alpha=0.4, label='IBTrACS_p (1.5°C)')
ax5.scatter(impact_knu25.at_event, payouts_future25[48].at_event, c='#d95f02', alpha=0.4,  label='IBTrACS_p (2.5°C)')
ax5 = fct.add_annotation_MOZ(impact_his, payouts_his, one_best, ax5)
plt.legend()
ax5.figure
plt.savefig(figures_path+'Fig5.pdf', format='pdf', bbox_inches='tight')


"""Figure 6"""
impact_colour = 'k'
payout_colour = 'k'
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
fig6.set_figwidth(fig_width*cm)
for idx in range(4):
    ax6.bar(bins_impact[1:-1], test[idx,0:len(bins_impact)-2], 0.08, label=str(idx*25)+'% Payout', 
            bottom=test[idx-1,0:len(bins_impact)-2], alpha=0.15*idx+0.15, color=payout_colour)

plt.xlim(0.05,0.85)
plt.ylim(0,1)
locs, xlabels = plt.xticks()

new_locs = bins_impact[1:-1]
newlabel = []
for idx, i in enumerate(bins_impact[:-2]):
    newlabel.append(str((i*100).astype(int))+'-'+str(((i+steps_impact)*100).astype(int)))
     
ax6.set_xticks(new_locs[:-1], newlabel[:-1])

ax6.set_xlabel('Modeled damage per affected unit and event (%)')
ax6.set_ylabel('Fraction of cases (%)')
locs, ylabels = plt.yticks()
ax6.set_yticks(locs, (locs*100).astype(int))

ax62 = ax6.twiny()
ax62.set_xlim(ax6.get_xlim())
ax62.set_xticks(new_locs[:-1],nr[:-2].tolist(),c=impact_colour)
ax62.set_xlabel('Number of cases', c=impact_colour)

ax6.scatter(bins_impact[1:-2], mean_impact, label='Mean damage', marker='D', c=impact_colour)
ax6.scatter(bins_impact[1:-2], mean_per_imp_bin[:-3]/100, c=payout_colour, label='Mean payout')

ax6.legend(loc='lower right')
plt.savefig(figures_path+'Fig6.pdf', format='pdf', bbox_inches='tight')



"""Case study 2: winter storms in France"""
country_iso_FRA = 'FRA'
country_name_FRA  = 'france'
basin_FRA  = None
impf_id_FRA  = 1
hazard_name_FRA = 'storm_europe'
bounds_FRA  = None
steps_FRA  = 10 #payout
categories_FRA  = np.arange(40,70,10) #intensity
radius_FRA = np.array([10, 20, 30, 40, 50])*1000 #meter

# Generate or load winter storms in France
winter_storms, ws_present_folder, hospitals_FRA, max_payout_FRA = fct.wrapper_haz_exp(hazard_name_FRA, country_iso_FRA, country_name_FRA, impf_id_FRA, 
                                             project_path, basin_FRA, bounds_FRA)
fct.wrapper_generate_haz_in_a_circle(winter_storms, hospitals_FRA, hazard_name_FRA, ws_present_folder, radius_FRA)

#Compute impacts
impf = ImpactFuncSet()
impf.append(ImpfStormEurope.from_welker())

#Synthetic impacts 
damages_syn_FRA = Impact()
damages_syn_FRA.calc(hospitals_FRA, impf, winter_storms.select(orig=False), save_mat=True)
#Historic damages
damages_his_FRA = Impact()
damages_his_FRA.calc(hospitals_FRA, impf, winter_storms.select(orig=True), save_mat=True)

#Parametric options
payf_FRA = fct.generate_payout_options(steps_FRA, categories_FRA)

payouts_syn, list_radius, list_payout_structure = fct.compute_payout(hazard_name_FRA, payf_FRA, categories_FRA, 
                                                                    radius, ws_present_folder, max_payout_FRA, impf_id_FRA, orig=False)

#Optimization
RMSE_all = fct.RMSE_payouts(payouts_syn, impact_syn)
idx_min_errors = np.where(RMSE_all == np.min(RMSE_all)) 
min_radius = np.asarray(list_radius)[idx_min_errors]
min_structure = np.asarray(list_payout_structure)[idx_min_errors]
one_best = idx_min_errors[0][0]

#Payouts for historic storms
payouts_his, _, _ = fct.compute_payout(hazard_name, payout_options, categories, 
                                      radius, hazard_present_folder, insured, impf_id, orig=True)

"""Fig7: Plotting parametric products and basis risk"""
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
fig7.set_figwidth(fig_width*cm)
categories_plot = np.arange(35,70,5)
ax7 = fct.plot_vul_his_FRA(hazard, categories_plot, impf, impf_id, min_structure[-1])
plt.savefig(figures_path+'Fig7.pdf', format='pdf', bbox_inches='tight')




