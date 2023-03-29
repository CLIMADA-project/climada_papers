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

Define functions to design parametric insurance products and create plots
"""

# from geopandas import GeoDataFrame
# from shapely.geometry import Point
# import pandas as pd
# import random
import geopandas
import copy
import numpy as np
import matplotlib.pyplot as plt
import contextily as ctx
from pathlib import Path
from scipy import sparse
import os
from sklearn.metrics import mean_squared_error
import math

from itertools import combinations_with_replacement
import cartopy.crs as ccrs

from climada.hazard.base import Hazard
from climada.hazard import TCTracks
from climada.hazard import Centroids, TropCyclone
from climada.entity import ImpactFunc, ImpactFuncSet
from climada.engine import Impact
from climada.entity import Exposures
from climada_petals.entity.exposures.openstreetmap.osm_dataloader import OSMRaw, OSMFileQuery
from climada import CONFIG
from climada.util.api_client import Client
client = Client()
from climada.entity.impact_funcs.storm_europe import ImpfStormEurope


def wrapper_haz_exp(hazard_name, country_iso, country_name, impf_id, project_path='data/', basin=None, bounds=None):
    
    Path(project_path).mkdir(parents=True, exist_ok=True)
    
    hazard_present_folder = project_path + country_iso + '/Hazard/'+ hazard_name+'/present_day'
    #hazard_future_folder = project_path + country_iso + '/Hazard/'+ hazard_name+'/future'
    hazard_file = hazard_name + '_'+country_iso+'.h5'
    exposure_folder = project_path + country_iso + '/Exposure'
    hospitals_file = 'hospitals_'+country_iso+'.h5'
    insured_file = 'insured_'+country_iso+'.h5'

    """Read or generate hazard and or exposure"""
    #read hazard file if it already exists, generate otherwise
    if os.path.exists((os.path.join(hazard_present_folder, hazard_file))):
        #read hazard
        hazard = Hazard.from_hdf5(os.path.join(hazard_present_folder, hazard_file))
            
    else:
        #generate output directories if they do not exist yet
        Path(hazard_present_folder).mkdir(parents=True, exist_ok=True)
        if hazard_name == 'tc':
            #generate TCs
            hazard = generate_TCs(basin, bounds)
        elif hazard_name == 'storm_europe':
            hazard = client.get_hazard('storm_europe', properties={'country_iso3alpha': country_iso})
        #save hazard
        hazard.write_hdf5(os.path.join(hazard_present_folder, hazard_file))

    #read exposure file if it already exists, generate otherwise
    if os.path.exists((os.path.join(exposure_folder, hospitals_file))):
        #read exposure
        hospitals = Exposures.from_hdf5(os.path.join(exposure_folder, hospitals_file))
        insured = Exposures.from_hdf5(os.path.join(exposure_folder, insured_file))
    else:
        # generate output directories if they do not exist yet
        Path(exposure_folder).mkdir(parents=True, exist_ok=True)
        #generate exposure
        hospitals = generate_osm_hospital(country_iso, country_name, hazard.tag.haz_type, impf_id)
        
        #hospitals = pf.preprocess_hospitals(hospitals)
        #exclude exposure points that are less than radius away from each other 
        #(consider them as one exposure point)
        # distance = 300 #meter
        # distance_idx_list = pf.centroids_in_a_circle(hospitals.gdf.geometry, hospitals.gdf, distance)
        # for idx in distance_idx_list:
        #     if len(idx) >= 2:
        #         hospitals.gdf = hospitals.gdf.drop(idx[1:], errors='ignore')
        
        hospitals.write_hdf5(os.path.join(exposure_folder, hospitals_file))
        
        if hazard_name == 'storm_europe':
            impf_storm = ImpfStormEurope.from_welker()
            impf_storm.plot()
            max_impact = impf_storm.calc_mdr(120) #maximal impact
        elif hazard_name == 'tc':
            max_impact = 1
        
        insured = copy.deepcopy(hospitals)
        insured.gdf.value = max_impact
        insured.write_hdf5(os.path.join(exposure_folder, insured_file))
        
        
    hospitals.assign_centroids(hazard)
    hospitals.gdf.index = np.arange(0, hospitals.gdf.index.shape[0])
    insured.assign_centroids(hazard)
    insured.gdf.index = np.arange(0, hospitals.gdf.index.shape[0])
    
    return hazard, hazard_present_folder, hospitals, insured


def wrapper_generate_haz_in_a_circle(hazard, hospitals, hazard_name, hazard_path, radius):
    """Generate cats in circles"""
    for rad in radius:
        circle_filename = hazard_name+'_circle_'+str(int(rad/1000))+'km.h5'
        circle_max_filename = hazard_name+'_circle_max_'+str(int(rad/1000))+'km.h5'
        if not os.path.exists((os.path.join(hazard_path, circle_filename))) or not os.path.exists(
                (os.path.join(hazard_path, circle_max_filename))):
            hazard.centroids.set_geometry_points()
            list_idx = centroids_in_a_circle(hazard.centroids.geometry, hospitals.gdf, rad)
            hazard_circle = cat_in_a_circle(hazard, list_idx)
            hazard_circle.write_hdf5(os.path.join(hazard_path, circle_filename))
            hazard_circle_max = max_in_circle(hazard_circle, hospitals.gdf, list_idx)
            hazard_circle_max.write_hdf5(os.path.join(hazard_path, circle_max_filename))


def generate_TCs(basin, bounds):
    """
    Read historic TC tracks, generate perturbed tracks and compute wind fields

    Parameters
    ----------
    basin : str
        The basin where a TC is formed is not defined in IBTrACS. However, this filter option
        allows to restrict to storms whose first valid eye position is in the specified basin,
        which simulates the genesis location. Note that the resulting genesis basin of a
        particular track may depend on the selected `provider` and on `estimate_missing`
        because only the first *valid* eye position is considered. Possible values are 'NA'
        (North Atlantic), 'SA' (South Atlantic), 'EP' (Eastern North Pacific, which includes
        the Central Pacific region), 'WP' (Western North Pacific), 'SP' (South Pacific),
        'SI' (South Indian), 'NI' (North Indian). If None, this filter is not applied.
        Default: None.
    bounds : tuple (lon_min, lat_min, lon_max, lat_max)

    Returns
    -------
    tc : TropCyclone

    """
    
    #Historic tracks
    tracks = TCTracks.from_ibtracs_netcdf(basin=basin, year_range=(1980,2022))
    #synthetic tracks
    tracks.equal_timestep()
    tracks.calc_perturbed_trajectories()
    #centroids 
    centrs = Centroids.from_pnt_bounds(bounds, res= 0.05)
    centrs.check()
    # Using the tracks, compute the windspeed at the location of the centroids
    tc = TropCyclone.from_tracks(tracks, centroids=centrs)
    tc.check()

    return tc

def generate_osm_hospital(country_iso, country_name, haz_type, impf_id):
    #hospitals
    DATA_DIR = CONFIG.exposures.openstreetmap.local_data.dir()
    OSMRaw().get_data_geofabrik(country_iso, file_format='pbf', save_path=DATA_DIR)
    HNDFileQuery = OSMFileQuery(Path(DATA_DIR,country_name+'-latest.osm.pbf'))
    gdf_health = HNDFileQuery.retrieve_cis('healthcare')
    #gdf_health = HNDFileQuery.retrieve_cis('healthcare')
    gdf_hospital = gdf_health[gdf_health.amenity=='hospital']
    hospitals = Exposures(gdf_hospital)
    #exp_nl_poly.gdf['impf_WS'] = 1
    hospitals.gdf.head()
    hospitals.gdf['value'] = 1
    hospitals.gdf['geometry'] = hospitals.gdf['geometry'].centroid
    hospitals.gdf['longitude'] = hospitals.gdf['geometry'].x
    hospitals.gdf['latitude'] = hospitals.gdf['geometry'].y
    hospitals.gdf['impf_'+haz_type]= impf_id #e.g. 5 for 508=MOZ
    
    return hospitals

def preprocess_hospitals(hospitals):
    hospitals.gdf.index = np.arange(0, hospitals.gdf.index.shape[0])


    names = hospitals.gdf['name']
    for idx_name, name in enumerate(names):
        print(name)
        if name is None:
            hospitals.gdf = hospitals.gdf.drop(idx_name)
        elif 'ospital' not in name:
            hospitals.gdf = hospitals.gdf.drop(idx_name)

    hospitals.gdf = hospitals.gdf.drop_duplicates(subset='name')
    
    return hospitals

def get_different_hospital_types(hospitals):
    names = hospitals.gdf['name']
    rural = []
    central = []
    military = []
    district = []
    general = []
    province = []
    other = []

    for name in names:
        if 'ural' in name:
            rural.append(name)
        elif 'ilitar' in name:
            military.append(name)
        elif 'istri' in name:
            district.append(name)
        elif 'alma' in name:
            district.append(name)
        elif 'acomia' in name:
            district.append(name)
        elif 'Namacurra' in name:
            district.append(name)
        elif 'eral' in name:
            general.append(name)
        elif 'rovinc' in name:
            province.append(name)
        elif 'entral' in name:
            central.append(name)
        else:
            other.append(name)
    
    return rural, central, military, district, general, province, other
    
    

def centroids_in_a_circle(geometry, exp_gdf, radius):
    """Cut out hazard in circles with specific radius"""
    crs_orig = exp_gdf.geometry.crs
    exp_gdf_metric = exp_gdf.geometry.to_crs(epsg=6933)
    circles_metric = exp_gdf_metric.buffer(radius, cap_style=1)
    circles  = circles_metric.to_crs(crs_orig)
    circles.index = np.arange(0, circles.index.shape[0])
    
    list_idx = []
    for circle in circles:
        clipped = geopandas.clip(geometry, circle)
        list_idx.append(clipped.index)
        
    return list_idx
 
   
def cat_in_a_circle(tc, list_idx):
    all_idx_in_one = []
    for indices in list_idx:
        all_idx_in_one.extend(indices.tolist())
        
    unique_idx= np.unique(np.asarray(all_idx_in_one))
    tc_circle = copy.deepcopy(tc)
    new_intensity = sparse.lil_matrix((tc.intensity.shape))
    new_intensity[:,unique_idx] = tc.intensity[:,unique_idx]
    tc_circle.intensity = new_intensity.tocsr()
    
    return tc_circle
    

def max_in_circle(tc_circle, exp_gdf, list_idx):    
    events = tc_circle.intensity.shape[0]
    max_per_hospital = np.zeros([events, exp_gdf.index.shape[0]])
    for idx_hospital, indices in enumerate(list_idx):
        max_per_hospital[:, idx_hospital] = np.asarray(((np.max(tc_circle.intensity[:, indices], axis=1)).todense())).flatten()        
    
    
    centr_hospitals = exp_gdf['centr_'+tc_circle.tag.haz_type].values
    
    tc_circle_max = copy.deepcopy(tc_circle)
    new_intensity = sparse.lil_matrix((tc_circle.intensity.shape))
    new_intensity[:,centr_hospitals] = max_per_hospital
    tc_circle_max.intensity = new_intensity.tocsr()
    
    return tc_circle_max


def payout_structure(steps, percentage, impf_id, haz_type):    
    impf = ImpactFunc()
    impf.id = impf_id
    impf.haz_type = haz_type
    impf.tag = ''
    impf.intensity = steps
    impf.paa = np.ones(len(steps))
    impf.mdd = percentage
    impf.check()
    impf_set = ImpactFuncSet()
    impf_set.append(impf)
    impf_set.check()
    
    return impf_set

def payout_from_TCcategory(percentages, impf_id, haz_type='TC'):
    #steps based on TC categories
    steps1 = np.array([0, 32.999, 33, 42.999, 43, 49.999, 50, 58.999, 59, 69.999, 70, 100])
    percentages1 = np.array([0, 0, 0, 0, percentages[0], percentages[0], percentages[1], percentages[1], 
                             percentages[2], percentages[2], percentages[3], percentages[3]])*0.01

    impf_payout = payout_structure(steps1, percentages1, impf_id, haz_type)
    
    return impf_payout

def create_impf_step_fct(categories, percentages, impf_id, haz_type):
    #steps based on TC categories
    steps = np.array([0])
    for category in categories:
        steps = np.append(steps, category)
        steps = np.append(steps, category)
    
    payout_fraction = np.array([0, 0])
    for percentage in percentages:
        payout_fraction = np.append(payout_fraction, percentage)
        payout_fraction = np.append(payout_fraction, percentage)
    
    payout_fraction = payout_fraction[:-1]*0.01

    impf_payout = payout_structure(steps, payout_fraction, impf_id, haz_type)
    
    return impf_payout

def generate_payout_options(steps, categories):
    nr_categories = len(categories)
    percentages = np.arange(0,100+steps,steps)
    payout_options = list(combinations_with_replacement(percentages, nr_categories))
    
    return payout_options
    

def compute_payout(haz_name, payout_options, categories, radius, output_folder, hospitals, impf_id, orig=False):
    list_payouts = []
    list_radius = []
    list_payout_structure = []
    
    for idx_rad, rad in enumerate(radius):
        circle_max_filename = haz_name+'_circle_max_'+str(int(rad/1000))+'km.h5'
        tc_circle_max = Hazard.from_hdf5(os.path.join(output_folder, circle_max_filename))
        
        for idx_opt, option in enumerate(payout_options):
            if haz_name == 'tc':
                impf_payout = payout_from_TCcategory(option, impf_id)
            elif haz_name == 'storm_europe':
                impf_payout = create_impf_step_fct(categories, option, impf_id, tc_circle_max.tag.haz_type)
            payout_hospitals = Impact()
            payout_hospitals.calc(hospitals, impf_payout, tc_circle_max.select(orig=orig), save_mat = True)
            list_payouts.append(payout_hospitals)
            list_radius.append(rad)
            list_payout_structure.append(option)
            
    return list_payouts, list_radius, list_payout_structure

def RMSE_payouts(list_payouts, impact_syn, coord_exp=None):
    MSE = np.zeros(len(list_payouts))
    RMSE = np.zeros(len(list_payouts))
    for idx_product, product in enumerate(list_payouts):
        MSE[idx_product] = mean_squared_error(impact_syn.select(coord_exp=coord_exp).at_event, product.select(coord_exp=coord_exp).at_event)
        RMSE[idx_product] = math.sqrt(MSE[idx_product])
    
    return RMSE



def regional_imp_pay(list_payouts, impact_syn, coord_exp=None):
    regional_imp = impact_syn.select(coord_exp=coord_exp)
    regional_pay = []
    for idx_product, product in enumerate(list_payouts):
        regional_pay.append(product.select(coord_exp=coord_exp))
        
    return regional_imp, regional_pay

def plot_schematic_figure(ax1):
    categories_dummy = np.arange(20,90,10)
    payout_categories_dummy = [25, 50, 75, 100]
    TS, cat1, cat2, cat3, cat4, cat5, upper_lim = categories_dummy
    steps1 = np.array([0, cat1-0.001, cat1, cat2-0.001, cat2, cat3-0.001, cat3, cat4-0.001, cat4, cat5-0.001, cat5, upper_lim])
    step_cat2, step_cat3, step_cat4, step_cat5 = payout_categories_dummy #step_TS, step_cat1, 
    percentages = np.array([0, 0, 0, 0, step_cat2, step_cat2, step_cat3, step_cat3, 
                                 step_cat4, step_cat4, step_cat5, step_cat5])
    linear_x = np.arange(30, 90, 10)    
    linear_y = np.array([0,0,  33.33, 66.66 , 100, 100])
    ax1.plot(steps1, percentages, label='Step payout \nfunction', c='k', ls='--', zorder=3)
    ax1.plot(linear_x, linear_y, label='Linear payout \nfunction', c='k', ls='-', zorder=2)
    plt.xlim(30, 80)
    plt.ylim(0, 110)
    locs, xlabels = plt.xticks()
    newlabel = ['thr 1', 'thr 2', 'thr 3', ' thr 4']
    plt.xticks(categories_dummy[2:-1], newlabel,horizontalalignment='left')
    ax1.set_xlabel('Parametric index')
    locs, ylabels = plt.yticks()
    newylabel = ['0%', 'Payout thr 1', 'Payout thr 2', 'Payout thr 3', 'Maximum payout = 100%']
    # ax.set_ylabel('Payout')
    plt.yticks([0]+payout_categories_dummy, newylabel)
    ax1.legend(loc='lower right')
    
    return ax1


def plot_hospitals_map(hospitals):
    # plot results
    figures_path = '/Users/carmenst/Library/CloudStorage/Dropbox/Aplicaciones/Overleaf/WP1_Environment_Systems_and_Decisions/sn-article-template/art/'
    ax = hospitals.gdf.set_crs(epsg=4326).to_crs(epsg=3857).plot(figsize=(15, 15), alpha=1, markersize=5, color='blue', 
                        edgecolor='blue', label='Hospitals MZB')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles, loc='upper left')
    #ax.set_title('Hospitals Mozambique', fontsize=25)
    ctx.add_basemap(ax)
    plt.savefig(figures_path+'hospital_locations_MOZ'+'.pdf', format='pdf', bbox_inches='tight')
    
    return ax

def plot_event_in_circle(rad, event, hazard_folder):    
    cmap = plt.cm.get_cmap("GnBu").copy()
    cmap.set_under(color='white') 
    circle_file = 'TC_circle_'+str(int(rad/1000))+'km.h5'
    tc_circle = TropCyclone.from_hdf5(os.path.join(hazard_folder, circle_file))
    projection= ccrs.PlateCarree()
    fig, axes = plt.subplots(subplot_kw={'projection':projection},figsize=(12,8))
    #ax = tc.plot_intensity(event=event, cmap=cmap, alpha=0.5, smooth=False, vmin=0.00001, vmax=80)
    tc_circle.plot_intensity(event=event, cmap=cmap, smooth=False, vmin=0.00001, vmax=80, axis=axes) #axis=ax
    
    return fig, axes

def plot_vul_his(tc, categories, impf_TC, impf_id, min_structure):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1 = plot_TC_hist(tc, categories, ax1)
    ax2 = ax1.twinx()
    ax2 = plot_vulnerability(impf_TC, impf_id, categories, min_structure, ax2)

    return ax1

def plot_vul_his_FRA(tc, categories, impf_TC, impf_id, min_structure):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1 = plot_FRA_hist(tc, categories, ax1)
    ax2 = ax1.twinx()
    ax2 = plot_vulnerability_FRA(impf_TC, impf_id, categories, min_structure, ax2)

    return ax1

def plot_vulnerability_FRA(impf_TC, fun_id, categories, payout_structures, ax):
    #plot impf and payout structure    
    ax = plot_impf(impf_TC, fun_id, ax)
    ax = plot_payout_structure_FRA(payout_structures, categories, ax)

    handles, labels = plt.gca().get_legend_handles_labels()
    # labels will be the keys of the dict, handles will be values
    temp = {k:v for k,v in zip(labels, handles)}    
    # ax.legend(temp.values(), temp.keys(), loc="upper center", bbox_to_anchor=(0.7, 1))
    ax.legend(temp.values(), temp.keys(), loc="upper right")
    ax.set_ylabel('Damage / payout (%)')
    plt.ylim(0, 100)

    return ax

def plot_vulnerability(impf_TC, fun_id, categories, payout_structures, ax):
    #plot impf and payout structure    
    ax = plot_impf(impf_TC, fun_id, ax)
    ax = plot_payout_structure(payout_structures, categories, ax)

    handles, labels = plt.gca().get_legend_handles_labels()
    # labels will be the keys of the dict, handles will be values
    temp = {k:v for k,v in zip(labels, handles)}    
    ax.legend(temp.values(), temp.keys(), loc="lower right")
    ax.set_ylabel('Damage / payout (%)')
    plt.ylim(0, 100)

    return ax

def plot_impf(impf_TC, fun_id, ax):
    mdd_impf = impf_TC.get_func(fun_id=fun_id)[0].mdd*100
    intensity_mdd = impf_TC.get_func(fun_id=fun_id)[0].intensity
    ax.plot(intensity_mdd, mdd_impf, label='Impact \nfunction', c='k', zorder=2)
    
    return ax 

def plot_payout_structure(payout_structures, categories, ax):
    
    TS, cat1, cat2, cat3, cat4, cat5, upper_lim = categories
    steps1 = np.array([0, cat1-0.001, cat1, cat2-0.001, cat2, cat3-0.001, cat3, cat4-0.001, cat4, cat5-0.001, cat5, upper_lim])

    # for idx_structure, structure in enumerate(payout_structures):
    #     step_cat2, step_cat3, step_cat4, step_cat5 = structure #step_TS, step_cat1, 
    #     percentages = np.array([0, 0, 0, 0, step_cat2, step_cat2, step_cat3, step_cat3, 
    #                              step_cat4, step_cat4, step_cat5, step_cat5])
        
    #     ax.plot(steps1, percentages, label='Payout \nstructures', c='k', ls='--', zorder=2)
    
    step_cat2, step_cat3, step_cat4, step_cat5 = payout_structures #step_TS, step_cat1, 
    percentages = np.array([0, 0, 0, 0, step_cat2, step_cat2, step_cat3, step_cat3, 
                                 step_cat4, step_cat4, step_cat5, step_cat5])
        
    ax.plot(steps1, percentages, label='Payout \nfunction', c='k', ls='--', zorder=2)

    # plt.ylim(0, 100)      
    
    
    return ax

def plot_payout_structure_FRA(payout_structures, categories, ax):
    
    TS, cat1, cat2, cat3, cat4, cat5, upper_lim = categories
    steps1 = np.array([0, cat1-0.001, cat1, cat2-0.001, cat2, cat4-0.001, cat4, cat5-0.001, cat5, upper_lim])
    
    # for idx_structure, structure in enumerate(payout_structures):
    #     step_cat2, step_cat3, step_cat4, step_cat5 = structure #step_TS, step_cat1, 
    #     percentages = np.array([0, 0, 0, 0, step_cat2, step_cat2, step_cat3, step_cat3, 
    #                              step_cat4, step_cat4, step_cat5, step_cat5])
        
    #     ax.plot(steps1, percentages, label='Payout \nstructures', c='k', ls='--', zorder=2)
    
    step_cat2, step_cat3, step_cat4  = payout_structures #step_TS, step_cat1, 
    percentages = np.array([0, 0, 0, 0, step_cat2, step_cat2, step_cat3, step_cat3, 
                                 step_cat4, step_cat4])
        
    # ax.plot(steps1, percentages, label='Payout \nfunction', c='k', ls='--', zorder=2)
    # plt.ylim(0, 100)      
    steps2 = np.array([0, 40-0.001, 40, 50-0.001, 50, 60-0.001, 60, 100])
    percentages2 = np.array([0,0, 10, 10, 10, 10, 100, 100])
    ax.plot(steps2, percentages2, label='Payout \nfunction', c='k', ls='--', zorder=2)
    
    
    return ax

def plot_TC_hist(tc, categories, ax):
    tc.centroids.set_on_land()
    tc.centroids.region_id = tc.centroids.on_land
    #should I limit it to on land centroids?
    dense_intensity_his = np.squeeze(np.asarray(tc.select(orig=True, reg_id=True).intensity.todense()).flatten())
    dense_intensity_syn = np.squeeze(np.asarray(tc.select(orig=False, reg_id=True).intensity.todense()).flatten())
    
    ax.hist(dense_intensity_his, bins=categories[1:], density=True, alpha=0.5, color='k', edgecolor='k', lw=2, zorder=3, label='IBTrACS \n(1980-2022)')
    ax.hist(dense_intensity_syn, bins=categories[1:], density=True, alpha=0.3, color='k', lw=0, label='IBTrACS_p \n(1980-2022)')
    # ax.hist(dense_intensity_his, bins=categories[1:], density=True, alpha=0.3, color='#7570b3', edgecolor='#7570b3', lw=2, zorder=3, label='IBTrACS \n(1980-2022)')
    # ax.hist(dense_intensity_syn, bins=categories[1:], density=True, alpha=0.3, color='#d95f02', edgecolor='#d95f02', lw=2, label='IBTrACS_p \n(1980-2022)')
    # ax.hist(dense_intensity_his, bins=categories[1:], density=True, fill=False, linewidth= 0, edgecolor=None, hatch='..', label='IBTrACS \n(1980-2022)')
    # ax.hist(dense_intensity_syn, bins=categories[1:], density=True, color='k', alpha=0.2, zorder=0 , label='IBTrACS_p \n(1980-2022)')
    # ax.set_ylabel('Spatial cumulative distribution of \non land storms $\geq$ CAT 1 (%)')
    # ax.set_ylabel('Land centroids affected by TC $\geq$ Cat. 1 (%)')
    ax.set_ylabel('Tropical cyclone intensity over land (%)')
    
    
    locs, xlabels = plt.xticks()
    newlabel = ['  Trop. storm', '   Cat. 1', ' Cat. 2', '   Cat. 3', '    Cat. 4', '    Cat. 5', '']
    ax.set_xlabel('Tropical cyclone category')
    plt.xticks(categories, newlabel,horizontalalignment='left')
    plt.xlim(10, 80)
    ax.legend(loc="upper left")
    
    return ax

def plot_FRA_hist(tc, categories, ax):
    tc.centroids.set_on_land()
    tc.centroids.region_id = tc.centroids.on_land
    #should I limit it to on land centroids?
    dense_intensity_his = np.squeeze(np.asarray(tc.select(orig=True, reg_id=True).intensity.todense()).flatten())
    dense_intensity_syn = np.squeeze(np.asarray(tc.select(orig=False, reg_id=True).intensity.todense()).flatten())
    
    ax.hist(dense_intensity_his, bins=categories[1:], density=True, alpha=0.5, color='k', edgecolor='k', lw=2, zorder=3, label='WISC \n(1940-2011)')
    ax.hist(dense_intensity_syn, bins=categories[1:], density=True, alpha=0.3, color='k', lw=0, label='WISC_p \n(1940-2011)')
    # ax.hist(dense_intensity_his, bins=categories[1:], density=True, fill=False, linewidth= 0, edgecolor=None, hatch='..', label='WISC \n(1940-2011)')
    # ax.hist(dense_intensity_syn, bins=categories[1:], density=True, color='k', alpha=0.2, zorder=0 , label='WISC_p \n(1940-2011)')
    # ax.hist(dense_intensity_his, bins=categories[1:], density=True, alpha=0.3, color='#7570b3', edgecolor='#7570b3', lw=2, zorder=3, label='WISC \n(1940-2011)')
    # ax.hist(dense_intensity_syn, bins=categories[1:], density=True, alpha=0.3, color='#d95f02', edgecolor='#d95f02',lw=2, label='WISC_p \n(1940-2011)')
    # ax.set_ylabel('Spatial cumulative distribution of \non land storms $\geq$ CAT 1 (%)')
    # ax.set_ylabel('Land centroids affected by wind speed (%)')
    ax.set_ylabel('Winter storm intensity over land (%)')
    
    # locs, xlabels = plt.xticks()
    # newlabel = ['  Trop. storm', '   Cat. 1', ' Cat. 2', '   Cat. 3', '    Cat. 4', '    Cat. 5', '']
    ax.set_xlabel('Wind speed (m/s)')
    # plt.xticks(categories, newlabel,horizontalalignment='left')
    plt.xlim(20, 76)
    ax.legend(loc="upper left")
    
    return ax

def plot_products(payout_syn, impact_syn, payout_his, impact_his, idx_min):
    #colours = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
    
    fig3 = plt.figure()
    ax1 = fig3.add_subplot(111)    

    
    product_idx= 0
    product = idx_min[0][product_idx]
    
    maximum = np.round(np.max(np.concatenate((impact_his.at_event, payout_his[product].at_event, 
                               impact_syn.at_event, payout_syn[product].at_event)))*1.1)

    # ax1.scatter(impact_his.at_event, payout_his[product].at_event, 
    #             label='IBTrACS (1980-2022)', c='#1b9e77', marker='D', s=50, zorder=2, edgecolor='k')  #colours[product_idx]
    # ax1.scatter(impact_syn.at_event, payout_syn[product].at_event, 
    #             label='IBTrACS_p (1980-2022)', marker='D', facecolors='none', edgecolor='#1b9e77')
    ax1.scatter(impact_his.at_event, payout_his[product].at_event, 
                label='IBTrACS (1980-2022)',facecolors='none', marker='D', s=40, zorder=2, edgecolor='k') #c='#1b9e77',  
    ax1.scatter(impact_syn.at_event, payout_syn[product].at_event, 
                label='IBTrACS_p (1980-2022)', c='#1b9e77', alpha=0.4)
    
    basis_risk_line = np.linspace(0, maximum*1.5)
    ax1.plot(basis_risk_line, basis_risk_line, c='k', ls='--', alpha=0.2)
    ax1.set_ylabel('Calculated payout (Hospital units)')
    ax1.set_xlabel('Modeled damage (Hospital units)')
    ax1.set_ylim(0, maximum)   
    ax1.set_xlim(0, maximum)
    plt.legend()
    
    return ax1

def plot_products_FRA(payout_syn, impact_syn, payout_his, impact_his, idx_min):
    #colours = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
    
    fig3 = plt.figure()
    ax1 = fig3.add_subplot(111)    

    
    product_idx= 0
    product = idx_min[0][product_idx]
    
    maximum = np.round(np.max(np.concatenate((impact_his.at_event, payout_his[product].at_event, 
                               impact_syn.at_event, payout_syn[product].at_event)))*1.1)

    
    ax1.scatter(impact_his.at_event, payout_his[product].at_event, 
                label='WISC (1980-2022)', c='#1b9e77', marker='D', s=40, zorder=2, edgecolor='k') #colours[product_idx]
    ax1.scatter(impact_syn.at_event, payout_syn[product].at_event, 
                label='WISC_p (1980-2022)', c='#1b9e77', alpha=0.4)
    
    basis_risk_line = np.linspace(0, maximum*1.5)
    ax1.plot(basis_risk_line, basis_risk_line, c='k', ls='--', alpha=0.2)
    ax1.set_ylabel('Modelled payout (Hospital units)')
    ax1.set_xlabel('Modelled impact (Hospital units)')
    ax1.set_ylim(0, maximum)   
    ax1.set_xlim(0, maximum)
    plt.legend()
    
    return ax1

def add_annotation_MOZ(impact_his, payouts_his, one_best, ax1):

    big_events = np.argsort(impact_his.at_event)[-3:]
    position_correction_x = [0.09, -0.15, -0.15] #-0.3, 0.3, -0.6, 
    position_correction_y = [-0.57, 0.04, 0.03] #0, -1.5, 0.2, 
    impact_big_events = impact_his.at_event[big_events]
    payouts_big_events = payouts_his[one_best].at_event[big_events]
    #event_ids = np.asarray(impact_his.event_name)[big_events]
    event_names = [ 'Favio (2007)', 'Kenneth (2019)', 'Idai (2019)']
    #'Eline:Leone (2000)', 'Japhet (2003)', 'Filao (1988)', 
    for idx_name, name in enumerate(event_names):
        ax1.annotate(name, (impact_big_events[idx_name]+position_correction_y[idx_name], 
                            payouts_big_events[idx_name]+position_correction_x[idx_name]))
    return ax1
    

# def construct_point_from_coord(lat, lon, crs):
#     point = [lat, lon]
#     #point = [-66.1, 18.4]
#     geometry = [Point(point)]
#     gdf = GeoDataFrame(crs=crs, geometry=geometry)
    
#     return gdf

# def construct_circle(gdf, radius, crs_orig):
#     # convert CRS to equal-area projection, the length unit is now `meter`
#     #enthält gdf sein crs?
#     gdf_meter = gdf.geometry.to_crs(epsg=6933)
#     circle_meter  = gdf_meter.buffer(radius, cap_style=1)
#     circle_degree  = circle_meter.to_crs(crs_orig)
    
#     return circle_degree

# def clip_intensity(tc, circles):
#     im_val = np.max(tc.intensity, axis=0).toarray().transpose()
#     values = pd.DataFrame(im_val, columns={'values'})
#     lon = tc.centroids.lon
#     lat = tc.centroids.lat
#     haz_gdf = geopandas.GeoDataFrame(
#         values, geometry=geopandas.points_from_xy(lon, lat))
#     haz_gdf.crs = tc.centroids.crs
    
#     idx_list = []
#     for circle in circles.geometry: 
#         clipped = geopandas.clip(haz_gdf, circle)
#         idx = clipped.index.values
#         idx_list.extend(idx.tolist())

#     unique_idx= np.unique(np.asarray(idx_list))
#     tc_circle2 = copy.deepcopy(tc)
#     tc_circle2.intensity[:] = 0
#     tc_circle2.intensity[:,unique_idx] = tc.intensity[:,unique_idx].tocsr()
    
#     return tc_circle2


# def EDI(list_payouts, impact_syn, coord_exp=None):
#     F = np.zeros(len(list_payouts))
#     H = np.zeros(len(list_payouts))
#     EDI = np.zeros(len(list_payouts))
#     F = np.full(exp_gdf.index.shape[0], False)
#     exp_points[random_idx] = True
#     for idx_product, product in enumerate(list_payouts):
#         false_positive = (np.where(impact_syn.imp_mat.todense() == 0) and np.where(product.imp_mat.todense() !=0))
        
#         F[idx_product] = mean_squared_error(impact_syn.select(coord_exp=coord_exp).at_event, product.select(coord_exp=coord_exp).at_event)
#         RMSE[idx_product] = math.sqrt(MSE[idx_product])
    
#     return RMSE    

# def calibrate_parametric_product(payout_options, radius, output_folder, hospitals, impact_syn):
#     MSE = np.zeros((len(radius), len(payout_options)))
#     RMSE = np.zeros((len(radius), len(payout_options)))
#     for idx_rad, rad in enumerate(radius):
#         circle_max_filename = 'TC_circle_max_'+str(int(rad/1000))+'km.h5'
#         tc_circle_max = TropCyclone.from_hdf5(os.path.join(output_folder, circle_max_filename))
        
#         for idx_opt, option in enumerate(payout_options):
#             impf_payout = payout_from_TCcategory(option)
#             payout_hospitals = Impact()
#             payout_hospitals.calc(hospitals, impf_payout, tc_circle_max.select(orig=False), save_mat = True)
            
#             MSE[idx_rad, idx_opt] = mean_squared_error(impact_syn.at_event, payout_hospitals.at_event)
#             RMSE[idx_rad, idx_opt] = math.sqrt(MSE[idx_rad, idx_opt])
    
#     idx_min_errors = np.where(RMSE == np.min(RMSE))
    
#     return idx_min_errors

# def calibrate_parametric_product_coast(payout_options, radius, output_folder, hospitals, impact_syn, inland_thr):
#     MSEcoast = np.zeros((len(radius), len(payout_options)))
#     RMSEcoast = np.zeros((len(radius), len(payout_options)))
    
#     MSEinland = np.zeros((len(radius), len(payout_options)))
#     RMSEinland = np.zeros((len(radius), len(payout_options)))
    
#     for idx_rad, rad in enumerate(radius):
#         circle_max_filename = 'TC_circle_max_'+str(int(rad/1000))+'km.h5'
#         tc_circle_max = TropCyclone.from_hdf5(os.path.join(output_folder, circle_max_filename))
        
#         tc_coast, tc_inland= split_by_coastdistance(tc_circle_max, inland_thr)
        
#         for idx_opt, option in enumerate(payout_options):
#             impf_payout = payout_from_TCcategory(option)
            
#             payout_coast = Impact()
#             payout_coast.calc(hospitals, impf_payout, tc_coast.select(orig=False))
#             MSEcoast[idx_rad, idx_opt] = mean_squared_error(impact_syn.at_event, payout_coast.at_event)
#             RMSEcoast[idx_rad, idx_opt] = math.sqrt(MSEcoast[idx_rad, idx_opt])
            
#             payout_inland = Impact()
#             payout_inland.calc(hospitals, impf_payout, tc_inland.select(orig=False))
#             MSEinland[idx_rad, idx_opt] = mean_squared_error(impact_syn.at_event, payout_inland.at_event)
#             RMSEinland[idx_rad, idx_opt] = math.sqrt(MSEinland[idx_rad, idx_opt])
    
#     idx_min_error_coast = np.where(RMSEcoast == np.min(RMSEcoast))
#     idx_min_error_inland = np.where(RMSEinland == np.min(RMSEinland))
    
#     return idx_min_error_coast, idx_min_error_inland

# def compute_sub_pools(pool_size, repetitions, exp_gdf, payouts_syn, impact_syn):
#     coord_exp = np.vstack((exp_gdf.latitude.values, exp_gdf.longitude.values)).T

    
#     min_RMSE = np.zeros((len(pool_size), repetitions))
#     # RMSE_his = np.zeros((len(pool_size), repetitions))
#     # min_radius = np.zeros((len(pool_size), repetitions))
#     # min_structure = np.zeros((len(pool_size), repetitions, 4))
#     for idx_pool, risk_pool in enumerate(pool_size):
        
#         for repetition in range(repetitions):
        
#             random_idx = np.random.permutation(exp_gdf.index)[:risk_pool]
#             exp_points = np.full(exp_gdf.index.shape[0], False)
#             exp_points[random_idx] = True
#             RMSE = RMSE_payouts(payouts_syn, impact_syn, coord_exp[exp_points])
#             min_RMSE[idx_pool, repetition] = np.min(RMSE)
#             #idx_min_errors = np.where(RMSE == np.min(RMSE)) 
#             # min_radius[idx_pool, repetition] = np.asarray(list_radius)[idx_min_errors][0]
#             # min_structure[idx_pool,repetition, :] = np.asarray(list_payout_structure)[idx_min_errors][0]
            
#             #RMSE_his[idx_pool, repetition] = np.mean(RMSE_payouts([payouts_his[idx_min_errors[0][0]]], impact_his, coord_exp[exp_points]))

#     return min_RMSE







    
# def plot_subpool_distribution(min_RMSE, pool_size):
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.boxplot(min_RMSE.T)
#     locs, xlabels = plt.xticks()
#     newlabel = pool_size.tolist()
#     plt.xticks(locs, newlabel)
#     ax.set_xlabel('Pool size (Hospitals)')
#     ax.set_ylabel('RMSE')
    
#     return ax

# def plot_subpool(RMSE_syn, RMSE_his, pool_size):
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.scatter(pool_size, RMSE_syn, c='k', label='Probablisitic')
#     ax.scatter(pool_size,RMSE_his, c='r', label='Historic')
#     plt.legend()
#     ax.set_xlabel('Pool size (Hospitals)')
#     ax.set_ylabel('RMSE')
#     ax.set_ylim(0,0.12)   
#     ax.set_xlim(0,160)
    
#     return ax

# def plot_impact_comparison(payout_his, impact_his):
#     fig3 = plt.figure()
#     ax1 = fig3.add_subplot(111)
#     maximum = np.max(payout_his.at_event)
#     basis_risk_line = np.linspace(0, maximum*1.5)
#     ax1.plot(basis_risk_line, basis_risk_line, c='k', ls='--', alpha=0.2)
#     #ax1.fill_between(impact_real, impact_real, maximum*1.001, color='#7570b3', alpha=0.1)
#     #ax1.fill_between(impact_real, impact_real, color='#f33', alpha=0.1)
#     ax1.scatter(impact_his.at_event, payout_his.at_event, c ='#e6550d', label='Historic events', marker = 'x')
#     # ax1.scatter(imp_hospitals.at_event, payout_hospitals.at_event, c ='k', label='Synthetic events', marker = 'x')
    
#     #ax1.scatter(impact_real,impact2, c = '#d95f02',  label='Product 2', marker = 'x')
#     #ax1.scatter(impact_real,impact3, c = '#f33', label='Product 3', marker = '+')
#     ax1.set_ylabel('Payout (Hospitals)')
#     ax1.set_xlabel('Modelled impact (Hospitals)')
#     ax1.set_ylim(0,6)   
#     ax1.set_xlim(0,6)
    
#     return ax1




# def split_by_coastdistance(haz, inland_thr):    

#     distance_coast = haz.centroids.dist_coast
#     [centr_coast] = np.where(distance_coast <= inland_thr)
    
#     tc_coast = copy.deepcopy(haz)
#     coast_intensity = sparse.lil_matrix((haz.intensity.shape))
#     coast_intensity[:, centr_coast] = haz.intensity[:, centr_coast]
#     tc_coast.intensity = coast_intensity.tocsr()
    
#     tc_inland = copy.deepcopy(haz)
#     tc_inland.intensity[:, centr_coast] = 0
#     tc_inland.intensity.tocsr()
    
#     return tc_coast, tc_inland









