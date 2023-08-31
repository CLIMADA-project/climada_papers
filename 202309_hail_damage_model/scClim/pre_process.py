""" 
Function to pre-process damage data
Pre-processing is based on POH-values with 5km Buffer
"""
import numpy as np
import pandas as pd 
import xarray as xr
import sys, os
import timeit
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.markers as mmarkers
import datetime as dt
from climada import util
import geopandas as gpd
from climada.entity import Entity, Exposures
from climada.hazard import Hazard,Centroids
from scipy import sparse
import timeit
import cartopy.crs as ccrs
import cartopy.feature as cf
from climada.entity import ImpactFunc, ImpactFuncSet
from climada.engine import Impact
import climada.util.coordinates as u_coord

def correct_impact_dates(
    impacts,
    haz_PH,
    haz_poh,
    haz_meshs,
    poh,
    dmg_file_template,
    param_dict,
    log_file,
    
):
    np.testing.assert_array_equal(haz_PH.event_id,haz_poh.event_id,
                                  "PH and poh events don't match")
    np.testing.assert_array_equal(haz_PH.event_id,haz_meshs.event_id,
                                  "PH and meshs events don't match")

    if not isinstance(impacts,gpd.GeoDataFrame): 
        raise UserWarning("impacts as dict are deprecated. Use GeoDataFrame instead.")

    #Remove nonZero damage entries
    zero_entries = sum(impacts['value']==0)
    if zero_entries>0:
        print("Removing non-zero damage entries")
        impacts = impacts.loc[impacts['value']>0,:]
        impacts.reset_index()

    assign_centroids_gdf(impacts,haz_poh)
    pre_2002_entries = sum(pd.to_datetime(impacts.date_dt).dt.year<2002)
    non_season_entries = sum(((pd.to_datetime(impacts.date_dt).dt.month<4) | 
                                (pd.to_datetime(impacts.date_dt).dt.month>9)) &
                                (pd.to_datetime(impacts.date_dt).dt.year>=2002))
    n_valid_entries = impacts.shape[0]-non_season_entries-pre_2002_entries
    dates = [pd.Timestamp(date) for date in 
            impacts.date_dt.unique() if pd.Timestamp(date).year>=2002
            and pd.Timestamp(date).month>=4
            and pd.Timestamp(date).month<=9]

    dates.sort()
    # date = dates[-79] #21.6.2021
    # date = dates[-58] #12.7.2021
    # date = dates[86] #12.8.2004

    #Initialize dataframe
    corrected_dmg = dmg_file_template
    if 'date_dt' not in corrected_dmg.columns:
        ValueError("'date_dt' is not in columns")
    corrected_dmg = corrected_dmg.assign(POH=[],MESHS=[])
    removed_entries_total = 0 #initiate counter for removed entries
    moved_entries_total = 0 
    moved_to_previous_day = 0
    #Loop over all dates
    for date in dates:
        print(date)
        # plot = True if date.year == 2021 else False
        plot=False
        # determine string_date/event_index(not id) and select today's damages
        str_date = dt.datetime.strftime(date, "%d%m%Y")
        event_index = np.where(haz_poh.date == date.toordinal())[0][0]
        
        #select entries of this date
        if isinstance(impacts,dict): 
            raise UserWarning("impacts as dict are deprecated. Use GeoDataFrame instead.")
            # dmg_gdf = exp_dmg_all[str_date].gdf
            # assign_centroids_gdf(dmg_gdf,haz_poh)
        elif isinstance(impacts,gpd.GeoDataFrame):
            dmg_gdf = impacts.loc[impacts.date_dt == date,:]


        # Assign a POH value based on climada centroids
        dmg_gdf['POH_orig']=haz_poh.intensity[event_index,dmg_gdf.centr_HL].toarray().squeeze()
        dmg_gdf['MESHS_orig']=haz_meshs.intensity[event_index,dmg_gdf.centr_HL].toarray().squeeze()
        
        #get day_adjustment array depending on min_day and max_day
        day_adj_arr = np.arange(param_dict['min_day'],param_dict['max_day']+1) #for +-2days: [-2,-1,0,1,2]
        win_length = len(day_adj_arr)

        #Check if day is within available hazard information (not at beginning/end of season)
        date_exists = np.array([date.toordinal()+day_adj in haz_poh.date for 
                                day_adj in day_adj_arr])
        max_ind = event_index+max(np.array(day_adj_arr)[date_exists])+1
        min_ind = event_index+min(np.array(day_adj_arr)[date_exists])
        index_day0 = np.where(np.array(day_adj_arr)[date_exists]==0)[0][0] #make sure to get correct day0 index
        # event_indices = np.arange(event_index-2,event_index+2+1)[date_exists]

        #Get Hazard info (PH,POH,MESHS) for neighbouring days too
        now_PH=haz_PH.intensity[min_ind:max_ind,dmg_gdf.centr_HL].toarray()
        now_POH=haz_poh.intensity[min_ind:max_ind,dmg_gdf.centr_HL].toarray()
        now_MESHS=haz_meshs.intensity[min_ind:max_ind,dmg_gdf.centr_HL].toarray()

        
        if param_dict["use_likelihood"]==True:
            assert(param_dict['min_POH_diff'] is None) #likelyhood cannot be combined with min_POH_diff
            assert(win_length==5) #only implemented for +-2 days
            #Set default weights. relative to day0 with a value of 1
            # e.g. for day-1 (which is second most likely due to common hail events at night
            # that are reported the next day) the weight is 0.5. Thus the POH-based likelihood 
            # on day-1 must be 50% higher for the date to change
            order_weight = np.array([0.1, 0.5, 1, 0.1, 0.1])

        else:
            # define default order weight (in case POH-difference is not analysed)
            # create an array with maximum weight for day0, and decreasing weight for furhter away days
            # between +X and -X days, the weight for minus is higher, because often hail is 
            # reported on the next day (thunderstorm at night) rather than a day too early
            order_weight = win_length-2*abs(day_adj_arr)+(day_adj_arr<0).astype(int)
            #default (for +-2d) is [2,4,5,3,1]
            
        #choose sub-selection of order_weight if not all dates exist
        order_weight = order_weight[date_exists]
        
        if param_dict['min_POH_diff']:
            #If POH difference is considered, get the POH difference between day0 and all other days
            poh_diff = now_POH-now_POH[index_day0,:][None,:]
            #set POH_diff of day0 to min_POH diff: Then only if POH_diff>min_POH_diff the date will be changed
            poh_diff[index_day0,:]=param_dict['min_POH_diff']
            
            # multiply PH with (1+poh_diff): +1 because if there is no POH for both days, 
            # we still want to shift to the PH value (within buffer)
            weighted_PH = now_PH*(1+poh_diff)
            #adjust by weight only if POH is equal (mostly if POH=0 for all)
            weighted_PH = weighted_PH*(1+0.001*order_weight[:,None]) 
        else:
            #default weighing, only considering the "possible hail (PH)" variable
            weighted_PH = now_PH*order_weight[:,None]
            print('default weighing')

        #assign new date based on PH value and add new POH/MESHS values
        dmg_gdf['new_date_dt']=np.apply_along_axis(get_date,axis=0,arr=weighted_PH,date=date,
                                                   index_day0 = index_day0,mode='date')
        dmg_gdf['POH']=np.apply_along_axis(get_date,axis=0,arr=np.concatenate((weighted_PH,now_POH),axis=0),
                                           date=date,index_day0 = index_day0,mode='poh')
        dmg_gdf['MESHS']=np.apply_along_axis(get_date,axis=0,arr=np.concatenate((weighted_PH,now_MESHS),axis=0),
                                             date=date,index_day0 = index_day0,mode='meshs')
        dmg_gdf['category']=np.apply_along_axis(get_date,axis=0,arr=weighted_PH,
                                                date=date,index_day0 = index_day0,
                                                mode='category')    
        
        moved_entries_total += sum(dmg_gdf['new_date_dt']!=dmg_gdf['date_dt'])
        #calculate how many entries are moved to previous day
        prev_day = [date-dt.timedelta(days=1) for date in dmg_gdf['date_dt']]
        moved_to_previous_day += sum(np.array(prev_day)==dmg_gdf['new_date_dt'])
        
        # # Remaining plot functions
        # if (plot and dmg_gdf.shape[0]>100):# or dmg_gdf.shape[0]>1000:
        #     plot_dBZ_dates = ['23062002','20062021','21062021','12072021']
        #     fig = sc.plot_funcs.plot_pre_process(exp_dmg,date,poh,haz_poh,ds_PH,
        #     event_index,plot_dBZ_dates,poh_level=param_dict['poh_level'],min_day=param_dict['min_day'], 
        #     max_day=param_dict['max_day'],extent=ZRH_EXTENT) 
        #     fig.savefig(out_dir + '/explorative/pre_process/v%d_%s.png'%(version_id,dt.datetime.strftime(date, "%Y-%m-%d")),
        #                 dpi=250, bbox_inches='tight') 
        #     plt.close(fig)
        
        #Add values to dataframe 
        if param_dict['delete_nonPOH_dmgs']:
            #select only damage entries where possible hail was detected in +-2 days
            ph_bol = dmg_gdf['category'] != 'none'
            new_df = dmg_gdf.loc[ph_bol,corrected_dmg.columns]
            #new_df['Schadendatum'] = dmg_gdf.loc[ph_bol,'new_date']
            new_df['date_dt'] = dmg_gdf.loc[ph_bol,'new_date_dt']
            if sum(ph_bol) != new_df.shape[0]:
                raise ValueError('Error: Number of removed entries does not match!')
            removed_entries_total += sum(~ph_bol)
        else:
            new_df = dmg_gdf[corrected_dmg.columns]
            #new_df['Schadendatum'] = dmg_gdf['new_date']
            new_df['date_dt'] = dmg_gdf['new_date_dt']
        corrected_dmg = pd.concat([corrected_dmg,new_df],axis = 0,ignore_index = True)
    
    #Remove Schadendatum column! (work with date_dt going forward)
    if 'Schadendatum' in corrected_dmg.columns:
        corrected_dmg = corrected_dmg.drop(columns=['Schadendatum'])

    if not corrected_dmg.shape[0]+removed_entries_total+pre_2002_entries+non_season_entries == impacts.shape[0]:
        raise UserWarning('New shape does not match old shape-removed entries')
    with open(log_file, 'a') as f:
        f.write('---------------------------------------------\n')
        f.write('File: XXX_Hail_Loss_date_corrected%d.csv\n'%param_dict['version_id'])
        f.write(f'Removed {pre_2002_entries} entries before 2002 and {non_season_entries} not within Apr-Sept\n')
        f.write('Removed %d (of %d) entries with zero damage\n'%(zero_entries,zero_entries+n_valid_entries))
        f.write('min_day: %d, max_day: %d, poh_level: %d, buffer: %dkm, delete_nonPOH_dmgs:%s, min_POH_diff:%s\n'%(
        param_dict['min_day'], param_dict['max_day'], param_dict['poh_level'], param_dict['buffer_km'],
        param_dict['delete_nonPOH_dmgs'],param_dict['min_POH_diff']))
        f.write('removed entries: %d from %d \n'%(removed_entries_total,n_valid_entries))
        f.write('moved entries: %d from %d. (Of which %d to previous day) \n'%(moved_entries_total,
                                                                            corrected_dmg.shape[0],
                                                                            moved_to_previous_day))
            
    return(corrected_dmg)


def get_date(in_array,date,index_day0,mode='date'):
    """new get_date function: gets date or POH/MESHS values for the correct date
    based on weighted_PH(possible_hail) values of an event +-some days

    Parameters
    ----------
    in_array : np.array (1D, but applied to 2D; i.e. second axis)
        either weigthed_PH array: [w1,w2,w3,w4,w5]
        or weighted_PH+MESHS/POH: [w1,w2,w3,w4,w5,POH1,POH2,POH3,POH4,POH5]
    date : datetime
        date in question
    index_day0 : int
        index of original date with array of +-XX days
    mode : str
        whether to return the new date, or the plotting category

    Returns
    -------
    either date or category or MESHS/POH value
    
    """         
    #first check nan values
    if in_array.sum() == np.nan:
        raise ValueError("weight_array contains NaN values")

    if mode=='poh' or mode=='meshs':
        # for POH and MESHS calcs the weigthed_PH array is structured 
        # as weighted_PH concatenated with POH/MESHS
        weighted_PH,meshs_poh =np.split(in_array,2)
    else:
        weighted_PH = in_array

    if max(weighted_PH) == 0 : #if all values are zero (no possible hail)
        max_index = index_day0
        new_date = date
        category = 'none'
    else:
        max_index = weighted_PH.argmax()
        new_date = date + pd.Timedelta(max_index-index_day0,'d')
        if max_index == index_day0:
            category = 'day0'
        else:
            category = 'shifted'
    
    if mode == 'date':
        #return dt.datetime.strftime(new_date, "%d%m%Y")
        return new_date
    elif mode == 'category':
        # print(category)
        return category#np.array(category,dtype='<U10')
    elif mode == 'poh' or mode=='meshs':
        return meshs_poh[max_index]    

def assign_centroids_gdf(gdf, hazard, distance='euclidean',
                        threshold=u_coord.NEAREST_NEIGHBOR_THRESHOLD):
    """Assign for each exposure coordinate closest hazard coordinate.
    -1 used for disatances > threshold in point distances. If raster hazard,
    -1 used for centroids outside raster.
    Parameters
    ----------
    hazard : Hazard
        Hazard to match (with raster or vector centroids).
    distance : str, optional
        Distance to use in case of vector centroids.
        Possible values are "euclidean", "haversine" and "approx".
        Default: "euclidean"
    threshold : float
        If the distance (in km) to the nearest neighbor exceeds `threshold`,
        the index `-1` is assigned.
        Set `threshold` to 0, to disable nearest neighbor matching.
        Default: 100 (km)
    See Also
    --------
    climada.entity.exposures.base for details: same function for exposure.gdf's
    climada.util.coordinates.assign_coordinates: method to associate centroids to
        exposure points
    """

    if not u_coord.equal_crs(gdf.crs, hazard.centroids.crs):
        raise ValueError('Set hazard and exposure to same CRS first!')
    if hazard.centroids.meta:
        assigned = u_coord.assign_grid_points(
            gdf.longitude.values, gdf.latitude.values,
            hazard.centroids.meta['width'], hazard.centroids.meta['height'],
            hazard.centroids.meta['transform'])
    else:
        assigned = u_coord.assign_coordinates(
            np.stack([gdf.latitude.values, gdf.longitude.values], axis=1),
            hazard.centroids.coord, distance=distance, threshold=threshold)
    gdf['centr_' + hazard.tag.haz_type] = assigned
        