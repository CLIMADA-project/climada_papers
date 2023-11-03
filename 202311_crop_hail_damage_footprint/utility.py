# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:49:23 2023

Utility functions for Paper: 
Portmann R., Schmid T., Villiger L., Bresch D., Calanca P. Modelling crop hail damage footprints
with single-polarization radar: The roles of spatial resolution, hail intensity,
and cropland density, submitted to Natural Hazards and Earth System Sciences.

@authors: Raphael Portmann, Timo Schmid, Leonie Villiger
"""
import pickle
import pandas as pd
import numpy as np
import matplotlib
import geopandas
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from cartopy.io import shapereader
import cartopy.feature as cf
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import sys
import xarray as xr
from climada.entity import Exposures
from climada.hazard import Hazard, Centroids
from climada.engine import Impact

import shapely.geometry as sg

import geopandas as gpd
from pyproj import CRS
from scipy import sparse
import copy
import warnings

#define translations of names
en_names={'Weizen_Mais_Gerste_Raps': 'field crops',
          'Weizen_Mais_Raps_Gerste': 'field crops',
          'wheat_maize_rapeseed_barley': 'field crops',
          'wheat_maize_barley_rapeseed': 'field crops',
          'Weizen': 'wheat',
          'Mais': 'maize',
          'Raps': 'rapeseed',
          'Gerste':'barley',
          'Reben': 'grapevine'}

#data directory
data_dir='C:/Users/F80840370/projects/scClim/climada/scClim/subproj_D/papers/NHESS/code_and_data/data/'

#%% read data
def read_at_centroid_data(datadir,croptypes,variable='MESHS',sample_id=None):
        """
        
    
        Parameters
        ----------
        datadir : str
            directory of data at centroid .p files
        croptypes : list of str
            croptypes to be read, if more than one the data is merged together
        variable : str, optional
            hazard variable name. The default is 'MESHS'.
        sample_id : str, optional
            ID of the random sample (only used if random samples are read. The default is None.
    
        Returns
        -------
        at_centroid_data_out : pandas.DataFrame
            Dataframe including exposure, damage, and hazard information at hazard centroids
        croptype : str
            croptype (if more than one, they are concatenated with a '_').
    
        """
        
        at_centroid_list=[]
        for croptype in croptypes:
            # read dictionary from pickle file
            if sample_id is not None:
                name = f'values_at_centroid_{variable}_1_2_4_8km_{croptype}_{sample_id}.p'
                if croptype == 'Reben':
                    name = f'values_at_centroid_{variable}_1_2_4_8km_{croptype}_{sample_id}.p'
            else:
                if variable=='E_kin':
                    name = f'data_at_centroid_{variable}_CC_{croptype}.p'
                else:
                    name = f'data_at_centroid_{variable}_{croptype}.p'
            with open(datadir+name, 'rb') as file:
                at_centroid_list.append(pickle.load(file))           
        if len(croptypes)>1:
            at_centroid_data_AGG={}
            at_centroid_data_MESHS={}
            for km in at_centroid_list[0].keys():
                at_centroid_data_MESHS[km] = pd.concat([d[km] for d in at_centroid_list])
                at_centroid_data_AGG[km] = at_centroid_data_MESHS[km].groupby(['centr_HL','date']).agg({'n_exp': 'sum', 'n_dmgs': 'sum', variable: 'mean'})
                at_centroid_data_AGG[km]['PAA'] = at_centroid_data_AGG[km]['n_dmgs']/at_centroid_data_AGG[km]['n_exp']
                at_centroid_data_AGG[km]['PAA'].replace(0, np.nan, inplace=True)
                at_centroid_data_AGG[km]['PAA'].loc[at_centroid_data_AGG[km]['PAA']>1] = 1
                at_centroid_data_AGG[km]['date'] = at_centroid_data_AGG[km].index.get_level_values(1)

            croptype='_'.join(croptypes)
            at_centroid_data_out=at_centroid_data_AGG
        else:
            at_centroid_data_out=at_centroid_list[0]
    
    
        return at_centroid_data_out, croptype

def exposure_from_gpd(filename, description, value_colname = 'area_ha', value_unit = 'ha'):

    """
    read exposure from a file that is readable as geopandas.GeoDataFrame.
    Here used to read crop-specific exposure data gridded on a 1x1km grid

    Parameters
    ----------
        filename:  string
            Name of the file to read as GeoDataFrame and create Exposure from
        description: string
            String describing the exposure. Is added to the Exposure.tag.
        value_column: string
            Name of the column in the GeoDataFrame that is interpreted as 'value' column of the Exposure (default: 'area_ha')
        value_unit: string
            unit of the exposure values (default: ha)
        return_cell_gdf: boolean
            whether or not to return the original GeoDataFrame with geometry of the gridcell saved as polygons (which is neat for plotting)
            default: False
    Returns
    -----------
        exp: climada.entity.exposures.base.Exposures
            Climada exposure
        cell_gdf: GeoDataFrame
            orignal GeoDataFrame with gridcells explicityl specified in geometry
    """

    filename = str(filename)
    value_colname = str(value_colname)
    description = str(description)

    #read file
    gdf = gpd.read_file(filename)

    #gdf=gdf.fillna(0)
    #create column for impact function
    gdf['impf_'] = 1

    #define the value column
    gdf=gdf.rename(columns = {value_colname : "value"})

    #set up exposure
    exp = Exposures(gdf, value_unit = value_unit)

    #set latitude longitude columns from geometry points
    exp.set_lat_lon()

    #tag exposure (not used after CLIMADA 4.0)
    #exp.tag.file = filename
    #exp.tag.description = description

    #check
    exp.check()

    #remove exposure points with value zero
    exp.gdf = exp.gdf[exp.gdf['value']>0]

    return exp

def read_gridded_damages_as_impact(filename_dmg,impact_metric='n_fields_dmg', binary_damage=False):
    """
    Read damage netcdf data as climada.impact

    Parameters
    ----------
    filename_dmg : str
        name of the damage file
    impact_metric : str, optional
        variable used as impact metric. The default is 'n_fields_dmg'.
    binary_damage : Boolean, optional
        Whether to return a binary field (damage yes/no). The default is False.

    Returns
    -------
    imp_out : climada.enginge.Impact
        Resulting climada impact object 

    """

    damages_data=xr.open_dataset(filename_dmg)

    #Group by coordinates and get df 
    df_all=damages_data.to_dataframe().reset_index().fillna(0)
    df_all_new=df_all[df_all['n_fields_exp']>0]
    coord_df=df_all_new.groupby(['chy','chx']).first()[['lat','lon']]

    #make coordinate array from array of tuples
    coord_exp=coord_df[['lat','lon']].values

    #Initialize impact matrix as dataframe with coordinates as index
    imp_df = pd.DataFrame(index=coord_df.index)

    #get all damage dates
    dates=list(df_all_new.groupby('time').groups.keys())

    #loop over dates and extend concatenate impacts for individual dates as new columns to impact dataframe (imp_df)
    for date in dates:        
        temp_df=df_all_new[df_all_new['time']==date].set_index(['chy','chx'])[impact_metric]
        if binary_damage == True:   #if binary damage is True set impact to 1 at each grid cell
            temp_df=temp_df/temp_df
        imp_df=pd.concat([imp_df,temp_df],axis=1)

    #fill Nan values and reset coordinate index to integer index
    imp_df = imp_df.fillna(0).reset_index().drop(['chx','chy'],axis=1)

    #create a sparse matrix from the impact dataframe
    imp_mat = sparse.csr_matrix(imp_df.T.values)

    #Prepare climada Impact attributes
    #Date and frequency
    ord_dates = np.array([d.toordinal() for d in dates])
    ev_names = [d.strftime('ev_%Y-%m-%d') for d in dates]
    n_years = np.floor((max(ord_dates)-min(ord_dates))/365)
    frequency = np.ones(imp_mat.shape[0])
    event_id = np.arange(imp_mat.shape[0])+1


    #Damages
    imp_at_event = imp_mat.sum(axis=1).getA1()
    imp_eai_exp = imp_mat.sum(axis=0).getA1()/n_years
    imp_aai_agg = imp_mat.sum()/n_years

    #initialize impact object
    imp_out = Impact(
                    imp_mat = imp_mat,
                    coord_exp = coord_exp,
                    crs = 'EPSG:4326',
                    date = ord_dates,
                    event_name = ev_names,
                    frequency = frequency,
                    event_id = event_id,
                    at_event = imp_at_event,
                    eai_exp = imp_eai_exp,
                    aai_agg= imp_aai_agg,
                    unit = '')
    return imp_out

# Hazard
def hazard_from_radar(files, varname='MESHS', time_dim='time', forecast_init=None, 
                      ensemble_dim=None, spatial_dims = None, country_code=None,
                      extent=None, subdaily = False, month=None, ignore_date=False, 
                      n_year_input=None, get_xarray=False):
    """Create a new Hail hazard from MeteoCH radar data 
    or COSMO HAILCAST ouput (single- or multi-member)

    Parameters
    ----------
    files : list of str or xarray Dataset
        list of netcdf filenames (string) or xarray Dataset object
    varname : string
        the netcdf variable name to be read from the file
    time_dim : str
        Name of time dimension, default: 'time'
    forecast_init : datetime object
        List with datetimes of forecast initializations,
        needs to have same length as time, default: None
    ensemble_dim : str
        Name of ensemble dimension, default: None
    spatial_dims : list of str
        Names of spatial dimensions
    country_code : int
        ISO 3166 country code to filter the data
    extent : list / array
        [lon_min, lon_max, lat_min, lat_max]
    ignore_date : boolean
        If True: ignores netcdf dates (e.g. for synthetic data).
    n_year_input : int
        Number of years: will only be used if ignore_date=True
    Returns
    -------
    haz : Hazard object
        Hazard object containing radar data with hail intensity (MESHS)
    """

    #Initialize default values
    if spatial_dims is None: spatial_dims = ['chy','chx']

    #read netcdf if it is given as a path
    if type(files) == xr.core.dataset.Dataset:
        netcdf = files
    else:
        netcdf = xr.open_mfdataset(files, concat_dim=time_dim, combine='nested', 
                                   coords='minimal')

    #select month of the year if given
    if month:
        grouped=netcdf.groupby("time.month")
        netcdf=grouped[int(month)]

    #Cut data to selected country/area only
    if extent:
        lon_min, lon_max, lat_min, lat_max = extent
        lon_cond = np.logical_and(netcdf.lon >= lon_min, netcdf.lon <= lon_max)
        lat_cond = np.logical_and(netcdf.lat >= lat_min, netcdf.lat <= lat_max)
        netcdf = netcdf.where(np.logical_and(lat_cond,lon_cond),drop=True)

    # #stack data
    # stacked = netcdf.stack(new_dim=spatial_dims)

    # if country_code:
    #     c_code = util.coordinates.get_country_code(stacked.lat,stacked.lon)
    #     stacked = stacked.assign_coords(country_code=("new_dim",c_code))
    #     stacked = stacked.where(stacked.country_code==country_code,drop=True)

    #Select variable and set units
    varname_xr = varname #by default the varname corresponds to the xr name
    if varname == 'MESHS' or varname == 'MESHS_4km':
        varname_xr = 'MZC'
        unit = 'mm'
    elif varname == 'MESHSdBZ' or varname == 'MESHSdBZ_p3':
        varname_xr = 'MESHSdBZ'
        unit = 'mm'
    elif varname == 'POH':
        varname_xr = 'BZC'
        unit = '%'
    elif 'DHAIL' in varname:
        unit = 'mm'
    elif varname == 'dBZ' or varname =='dBZfiltered':
        varname_xr = 'CZC'
        unit = 'dBZ'
        # Filter values for efficient calculation. dBZ<40 are set to zero
        netcdf = netcdf.where(netcdf[varname_xr]>40,0)
    elif varname == 'possible_hail':
        unit = '[ ](boolean)'
    elif varname == 'durPOH':
        varname_xr = 'BZC80_dur'
        netcdf[varname_xr] = netcdf[varname_xr]*5 #times 5 to get minutes
        unit = '[min]'
    elif varname == 'MESHSweigh':
        unit = 'mm (scaled by duration)'
    elif varname == 'HKE':
        unit = 'Jm-2'
    elif varname == 'crowd' or varname=='crowdFiltered':
        warnings.warn('use smoothed data for crowd-sourced data')
        varname_xr = 'h_smooth'
        unit = 'mm'
    elif varname == 'E_kin' or varname=='E_kinCC': #E_kin from Waldvogel 1978, or Cecchini 2022
        varname_xr = 'E_kin'
        unit = 'Jm-2'
    elif varname == 'VIL':
        unit = 'g/m2'
        varname_xr = 'dLZC'
        # Filter values for efficient calculation. VIL<10g/m2 are set to zero
        netcdf = netcdf.where(netcdf[varname_xr]>10,0).round(0)
    else:
        raise ValueError(f'varname "{varname}" is not implemented at the moment')

    #prepare xarray with ensemble dimension to be read as climada Hazard
    if ensemble_dim:
        # omit extent if ensemble_dim is given
        if extent:
            warnings.warn("Do not use keyword extent in combination with "
                          "ensemble_dim. Plotting will not work.")
        # omit igonore_date if ensemble_dim is given
        if ignore_date:
            warnings.warn('Do not use keyword ignore_date in combination with '
                          'ensemble_dim. Event names are set differently.')
        # stack ensembles along new dimension
        netcdf = netcdf.stack(time_ensemble=(time_dim, ensemble_dim))
        # event names
        if forecast_init: #event_name = ev_YYMMDD_ensXX_init_YYMMDD_HH
            n_member, = np.unique(netcdf[ensemble_dim]).shape
            forecast_init = np.repeat(forecast_init, n_member)
            if netcdf[time_dim].size != len(forecast_init):
                warnings.warn("Length of forecast_init doesn't match time.")
            event_name = np.array([f"{pd.to_datetime(ts).strftime('ev_%y%m%d')}_ens{ens:02d}_{init.strftime('init_%y%m%d_%H')}"
                                   for (ts,ens),init in zip(netcdf.time_ensemble.values, forecast_init)])
        else: #event_name = ev_YYMMDD_ensXX
            event_name = np.array([f"{pd.to_datetime(ts).strftime('ev_%y%m%d')}_ens{ens:02d}"
                                   for ts,ens in netcdf.time_ensemble.values])
        #convert MultiIndex to SingleIndex
        netcdf = netcdf.reset_index('time_ensemble')
        netcdf = netcdf.assign_coords({'time_ensemble':netcdf.time_ensemble.values})
        # remove duplicates along new dimension for variables that are identical across members
        netcdf['lon'] = netcdf['lon'].sel(time_ensemble=0, drop=True)
        netcdf['lat'] = netcdf['lat'].sel(time_ensemble=0, drop=True)

    # get number of events and create event ids
    n_ev = netcdf[time_dim].size
    event_id = np.arange(1, n_ev+1, dtype=int)

    if ignore_date:
        n_years = n_year_input
        if 'year' in netcdf.coords:
            event_name = np.array(['ev_%d_y%d'%i for i in zip(event_id,netcdf.year.values)])
        else:
            event_name = np.array(['ev_%d'%i for i in event_id])
    elif ensemble_dim:
        n_years = netcdf[time_dim].dt.year.max().values-netcdf[time_dim].dt.year.min().values + 1
    else:
        n_years = netcdf[time_dim].dt.year.max().values-netcdf[time_dim].dt.year.min().values + 1
        if subdaily:
            event_name = netcdf[time_dim].dt.strftime("ev_%Y-%m-%d_%H:%M").values
        else:
            event_name = netcdf[time_dim].dt.strftime("ev_%Y-%m-%d").values

    #Create Hazard object
    if ensemble_dim:
        event_dim = 'time_ensemble'
    else:
        event_dim = time_dim
    coord_vars = dict(event=event_dim,longitude='lon',latitude='lat')
    haz = Hazard.from_xarray_raster(netcdf,'HL',unit,intensity=varname_xr,
                                    coordinate_vars=coord_vars)
    #set correct event_name, frequency, date
    haz.event_name = event_name
    haz.frequency = np.ones(n_ev)/n_years
    if ignore_date: 
        haz.date = np.array([], int)
    if ensemble_dim: 
        haz.date = np.array([pd.to_datetime(ts).toordinal() for ts in netcdf[time_dim].values])

    netcdf.close()
    haz.check()

    if get_xarray:
        return haz,netcdf
    else:
        return haz
    
def compute_affected_area_from_dmgs(exposure_damages_object, exposure):

    """ compute percent area affected from damages and exposure
    Parameters
    ---------

    exposure_damages_object: climada.entity.exposure
                damage claims of a single event stored as exposure object
    exposure:  climada.entity.exposure
                exposure data for a given croptype

    Returns
    -------

    exposure_damages_object_new: climada.entity.exposure
                a copy of exposure_damages_object but with values for PAA and total area affected"""


    exp_dmg=exposure_damages_object

    #get field area
    exposure.gdf['field_area']=exposure.gdf['value']/exposure.gdf['n_fields']
    #compute Swiss average field area to use if no field area can be defined
    field_area_mean=exposure.gdf['value'].sum()/exposure.gdf['n_fields'].sum()

    #add field area column
    exp_dmg.gdf['field_area']=np.nan

    #loop over all damages and add field area
    for index in exp_dmg.gdf.index:

           #get centroids
           centroid=exp_dmg.gdf.loc[index,'centr_HL']
           #read field area from exposure
           field_area=exposure.gdf.loc[exposure.gdf[exposure.gdf['centr_HL']==centroid].index,'field_area'].values
           if field_area.size == 1:
               exp_dmg.gdf.loc[index,'field_area']=field_area[0]
           elif field_area.size == 0:
               #issue warning if damage is recorded where no exposure is present
               exp_dmg.gdf.loc[index,'field_area']=field_area_mean
               warnings.warn('damage reported in region without exposure. swiss average field size used.')
           else:

               raise ValueError('more than one field area found.')



    return exp_dmg

def read_shv_dmg(excel_file, croptype, return_type = 'exp_dict', exposure = None, hazard = None):
    """Read Excel file

    Parameters
    ----------
    csv_file : csv
        file with GVZ exposure (or damage) data
    croptype : string
        string denoting croptype, for several croptypes use underscore (_) to separate between names
    on_grid : bool
        whether or not data should be interpolated on regular grid (not implemented yet)
    w_rel_dmg : bool
        if relative damage should be saved in gdf too
    return_type: str
        if 'imp': return imp object
        if 'exp_dict': returns dict of exposures
    haz: climada.hazard
        climada.hazard object to get event_ids. Only needed if return_type='imp'
    Returns
    -------
    exp_dict : dict
        Dictionary with Exposure objects for each event
    """

    #### start helper functions ####
    def get_damage_dict(shv_dmg, return_type):

           out_dict = {}
           dates=shv_dmg.groupby('Schadendatum').groups.keys()

           #Loop over all dates
           for date in dates:
               df_sel =  shv_dmg.loc[shv_dmg.Schadendatum==date]

               if return_type == 'exp_dict':
                   gdf_dmg = gpd.GeoDataFrame(df_sel,geometry = gpd.points_from_xy(df_sel['lon'],
                                          df_sel['lat'],crs = 'EPSG:4326'))

                   exp_dmg = Exposures(gdf_dmg,value_unit = '%')
                   exp_dmg.set_lat_lon()
                   exp_dmg.check()
                   out_dict.update({date:exp_dmg})
               elif return_type == 'gdf_dict':
                   gdf_dmg = gpd.GeoDataFrame(df_sel,geometry = gpd.points_from_xy(df_sel['X'],
                                          df_sel['Y'],crs = 'EPSG:2056'))
                   out_dict.update({date: gdf_dmg})

           return out_dict

    def adjust_croptype(croptype):
        if croptype == "Weizen":
            croptype = "Winterweizen"
        if croptype == "Gerste":
            croptype = "Wintergerste"
        if croptype == "Reben":
            croptype = "Wein"
        if croptype == "Aepfel":
            croptype = "Tafeläpfel"
        if croptype == "Mais":
            croptype = "Mais|mais"
        if croptype == "Raps":
            croptype = "Raps|raps"
        return croptype

    def read_damage_dataframe(excel_file, croptype):
        shv_dmg = pd.read_excel(excel_file) #,parse_dates=['Schadendatum'])
        # select croptype(s)
        #split if several croptypes
        croptypes=croptype.split("_")
        #adjust names to make sure all of them are found
        croptypes_corrected=[adjust_croptype(c) for c in croptypes]
        shv_dmg = shv_dmg.loc[shv_dmg['Kultur'].str.contains('|'.join(croptypes_corrected))]
        # remove zero values
        shv_dmg=shv_dmg.loc[shv_dmg['Ertragsverlust (%)']!=0]
        # set value
        shv_dmg["value"] = shv_dmg["Ertragsverlust (%)"].values

        return shv_dmg

    def get_impact_mats_for_imp_dict(shv_dmg,exposure,hazard):

            #Group by coordinates and create copys for output
            imp_df=shv_dmg.groupby(['lat','lon']).size().rename('count')
            nfields_df=imp_df.copy(deep=True)
            area_df=imp_df.copy(deep=True)
            ha_loss_df=imp_df.copy(deep=True)

            #make coordinate array from array of tuples
            coord_map=map(np.array,imp_df.index.values)
            coord_exp=np.array(list(coord_map))

            dates=shv_dmg.groupby('Schadendatum').groups.keys()

            if exposure:
                exp_dmgs=get_damage_dict(shv_dmg, return_type = 'exp_dict')

            for date in dates:

                if exposure and hazard:
                    #get dataframe with field area for each coordinate
                    exp_dmg = exp_dmgs[date]
                    exp_dmg.assign_centroids(hazard)
                    exposure.assign_centroids(hazard)
                    exp_dmg = compute_affected_area_from_dmgs(exp_dmg, exposure)
                    area_temp = exp_dmg.gdf.groupby(['lat','lon'])['field_area'].mean().rename(date)

                #get damages (loss in percent) dataframe at date
                df_sel =  shv_dmg.loc[shv_dmg.Schadendatum==date]

                #get damages in ha by multiplying with field area
                #get all damage values and set coordinates as index
                temp = df_sel.set_index(['lat','lon'])['Ertragsverlust (%)'].rename(date)


                if exposure and hazard:
                   #get total area lost (multiply harvest loss in % with area exposed)
                    col=[]
                    for i in temp.index:
                        col.append(area_temp.loc[area_temp.index==i].values[0])
                    temp_ha=temp.copy(deep=True)*np.array(col)/100

                #count number of damage claims at each coordinate for this date (if larger than 1 this indicates a coordinate related to community centroid rather than a specific field)
                nfields_temp=temp.groupby(['lat','lon']).count()

                #If several damages for same coordinate, then average them together
                if len(temp)!=len(np.unique(temp.index)):
                    temp = temp.groupby(['lat','lon']).mean()
                    if exposure and hazard:
                        temp_ha = temp_ha.groupby(['lat','lon']).sum().rename(date)

                #add the damages of 'date' to the impact matrix
                imp_df=pd.concat([imp_df,temp],axis=1)
                nfields_df=pd.concat([nfields_df,nfields_temp],axis=1)

                if exposure and hazard:
                    area_df=pd.concat([area_df,area_temp],axis=1)
                    ha_loss_df=pd.concat([ha_loss_df,temp_ha],axis=1)



            imp_df = imp_df.fillna(0)
            imp_df=imp_df.drop(columns='count')
            nfields_df = nfields_df.fillna(0)
            nfields_df=nfields_df.drop(columns='count')

            imp_mat = sparse.csr_matrix(imp_df.T.values)
            nfields_mat = sparse.csr_matrix(nfields_df.T.values)

            if exposure and hazard:
                area_df = area_df.fillna(0)
                area_df=area_df.drop(columns='count')
                area_mat = sparse.csr_matrix(area_df.T.values)
                ha_loss_df = ha_loss_df.fillna(0)
                ha_loss_df = ha_loss_df.drop(columns='count')
                ha_loss_mat = sparse.csr_matrix(ha_loss_df.T.values)
            else:
                ha_loss_mat=sparse.csr_matrix(imp_mat.shape())
                area_mat=sparse.csr_matrix(imp_mat.shape())

            return imp_mat, nfields_mat, ha_loss_mat, area_mat, coord_exp

    def get_impact_mats_for_imp(shv_dmg,exposure,hazard):

            #Group by ID and create copys for output
            coord_df = shv_dmg.groupby('ID').first()[['lat','lon']]

            #Initialize impact matrix as dataframe
            imp_df = pd.DataFrame(index=coord_df.index)
            #Initialize matrix to store exposed asset values
            aff_value_df = pd.DataFrame(index=coord_df.index)

            #make coordinate array from array of tuples
            coord_exp=coord_df[['lat','lon']].values

            #get list of dates
            dates=list(shv_dmg.groupby('Schadendatum').groups.keys())

            #get damage exposure dict (to get field area)
            exp_dmgs=get_damage_dict(shv_dmg, return_type = 'exp_dict')


            for date in dates:

                #get dataframe with field area for each exposure point
                exp_dmg = exp_dmgs[date]
                exp_dmg.assign_centroids(hazard)
                exposure.assign_centroids(hazard)
                exp_dmg = compute_affected_area_from_dmgs(exp_dmg, exposure)
                area_temp = exp_dmg.gdf.groupby(['ID'])['field_area'].mean().rename(date)

                #get damages (loss in percent) dataframe at date
                df_sel =  shv_dmg.loc[shv_dmg.Schadendatum==date]
                #get all damage values and set ID as index
                temp = df_sel.set_index(['ID'])['Ertragsverlust (%)'].rename(date)

                #get total area lost (multiply harvest loss in % with area exposed)
                col=[]
                for i in temp.index:
                        col.append(area_temp.loc[area_temp.index==i].values[0])
                temp_ha=temp.copy(deep=True)*np.array(col)/100

                #add the damages of 'date' to the impact matrix
                imp_df=pd.concat([imp_df,temp_ha],axis=1)

                #add affected area to impact matrix for affected area
                aff_value_df=pd.concat([aff_value_df,area_temp],axis=1)


            imp_df = imp_df.fillna(0)
            imp_mat = sparse.csr_matrix(imp_df.T.values)

            #create matrix of exposed asset values
            aff_value_df = aff_value_df.fillna(0)
            aff_mat = sparse.csr_matrix(aff_value_df.T.values)

            return imp_mat, aff_mat, coord_exp

    #### end helper functions ####

    if (return_type in ['imp','imp_dict']) and (not exposure):
        raise ValueError('For return_type {} an exposure need to be passed.'.format(return_type))

    #adjust croptype names if necessary
    croptype = adjust_croptype(croptype)

    #read damage data
    shv_dmg = read_damage_dataframe(excel_file,croptype)
    dates=shv_dmg.groupby('Schadendatum').groups.keys()

    if return_type in ['imp','imp_dict']:

        #get impact matrices
        if return_type == 'imp_dict':
            imp_mat, nfields_mat, ha_loss_mat, area_mat, coord_exp = \
               get_impact_mats_for_imp_dict(shv_dmg,exposure,hazard)
        elif return_type ==  'imp':
            imp_mat, aff_mat, coord_exp = \
                    get_impact_mats_for_imp(shv_dmg,exposure,hazard)

        #Date and frequency
        ord_dates = np.array([d.toordinal() for d in dates])
        ev_names = np.array([d.strftime('ev_%Y-%m-%d') for d in dates])
        n_years = np.floor((max(ord_dates)-min(ord_dates))/365)
        frequency = np.ones(imp_mat.shape[0])
        event_id = np.arange(imp_mat.shape[0])+1

        #Others
        tot_value = np.nan

        #Damages
        imp_at_event = imp_mat.sum(axis=1).getA1()
        imp_eai_exp = imp_mat.sum(axis=0).getA1()/n_years
        imp_aai_agg = imp_mat.sum()/n_years

        if return_type == 'imp_dict':
            nfields_at_event = nfields_mat.sum(axis=1).getA1()
            nfields_eai_exp = nfields_mat.sum(axis=0).getA1()/n_years
            nfields_aai_agg = nfields_mat.sum()/n_years

            #aggregated damages for area measures
            ha_loss_at_event = ha_loss_mat.sum(axis=1).getA1()
            ha_loss_eai_exp = ha_loss_mat.sum(axis=0).getA1()/n_years
            ha_loss_aai_agg = ha_loss_mat.sum()/n_years

            area_at_event = area_mat.sum(axis=1).getA1()
            area_eai_exp = area_mat.sum(axis=0).getA1()/n_years
            area_aai_agg = area_mat.sum()/n_years



            #initialize impact objects
            imp_out = Impact(
                imp_mat = imp_mat,
                coord_exp = coord_exp,
                crs = 'EPSG:4326',
                date = ord_dates,
                event_name = ev_names,
                frequency = frequency,
                event_id = event_id,
                at_event = imp_at_event,
                eai_exp = imp_eai_exp,
                aai_agg= imp_aai_agg,
                unit ='%',
                tot_value=tot_value)

            nfields_out = Impact(
                imp_mat = nfields_mat,
                coord_exp = coord_exp,
                crs = 'EPSG:4326',
                date = ord_dates,
                event_name = ev_names,
                frequency = frequency,
                event_id = event_id,
                at_event = nfields_at_event,
                eai_exp = nfields_eai_exp,
                aai_agg = nfields_aai_agg,
                unit='',
                tot_value=tot_value)

            ha_loss_out = Impact(
                    imp_mat = ha_loss_mat,
                    coord_exp = coord_exp,
                    crs = 'EPSG:4326',
                    date = ord_dates,
                    event_name = ev_names,
                    frequency = frequency,
                    event_id = event_id,
                    at_event = ha_loss_at_event,
                    eai_exp = ha_loss_eai_exp,
                    aai_agg= ha_loss_aai_agg,
                    unit='ha',
                    tot_value=tot_value)

            area_out = Impact(
                    imp_mat = area_mat,
                    coord_exp = coord_exp,
                    crs = 'EPSG:4326',
                    date = ord_dates,
                    event_name = ev_names,
                    frequency = frequency,
                    event_id = event_id,
                    at_event = area_at_event,
                    eai_exp = area_eai_exp,
                    aai_agg= area_aai_agg,
                    unit='ha',
                    tot_value=tot_value)

            #create output dict
            imp={'loss (%)': imp_out, 'fields affected': nfields_out, 'loss (ha)': ha_loss_out, 'area affected (ha)': area_out}

        elif return_type == 'imp':
                imp = Impact(
                    imp_mat = imp_mat,
                    coord_exp = coord_exp,
                    crs = 'EPSG:4326',
                    date = ord_dates,
                    event_name = ev_names,
                    frequency = frequency,
                    event_id = event_id,
                    at_event = imp_at_event,
                    eai_exp = imp_eai_exp,
                    aai_agg= imp_aai_agg,
                    unit ='ha',
                    tot_value=tot_value)
                #ad matrix of affected area
                imp.aff_mat=aff_mat

        return imp

    elif return_type in ['exp_dict', 'gdf_dict']:

        out_dict=get_damage_dict(shv_dmg, return_type = return_type)

        return out_dict

#%% Regridding
def add_zero_values(merged):
    """ Add for each grid cell of output grid, specified by 'index_right'
    in merged, rows with intensity zero such that the total number
    values for the aggregation is equal to n=gridsize_output**2/gridsize_input**2


    Parameters
    ----------
    merged : geopandas.GeoDataframe
        merged dataframe of input data and output grid, with m values per grid cell

    Returns
    -------
    new_merged : geopandas.GEoDataframe
        merged datafarame with n-m zero values added for each grid cell
    """
    dates = []
    indices_right = []
    x = []
    y = []
    intensity = []

    nmax=np.nanmax(merged.groupby(['date','index_right']).size())
    print(nmax)
    for (date, index), count in merged.groupby(['date','index_right']).size().loc[merged.groupby(['date','index_right']).size()<16].iteritems():

        dates = dates + list(np.ones(nmax-count)*date)
        indices_right = indices_right + list(np.ones(nmax-count)*index)
        intensity = intensity + list(np.zeros(nmax-count))

        #geometry
        geometry = merged.loc[merged['index_right'] == index]['geometry']
        x = x + list(np.ones(nmax-count)*geometry.x.values[0])
        y = y + list(np.ones(nmax-count)*geometry.y.values[0])

    gdf=gpd.GeoDataFrame({'intensity': intensity,
                        'date': dates, 'index_right': indices_right},
                        geometry = gpd.points_from_xy(x,y), crs=merged.crs)

    new_merged = pd.concat([merged, gdf])

    return new_merged

def aggregate_hazard(hazard_in, original_grid_epsg = 2056, extent_new = None, cell_size_new = 2000, 
                     projection_new_epsg = 4326, aggfunc = 'max', treat_zeros_as_nans = True, 
                     return_xr=False,dates_xr=None):
    """Aggregating a hazard object to a coarser grid


    Parameters
    ----------
    hazard_in : climada.hazard
        Climada hazard to aggregate to coarser grid
    original_grid_epsg : int, optional
        EPSG number of the original coordinate reference sytem
        of the hazard. The default is 2056.
    extent_new : list of ints, optional if original_grid_epsg is 2056.
        Extent of the new grid (xmin,ymin,xmax,ymax). The default is None.
    cell_size_new : float, optional
        Cell size of the new grid (in units of the original CRS). The default is 2000.
    projection_new_epsg : int, optional
        projection of the centroid coordinates of the new hazard. The default is 4326.
    aggfunc : str, optional
        Function to be used for aggregation. The default is 'mean'.
    treat_zeros_as_nans : boolean, optional
        If True, treat zero values as nans and neglect for aggregation
        If False, treat zero values as zeros and include for aggregation (time consuming).
        The default is True.
    return_xr: boolean, optional
        If True, return the aggregated hazard also as xarray.Dataset
        If False, only return hazard

    Raises
    ------
    ValueError
        If aggregation function is 'max' and treat_zeros_as_nans is False,
        an error is raised to avoid unnecessary time consuming computation.

    Returns
    -------
    hazard_out : climada.hazazrd
        Climada hazard on an aggregated grid

    """


    original_grid_epsg = int(original_grid_epsg)
    original_crs = CRS.from_epsg(original_grid_epsg)
    projection_new_epsg = int(projection_new_epsg)
    proj_new_crs=CRS.from_epsg(projection_new_epsg)
    
    intensities=[] #list to store new intensity values
    xr_out_list=[] #list to store xarray; only needed if return_xr=True
    if dates_xr is None:
        dates_xr=hazard_in.date
        
    if aggfunc == 'max' and treat_zeros_as_nans == False:
        raise ValueError('Set treat_zeros_as_nans for aggfunc = "max" to avoid uneccesary slow down of the code')

    if treat_zeros_as_nans == False:
        warnings.warn('treat_zeros_as_nans = False is currently very slow.')

    # create empty output grid
    print(f'create empty grid with extent {extent_new}, cell size: {cell_size_new}')
    cell, crs, extent = create_empty_grid(epsg = original_grid_epsg, 
                                             cell_size = cell_size_new, extent = extent_new)

    # create a geodataframe with all nonzero values of the hazard intensity over all dates
    # project grid back to original data if not the same
    gdf = gdf_from_hazard(hazard_in)

    # reproject geometry to original grid if hazard grid and original grid are not of identical projection
    if hazard_in.centroids.crs.to_epsg() != original_grid_epsg:
            gdf = gdf.to_crs(crs=original_crs)

    print('Merge hazard with output grid...')
    # merge hazard intensity data with output grid
    merged = gpd.sjoin(gdf, cell, how='left')

    # if missing values have to be treated as zeros add required columns to the dataframe
    if treat_zeros_as_nans == False:
        merged = add_zero_values(merged)

    print(f'Dissolve hazard intensity in output grid using the following aggregation function: {aggfunc}')
    # dissolve orignal intensity data in the new grid using the user specified aggfunc
    dissolve = merged.dissolve(by=['event_id','index_right'], aggfunc={
                "intensity": aggfunc})

    # loop over events and create sparse matrix of hazard intensities for the new grid
    year_now = 0 # counter to catch year changes
    if hazard_in.event_id.shape != hazard_in.date.shape:
        sys.exit("Number of events and dates in hazard are not equal. Abort hazard aggregation.")
    for event, date in zip(hazard_in.event_id, hazard_in.date):
        # print year first time it appears
        year = pd.Timestamp.fromordinal(date).year
        if year > year_now:
            year_now = year
            print(f"Aggregating hazard data from {year}.")

        # since gdf is produced from hazard_in it should always contain the event_id
        # instead, check if event still exists after dissolve (events with zero intensity
        # in the domain of extent_new are lost through the aggregation)
        if event in dissolve.index.get_level_values(0):

            #create deep copy of output grid
            cell_now=cell.copy(deep = True)

            #select the event subset of the dissolved gdf
            dissolve_now=dissolve.loc[event, :] # multi-index, thus two indices given

            #fill output cell with new intensities
            cell_now.loc[dissolve_now.index, 'intensity'] = dissolve_now['intensity'].values

            #create intensity sparse matrix
            cell_now=cell_now.fillna(0)
            intensities.append(sparse.csr_matrix(cell_now['intensity'].values))

            if return_xr == True:
                
                if date in dates_xr:
                    print(pd.Timestamp.fromordinal(date))
                    #create deep copy of data
                    cell_xr=cell_now.copy(deep = True)
                    cell_xr['geometry'] = cell_now.geometry.centroid
                    cell_xr["chx"] = cell_xr.geometry.x
                    cell_xr["chy"] = cell_xr.geometry.y
                    cell_xr = cell_xr.round({'chx': 0, 'chy': 0})
                    cell_xr['geometry'] = gpd.points_from_xy(cell_xr.chx, cell_xr.chy)
                 
                    #get lat lon values
                    geometry_latlon=cell_xr['geometry'].to_crs(crs=CRS.from_epsg(int(4326)))
                    cell_xr["lon"]=geometry_latlon.geometry.x
                    cell_xr["lat"]=geometry_latlon.geometry.y
            
                    cell_multiindex = cell_xr.drop(columns='geometry').set_index(['chy','chx'])
                    cell_xarray=cell_multiindex.to_xarray().set_coords(('chy','chx'))
                    cell_xarray = cell_xarray.expand_dims(time=[pd.Timestamp.fromordinal(date)])
                    xr_out_list.append(cell_xarray)
        else:
            intensities.append(sparse.csr_matrix(np.zeros(len(cell))))
                
    #stack sparse matrices together
    intensities_all = sparse.vstack(intensities)

    #compute hazard centroids
    cell.geometry = cell.geometry.centroid
    
    # if projection of the new centroids is not the same as coordinate reference 
    # system of the original grid, project centroids to new crs
    if projection_new_epsg != original_grid_epsg:
        cell = cell.to_crs(crs=proj_new_crs)
    centroids = Centroids(lat=cell.geometry.y.to_numpy(copy=True),
                          lon=cell.geometry.x.to_numpy(copy=True),
                          geometry=cell.geometry)
    #centroids=Centroids.from_geodataframe(cell)

    # get new hazard with aggregated intensity and new centroids
    hazard_out = copy.deepcopy(hazard_in)
    hazard_out.centroids = centroids
    hazard_out.intensity = intensities_all

    # adjust shape of hazard.fraction (to make hazard.check() pass) in case it contains no data
    if hazard_in.fraction.data.shape[0] == 0:
        hazard_out.fraction = sparse.csr_matrix(np.zeros(hazard_out.intensity.shape))
    else:
        print("Shape of aggregated hazard's intensity and fraction disagree. Hazard.check() will fail.")

    if return_xr==True:
        xr_out = xr.concat(xr_out_list, dim = 'time')
    else:
        xr_out = None
    return hazard_out, xr_out

def gdf_from_hazard(hazard):
    """Create GeoDataFrame from hazard object with columns 'intensity', 'date', 'event_id', and hazard centroids as geometry

    Parameters
    ----------
    hazard : climada.hazard
        Climada hazard
    Returns
    -------
    gdf : pandas.GeoDataFrame
        geodataframe with hazard information

    """

    #select date and centroid indices
    i_dates,i_centr=hazard.intensity.nonzero()

    #create gdf of geometry
    geometry = gpd.points_from_xy(hazard.centroids.lon[i_centr], hazard.centroids.lat[i_centr])


    gdf=gpd.GeoDataFrame({'intensity': np.squeeze(np.asarray(hazard.intensity[i_dates,i_centr])),
                              'date': hazard.date[i_dates],
                              'event_id': hazard.event_id[i_dates]},
                             geometry = geometry, crs=hazard.centroids.crs)

    return gdf


def create_empty_grid(epsg=2056, cell_size = 1000, extent = None):
   """ get an empty grid as geopandas dataframe based on passed epsg, resolution, and grid size

   Parameters
   ----------
       epsg: int
           String denoting the required coordinate system (EPSG format)
       cell_size: float or int
           size of the grid cells either in the unit of the specified coordinate system
       extent: list or tuple
           extent of the grid (xmin, ymin, xmax, ymax)
   Returns
   ----------
       cell: geopandas.GeoDataFrame
           a GeoDataFrame with each grid cell of the 1x1km regular grid represented as polygon.
   """

   # specify extent of the grid (here: LV95 bounds)
   if epsg == 2056:
       xmin, ymin, xmax, ymax = (2255000, 840000, 2964000, 1479000)
       if extent is not None:
           warnings.warn('Extent of grid EPSG 2056 is predefined. Argument extent is ignored.')
   else:
       xmin, ymin, xmax, ymax = extent

   # create the cells in a loop
   grid_cells = []
   for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
       for y0 in np.arange(ymin, ymax+cell_size, cell_size):
           #bounds
           x1 = x0+cell_size
           y1 = y0+cell_size
           grid_cells.append(sg.box(x0, y0, x1, y1)  )


   crs=CRS.from_epsg(int(epsg))

   #create geopandas.GeoDataFrame
   cell = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                                    crs = crs) #"EPSG:2056")

   return cell, crs, extent
#%%verification
def split_events_non_events(df, variable):
    """
    

    Parameters
    ----------
    df : pandas.Dataframe
        At centroid dataframe
    variable : str
        hazard variable.

    Returns
    -------
    events : pandas.Dataframe
        dataframe including all events (prediction or damages).
    no_damages : pandas.Dataframe
        dataframe including all centroids without damage

    """
    
    at_centroid_df=df.copy()
    at_centroid_df['n_dmgs']=at_centroid_df['n_dmgs'].fillna(0)
    #at_centroid_df['dmg_val']=at_centroid_df['dmg_val'].fillna(0)
    events=at_centroid_df.loc[(at_centroid_df[variable]!=0) | (at_centroid_df['n_dmgs']!=0)] 
    no_damages = at_centroid_df.loc[at_centroid_df['n_dmgs']==0]
    return events,no_damages

def compute_verification_stats(at_centroid_data, variable = 'MESHS',exposure_thresh=1):
    """
    

    Parameters
    ----------
    at_centroid_data : pandas.Dataframe
        Dataframe with verification data at centroids
    variable : str, optional
       hazard variable. The default is 'MESHS'.
    exposure_thresh : int, optional
       minimum number of fields for a grid point to count as exposure.
       The default is 1.

    Returns
    -------
    df : pandas.DataFrame
        dataframe with hazard values as indices and verification metrics as columns
    npred : int
        total number of predictions

    """
    
    #select only the crop area for verification (damages outside are disregarded)
    at_centroid_data['n_exp']=at_centroid_data['n_exp'].fillna(0)
    data=at_centroid_data[at_centroid_data['n_exp'] >= exposure_thresh]
 
    #make subsets events (damage or prediction), and no_damages (no damage, either false alarms or non events)
    events,no_damages=split_events_non_events(data,variable)

    #HITS: number of predictions with at least one observed damage (per centroid)
    A=events.groupby(variable).agg({'n_dmgs': lambda x: x[x != 0].count()})
    #number of prediction with no observed damage (per centroid)
    # FALSE ALARMS: 
    B=events.groupby(variable).agg({'n_dmgs': lambda x: x[x == 0].count()})
    # non
    D=no_damages.groupby(variable)['n_dmgs'].count()
    
    npred=np.int64(np.sum(A[1::]+B[1::]).values[0])
    Asum=np.ones(len(A.index[1::]))*np.nan
    Bsum=np.ones(len(A.index[1::]))*np.nan
    Csum=np.ones(len(A.index[1::]))*np.nan
    Dsum=np.ones(len(A.index[1::]))*np.nan
    
    index=range(len(A.index[1::]))
    for i in index:
        #indexing of original dataframe (to exclude values at haz intensity 0)
        l=i+1
        # count number hits above threshold
        Asum[i]=A['n_dmgs'].values[l::].sum()
        # count number false alarms above threshold
        Bsum[i]=B['n_dmgs'].values[l::].sum()
        # count number of misses
        Csum[i]=A['n_dmgs'].values[0:l].sum()
        # count number of non_events
        Dsum[i]=D.values[0:l].sum()

    
    #check here: https://www.cawcr.gov.au/projects/verification/#Types_of_forecasts_and_verifications
    FAR=Bsum/(Asum+Bsum)
    POD=Asum/(Asum+Csum)
    CSI=Asum/(Asum+Bsum+Csum)
    NER=Dsum/(Bsum+Dsum)
    POFD=Bsum/(Bsum+Dsum)
    PC=(Asum+Dsum)/(Asum+Bsum+Csum+Dsum) #proportion correct
    fyes=(Asum+Csum)/(Bsum+Dsum) #fraction of yes observations on all observations
    SEDI=np.nan#(np.log(POFD)-np.log(POD)+np.log(1-POD)-np.log(1-POFD))/(np.log(POFD)+np.log(POD)+np.log(1-POD)+np.log(1-POFD))
    HSS=2*(Asum*Dsum-Bsum*Csum)/((Asum + Csum)*(Csum + Dsum) + (Asum + Bsum)*(Bsum + Dsum))
    base_rate=(Asum+Csum)/(Asum+Bsum+Csum+Dsum)
    #Hansen Kuyper Discriminant
    HK = POD - POFD
    #Bias
    B=(Asum+Bsum)/(Asum+Csum)
    #total number of centroids with hail predictions
    Pred=Asum+Bsum    
    #total number of centroids with hail damage
    Obs=Asum+Csum
    
    df=pd.DataFrame({'FAR': FAR,
                     '1-FAR': 1-FAR,
                     'POD': POD,
                     'NER': NER,
                     'CSI': CSI,
                     'HSS': HSS,
                     'SEDI': SEDI,
                     'PC': PC,
                     'POFD': POFD,
                     'HK': HK,
                     'B':B,
                     'PRED': Pred,
                     'Fyes': fyes,
                     'OBS': Obs,
                     's': base_rate,
                     'hits': Asum,
                     'false alarms': Bsum,
                     'misses': Csum},
                    index = A.index[1::])
    return df, npred


def get_scores_from_exposure_threshs(at_centroid_df,variable,s,exposure_threshs,scores):
    """
    Get the model skill metrics for different exposure density thresholds from the at_centroid data frame
    
    Parameters
    ----------
    at_centroid_df: pandas.Dataframe
        Dataframe containing the observation and prediction data at centroids
    s : int,float
        hazard intensity threshold (needs to match hazard values in at_centroid_df)
    variable : str
        hazard variable
    exposure_threshs : np.array
        Array of exposure thresholds
    scores : list of strings
        List with scores to be returned for different exposure thresholds

    Returns
    -------
    df : pandas.Dataframe
        Dataframe with exposure thresholds as index and scores as columns

    """
    score_dict={sc: [] for sc in scores}
    for th in exposure_threshs:
        df, npred = compute_verification_stats(at_centroid_df,
                                                variable = variable,
                                                exposure_thresh=th)
        for sc in scores:
            score_dict[sc].append(df.loc[s][sc])

    df=pd.DataFrame(score_dict, index=exposure_threshs)
    return df

########  PLOTTING

def get_data_per_date(at_centroid_data,variable,thresh,measures):    
    """
    get verification statistics per date

    Parameters
    ----------
    at_centroid_data :  andas.Dataframe
        Dataframe containing hazard,exposure and impact data at centroids
    variable : str
        hazard variable
    thresh : int
        exposure threshold (n_thresh)
    measures : list of str
        verification scores


    Returns
    -------
    df_all_dates : pandas.DataFrame
        verification data for each date
    npred_per_date : list
        number of predictions per date
    dates_all : list
        dates

    """
    dates_all=sorted(list(set(at_centroid_data.date)))
    data_per_date=[]
    npred_per_date=[]     
    for date in dates_all:
                df, npred = compute_verification_stats(at_centroid_data.loc[at_centroid_data['date']==date],
                                                variable = variable,
                                                exposure_thresh=thresh)
                #count number of predictions for each date
                npred_per_date.append(npred)
                if variable == 'MESHS':
                    data_per_date.append(df.iloc[0][measures].rename(date))
                if variable == 'E_kin':
                    data_per_date.append(df.loc[100][measures].rename(date))
                if variable == 'POH':
                    if 80 in df.index:
                        data_per_date.append(df.loc[80][measures].rename(date))
                    elif 81 in df.index:
                        data_per_date.append(df.loc[81][measures].rename(date))
    #create pandas Dataframe from data_per_date
    df_all_dates=pd.concat(data_per_date,axis=1)
    return df_all_dates, npred_per_date, dates_all


def plot_measures(df,x,ax,npred=None):
    """
    
    Plot verification measures per hazard value
    
    Parameters
    ----------
    df : pandas.Dataframe
        dataframe with verification scores per hazard value
    x : np.array
        hazard values
    ax : matplotlib.axes
        axes to plot on
    npred : np.array or list, optional
        number of predictions per hazard valur. The default is None.

    Returns
    -------
    ax : matplotlib.axes
        axes with plotted measures.

    """
    
    lw=5
    lw2=3
    
    h1,=ax.plot(x,df.FAR[x],color='r',linewidth=lw2, label = 'FAR')
    #ax.fill_between(x,q05[resolution]['FAR'],q95[resolution]['FAR'],color='r',alpha=0.2)
    
    h2,=ax.plot(x,df.POD[x],color='blue',linewidth=lw2, linestyle = 'solid', label = 'POD')
    #ax.fill_between(x,q05[resolution]['POD'],q95[resolution]['POD'],color='blue',alpha=0.2)

    #h2,=ax.plot(x,df.POFD[x],color='orange',linewidth=lw2, linestyle = 'solid', label = 'POFD')

    #h2,=ax.plot(x,df.PC[x],color='green',linewidth=lw2, linestyle = 'solid', label = 'PC')

    #h3,=ax.plot(x,df.CSI[x],color='k',linewidth=lw2, linestyle = 'dashed', label = 'CSI')
    #ax.fill_between(x,q05[resolution]['CSI'],q95[resolution]['CSI'],color='grey',alpha=0.2)

    h4,=ax.plot(x,df.HSS[x],color='k',linewidth=lw, linestyle = 'dotted', label = 'HSS')
    #ax.fill_between(x,q05[resolution]['HSS'],q95[resolution]['HSS'],color='grey',alpha=0.2)

    #h6,=ax.plot(x,df.HK[x],color='k',linewidth=lw, linestyle = 'dashed', label = 'HK')

    h7,=ax.plot(x,df.CSI[x],color='k',linewidth=lw, linestyle = 'dashed', label = 'CSI')

    if npred is not None:
        h5,=ax.plot(x,df.PRED[x]/npred,color='k',linestyle='dotted',linewidth=lw2, label='Ratio of total predictions')
    #ax.fill_between(x,q05[resolution]['PRED']/npred,q95[resolution]['PRED']/npred,color='grey',alpha=0.2)
    
    ax.set_yticks(np.arange(0,1.1,0.2))
    ax.tick_params(axis='both', labelsize=25)
    return  ax

def plot_measures_per_event(at_centroid_data_MESHS, variable, ax, measures=['FAR','POD','HSS','HK'], thresh=0,fontsize=25):
    """
    plot measures per event

    Parameters
    ----------
    at_centroid_data_MESHS : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.
    measures : TYPE, optional
        DESCRIPTION. The default is ['FAR','POD','HSS','HK'].
    thresh : TYPE, optional
        DESCRIPTION. The default is 0.
    fontsize : TYPE, optional
        DESCRIPTION. The default is 25.

    Returns
    -------
    handles : TYPE
        DESCRIPTION.
    labs : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.
    dates_all : TYPE
        DESCRIPTION.
    npred_per_date : TYPE
        DESCRIPTION.
    all_data : TYPE
        DESCRIPTION.

    """
    
    #define the two resolutions to show
    resolutions=['1km','8km']

    #Plot properties
    markers=['X','o','d','v']
    colors=['red','blue','k','dimgrey']
    linestyles=['solid','solid','solid','solid']
    markersizes=[4,4,4,4]
    markersizes2=[6,6,6,6]
    font_l=fontsize
    matplotlib.rcParams.update({'font.size': font_l, 'axes.labelsize':font_l})

    handles=[]
    labs=[]

    #PLOT MESHS
    #get stats
    npred_per_date={r: [] for r in resolutions}
    all_data={}
    for i,resolution in enumerate(resolutions):
        
        #gat skill metrics per date
        df_all_dates, npreds, dates_all=get_data_per_date(at_centroid_data_MESHS[resolution],variable,thresh,measures)
        all_data[resolution]=df_all_dates

        npred_per_date[resolution]=npreds

        #set x-axis ticks    
        x=np.arange(0.5,len(dates_all),1)
    
        #plot measures 
        for i,measure in enumerate(measures):
            if resolution=='1km':
                h,=ax.plot(x,df_all_dates.loc[measure], marker=markers[i],markersize=markersizes[i],linewidth=0.5, color=colors[i],linestyle=linestyles[i],label=measure)
            else:
                h,=ax.plot(x,df_all_dates.loc[measure], marker=markers[i],markersize=markersizes2[i],markeredgewidth=0.5,linewidth=0.5, color=colors[i],linestyle='dashed',markeredgecolor=colors[i],markerfacecolor='none',label=measure)
            handles.append(h)
            labs.append(measure)

    ax.set_xticks(np.arange(0.5, len(df_all_dates.columns), 1), df_all_dates.columns.strftime('%Y-%m-%d'),rotation=80)
    ax.set_ylim([0,1])

    return handles, labs, ax, dates_all, npred_per_date, all_data

def plot_skill_haz_threshold(at_centroid_data_wheat,at_centroid_data_grapevine,variable,unit, figdir, thresh=0,resolution='1km'):
    """
    

    Parameters
    ----------
    at_centroid_data_wheat : TYPE
        DESCRIPTION.
    at_centroid_data_grapevine : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.
    unit : TYPE
        DESCRIPTION.
    figdir : TYPE
        DESCRIPTION.
    thresh : TYPE, optional
        DESCRIPTION. The default is 0.
    resolution : TYPE, optional
        DESCRIPTION. The default is '1km'.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.

    """
    
     #, '4km', '8km']
    labels=['a)','b)','c)','d)','e)','f)','g)','h)','i)']
    
    font_s=18
    font_l=18
    
    matplotlib.rcParams.update({'font.size': font_l, 'axes.labelsize':font_l})
    
    fig, axes = plt.subplots(figsize=(12.0, 5.0) , nrows=1, ncols=2, sharey='row') 
    axs=axes.flatten()
    for i,ax in enumerate(axs):
        ax.text(0,1.08,labels[i],
                transform=ax.transAxes,fontweight='bold')
      
    plt.subplots_adjust(wspace=0.15, hspace=0.45)
        
    #### SUBPLOT A: field crops 
    ax=axes[0]
    #get stats
    df, npred = compute_verification_stats(at_centroid_data_wheat[resolution],
                                               variable = variable,
                                               exposure_thresh=thresh)
                

    imax=51 #maximum index in df to be plotted
    x=df.index[0:imax]

    #plot lines
    ax=plot_measures(df,x,ax)

    #setup axis
    ax.set_title('field crops',fontweight='bold',fontsize=font_l,pad=40)  
    ax.set_xlim([20,60])
    ax.set_xticks(np.arange(20,61,10))
    ax.text(0.5,1.01,'n = {}'.format(npred),
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=ax.transAxes, fontsize=font_s)
    ax.set_ylim([0,1])

    #plot bias
    col='green'
    ax.axvline(44,color=col,linewidth=1)
    ax.text(44,0.99,'B=1',color=col,ha='left',va='top')
    ax.axvline(30,color=col,linewidth=1,linestyle='dashed')
    ax.text(30,0.99,'B=2',color=col,ha='left',va='top')
    ax.set_xlabel(unit,fontsize=font_l)

    #set up secondary axis
    ax2=ax.twinx()
    ax2.bar(x,df['PRED'][x],width=1,color='lightgrey',zorder=0)
    ax2.set_yticks(np.arange(0,16000,4000),[],color='grey')
    ax2.yaxis.label.set_color('lightgrey')
    ax.spines['right'].set_color('grey')
    ax2.spines['right'].set_color('grey')
    ax2.tick_params(colors='grey')  # 'both' refers to minor and major axes

    #set ax1 above ax2
    ax.set_zorder(ax2.get_zorder()+1)
    ax.set_frame_on(False)

    ##### SUBPLOT B: grapevine
    ax=axes[1] #fig.add_subplot(3,3,i+4)
    #get stats
    df, npred = compute_verification_stats(at_centroid_data_grapevine[resolution],
                                               variable = variable,
                                               exposure_thresh=thresh)
    imax=51 #maximum index in df to be plotted
    x=df.index[0:imax]
    #plot lines
    ax=plot_measures(df,x,ax)

    #setup and customize axis
    ax.set_title('grapevine',fontweight='bold',fontsize=font_l,pad=40)  
    ax.set_xlim([20,60])
    ax.set_xticks(np.arange(20,61,10))
    ax.text(0.5,1.01,'n = {}'.format(npred),
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=ax.transAxes, fontsize=font_s)
    
    #plot vertical lines for bias
    ax.axvline(46,color=col,linewidth=1)
    ax.text(46,0.99,'B=1',color=col,ha='left',va='top')
    ax.axvline(34,color=col,linewidth=1,linestyle='dashed')
    ax.text(34,0.99,'B=2',color=col,ha='left',va='top')

    ax.set_ylim([0,1])
    ax.set_xlabel(unit,fontsize=font_l)

    #plot legend 
    ax.legend(ncol=6,loc='center',bbox_to_anchor=(0 , -0.3))

    #plot secondary vertical axis
    ax2=ax.twinx()
    ax2.bar(x,df['PRED'][x],width=1,color='lightgrey',zorder=0)
    ax2.set_yticks(np.arange(0,16000,4000),np.arange(0,16000,4000),color='lightgrey')
    ax2.yaxis.label.set_color('lightgrey')
    ax.spines['right'].set_color('grey')
    ax2.spines['right'].set_color('grey')
    ax2.tick_params(colors='grey')  # 'both' refers to minor and major axes
   
    #make sure ax1 is above ax2
    ax.set_zorder(ax2.get_zorder()+1)
    ax.set_frame_on(False)
            
    return fig, ax

def plot_CH_map(ncols=2,nrows=2, proj=ccrs.epsg(2056),stamen_map=True, 
                figsize=(15,10),lakes=True,cantons=False,pads=0.1,edgecolor='black'):
    
    """Plot a Figure of a map of Switzerland optionally using the Stamen terrain background map (http://maps.stamen.com/#toner/12/37.7706/-122.3782)
    Parameters
    ----------  
    
    stamen_map: boolean
            if True, plot Stamen terrain map, otherwise no background (default is True)
    figsize: tuple
            the figsize keyword argument for matplotlib.figure

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot

    """
    
    ch_shp_path = 'C:/Users/F80840370/projects/scClim/climada/data/TLM/'

    def rect_from_bound(xmin, xmax, ymin, ymax):
        """Returns list of (x,y)'s for a rectangle"""
        xs = [xmax, xmin, xmin, xmax, xmax]
        ys = [ymax, ymax, ymin, ymin, ymax]
        return [(x, y) for x, y in zip(xs, ys)]
    
    # request data for use by geopandas
    resolution = '10m'
    category = 'cultural'
    name = 'admin_0_countries'
    
    shpfilename = shapereader.natural_earth(resolution, category, name)
    
    df = geopandas.read_file(shpfilename)
    
    # get geometry of a country
    poly = [df.loc[df['ADMIN'] == 'Switzerland']['geometry'].values[0]]
    
    stamen_terrain = cimgt.Stamen('terrain-background')
    
    # projections that involved
    if stamen_map == True:
        st_proj = stamen_terrain.crs  #projection used by Stamen images
        ll_proj = ccrs.PlateCarree()
    else:
        ll_proj = ccrs.AlbersEqualArea(central_longitude=8) #proj #ccrs.PlateCarree() #ccrs.AlbersEqualArea(central_longitude=8) #ccrs.PlateCarree() #  #CRS for raw long/lat
        st_proj = ll_proj

    if isinstance(pads,list):
        if len(pads)==4:
            exts = [poly[0].bounds[0] - pads[0], poly[0].bounds[2] + pads[1], 
                    poly[0].bounds[1] - pads[2], poly[0].bounds[3] + pads[3]];
        else:
            raise ValueError('Argument pads needs to be sequence of length 4 or float.')
    else:
        pad1=pads
        exts = [poly[0].bounds[0] - pad1, poly[0].bounds[2] + pad1, 
                poly[0].bounds[1] - pad1, poly[0].bounds[3] + pad1];
                
    # create fig and axes using intended projection
    fig = plt.figure(figsize = figsize)
    
    axes=[]

    subplot=0
    
    for row in range(nrows):
        for col in range(ncols):
            subplot+=1
            ax = fig.add_subplot(nrows,ncols,subplot, projection=st_proj) #st_proj)
            ax.add_geometries(poly, crs=ll_proj, facecolor='none', edgecolor=edgecolor)
            
            #pad1 = .1  #padding, degrees unit
            #exts = [poly[0].bounds[0] - pad1, poly[0].bounds[2] + pad1, poly[0].bounds[1] - pad1, poly[0].bounds[3] + pad1];
            ax.set_extent(exts, crs=proj) 
            ax.spines['geo'].set(edgecolor=edgecolor)
    
            # make a mask polygon by polygon's difference operation
            # base polygon is a rectangle, another polygon is simplified switzerland
            msk = Polygon(rect_from_bound(*exts)).difference( poly[0].simplify(0.01) )
            msk_stm  = st_proj.project_geometry (msk, ll_proj)  # project geometry to the projection used by stamen
            
            if stamen_map == True:
                # get and plot Stamen images
                ax.add_image(stamen_terrain, 8) # this requests image, and plot
            
                # plot the mask using semi-transparency (alpha=0.65) on the masked-out portion
                ax.add_geometries( msk_stm, st_proj, zorder=12, facecolor='white', edgecolor='none', alpha=0.65)
        
            #ax.gridlines(draw_labels = False)
            #ax.add_feature(cf.LAKES)
            ax.add_feature(cf.BORDERS)
           
            if lakes == True:
                reader2 = shapereader.Reader("%s/Hydrography/swissTLMRegio_Lake.shp"%ch_shp_path)
                geometry = reader2.geometries()
                geometry = np.array([g for g in geometry])
                lakesize = np.array([a.area for a in reader2.geometries()])
                geometry = geometry[lakesize>2e7]
                shape_feature2 = cf.ShapelyFeature(geometry,
                                                   ccrs.epsg(2056), edgecolor='none',facecolor = "lightblue")   
                ax.add_feature(shape_feature2,zorder=0.5)#0.5,2
            if cantons == True:
                reader = shapereader.Reader("%s/swissTLMRegio_KANTONSGEBIET_LV95.shp"%ch_shp_path)
                sel_cantons = [place for place in reader.records() if place.attributes["OBJEKTART"]=='Kanton' and place.attributes["ICC"]=='CH']
                for sel_canton in sel_cantons:
                    shape_feature = cf.ShapelyFeature([sel_canton.geometry],
                                                           ccrs.epsg(2056), edgecolor = 'dimgrey',facecolor = "none")   
                    ax.add_feature(shape_feature,linewidth=0.5,zorder=0.5)
            
            axes.append(ax)
            

    return fig, axes,exts

### Plot performance diagram
def get_uneven_cmap(cmap_rb,uneven_levels):
    """
    get an uneven colormap

    Parameters
    ----------
    cmap_rb : TYPE
        DESCRIPTION.
    uneven_levels : TYPE
        DESCRIPTION.

    Returns
    -------
    cmap : TYPE
        DESCRIPTION.
    norm : TYPE
        DESCRIPTION.

    """
    import matplotlib.colors as mcolors    
    
    cols = cmap_rb(np.linspace(0, 1, len(uneven_levels) - 1))
    cmap, norm = mcolors.from_levels_and_colors(uneven_levels, cols)
    
    return cmap, norm

def plot_performance_diagram(ax):
    """
    plot a performance diagram
 
    parameters
    ----------
    ax : matplotlib.axes
        axes to plot on
 
    returns
    -------
    ax : matplotlib.axes
        axes with performance  diagram
    h : handle
        contour handle for csi
    hl : handle
        contour handle for bias
 
    """
    fars=np.arange(0,1.001,0.001)
    pods=np.arange(0,1.001,0.001)
    x1,y1=np.meshgrid(fars,pods)
    B=y1/x1
    CSI=1/(1/x1+1/y1-1)

    B_levs=[0.25,0.5,0.75,1,1.5,2,4]
    CSI_levs=np.arange(0,1.1,0.1)
    hl=ax.contour(x1,y1,B,levels=B_levs,colors='k',alpha=0.8,linestyles='dashed')
    h=ax.contourf(x1,y1,CSI,cmap=cc.cm.CET_CBL1_r,levels=CSI_levs)
    ax.set_ylabel('POD')
    return ax,h,hl

def plot_roc_diagram(ax):
    """
    Plot a ROC diagram

    Parameters
    ----------
    ax : matplotlib.axes
        axes to plot on

    Returns
    -------
    ax : matplotlib.axes
        axes with ROC diagram
    h : handle
        contour handle for HK

    """
    
    pofds=np.arange(0,0.5,0.001)
    pods=np.arange(0,1.001,0.001)
    x2,y2=np.meshgrid(pofds,pods)
    hk=y2-x2
    HK_levs=np.arange(-1,1.1,0.2)

    cmap,norm=get_uneven_cmap(cc.cm.bwy_r, HK_levs)
    h=ax.contourf(x2,y2,hk,cmap=cmap,norm=norm,levels=HK_levs)

    return ax,h

def plot_diagram_points(at_centroid_data_var,variable,opt_resolutions,exp_thresh,opt_threshs,unit,score,colors,markers,fontsize,ax,labeled=True,label_resolutions=True,loc_thresh_labels='l'):
    """
    Plot points in Performance diagram

    Parameters
    ----------
    at_centroid_data_var : TYPE
        DESCRIPTION.
    variable : TYPE
        DESCRIPTION.
    opt_resolutions : TYPE
        DESCRIPTION.
    exp_thresh : TYPE
        DESCRIPTION.
    opt_threshs : TYPE
        DESCRIPTION.
    unit : TYPE
        DESCRIPTION.
    score : TYPE
        DESCRIPTION.
    colors : TYPE
        DESCRIPTION.
    markers : TYPE
        DESCRIPTION.
    fontsize : TYPE
        DESCRIPTION.
    ax : TYPE
        DESCRIPTION.
    labeled : TYPE, optional
        DESCRIPTION. The default is True.
    label_resolutions : TYPE, optional
        DESCRIPTION. The default is True.
    loc_thresh_labels : TYPE, optional
        DESCRIPTION. The default is 'l'.

    Returns
    -------
    ax : TYPE
        DESCRIPTION.

    """
    
        x_threshs={thresh: [] for thresh in opt_threshs}
        y_threshs={thresh: [] for thresh in opt_threshs}
        for i,res in enumerate(opt_resolutions):
        
            df, npred = compute_verification_stats(at_centroid_data_var[res],
                                               variable = variable,
                                               exposure_thresh=exp_thresh)
     
            
            #df_sel=df[(df['B']>0.5) & (df['B']<2) & (df['HK']>0.4) & (df['HSS']>0.3)]
            #h=ax.scatter(df_sel[score],df_sel.POD,color=colors[i],marker=markers[i],label=res)
            
            #plot POD vs Score
            x,y=(df[score][opt_threshs],df.POD.loc[opt_threshs])
            for thresh in opt_threshs:
                x_threshs[thresh].append(x[thresh])
                y_threshs[thresh].append(y[thresh])

            if labeled==True:
                ax.plot(x,y,color='k',linewidth=1)
                ax.scatter(df[score][opt_threshs],df.POD.loc[opt_threshs],color=colors[i],marker=markers[i],s=100,edgecolor='k',label=res,zorder=10)
                x,y=df[score][opt_threshs],df.POD.loc[opt_threshs]  

            else:
                ax.plot(x,y,color='r',linewidth=2,alpha=1)
                #h=ax.scatter(df[score][opt_threshs],df.POD.loc[opt_threshs],color='dimgrey',marker=markers[i],s=100,edgecolor='k',label=res,zorder=10,alpha=0.5)
   
            if labeled==True:
                #Annotate resolutions
                if label_resolutions==True:
                    ax.annotate(f'{res}',(x.values[-1],y.values[-1]),
                                xytext=(x.values[-1]+0.005,y.values[-1]-0.005),
                                xycoords='data',
                               fontsize=fontsize,va='top',bbox=dict(facecolor='white', edgecolor='none',alpha=0.5))   
                    
                #annotate thresholds
                x,y=df[score][opt_threshs],df.POD.loc[opt_threshs]  
                if res=='1km':
                    for k,thresh in enumerate(opt_threshs):
                        if loc_thresh_labels=='l':
                            ax.annotate(f'{thresh}{unit}',(x.values[k],y.values[k]),
                                        xytext=(x.values[k]-0.005,y.values[k]),
                                        xycoords='data',
                                       fontsize=fontsize,va='center',ha='right',bbox=dict(facecolor='white', edgecolor='none',alpha=0.5))   
                        elif loc_thresh_labels=='r':
                            ax.annotate(f'{thresh}{unit}',(x.values[k],y.values[k]),
                                        xytext=(x.values[k]+0.005,y.values[k]),
                                        xycoords='data',
                                       fontsize=fontsize,va='center',ha='left',bbox=dict(facecolor='white', edgecolor='none',alpha=0.5))   
                    
        for thresh in opt_threshs:
            if labeled==True:
                ax.plot(x_threshs[thresh],y_threshs[thresh],color='k',linestyle='dashed',linewidth=0.5)
            else:
                ax.plot(x_threshs[thresh],y_threshs[thresh],color='r',linestyle='dashed',linewidth=2)

        return ax