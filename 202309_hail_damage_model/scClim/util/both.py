# -*- coding: utf-8 -*-
"""
Utility functions for both/all Subprojects
"""
import datetime as dt
import os

import climada.util.coordinates as u_coord
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from climada import CONFIG
from scipy.signal import convolve
from geopy.geocoders import Nominatim, Bing
from geopy.extra.rate_limiter import RateLimiter
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import shapely
from pyproj import CRS



import scClim.hail_climada as fct

def interp2d_na(zs,method='linear'):
    """interpolate 2d array with nan values 

    Args:
        zs (np.array): 2d array
        method (str, optional): interpolation method. Defaults to 'linear'.

    Returns:
        np.array: filled array
    """
    x_indx, y_indx = np.meshgrid(np.arange(0, zs.shape[1]),
                                np.arange(0, zs.shape[0]))

    # mask all invalid values
    zs_masked = np.ma.masked_invalid(zs)

    # retrieve the valid, non-Nan, defined values
    valid_xs = x_indx[~zs_masked.mask]
    valid_ys = y_indx[~zs_masked.mask]
    valid_zs = zs_masked[~zs_masked.mask]

    # generate interpolated array of z-values
    zs_interp = griddata((valid_xs, valid_ys), valid_zs.ravel(),
                                    (x_indx, y_indx), method=method)
    return zs_interp


def geocode_df(df,address_column='address',geoloc_type='Nominatim'):
    
    df_out = df.copy(deep=True)
    if any(np.isin(['location','longitude','latitude'],df_out.columns)):
        raise ValueError('columns location/longitude/latitude already exists')
    #Creating an instance of Nominatim Class
    if geoloc_type == 'Nominatim':
        print(f'Estimated time: {df_out.shape[0]/60:.1f} min')
        geolocator = Nominatim(user_agent="my_request")
        #applying the rate limiter wrapper
        geocode = RateLimiter(geolocator.geocode, min_delay_seconds=1)
    elif geoloc_type == 'Bing':
        if df_out.shape[0]>20000:
            #so far used: 16k GVL, 16k AGV,
            raise TimeoutError('Careful, only 1250000 requests for Bing allowed')
        geolocator = Bing(api_key='***')
        geocode = RateLimiter(geolocator.geocode, min_delay_seconds=0.0)

    #Applying the method to pandas DataFrame
    df_out['location'] = df_out[address_column].apply(geocode)
    df_out['latitude'] = df_out['location'].apply(lambda x: x.latitude if x else None)
    df_out['longitude'] = df_out['location'].apply(lambda x: x.longitude if x else None)
    return df_out

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


def assign_centroids_imp(imp,hazard,distance='euclidean',
                        threshold=u_coord.NEAREST_NEIGHBOR_THRESHOLD):
    """
    same as assign_centroids_gdf, but for impact object 

    """
    #make sure the imp crs is epsg:4326, as it is required by the u_coord methods
    assert(imp.crs == 'EPSG:4326')
    coord_df = pd.DataFrame(imp.coord_exp,columns=['latitude', 'longitude'])
    gdf = gpd.GeoDataFrame(coord_df,geometry = gpd.points_from_xy(imp.coord_exp[:,1],imp.coord_exp[:,0],crs = 'EPSG:4326'))
    assign_centroids_gdf(gdf,hazard,distance,threshold)

    setattr(imp,f'centr_{str(hazard.tag.haz_type)}' , gdf['centr_' + hazard.tag.haz_type] )
    

def npy_to_netcdf(npy_file,var):

    #load example netcdf file with correct coords
    nc_filepath = str(CONFIG.mch_nc_example)
    ncfile = xr.open_dataset(nc_filepath)
    
    if npy_file.endswith('.npy'):
        arr = np.flip(np.load(npy_file),axis=[0])
        #convert MESHS from cm to mm
        if var in ['MZC','meshs']:
            arr = arr*10
        # arr = np.load(npy)
        ds_out = xr.Dataset({var: (("chy","chx"),arr)},
                        coords = ncfile.coords).drop('time')   
    elif os.path.isdir(npy_file):
        npy_files = os.listdir(npy_file)
        for file in npy_files:
            timestamp = file.replace('%s_'%var,'').replace('.npy','')
            arr = np.flip(np.load(os.path.join(npy_file,file)),axis=[0])
            #convert MESHS from cm to mm
            if var in ['MZC','meshs']:
                arr = arr*10
            # arr = np.load(npz)
            ds = xr.Dataset({var: (("chy","chx"),arr)},
                            coords = ncfile.coords).isel(time=0)
            ds['time'] =dt.datetime.strptime(timestamp,"%Y%m%d%H%M%S")
            if file == npy_files[0]:
                ds_out = ds
            else:
                ds_out = xr.concat([ds_out,ds],dim='time')  
    else:
        TypeError('npy_file is neither .npy file nor directory')
    return ds_out

def get_possible_hail(date, poh, haz_poh, extent, poh_thresh=80, buffer_km = 4,
                      return_type='oneDay',get_likelihood=False):
    """get all the location where hail is expected with, based on POH

    Parameters
    ----------
    date : datetime
        date in question
    poh : xarray.DataArray
        POH values 
    haz_poh : climada.hazard
        hazard object based on the POH values from the poh DataArray
    extent : list / array
        [lon_min, lon_max, lat_min, lat_max]
    poh_thresh : int
        poh threshold
    buffer_km : int
        added buffer in km
    return_type : str
        'oneDay': returns hazard + DataArray of selected date
        'all' : ignored 'date' and returns DataArray of all timesteps
    get_likelihood : bool
        if True, returns likelyhood of hail, based on POH and buffer_km.
        To be used for MFZrandom only!
        if False, return boolean for possible hail yes/no. Used for buildings

    Returns
    -------
    da : xarray.DataArray, xarray.Dataset
    
    """    
    k_size = buffer_km*2+1
    center_idx = buffer_km
    kernel = np.fromfunction(lambda i,j: ((center_idx - i) ** 2 + (center_idx - j) ** 2)<=buffer_km**2,(k_size,k_size)).astype(int)
    if return_type == 'oneDay':
        #calculate convolution of POH>threshold
        # over_thresh = poh.sel(time=date) > poh_thresh
        poh_sel = poh.sel(time=date) 
        if get_likelihood:
            #here possible hail will return a likelihood based on the number of pixels with POH>threshold
            possible_hail = convolve(poh_sel> poh_thresh,kernel,mode='same') / np.sum(kernel)
        else: #here possible hail return a yes/no field of locations with possible hail
            possible_hail = convolve(poh_sel> poh_thresh,kernel,mode='same') > 0
        ds = poh_sel.to_dataset()
        ds=ds.assign(possible_hail=(('chy','chx'),possible_hail))
        ds=ds.expand_dims('time')
        haz_sel = fct.hazard_from_radar(ds,varname='possible_hail', extent = extent)
        #assert that shape is the same as haz_poh!
        np.testing.assert_array_equal(haz_poh.centroids.lat,haz_sel.centroids.lat)
        np.testing.assert_array_equal(haz_poh.centroids.lat,haz_sel.centroids.lat)
        return haz_sel,ds.possible_hail
    elif return_type == 'all':
        for date_now in poh.time:
            print(date_now.values)
            poh_sel = poh.sel(time=date_now) 
            if get_likelihood:
                #here possible hail will return a likelyhood based on the number of pixels with POH>threshold
                possible_hail = convolve(poh_sel> poh_thresh,kernel,mode='same') / np.sum(kernel)
            else: #here possible hail return a yes/no field of locations with possible hail
                possible_hail = convolve(poh_sel> poh_thresh,kernel,mode='same') > 0
            possible_hail = np.expand_dims(possible_hail,axis=0)
            if date_now == poh.time[0]:
                ph_all = possible_hail
            else:
                ph_all = np.concatenate([ph_all,possible_hail],axis=0)
                
            # ds = poh_sel.to_dataset()
            # ds=ds.assign(possible_hail=(('chy','chx'),possible_hail))
            # ds=ds.expand_dims('time')    
        ds_PH = poh.copy(deep=True).to_dataset()
        ds_PH=ds_PH.assign(possible_hail=(('time','chy','chx'),ph_all))
        ds_PH=ds_PH.drop_vars('BZC')       
        np.testing.assert_array_equal(poh.lat,ds_PH.lat)
        np.testing.assert_array_equal(poh.time,ds_PH.time)
        return ds_PH
    
def get_date_or_category(date,min_day,max_day,m2,m1,now,p1,p2,mode = 'date'):
    raise Warning("Function is deprecated, use get_date instead")
    
def get_date(in_array,date,index_day0,mode='date'):
    raise Warning("Function is deprecated, use get_date from pre_process.py instead")

def smooth_monotonic(x,y,plot=False):
    """
    monotonic smoother based on https://stats.stackexchange.com/questions/467126/monotonic-splines-in-python
    x must be ordered increasing!
    
    """
    assert(len(y)==len(x))
    N=len(y)
    # Prepare bases (Imat) and penalty
    dd = 3
    E  = np.eye(N)
    D3 = np.diff(E, n = dd, axis=0)
    D1 = np.diff(E, n = 1, axis=0)
    la = 100
    kp = 10000000

    # Monotone smoothing
    ws = np.zeros(N - 1)

    for it in range(30):
        Ws      = np.diag(ws * kp)
        mon_cof = np.linalg.solve(E + la * D3.T @ D3 + D1.T @ Ws @ D1, y)
        ws_new  = (D1 @ mon_cof < 0.0) * 1
        dw      = np.sum(ws != ws_new)
        ws      = ws_new
        if(dw == 0): break  
        
    if plot:
        # Monotonic and non monotonic fits
        z  = mon_cof
        z2 = np.linalg.solve(E + la * D3.T @ D3, y)
        # Plots
        plt.scatter(x, y, linestyle = 'None', color = 'gray', s = 0.5, label = 'raw data')
        plt.plot(x, z, color = 'red', label = 'monotonic smooth')
        plt.plot(x, z2, color = 'blue', linestyle = '--', label = 'unconstrained smooth')
        plt.legend(loc="lower right")
        plt.show()
    return mon_cof


