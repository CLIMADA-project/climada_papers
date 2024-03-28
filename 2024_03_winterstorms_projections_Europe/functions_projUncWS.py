#import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import copy as cp
import pandas as pd
from scipy import sparse
import cartopy.crs as ccrs
##from timeit import default_timer as timer
#from os import# mkdir, remove, rmdir

#from climada.engine import Impact, ImpactCalc
#from climada.entity import ImpactFunc,ImpactFuncSet
from climada.hazard import Hazard,Centroids

## Define functions
def set_centroids(da,stack=True,haztype='WS',timeres="day",plot=False, printout=False):
    '''Function which takes a xarray.DataArray as an input and returns a climada.hazard object, with centroids corresponding to
        latitude and longitude of the DataArray.'''

    period = da.name
    lat = da.lat
    lon = da.lon

    # number of points
    n_lat = len(lat)
    n_lon = len(lon)

    #stack first members and time resolution, and then latitutde and longitude to obtain the event matrix
    if stack:
        events = da.stack(events=("member",timeres)) #need to be of size n_ev * ncent = nmem*180 * ncent
        nmem = len(da["member"])
    else:
        events = da.rename({timeres:"events"})
        nmem = 1
    evmat = events.stack(cent=("lat","lon"))
    #print('Evmat: '+str(evmat.shape))
    n_ev = len(evmat["events"])

    #intiate Hazard object, using Centroids.from_pix_bounds
    haz = Hazard(haztype)

    #initiate lat and lon array for Centroids.from_lat_lon
    lat_ar = lat.values.repeat(len(lon))
    lon_ar = lon.values.reshape(1,n_lon).repeat(n_lat,axis=0).reshape(n_lon*n_lat)

    haz.centroids = Centroids.from_lat_lon(lat_ar,lon_ar)
    haz.intensity = sparse.csr_matrix(evmat)
    haz.units = 'm/s'
    ev_id = np.arange(n_ev, dtype=int)
    haz.event_id = ev_id
    ev_names = events.coords["events"].to_numpy().tolist()
    #ev_names = pd.to_datetime(evmat.events,format='YYYY-MM-DD')
    haz.event_name = ev_names
    haz.orig = np.zeros(n_ev, bool)
    haz.frequency = np.ones(n_ev)/(nmem*30)
    haz.fraction = haz.intensity.copy()
    haz.fraction.data.fill(1)
    haz.centroids.set_meta_to_lat_lon()
    haz.centroids.set_geometry_points()
    haz.check()
    if printout:
        print('Lat resolution original: '+str(-latreso)+' ,rounded: '+str(-latres)+
              '\nLon resolution original: '+str(lonreso)+' ,rounded: '+str(lonres))
        print('lat: '+str(n_lat)+', lon: '+str(n_lon))
        print('Check centroids borders:', haz.centroids.total_bounds)
    if plot:
        haz.centroids.plot()

    return haz

def sel_reg_exp(reg_ids,exp):
    """function that takes a list of country ISO codes reg_ids and a CLIMADA LitPop exposure
    object exp as inputs and returns a copy of the exposure object with only the regions selected in reg_ids"""
    sel_exp = cp.deepcopy(exp)
    sel_exp.gdf = sel_exp.gdf.where(sel_exp.gdf['region_id'].isin(reg_ids)).dropna()
    return sel_exp

def get_lat_lon_res(ds):
    '''Function to obtain the average lat and lon gridspacing from a dataset of a non regular model grid. '''
    lat = ds.coords['lat']
    lon = ds.coords['lon']
    difflat = lat - lat.shift(lat=1)
    latres = difflat.mean().to_numpy()
    difflon = lon - lon.shift(lon=1)
    lonres = difflon.mean().to_numpy()
    return latres, lonres

def make_fn(addlist,basename="",sep="_",filetype=''):
    """Function to facilitate the creation of file names that takes a list of
    strings addlist, a string basename, a string sep and a string filetype
    as input and returns the string consisting of the keywords from addlist separated
    by the sep separator and concatenated with the string basename and the string filetype"""
    return sep.join(addlist)+sep+basename+filetype

#each impact function gets a preprocessing function to prepare the hazard data
def mask_qt(ds,q,mask_abs=None,timeres='day',stack=True,pastname='historical',futname='ssp585',cutarea=1000000):
    '''Function taking a dateset as an input, and returning it with the values below the quantile q at each grid cell
       masked. The mask is computed for each gridcell for the past period. Fields for which less than cutarea is above the
       quantile are dropped'''

    Upast = ds[pastname]
    Ufut = ds[futname]
    if stack:
        Upast = Upast.stack(real=("member",timeres))
        Ufut = Ufut.stack(real=("member",timeres))
        dim="real"
    else:
        dim = timeres

    Upast_qt = Upast.quantile(q,dim=dim)
    U_mask_past = Upast.where(Upast>Upast_qt)
    U_mask_fut = Ufut.where(Ufut>Upast_qt)

    #mask values below threshold
    if mask_abs:
        U_mask_past = U_mask_past.where(U_mask_past>=mask_abs)
        U_mask_fut = U_mask_fut.where(U_mask_fut>=mask_abs)

    latres, lonres = get_lat_lon_res(ds)
    gcarea = latres*lonres*100*100 #gridcell area approximated: 1 deg corresponds to 100km
    threshold = round(cutarea/gcarea)

    U_mask_past = U_mask_past.dropna(dim=dim,thresh=threshold) #keep fields for which at least X values are not NaN
    U_mask_fut = U_mask_fut.dropna(dim=dim,thresh=threshold) #keep fields for which at least X values are not NaN

    #unstack and assemble
    U_mask_past = U_mask_past.unstack().fillna(0)
    U_mask_fut = U_mask_fut.unstack().fillna(0)
    if futname != pastname:
        U_mask_fut.name = futname
        U_mask = xr.combine_by_coords([U_mask_past, U_mask_fut]).fillna(0)
    else:
        U_mask = U_mask_past

    return U_mask

def scale_qt(ds,q,mask_abs=None,timeres='day',stack=True,pastname='historical',futname='ssp585',cutarea=1000000):
    '''Same as mask_qt but scaling by the q quantile.'''
    Upast = ds[pastname]
    Ufut = ds[futname]

    if stack:
        Upast = Upast.stack(real=("member",timeres))
        Ufut = Ufut.stack(real=("member",timeres))
        dim = "real"
    else:
        dim = timeres

    Upast_qt = Upast.quantile(q,dim=dim)
    U_mask_past = Upast.where(Upast>Upast_qt)
    U_mask_fut = Ufut.where(Ufut>Upast_qt)
    #mask values below threshold
    if mask_abs:
        U_mask_past = U_mask_past.where(U_mask_past>=mask_abs)
        U_mask_fut = U_mask_fut.where(U_mask_fut>=mask_abs)

    U_scaled_past = (U_mask_past-Upast_qt)/Upast_qt
    U_scaled_fut = (U_mask_fut-Upast_qt)/Upast_qt

    latres, lonres = get_lat_lon_res(ds)
    gcarea = latres*lonres*100*100 #gridcell area approximated: 1 deg corresponds to 100km
    threshold = round(cutarea/gcarea)

    U_scaled_past = U_scaled_past.dropna(dim=dim,thresh=threshold) #keep fields for which at least X values are not NaN
    U_scaled_fut = U_scaled_fut.dropna(dim=dim,thresh=threshold) #keep fields for which at least X values are not NaN

    #unstack and assemble
    U_scaled_past = U_scaled_past.unstack().fillna(0)
    U_scaled_fut = U_scaled_fut.unstack().fillna(0)
    if futname != pastname:
        U_scaled_fut.name = futname
        U_scaled = xr.combine_by_coords([U_scaled_past, U_scaled_fut]).fillna(0)
    else:
        U_scaled = U_scaled_past
    return U_scaled

#put processing funcs in a dict
pp_func_dic ={}
pp_func_dic['Sw2010'] = mask_qt
pp_func_dic['CubEOT'] = scale_qt
