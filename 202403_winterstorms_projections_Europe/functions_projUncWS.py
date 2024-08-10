#import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import copy as cp
import pandas as pd
from scipy import sparse
import cartopy.crs as ccrs
from climada.hazard import Hazard,Centroids

def read_sfc(filepath, datestart, dateend):
    """Read in a time slice of a surface field from datestart to dateend.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = directory where files are located
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    #First try opening and doing the select assuming everything is working ok with the time axis
    try:
        dat = \
        xr.open_mfdataset\
        (filepath, coords="minimal", join="override", decode_times=True, use_cftime=True).\
        sel(time=slice(datestart, dateend))
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
            print("changing longitude --> lon, latitude --> lat")
        except: pass

    except:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override", decode_times = False)
        try:
            dat=dat.rename({"longitude":"lon", "latitude":"lat"}) #problematic coord names
        except: pass

        dat = xr.decode_cf(dat, use_cftime = True)
        dat = dat.sel(time=slice(datestart, dateend))
        datetimeindex = dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex
        print("Something's wierd about the time axis, decoding manually")

    return dat

def read_field(filepath, datestart, dateend,latmin,latmax,lonmin,lonmax,plev):
    """Read in a time slice from datestart to dateend and calculate the zonal mean.
    Try using datetime64 and if that doesn't work decode times manually.
    Args:
        filepath (string) = path to files e.g., "/path/to/files/*.nc"
        datestart (string) = start date for time slice
        dateend (string) = end date for time slice
    """

    try:
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override",
                 decode_times=True, use_cftime=True).\
                 sel(time=slice(datestart, dateend),lat=slice(latmin,latmax),lon=slice(lonmin,lonmax))

        if len(plev) == 1:
            dat = dat.sel(plev=plev,method="nearest", tolerance=1) #avoid issue for models with inaccurate plevs
        else:
            dat = dat.sel(plev=slice(plev[0]+1,plev[1]-1)) #manual tolerance of 1 because method="nearest" is not implemented for slices

    except:
        print("Something's wierd about the time axis, decoding manually")
        dat = xr.open_mfdataset(filepath, coords="minimal", join="override",
                   decode_times=False)

        dat=xr.decode_cf(dat, use_cftime=True)

        dat=dat.sel(time=slice(datestart, dateend),lat=slice(latmin,latmax),lon=slice(lonmin,lonmax))
        if len(plev) == 1:
            dat = dat.sel(plev=plev,method="nearest", tolerance=1) #avoid issue for models with inaccurate plevs
        else:
            dat = dat.sel(plev=slice(plev[0]+1,plev[1]-1))  #manual tolerance of 1 because method="nearest" is not implemented for slices

        datetimeindex=dat.indexes['time'].to_datetimeindex()
        dat['time'] = datetimeindex

    return dat

def get_lat_lon_res(ds):
    """Function to obtain the average lat and lon gridspacing from a dataset of a non regular model grid.
    Args:
        ds (xr.dataset or xr.dataarry) = input dataset from which average lat and lon resolutions must be calculated
    Ouputs:
        latres, lonres = average latitudinal and longitudinal resolutions
    """
    lat = ds.coords['lat']
    lon = ds.coords['lon']
    difflat = lat - lat.shift(lat=1)
    latres = difflat.mean().to_numpy()
    difflon = lon - lon.shift(lon=1)
    lonres = difflon.mean().to_numpy()
    return latres, lonres

def def_domain(ds,min_lat,max_lat,min_lon,max_lon):
    """Function that takes xr.dataset or xr.dataarry and box coordinates from a domain and crops the dataset or datarray to
    the domain defined by the box coordinates.
    Args:
        ds (xr.dataset or xr.dataarry) = input dataset which lat and lon needs to be cropped
        min_lat, max_lat = minimum and maximum latitudinal coordinates
        min_lon, max_lon = minimum and maximum longitudinal coordinates
    Ouputs:
        ds = xr.dataset or xr.dataarry cropped to the input domain
    """
    LatIndexer, LonIndexer = 'lat', 'lon'
    ds = ds.loc[{LatIndexer: slice(min_lat, max_lat),
                      LonIndexer: slice(min_lon, max_lon)}]
    return ds

def norm_lon(ds):
    """Function that takes xr.dataset or xr.dataarry and normalizes its longitude coordinate
    Args:
        ds (xr.dataset or xr.dataarry) = input dataset which lon coordinates needs to be normalized

    Ouputs:
        ds = xr.dataset or xr.dataarry with normalized lon coordinates
    """
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    return ds.sortby(ds.lon)

def get_ONDJFM_day(ds, months=[1,2,3,10,11,12],timedim="day"):
    """Function that takes xr.dataset or xr.dataarry and months of the year and returns a dataset corresponding to
    the specified months
    Args:
        ds (xr.dataset or xr.dataarry) = input dataset from which monthly data needs to be selected
        months = months which are requested, in ordinal format (1=January, 2=February, ...,12=December)
        timedim = name of the temporal coordinate of the dataset

    Ouputs:
        ds = xr.dataset or xr.dataarry corresponding to the months selected
    """
    return ds.isel({timedim:ds[timedim].dt.month.isin(months)})

## Define functions
def set_centroids(da,stack=True,haztype='WS',timeres="day",plot=False):
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
    #haz.centroids.set_geometry_points()
    haz.check()
    if plot:
        haz.centroids.plot()

    return haz

def sel_reg_exp(reg_ids,exp):
    """function that takes a list of country ISO codes reg_ids and a CLIMADA LitPop exposure
    object exp as inputs and returns a copy of the exposure object with only the regions selected in reg_ids"""
    sel_exp = cp.deepcopy(exp)
    sel_exp.gdf = sel_exp.gdf.where(sel_exp.gdf['region_id'].isin(reg_ids)).dropna()
    return sel_exp

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
