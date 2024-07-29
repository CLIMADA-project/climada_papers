# -*- coding: utf-8 -*-
"""
Functions to do regridding of climada hazards
currently includes
    - aggregation of climada hazard to a regular grid with lower resolution
    - gridding of Crowd-sourced reports
@author: Raphael Portmann, Timo Schmid, Leonie Villiger
"""
import numpy as np
import sys
import xarray as xr
from pathlib import Path
import datetime as dt
import shapely.geometry as sg
from pyproj import Proj
import cartopy.crs as ccrs
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.signal import convolve
from climada import CONFIG
from climada.hazard import Centroids
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
import geopandas as gpd
from pyproj import CRS
from scipy import sparse
import pandas as pd
import copy
import warnings

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
    for (date, index), count in merged.groupby(['date','index_right']).size().loc[
        merged.groupby(['date','index_right']).size()<16].iteritems():

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

def gdf_from_hazard(hazard):
    """Create GeoDataFrame from hazard object with columns 'intensity',
    'date', 'event_id', and hazard centroids as geometry

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
    cell, crs, extent = sc.create_empty_grid(epsg = original_grid_epsg, 
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

def grid_crowd_source(crowd_source_data,pop_path,date_sel='20210628',
                      kernel_size_gap_fill=7,min_crowd_count=100,MESHS_nc=None):
    """function to grid crowd-sourced data

    Args:
        crowd_source_data (Path,str,pd.DataFrame): dataframe of crowd-source reports
        pop_path (Path,str): path for population data source
        date_sel (str, optional): 'all' for all dates, otherwise '%Y%m%d'. Defaults to '20210628'.
        kernel_size_gap_fill (int, optional): size of kernel (must be odd). Defaults to 7.
        min_crowd_count (int, optional): min count of crowd-source reports. 
            Only used if date_sel='all'. Defaults to 100.

    Returns:
        ds: xr.Dataset of gridded crowd-sourced data
    """

    assert(kernel_size_gap_fill%2 == 1), 'kernel_size_gap_fill must be an odd number'

    if isinstance(crowd_source_data,(str,Path)):
        path = crowd_source_data
        crowd = pd.read_csv(path,sep=',')
    elif isinstance(crowd_source_data,pd.DataFrame):
        crowd = crowd_source_data

    #filter out filtered_out points
    crowd = crowd.loc[crowd['FILTEREDOUT'] == 0,:]

    #population pre-processing
    population = pd.read_csv(pop_path,sep=';')
    population['chx'] = np.round(population['E_KOORD']/1000)
    population['chy'] = np.round(population['N_KOORD']/1000)
    pop_sum = population.groupby(['chx','chy']).agg({'B21BTOT': 'sum'})
    chx_range = [2484, 2838]
    chy_range = [1073, 1299]

    chx = np.arange(chx_range[0],chx_range[1])
    chy = np.arange(chy_range[0],chy_range[1])

    pop = np.zeros((len(chy),len(chx)))
    for i in range(len(pop_sum)): #tqdm()
        x_coord = chx==pop_sum.index.get_level_values('chx')[i]
        y_coord = chy==pop_sum.index.get_level_values('chy')[i]
        pop[y_coord, x_coord] = pop_sum.iloc[i]['B21BTOT']


    #crowd source pre processing
    crowd['Time'] = pd.to_datetime(crowd['Time'], format='%Y-%m-%d %H:%M:%S')
    # round to nearest minute
    crowd['Time'] = crowd['Time'].dt.round('min')
    # # select times on 2021-06-28
    # crowd = crowd.loc[crowd['Time'].dt.date == pd.to_datetime('2021-06-28').date()]
    # select times in range (1 hailday: 6UTC-6UTC)

    #create copy of crowd dataframe with all values
    crowd_all = crowd.copy(deep=True)

    if date_sel == 'all':
        # use -6h, to get haildays correctly! (1 hailday: 6UTC-6UTC)
        crowd['hail_day'] = (crowd['Time'] - pd.Timedelta('6h')).dt.date
        grpy = crowd.groupby('hail_day').ID.count()
        sel_dates = grpy.index[grpy>min_crowd_count] #min *min_crowd_count* reports per day (CH wide)
        dates = [d.strftime('%Y%m%d') for d in sel_dates if d.year >=2017] #only select dates after 2017
        if MESHS_nc is not None:
            #only select dates after 2017 and up to 2021 and only Hail season (April-September)
            dates = [d.strftime('%Y%m%d') for d in sel_dates if 
                     d.year >=2017 and d.year<=2021 and d.month>=4 and d.month<=9]
    else:
        dates = [date_sel]

    total_removed_reports = 0 #for MESHS removal
    for date in dates:
        print(date)
        date_plus_one = (pd.Timestamp(date) + pd.Timedelta(days=1)).strftime('%Y%m%d')
        datelist_range = [f'{date}0600',f'{date_plus_one}0600']
        crowd = crowd_all.loc[(crowd_all['Time'] >= datelist_range[0]) & 
                              (crowd_all['Time'] <= datelist_range[1])]

        # convert crowd size classes to mm todo: check this & find more robust way
        crowd['size_mm'] = np.nan
        crowd['size_mm'][crowd['size']==(6+10)] = 68
        crowd['size_mm'][crowd['size']==(5+10)] = 43
        crowd['size_mm'][crowd['size']==(4+10)] = 32
        crowd['size_mm'][crowd['size']==(3+10)] = 23
        crowd['size_mm'][crowd['size']==(2+10)] = 8
        crowd['size_mm'][crowd['size']==(1+10)] = 5
        crowd['size_mm'][crowd['size_mm']==0] = 0

        # round to nearest km (for later groupby operation)
        crowd['chx'] = np.round((crowd['x'])/1000-600)
        crowd['chy'] = np.round((crowd['y'])/1000-200)


        # group by chx, chy to speed up array loop
        crowd_mean = crowd.groupby(['chx','chy']).agg({'size_mm': 'mean'})
        crowd_count = crowd.groupby(['chx','chy']).agg({'size_mm': 'count'})


        crowd_sizes = np.zeros((len(chy),len(chx)))
        crowd_sizes[:] = np.nan
        crowd_counts = np.zeros((len(chy),len(chx)))

        # assign grouped point values to array
        for i in range(len(crowd_mean)):
            x_coord = chx-2600==crowd_mean.index.get_level_values('chx')[i]
            y_coord = chy-1200==crowd_mean.index.get_level_values('chy')[i]
            crowd_sizes[y_coord, x_coord] = crowd_mean.iloc[i]['size_mm']
            crowd_counts[y_coord, x_coord] = crowd_count.iloc[i]['size_mm']

        # reporting_fraction = crowd_counts/pop
        # reporting_fraction = np.nan_to_num(reporting_fraction, nan=0.0, posinf=1)
        # reporting_fraction[reporting_fraction<=0] = np.nan
        # print(reporting_fraction.shape)

        crowd_estimate = {}
        crowd_estimate['size'] = crowd_sizes.copy()
        crowd_estimate['count'] = crowd_counts.copy()
        crowd_estimate['chx'] = chx
        crowd_estimate['chy'] = chy

        ########################
        #filter with MESHS
        if MESHS_nc is not None:
            meshs_here = MESHS_nc.sel(chx=slice(chx[0]*1000+500,chx[-1]*1000+500),
                                      chy=slice(chy[0]*1000+500,chy[-1]*1000+500),time=date)
            meshs_arr = meshs_here.MZC.values
            assert(meshs_arr.shape == crowd_estimate['size'].shape)
            meshs_bol = meshs_arr>0
            #
            buffer_km = 2 #5km kernel
            k_size = buffer_km*2+1
            center_idx = buffer_km
            kernel = np.fromfunction(lambda i,j: ((center_idx - i) ** 2 + (center_idx - j) ** 2)<=buffer_km**2,(k_size,k_size)).astype(int)
            possible_extreme_hail = convolve(meshs_bol,kernel,mode='same') > 0
            false_extreme_reports = (~possible_extreme_hail) & (crowd_estimate['size']==68) #tennis ball size report
            total_removed_reports += np.sum(false_extreme_reports)
            crowd_estimate['size'][false_extreme_reports] = np.nan

        ########################
        # don't trust gridpoints with only 1 report, big population and big hail
        # crowd_estimate['size'][(crowd_estimate['count']==1) & (pop>500) & (crowd_estimate['size']>30)] = np.nan

        # gridpoints with high enough population and no reports are assumed to have no hail
        crowd_estimate['size'][(pop>2000) & (np.isnan(crowd_sizes))] = 0

        crowd_estimate['size_interpolated'] = crowd_estimate['size'].copy()

        # kernel marching for filling in the gaps todo: convert to disk kernel with gaussian weighting
        kernel_size = kernel_size_gap_fill # odd number
        kernel_size = kernel_size//2 # kernel is implemented as indexing with +- kernel_size
        for i in range(crowd_estimate['size'].shape[0]): # loop over x
            for j in range(crowd_estimate['size'].shape[1]): # loop over y
                if np.isnan(crowd_estimate['size'][i,j]): # only fill gridpoints with nan
                    if np.sum(crowd_estimate['size'][i-kernel_size:i+kernel_size,j-kernel_size:j+kernel_size]>=0) > 3: # only fill if there are at least 4 valid gridpoints whitinin the kernel
                        crowd_estimate['size_interpolated'][i,j] = np.nanmean(crowd_estimate['size'][i-kernel_size:i+kernel_size,j-kernel_size:j+kernel_size]) # fill with mean of kernel

        # filter out high frequency noise with small kernel
        noisy = crowd_estimate['size_interpolated'].copy()
        denoised = np.zeros((len(chy),len(chx)))
        denoised[:] = np.nan

        kernel_size = 3//2
        for i in range(denoised.shape[0]):  # loop over x
            for j in range(denoised.shape[1]):  # loop over y
                if np.sum(noisy[i-kernel_size:i+kernel_size,j-kernel_size:j+kernel_size]>=0) > 1: # only alter gridpoints with at least 2 valid gridpoints whitinin the kernel
                    denoised[i,j] = np.nanmean(noisy[i-kernel_size:i+kernel_size,j-kernel_size:j+kernel_size]) # fill with mean of kernel

        crowd_estimate['denoised'] = denoised.copy()

        #convert to xr.DataSet
        ds = xr.Dataset(data_vars={'h_raw': (['chy', 'chx'], crowd_estimate['size']),
                            'h_grid': (['chy', 'chx'], crowd_estimate['size_interpolated']),
                            'h_smooth': (['chy', 'chx'], crowd_estimate['denoised']),
                            'n_rep': (['chy', 'chx'], crowd_estimate['count'])},
                    coords={'chx': (['chx'], crowd_estimate['chx']*1000),'chy': (['chy'], crowd_estimate['chy']*1000)})

        ds = ds.expand_dims({'time': [dt.datetime.strptime(date, '%Y%m%d')]})
        if date == dates[0]:
            ds_all = ds.copy(deep=True)
        else:
            ds_all = xr.concat([ds_all, ds], dim='time')

    ##get lat lon (testing)
    projdef = ccrs.epsg(2056)#.proj4_init
    projdef
    #create meshgrid
    meshX,meshY = np.meshgrid(ds_all.chx, ds_all.chy)
    p = Proj(projdef)
    lon, lat = p(meshX,meshY, inverse=True)
    ds_all=ds_all.assign_coords({'lon':(('chy','chx'),lon)})
    ds_all=ds_all.assign_coords({'lat':(('chy','chx'),lat)})


    if MESHS_nc is not None:
        print(f'Total removed reports (Tennis ball with no MESHS within {kernel_size}km kernel): '+str(total_removed_reports))
    return ds_all


#general regriding function inspired by https://james-brennan.github.io/posts/fast_gridding_geopandas/
def create_cell(extent2056,cell_size,epsg=2056):
    """create regular grid of cells (e.g. 1kmx1km grid)

    Args:
        extent2056 (tuple): (xmin, xmax, ymin, ymax) values in swiss coordinates (epsg:2056)
        cell_size (int): cell size in coordinate units (with epsg:2056, 1 equals 1km)
        epsg (int, optional): EPSG. Defaults to 2056.

    Returns:
        GeoDataFrame: gdf with one row for each cell
    """
    xmin, xmax, ymin, ymax = extent2056
    # projection of the grid
    if epsg==2056:
        crs = ccrs.epsg(epsg)
    else:
        raise UserWarning('Function is only tested with epsg=2056 so far!')
    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            # bounds
            x1 = x0-cell_size
            y1 = y0+cell_size
            grid_cells.append(sg.box(x0, y0, x1, y1)  )
    cell = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                                     crs=crs)

    #get Coordinates from geometry
    cell['KoordinateOst']=cell.geometry.centroid.x
    cell['KoordinateNord']=cell.geometry.centroid.y
    cell['lon']=cell.geometry.centroid.to_crs(epsg=4326).x
    cell['lat']=cell.geometry.centroid.to_crs(epsg=4326).y
    return cell

def group_to_cell(cell,gdf_exp,aggfunc,count_var,plot=False):
    """Group values to pre-defined cells (regular grid)

    Args:
        cell (GeoDataFrame): gdf of cells (regular grid of square cells)
        gdf_exp (GeoDataFrame): gdf of exposure points
        aggfunc (dict): Dictionary with signature {var_1:*aggregation type*,var2:...},
            aggregation type can be e.g.: 'sum', 'mean', 'first'.
        count_var (str): Column name used to count the entries. Can be any column
        plot (bool, optional): Plot the grid. Defaults to False.

    Returns:
        GeoDataFrame: gdf with aggregated values
    """
    cell_out = cell.copy(deep=True)
    columns = list(aggfunc.keys())+['n_count']
    aggfunc_w_count = dict(aggfunc,**{count_var:'count'})
    # Merge the point data to grid
    merged = gpd.sjoin(gdf_exp, cell_out, how='left', op='within')
    # Compute stats per grid cell -- aggregate fires to grid cells with dissolve
    dissolve = merged.dissolve(by="index_right", aggfunc=aggfunc_w_count).rename(columns={count_var:'n_count'})

    for column in columns:
        cell_out.loc[dissolve.index, column] = dissolve[column].values
        if plot:
            #Plot example
            if column in ['n_count','Baujahr','Volumen']:
                norm = None
            else:
                norm = colors.LogNorm(max(10,cell_out[column].min()),cell_out[column].max())
            fig,ax= plt.subplots(1,1,subplot_kw={'projection':ccrs.epsg(2056)},figsize=(6,6))
            cell_out.plot(column=column, figsize=(12, 8), cmap='viridis', edgecolor="grey",ax=ax,zorder=0.1,
                        legend=True,legend_kwds={'label':column},norm=norm)#locator=ticker.LogLocator()
            sc.plot_canton(ax,canton=['Zürich','Bern','Luzern','Aargau'])
            # if save_path:
            #     plt.savefig(save_path,bbox_inches='tight',dpi=200) #out_dir + '/explorative/gvz/grid_%s.png'%(column)

    return cell_out

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