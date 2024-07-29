import datetime
from dateutil.parser import parse
import numpy as np
import sys
import pandas as pd
import xarray as xr
import geopandas as gpd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from pyproj import Transformer
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
data_dir = str(CONFIG.local_data.data_dir)
crowd_url = str(CONFIG.crowd_url)

def plot_crowd(meas_crowd,dpi=200,relief=True,figsize=(11,10)):
    """Plotting function derived from Jerome Kopp

    Args:
        meas_crowd (pd.Dataframe): DF with crowd sourced data
        dpi (int, optional): DPI of figure. Defaults to 200.
        relief (bool, optional): Whether to plot relief. Defaults to True.
    Returns:
        fig,ax: Figure and axes object
    """

    #load data
    if relief:
        # read relief and data
        da_relief = xr.open_rasterio(data_dir+'/ch_shapefile/relief_georef_clipped_swiss.tif')

        # Compute the lon/lat coordinates with rasterio.warp.transform
        ny, nx = len(da_relief['y']), len(da_relief['x'])
        x, y = np.meshgrid(da_relief['x'], da_relief['y'])
        # Rasterio works with 1D arrays
        outProj = 'epsg:4326' # WGS84, see https://epsg.io/4326
        transformer = Transformer.from_crs(da_relief.crs, outProj)
        lat, lon = transformer.transform(x.flatten(), y.flatten())
        lon = np.asarray(lon).reshape((ny, nx))-0.01
        lat = np.asarray(lat).reshape((ny, nx))
        da_relief.coords['lon'] = (('y', 'x'), lon)
        da_relief.coords['lat'] = (('y', 'x'), lat)

        # get band
        da_relief = da_relief.isel(band=0, drop=True)
        da_relief = da_relief.where(da_relief > 1, drop=True)

        # let's load Lakes from Swisstopo
        gdf_lakes = gpd.read_file(data_dir+"/ch_shapefile/swiss-lakes-maps.json")

    #plot data
    fig = plt.figure(figsize=figsize, dpi=dpi)
    prj = ccrs.AlbersEqualArea(8.222665776, 46.800663464)
    prj2 = ccrs.PlateCarree()
    ax = plt.axes(projection=prj)
    Ticino = [8.55, 9.4, 45.8, 46.3]
    Central = [7.6, 8.5, 46.5, 47.2]
    ZRH = [8.5,8.6,47.3,47.45]
    CH = [5.8, 10.7, 45.7, 47.95]

    extent = CH
    ax.set_extent(extent)

    fs = 14
    nudge = 0.05

    # add Relief
    if relief:
        da_relief.plot(ax=ax, x='lon', y='lat', cmap="Greys_r",
                    norm=colors.Normalize(vmin=110, vmax=255),
                    add_colorbar=False, transform=prj2)

    # Crowdsourced data
    cmap_crowd = cm.get_cmap('Purples')

    scat_kw = dict(marker='o', zorder=4, s=16, linewidth=0.4, c=meas_crowd['size'],
                   cmap=cmap_crowd ,edgecolors="purple", alpha=1, transform=ccrs.Geodetic())
    crowd = ax.scatter(meas_crowd['Lon'],meas_crowd['Lat'],**scat_kw)

    handles, labels = crowd.legend_elements()
    labels = ['< coffee bean', 'Coffee bean', 'One franc coin', 'Five franc coin', 'Golf ball', 'Tennis ball']
    for ha in handles:
        ha.set_markeredgecolor("purple")
    legend1 = ax.legend(handles, labels, loc="upper left", title="Report size",
                        prop={'size': fs-2}, framealpha=1, markerscale=2)
    legend1.get_title().set_fontsize(fs)
    ax.add_artist(legend1)

    if relief:
        # add lakes
        gdf_lakes.plot(ax=ax, edgecolor='none', color="cornflowerblue", transform=prj2)
        # add country border
        ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1.5)
        # add urban areas
        #gdf_urban.boundary.plot(ax=ax,transform=prj2,zorder=1.4,linewidth=1,alpha=0.4, color='black')

    # format gridlines
    gl = ax.gridlines(crs=prj2, draw_labels=True, linewidth=0.5, color='gray',
                      alpha=0.4, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': fs}
    gl.ylabel_style = {'size': fs}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.title('')
    return fig, ax

def plot_crowd_date(date:str,processed:bool=False,dpi:int=200,relief:bool=True,
                    figsize:tuple=(11,10)):
    """plot corwd sourced data for a given date

    Args:
        date (str): date
        processed (bool, optional): Wheter to use locally saved processed data,
        or download raw data. Defaults to False.
        dpi (int, optional): Figure DPI. Defaults to 200.
        relief (bool, optional): Whether to plot a relief. Defaults to True.

    Returns:
        fig, ax (tuple): Figure and axes object 
    """    
    sel_day = datetime.datetime.strptime(date, '%Y-%m-%d')

    if processed:
        assert(sel_day<=datetime.datetime(2021,12,31))
        data_dir = 'C:/Users/timo_/Documents/PhD/data/'
        file_crowd = 'crowd-source/Reports_20150512_20220901_filtered.csv'
        crowd_data = pd.read_csv(data_dir + file_crowd, 
                                 dtype={'ID':str, 'x':int, 'y':int, 'size':int, 'CustomTime':str,
                                    'CustomLocation':str, 'OsVersion':str, 'AppVersion':str,
                                    'Language':str, 'EventDate':int, 'maxCZC':float,
                                    'CZCfilteredout':int, 'TDiff':float, 'SubvsEventbad':int,
                                    'Multi1':int, 'Multi2':int, 'Multi3':int,
                                    'blacklisted':int, 'xr':int, 'yr':int,
                                    'FILTEREDOUT':int})
        
        crowd_data=process_crowd_data(crowd_data)
        crowd_filtered = crowd_data[(crowd_data['size'].isin(np.arange(11,17))) &
                                    (crowd_data['FILTEREDOUT'] == 0)]
        #select 21.6.2021
        date_plot = sel_day.strftime("%Y-%m-%d")
        date_crowd = pd.to_datetime(date_plot).date()
        crowd_sel_day = crowd_filtered.loc[crowd_filtered['hailday']==date_crowd]
    else:
        #download crowd-sourced data from MCH server
        tmr = sel_day + datetime.timedelta(days=2)
        sel_day_str = sel_day.strftime("%Y-%m-%d")
        tmr_str = tmr.strftime("%Y-%m-%d")
        crowd_sel_day = pd.read_csv(f'{crowd_url}?start={sel_day_str}&end={tmr_str}')
        crowd_sel_day=sc.process_crowd_data(crowd_sel_day,processed=False)

    meas_crowd = crowd_sel_day.loc[crowd_sel_day['hailday']==sel_day.date()]
    meas_crowd = meas_crowd.sort_values('size')
    fig,ax=plot_crowd(meas_crowd,dpi=dpi,relief=relief,figsize=figsize)
    return fig, ax

def process_crowd_data(crowd_data,processed=True):
    """Processing raw crowd sourced data

    Args:
        crowd_data (pd.DataFrame): Dataframe of crowd sourced data
        processed (bool, optional): whether data is already pre processed. Defaults to True.

    Returns:
        crowd_data: pandas.DataFrame with crowd sourced data
    """

    #remove points with default locations
    crowd_data = crowd_data.loc[~((crowd_data['x']==717144) & (crowd_data['y']==95751))]
    crowd_data = crowd_data.loc[~((crowd_data['x']==537969) & (crowd_data['y']==152459))]
    crowd_data = crowd_data.loc[~((crowd_data['x']==683472) & (crowd_data['y']==247852))]

    #Convert time columns to datetime formate
    if processed:
        crowd_data['Time'] = pd.to_datetime(crowd_data['Time'])
        crowd_data['SubmissionTime'] = pd.to_datetime(crowd_data['SubmissionTime'])
        crowd_data['Timer'] = pd.to_datetime(crowd_data['Timer'])
    else:
        
        crowd_data = crowd_data.rename(columns={'Type':'size'})
        crowd_data = crowd_data.loc[((crowd_data['size']>=10) & (crowd_data['size']<=16))]
        crowd_data['Time'] = pd.to_datetime(crowd_data['Time'],unit='ms')
        crowd_data['SubmissionTime'] = pd.to_datetime(crowd_data['SubmissionTime'],unit='ms')

    #add columns for date, hailday, month and year
    crowd_data['days'] = crowd_data['Time'].dt.date
    crowd_data['hailday'] = (crowd_data['Time']-pd.Timedelta(hours=6)).dt.date #use -6h, to get haildays correctly! (1 hailday: 6UTC-6UTC)
    crowd_data['months'] = crowd_data['Time'].dt.to_period('M')
    crowd_data['years'] = crowd_data['Time'].dt.to_period('Y')

    conditions = [
        (crowd_data['size'] == 10),(crowd_data['size'] == 11),
        (crowd_data['size'] == 12),(crowd_data['size'] == 13),
        (crowd_data['size'] == 14),(crowd_data['size'] == 15),
        (crowd_data['size'] == 16)]

    # create a list of the values we want to assign for each condition
    values = ['No hail', '< coffee bean', 'Coffee bean', 'One franc coin',
              'Five franc coin', 'Golf ball', 'Tennis ball']
    values2 = ['No hail', '>0-5 [mm]', '5-8 [mm]', '23 [mm]',
               '32 [mm]', '43 [mm]', '68 [mm]']


    # create a new column and use np.select to assign values to it using our lists as arguments
    crowd_data['size_text'] = np.select(conditions, values)
    crowd_data['size_mm'] = np.select(conditions, values2)
    #inProj = 'epsg:2056' # CH1903+ / LV95, see https://epsg.io/2056
    inProj = 'epsg:21781' # CH1903 / LV03, see https://epsg.io/21781
    outProj = 'epsg:4326' # WGS84, see https://epsg.io/4326

    transformer = Transformer.from_crs(inProj, outProj)
    coords_crowd = transformer.transform(crowd_data['x'],crowd_data['y'])

    crowd_data['Lat'] = coords_crowd[0].tolist()
    crowd_data['Lon'] = coords_crowd[1].tolist()
    return crowd_data