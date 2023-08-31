# -*- coding: utf-8 -*-
"""
Collection of functions for reading and processing data for the use in Climada 
in the framework of the scClim project (scclim.ethz.ch)
 - reading hazard data from radar
 - reading exposure data from .csv/.gpkg files

"""

import datetime as dt
import geopandas as gpd
import numpy as np
import pandas as pd
import sys
import xarray as xr
from climada.entity import Exposures, ImpactFunc, ImpactFuncSet
from climada.hazard import Hazard
from climada.hazard.tag import Tag
from climada.entity.tag import Tag as E_tag
from climada import CONFIG
from scipy import sparse
import warnings
import climada.util.lines_polys_handler as u_lp
from climada.engine import Impact, ImpactCalc


import matplotlib.pyplot as plt
sys.path.append(str(CONFIG.local_data.func_dir))
from scClim.constants import BAUINDEX,BAUJAHR_DICT

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


def read_MFZ_exposure(shapefile,input_type = 'random_points',
                         return_type = 'exp_pnt'):
    """Read MFZ exposure shapefile (.gpkg) to exposure object

    Parameters
    ----------
    shapefile : csv
        file with GVZ exposure (or damage) data
    input_type : string
        'random_points' for point gdf, 'poly' for polygon gdf
    return_type : string
        'exp_pnt' for point exposure, 'exp_shp' for polygon exposure
    Returns
    -------
    exp : Exposure object
        Exposure object
    """

    exp_gdf = gpd.read_file(shapefile)
    exp = Exposures(exp_gdf.to_crs(epsg=4326),unit='CHF')


    if return_type == 'exp_pnt' and input_type=='poly':
        exp_pnt = u_lp.exp_geom_to_pnt(
            exp, res=1000, to_meters=True,
            disagg_met=u_lp.DisaggMethod.DIV, disagg_val=None
        )
        exp_pnt.check()
        return(exp_pnt)
    elif return_type == 'exp_shp' and input_type=='random_points':
        ValueError('Cannot return shapefile with points as input')
    elif return_type == 'exp_shp' and input_type=='poly':
        return(exp)
    elif return_type == 'exp_pnt' and input_type=='random_points':
        return(exp)

def read_MFZ_dmg(gdf_path,years=None):
    gdf = gpd.read_file(gdf_path)

    #remove rows without geometries
    no_geom = gdf.geometry.isna()
    if sum(no_geom)>0.01*len(no_geom):
        ValueError('More than 1% of geometries are invalid. Please check!')
    gdf = gdf[~no_geom]

    #add date
    gdf['date_dt']=pd.to_datetime(gdf.CLAIMS_DATE)

    #rename damage as 'value'
    gdf = gdf.rename(columns={'PAID':'value'})

    #filter years
    if years:
        year_dt = gdf.date_dt.dt.year
        gdf = gdf.loc[(year_dt>=years[0]) & (year_dt<=years[1]),:]

    #filter months (only hail season April-September)
    month_sel = (gdf.date_dt.dt.month>=4) & (gdf.date_dt.dt.month<=9)
    if (sum(gdf.loc[~month_sel,'value'])/sum(gdf['value']))>0.05:
        ValueError('More than 5% of reported damages are not in the \
                    hail season. Please check input!')
    gdf = gdf.loc[month_sel,:]

    assert(gdf.crs =='EPSG:4326')
    gdf['lon'] = gdf.geometry.x
    gdf['lat'] = gdf.geometry.y

    imp_out = gdf_to_imp(gdf,id_col='POLNR',exp_val_col='VEHVALUE')
    return(imp_out)

def read_gvz_exposure(csv_file,on_grid = False,exp_n_assets=False, 
                      exterior_only=False,crs = 'EPSG:2056'):
    """Read CSV file

    Parameters
    ----------
    csv_file : csv
        file with GVZ exposure data
    exp_n_assets : bool
        if True: read number of assets as exposure (rather than their value)
    on_grid : bool
        whether or not data should be interpolated on regular grid (not implemented yet)
    exterior_only : bood
        if yes: estimate value for building exterior (which is easily damaged) from total value
    Returns
    -------
    exp : Exposure object
        Exposure object
    """
    if crs.upper() == 'EPSG:2056':
        x_coord='KoordinateOst'
        y_coord= 'KoordinateNord'
    elif crs.upper() == 'EPSG:4326':
        x_coord='longitude'
        y_coord= 'latitude'

    gvz_df = pd.read_csv(csv_file,sep=";")

    if exp_n_assets:
        unit = ''
        gvz_df['value'] = 1 #set value for each exposure to 1
    else:
        unit = 'CHF'
        gvz_df=gvz_df.rename({'Versicherungssumme':'value'},axis =1)

        if exterior_only:
            # values for all cantons (see grid_cantonal_data.py (5.4.2023))
            # 0.27 is empirical powerlaw, 1560 corresponds to 99% of damages below new estimate
            gvz_df['value'] = np.minimum(gvz_df['value'],1560*gvz_df['value']**0.27) 
            # Note that in theory the scaling from volume to surface would be 
            # a powerlaw with 0.66 as exponent (2/3)
    gdf = gpd.GeoDataFrame(gvz_df, geometry = gpd.points_from_xy(gvz_df[x_coord],
                                                                 gvz_df[y_coord],
                                                                 crs = crs))

    #Create Exposure
    exp = Exposures(gdf,value_unit = unit)
    exp.set_lat_lon()
    exp = exp.to_crs(epsg=4326)
    exp.check()

    return exp

def read_gvz_dmg(csv_file, cant_insurance ,min_dmg = 0,w_rel_dmg=False,
                return_type='imp',exp_n_assets=False,haz=None,years: tuple=None,
                index_dmgs=True,crs = 'EPSG:2056',id_col='VersicherungsID',
                baujahr_filter=''):
    """Read CSV file

    Parameters
    ----------
    csv_file : csv
        file with GVZ exposure (or damage) data
    cant_insurance : str
        Name of cantonal building insurance company (e.g. 'GVZ')
    on_grid : bool
        whether or not data should be interpolated on regular grid 
        (not implemented yet)
    w_rel_dmg : bool
        if relative damage should be saved in gdf too
    return_type: str
        if 'imp': return imp object
        if 'gdf': return gpd.GeoDataFrame
        if 'imp_df': return pd.Dataframe to be used for calib_optimize
        if 'imp_df_yearly': same, but with yearly damages
    exp_n_assets : bool
        if True: read number of assets as exposure (rather than their value)
    haz: climada.hazard
        climada.hazard object to get event_ids. Only needed if return_type='imp_df'
    years: tuple
        tuple of (yearMin,yearMax), to only select the years in this range
    baujahr_filter: str
        if not '', filters damages only for assets with certain year of construction
    Returns
    -------
    xxx : *see "return type" input
    """
    if crs.upper() == 'EPSG:2056':
        x_coord='KoordinateOst'
        y_coord= 'KoordinateNord'
    elif crs.upper() == 'EPSG:4326':
        x_coord='longitude'
        y_coord= 'latitude'

    gvz_dmg = pd.read_csv(csv_file,sep=";") #,parse_dates=['Schadendatum'])
    if 'date_dt' in gvz_dmg.columns:
        gvz_dmg['date_dt'] = pd.to_datetime(gvz_dmg.date_dt)
    else:
        gvz_dmg['Schadendatum'] = gvz_dmg['Schadendatum'].astype(str).str.zfill(8)
        gvz_dmg['date_dt']=pd.to_datetime(gvz_dmg.Schadendatum,format='%d%m%Y')

    #drop 'Schadendatum' column to avoid further confusion
    if 'Schadendatum' in gvz_dmg.columns:
        gvz_dmg.drop(columns=['Schadendatum'])

    #filter years
    if years:
        year_dt = gvz_dmg.date_dt.dt.year
        gvz_dmg = gvz_dmg.loc[(year_dt>=years[0]) & (year_dt<=years[1]),:]

    #filter by year of construction
    if not baujahr_filter=='':
        yearMin, yearMax = BAUJAHR_DICT[baujahr_filter]
        if sum(gvz_dmg.Baujahr.isna())/len(gvz_dmg.Baujahr)>0.01:
            raise ValueError('Over 1% of exposure points without valid Baujahr')
        # Missing entries  for 'Baujahr' are mainly from Bern, when buildings 
        # were renovated! i.e. they were already standing before 2002, 
        # when the radar data begins ->thus set to 2000
        gvz_dmg.Baujahr[gvz_dmg.Baujahr.isna()] = 2000
        gvz_dmg = gvz_dmg.loc[(gvz_dmg.Baujahr>=yearMin) & (gvz_dmg.Baujahr<=yearMax),:]

    if exp_n_assets:
        unit = ''
        gvz_dmg['value'] = 1 #set value for each impact to 1
    else:
        unit = 'CHF'

        if index_dmgs:#Index damages to year 2021
            gvz_dmg['Bauindex'] = gvz_dmg.date_dt.dt.year.map(BAUINDEX[cant_insurance])
            gvz_dmg['value'] = BAUINDEX[cant_insurance].loc[2021]/gvz_dmg['Bauindex']*gvz_dmg['Schadensumme']
        else:
            print('Warning: No indexing of damages')
            gvz_dmg['value'] = gvz_dmg['Schadensumme']

    grpy = gvz_dmg.groupby('date_dt').sum(numeric_only=True)
    grpy = grpy.reset_index()
    dates= grpy['date_dt'][grpy.value>min_dmg]
    #dates = grpy.index[grpy.value>min_dmg] #.strftime('%Y-%m-%d')
    if w_rel_dmg:
        gvz_dmg['rel_dmg'] = gvz_dmg['value'] / gvz_dmg['Versicherungssumme']

    if return_type == 'gdf':
        gdf_dmg = gpd.GeoDataFrame(gvz_dmg,geometry = gpd.points_from_xy(gvz_dmg[x_coord],
                               gvz_dmg[y_coord],crs = crs))
        gdf_dmg = gdf_dmg.to_crs(4326)
        gdf_dmg['latitude'] = gdf_dmg.geometry.y
        gdf_dmg['longitude'] = gdf_dmg.geometry.x
        gdf_dmg['id_col'] = gdf_dmg[id_col]
        return gdf_dmg

    elif return_type=='imp':
        imp_out = gdf_to_imp(gvz_dmg,id_col=id_col,dates=dates,exp_val_col='Versicherungssumme',
                            x_coord=x_coord,y_coord=y_coord,crs=crs,unit=unit)
        return(imp_out)

    elif return_type=='id_col':
        coord_df = gvz_dmg.groupby(id_col).first()[[y_coord,x_coord]]
        return(np.array(coord_df.index))

    elif return_type == 'imp_df_yearly':
        years = np.array(dates.dt.year)
        grpy['year']=years
        grpy_year = grpy.groupby('year').sum()
        impact_data_source = grpy_year.reset_index()[['year','value','Schadensumme']].rename(
            columns={'value':'impact_scaled','Schadensumme':'impact'})
        for col,val in zip(['ISO','regiond_id','reference_year'],['CHE',756,2021]):
            impact_data_source[col]=val
        return impact_data_source

    elif return_type == 'imp_df':
        #get ordinal date and year to add to grpy_df
        ord_dates = np.array(dates.map(dt.datetime.toordinal))
        grpy['date'] = ord_dates
        grpy_sel = grpy[['value','Schadensumme','date']].rename(
            columns={'value':'impact_scaled','Schadensumme':'impact'})

        #initialize dataframe with rows for each event from haz object
        impact_data_source=pd.DataFrame({'event_id':haz.event_id,
                    'event_name':haz.event_name,
                    'year':[dt.datetime.fromordinal(d).year for d in haz.date],
                    'date':haz.date})
        impact_data_source = impact_data_source.merge(grpy_sel,how='outer',on='date')
        return impact_data_source

def gdf_to_imp(gdf,id_col,dates=None,exp_val_col=None,x_coord='lon',
               y_coord='lat',crs = 'EPSG:4326',unit='',haz_tag='HL'):

        #assign dates if not given
        if dates is None:
            dates = pd.Series(gdf.date_dt.unique())

        #create coordinate df
        coord_df = gdf.groupby(id_col).first()[[y_coord,x_coord]]
        coord_df = gpd.GeoDataFrame(coord_df,geometry = gpd.points_from_xy(coord_df[x_coord],
                               coord_df[y_coord],crs = crs))
        if not coord_df.crs == 'EPSG:4326':
            coord_df = coord_df.to_crs(epsg=4326)
        coord_df['latitude'] = coord_df.geometry.y
        coord_df['longitude'] = coord_df.geometry.x

        #Initialize impact matrix as dataframe
        imp_df = pd.DataFrame(index=coord_df.index)
        #Initialize matrix to store affected asset values
        aff_value_df = pd.DataFrame(index=coord_df.index)
        for date in dates:
            df_sel =  gdf.loc[gdf.date_dt==date,:]
            temp = df_sel.set_index(id_col)['value'].rename(date)
            #In rare cases there are 2 damages for the same house and date, 
            # then add them together
            if len(temp)!=len(np.unique(temp.index.astype(str))):
                temp = temp.groupby(level=0).sum()
            #add the damages of 'date' to the impact matrix
            imp_df=pd.concat([imp_df,temp],axis=1)

            #Add values to affected_asset_value matrix:
            if exp_val_col != None:
                temp = df_sel.set_index(id_col)[exp_val_col].rename(date)
                #In rare cases there are 2 damages for the same house and date,
                #just select the Versicherungsumme *1
                if len(temp)!=len(np.unique(temp.index.astype(str))):
                    temp = temp.groupby(level=0).first()
                aff_value_df=pd.concat([aff_value_df,temp],axis=1)

        imp_df = imp_df.fillna(0)
        imp_mat = sparse.csr_matrix(imp_df.T.values)


        imp_out = Impact()
        imp_out.coord_exp = coord_df[['latitude','longitude']].values
        imp_out.crs = coord_df.crs

        #Date and frequency
        ord_dates = np.array(dates.map(dt.datetime.toordinal))
        ev_names = np.array([date.strftime('ev_%Y-%m-%d') for date in dates])
        imp_out.date = ord_dates
        imp_out.event_name = ev_names
        n_years = np.ceil((max(ord_dates)-min(ord_dates))/365) #alternative: count unique years
        imp_out.frequency = np.ones(imp_df.shape[1])/n_years
        imp_out.event_id = np.arange(imp_df.shape[1])+1

        #Damages
        imp_out.imp_mat = imp_mat
        imp_out.at_event = imp_mat.sum(axis=1).getA1()
        imp_out.eai_exp = imp_mat.sum(axis=0).getA1()/n_years
        imp_out.aai_agg = imp_mat.sum()/n_years

        #Haz Tag information
        if haz_tag=='HL':
            imp_out.tag = {'haz':Tag('HL'), 'exp':E_tag(),'impf_set':E_tag()}

        #Exposed asset value
        if exp_val_col != None:
            aff_value_df = aff_value_df.fillna(0)
            aff_mat = sparse.csr_matrix(aff_value_df.T.values)
            imp_out.aff_mat = aff_mat
            #NOTE: aff_mat (affected exposure) is a matrix with the total value of affected assets
        #Others
        imp_out.unit = unit
        imp_out.affected_total_value = np.nan
        return imp_out

def read_xr_exposure(nc_file,var_name,val_unit='CHF'):
    """Read exposure from netCDF file"""
    if isinstance(nc_file,str):
        nc = xr.open_dataset(nc_file)[var_name]
    elif isinstance(nc_file,xr.Dataset):
        nc = nc_file[var_name]
    nc = nc.rename('value')
    df_exp = nc.to_dataframe()

    if 'lon' in df_exp.columns:
        df_exp = df_exp.rename(columns={'lon':'longitude','lat':'latitude'})
    gdf_exp = gpd.GeoDataFrame(df_exp,geometry = gpd.points_from_xy(df_exp.longitude,
                               df_exp.latitude,crs = 'EPSG:4326'))
    exp = Exposures(gdf_exp,value_unit = val_unit)
    exp.check()
    return exp

def read_xr_impact(nc_file,var_name,spatial_dims=['chx','chy'],time_dim='date',
                   unit='CHF',years=None):
    """Read impact from netCDF file"""

    if isinstance(nc_file,str):
        nc = xr.open_dataset(nc_file)[var_name]
    elif isinstance(nc_file,xr.Dataset):
        nc = nc_file[var_name]
    nc = nc.rename('value')
    stacked = nc.stack(new_dim=spatial_dims).fillna(0)

    # filter by years
    if years is not None:
        stacked = stacked.sel({time_dim:slice(str(years[0]),str(years[-1]))})

    # df_imp = nc.to_dataframe()
    n_ev = len(stacked[time_dim])
    n_years = len(np.unique(stacked[time_dim].dt.year))

    if 'lon' in stacked.coords:
        coord_exp = np.array([stacked.lat.values,stacked.lon.values]).T
        crs = 'EPSG:4326'
    else:
        raise NotImplementedError('Only lat/lon coordinates are supported')

    imp = Impact(
        event_id = np.arange(n_ev)+1,
        event_name = np.array([d.strftime('ev_%Y-%m-%d') for d in stacked[time_dim].dt.date.values]),
        date = np.array([d.toordinal() for d in stacked[time_dim].dt.date.values]),
        coord_exp = coord_exp,
        imp_mat = sparse.csr_matrix(stacked),
        crs = crs,
        eai_exp = stacked.sum(dim=[time_dim]).values/n_years,
        at_event = stacked.sum(dim='new_dim').values,
        frequency = np.ones(n_ev)/n_years,
        aai_agg = float(stacked.sum(dim=["new_dim",time_dim]).values)/n_years,
        unit = unit
    )

    return imp


def impf_from_csv(csv,smooth,emanuel_fit=False,PAA_only=False,plot=False,
                  return_impf=False):

    if smooth:
        mdr,mdd,paa= 'MDR_smooth','MDD_smooth','PAA_smooth'
    else:
        mdr,mdd,paa= 'MDR','MDD','PAA'
    if emanuel_fit:
        if smooth: raise ValueError('Cannot have smooth and emanuel_fit=True')
        mdr,mdd,paa='MDR_emanuel','MDD_emanuel','PAA_emanuel'

    df_impf = pd.read_csv(csv,index_col = 0)
    if PAA_only: #Impact function for PAA only (MDD=1)
        imp_fun1 = ImpactFunc()
        imp_fun1.id = 1
        imp_fun1.intensity = df_impf.index.values
        imp_fun1.mdd = np.ones_like(df_impf[paa].values)
        imp_fun1.paa = df_impf[paa]
        imp_fun1.haz_type = 'HL'
        imp_fun1.check()

        imp_fun_set = ImpactFuncSet()
        imp_fun_set.append(imp_fun1)
        if plot:
            ax = imp_fun1.plot()
            ax.set_ylim([0,max(df_impf[paa].values)*100*1.1])
    else:

        imp_fun1 = ImpactFunc()
        imp_fun1.id = 1
        imp_fun1.intensity = df_impf.index.values
        imp_fun1.mdd = df_impf[mdr].values
        imp_fun1.paa = np.ones_like(df_impf[mdr].values)
        imp_fun1.haz_type = 'HL'
        imp_fun1.check()

        # second impact function based on PAA*MDD Note that it has a bias 
        # because of expensive buildings are damaged more often!
        imp_fun2 = ImpactFunc()
        imp_fun2.id = 2
        imp_fun2.intensity = df_impf.index.values
        imp_fun2.mdd = df_impf[mdd].values
        imp_fun2.paa = df_impf[paa]
        imp_fun2.haz_type = 'HL'
        imp_fun2.check()
        if plot:
            ax = imp_fun1.plot()
            ax.set_ylim([0,max(df_impf[mdr].values)*100*1.1])
            imp_fun2.plot()

        imp_fun_set = ImpactFuncSet()
        imp_fun_set.append(imp_fun1)
        imp_fun_set.append(imp_fun2)
    return imp_fun_set