# -*- coding: utf-8 -*-
"""
Ploting functions for climada and non-climada objects used in scClim

"""

import numpy as np
import pandas as pd
import xarray as xr
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize,ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib import gridspec
import matplotlib.patheffects as PathEffects
from scipy.ndimage import gaussian_filter as g_filter
import datetime as dt
import cartopy.crs as ccrs
import cartopy.feature as cf
import cartopy.io.shapereader as shpreader
from climada.engine import Impact
from climada import CONFIG
from climada.util.constants import CMAP_IMPACT

from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Polygon
from cartopy.io import shapereader
import cartopy.io.img_tiles as cimgt
import cartopy.geodesic as cgeo
import geopandas
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
import colorcet as cc
cmap_imp = CMAP_IMPACT
cmap_imp.set_under('white',alpha=0)
from scClim.plot_cmaps import CMAP_IMPACT_CC, CMAP_IMPACT_CHF
from scClim.constants import CH_EXTENT_EPSG2056, CMAP_VIR

# %% General plotting functions
def plot_canton(ax,canton='Zürich',edgecolor = 'black',facecolor = "none",
                lakes=True,ch_border=None,zorder=0.5,lakeEdgeColor='none'):

    ch_shp_path = str(CONFIG.ch_shp_path)
    reader = shpreader.Reader("%s/swissTLMRegio_KANTONSGEBIET_LV95.shp"%ch_shp_path)

    if canton=='all':
        sel_cantons = [place for place in reader.records() if
                       place.attributes["OBJEKTART"]=='Kanton' and
                       place.attributes["ICC"]=='CH']
        for sel_canton in sel_cantons:
            shape_feature = cf.ShapelyFeature([sel_canton.geometry],
                                               ccrs.epsg(2056), edgecolor=edgecolor,
                                               facecolor = facecolor)
            ax.add_feature(shape_feature,zorder=zorder)
    elif type(canton)==list:
        sel_cantons = [place for place in reader.records() if
                       place.attributes["NAME"] in canton]
        for sel_canton in sel_cantons:
            shape_feature = cf.ShapelyFeature([sel_canton.geometry],
                                               ccrs.epsg(2056), edgecolor=edgecolor,
                                               facecolor = facecolor)
            ax.add_feature(shape_feature,zorder=zorder)
    elif type(canton)==str:
        sel_canton = [place for place in reader.records() if
                      place.attributes["NAME"]==canton][0]
        shape_feature = cf.ShapelyFeature([sel_canton.geometry],
                                           ccrs.epsg(2056), edgecolor=edgecolor,
                                           facecolor = facecolor)
        ax.add_feature(shape_feature)
    # add lake
    if lakes:
        reader2 = shpreader.Reader("%s/Hydrography/swissTLMRegio_Lake.shp"%ch_shp_path)
        geometry = reader2.geometries()
        geometry = np.array([g for g in geometry])
        lakesize = np.array([a.area for a in reader2.geometries()])
        geometry = geometry[lakesize>1e7] #default: 2e7
        shape_feature2 = cf.ShapelyFeature(geometry,
                                           ccrs.epsg(2056), edgecolor=lakeEdgeColor,
                                           facecolor = "lightblue")
        ax.add_feature(shape_feature2,zorder=zorder)#0.5,2
    if ch_border:

        reader3 = shpreader.Reader("%s/swissTLMRegio_LANDESGEBIET_LV95.shp"%ch_shp_path)
        sel_ctry = [place for place in reader3.records() if
                    place.attributes["ICC"]=='CH'][0]
        shape_feature3 = cf.ShapelyFeature([sel_ctry.geometry],
                                           ccrs.epsg(2056), edgecolor=ch_border,
                                           facecolor = "none")
        ax.add_feature(shape_feature3)
        shape_feature4 = cf.ShapelyFeature([sel_ctry.geometry],
                                           ccrs.epsg(2056), edgecolor='None',
                                           facecolor = 'white')
        ax.add_feature(shape_feature4,zorder=0)


def init_fig(projection=ccrs.PlateCarree,figsize=None,ch_border = True,title=''):

    if not figsize: figsize = (9,6)
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':projection})
    ax.set(title=title)
    plot_canton(ax,edgecolor="none",ch_border='grey')
    return fig, ax

def plot_nc(nc,ax=None,fillna='none',extent = None,title='',discrete='none',vmax=None,
            vmin=None,borders=True,cbar_lbl='',pl_type='',crs='default',canton=None,
            logScale=False,cbar_horizontal=False,extend=None,cbar=True,return_plot=False,
            border_color='black',**kwargs):
    """
     Parameters
    ----------
    nc : xr.DataArray
        with 1 timestep only
    """


    if crs == 'default': #PlateCaree()
        transform = ccrs.PlateCarree()
        x_dim = 'lon'
        y_dim = 'lat'
    elif crs=='EPSG:2056':
        transform = ccrs.epsg(2056)
        x_dim = 'chx'
        y_dim = 'chy'
    elif crs=='EPSG:21781':
        transform = ccrs.epsg(21781)
        x_dim = 'chx'
        y_dim = 'chy'
    elif type(crs) == ccrs.AlbersEqualArea: #crs given by rxr object
        transform = crs #nc.rio.crs
        x_dim = 'x'
        y_dim = 'y'

    if not fillna=='none':
        nc = nc.fillna(fillna)
    if not discrete=='none':
        nc = nc>discrete
    if not ax:
        fig,ax = plt.subplots(1,1,subplot_kw={'projection':transform})
    if extent: #
        ax.set_extent(extent,crs=transform)

    if logScale:
        vmin = vmin if vmin else nc.min().values
        vmax = vmax if vmax else nc.max().values
        vmin = max(1,vmin) #avoid log(0)
        vmax = max(1,vmax) #avoid empty norm
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        kwargs['norm'] = norm
        vmin = None; vmax = None
    else:
        norm = None

    if pl_type=='contour':
        plot = ax.contour(nc[x_dim],nc[y_dim],nc,transform=transform,
                          extend=extend,**kwargs)
    elif pl_type=='contourf':
        plot = ax.contourf(nc[x_dim],nc[y_dim],nc,transform=transform,
                           extend=extend,**kwargs)
    elif pl_type=='bool_field':
        #create bool colormap
        if 'cmap' in kwargs:
            cmap = ListedColormap(['white',kwargs['cmap']])
            # cmap = plt.cm.get_cmap(kwargs['cmap']).copy()
            cmap.set_under('white',alpha=0)
            kwargs['cmap'] = cmap
        if not nc.max()==0:
            plot = ax.pcolormesh(nc[x_dim],nc[y_dim],nc,vmax=vmax,
                                 transform=transform,vmin=0.5,**kwargs)
    else:
        plot = ax.pcolormesh(nc[x_dim],nc[y_dim],nc,vmax=vmax,vmin=vmin,
                             transform=transform,**kwargs)
        if cbar:
            add_cbar(ax,plot,cbar_lbl,extend=extend,horizontal=cbar_horizontal)

    if borders:
        ax.add_feature(cf.BORDERS,edgecolor=border_color)
    if title:
        try:
            ax.set_title(title+nc.name+nc.time)
        except:
            ax.set_title(title+nc.name)
    if canton:
        plot_canton(ax,canton=canton)

    if return_plot:
        return ax, plot
    else:
        return ax

def add_cbar(ax,plot,cbar_lbl='',extend=None,horizontal=False,pad=0.1):
        if horizontal:
            cbax = make_axes_locatable(ax).append_axes(
                'bottom', size="6.5%", pad=pad, axes_class=plt.Axes)
            cbar = plt.colorbar(plot, cax=cbax, orientation='horizontal',
                                extend=extend)
        else:
            cbax = make_axes_locatable(ax).append_axes(
                'right', size="6.5%", pad=pad, axes_class=plt.Axes)
            cbar = plt.colorbar(plot, cax=cbax, orientation='vertical',
                                extend=extend)
        cbar.set_label(cbar_lbl)

def cut_cmap(cmap, minval=0.0, maxval=1.0, n=-1,return_weak=False):

    if n == -1:
        n = cmap.N
    new_cmap = LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n))**1.1)
    new_cmap_light = LinearSegmentedColormap.from_list(
         'truncLight({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(np.linspace(minval, maxval, n))**0.3)
    if return_weak:
        return new_cmap,new_cmap_light
    else:
        return new_cmap

def scale_bar(ax,point,length,crs='EPSG:2056',dy_label=4e3,fontsize=10,line_kwargs=None):
    """plot scale bar on map

    Args:
        ax (plt.Axes): figure axes
        point (tuple): starting point of scale (in fig coordinates)
        length (int): length in coord units (m by default)
        crs (str, optional): CRS. Defaults to 'EPSG:2056'.
    """
    #get the geographic coordinates from the figure coordinates
    if not crs=='EPSG:2056':
        raise NotImplementedError('Only implemented for EPSG:2056')

    p_a_disp = ax.transAxes.transform(point)
    xy = ax.transData.inverted().transform(p_a_disp)

    scaleBar = Line2D((xy[0],xy[0]+length),(xy[1],xy[1]),
                      transform = ccrs.epsg(2056),color='black',**line_kwargs)
    ax.add_line(scaleBar)
    ax.text(xy[0]+length, xy[1]+dy_label, f'{length/1000:.0f} km',
            transform=ccrs.epsg(2056),fontsize=fontsize,ha='right',va='bottom')

def set_up_broken_axis(ylim1,ylim2,hspace=0.05,height_ratio=(1,3),figsize=None,
                       figAxAx=None):

    if figAxAx is None:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize,
                                    gridspec_kw={'height_ratios': height_ratio,'hspace':hspace})
    else:
        fig,ax1,ax2 = figAxAx
    # fig.subplots_adjust(hspace=hspace)  # adjust space between axes

    # plot the same data on both axes
    # ax1.plot(pts)
    # ax2.plot(pts)

    # zoom-in / limit the view to different portions of the data
    ax1.set_ylim(ylim2)  # outliers only
    ax2.set_ylim(ylim1)  # most of the data

    # hide the spines between ax and ax2
    ax1.spines.bottom.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

    # Now, let's turn towards the cut-out slanted lines.
    # We create line objects in axes coordinates, in which (0,0), (0,1),
    # (1,0), and (1,1) are the four corners of the axes.
    # The slanted lines themselves are markers at those locations, such that the
    # lines keep their angle and position, independent of the axes size or scale
    # Finally, we need to disable clipping.

    if ylim1[1] == ylim2[0] and hspace==0: #do not draw break if axis is continous
        ax2.axhline(ylim1[1],color='k',lw=1,linestyle='--')
    else:
        d = .5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                    linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
        ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
    return fig, (ax1, ax2)

# %% Specific plotting functions

########################### Subproject E #####################################
def plot_poh_sel(date,haz_poh,poh,exp_dmg,extent=None,poh_level=80):
    """
    Plot POH (+-2d) and reported damages per
    """
    raise DeprecationWarning('This function is deprecated. See event_definition.py for pre-processing plots')


def plot_event(date,imp,haz,act_dmgs,exp=None,dur_path=None,projection=None,extent=None,
    pl_type='default',canton=None,ch2056_filter = None,draw_grid=False,haz_var=None,
    vmax=0,vmax_haz=None,cmap_log=False,sel_dates=None,**kwargs):
    """Plot various plots per event

    Parameters
    ----------
    date : datetime
        date
    imp : Impact
        CLIMADA impact
    haz : Hazard
        CLIMADA hazard (MESHS)
    act_dmgs : dict of CLIMADA Exposure OR CLIMADA impact object
        Dictionary of actual damages read into a CLIMADA Exposure object
        The keys are the dates, formated as strings "ddmmYYYY"
    exp : CLIMADA Exposure
        Exposure
    dur_path : str
        path to duration files with YYYY=year and XXX = variable
    projection : cartopy.crs
        default = ccrs.PlateeCaree()
    extent : list
        [lon_min,lon_max,lat_min, lat_max]
    pl_type : str
        'default' to plot 1d, including duration
        'temp_ev' to ploat also +- 1 day, but no duration
        'sel_dates' to plot only the selected dates
    canton : str
        which (if any) canton borders to plot
    vmax: float
        maximum value for colorbar
        if vmax=0: set vmax as maximum of observed and modelled damages (default)
    Returns
    -------
    plot : matplotlib.figure.Figure
        Figure with plots
    """

    #Determine event IDs
    event_id = haz.event_id[haz.date==date.toordinal()]
    if not imp.event_id[haz.date==date.toordinal()] == event_id:
        ValueError("Event id is not consistent between Hazard and Impact")

    #set Projection
    if not projection:
        projection = ccrs.PlateCarree()

    #turn off borders if canton is plotted explicitly
    shapes = True # shapes = False if canton else True

    #determine unit
    assert(imp.unit==act_dmgs.unit)

    #Set vmin_for logatithmic impact scale
    if cmap_log:
        if imp.unit=='CHF' or imp.unit=='USD':
            vmin_log = 1e2
        elif imp.unit=='':
            vmin_log = 1
        else:
            raise ValueError('Impact Unit not recognized')

    #Initialize figure
    if pl_type == 'default':
        fig,axes = plt.subplots(2,3,gridspec_kw={'hspace':-0.1,'wspace':0.35},
                                subplot_kw={'projection':projection},figsize=(12,8))
        #plot duration
        if dur_path:
            for var,v_name, col in zip(['MZC','BZC'],['Duration of MESHS>20mm','Duration of POH>80%'],[1,2]):
                duration = dur_path.replace('YEAR',str(date.year)).replace('XXX',var)
                if os.path.exists(duration):
                    dur_da = xr.open_dataarray(duration).sel(time=date)*5
                    if ch2056_filter:
                        dur_da=dur_da.sel(chx=slice(ch2056_filter[0],ch2056_filter[1]),
                                          chy=slice(ch2056_filter[2],ch2056_filter[3]))
                    plot = axes[0,col].pcolormesh(dur_da.lon,dur_da.lat,dur_da)
                    axes[0,col].add_feature(cf.BORDERS)
                    axes[0,col].set_title(v_name)
                    add_cbar(axes[0,col],plot,'Duration (min)')
                    #plt.colorbar(plot,ax=axes[0,col],shrink = 0.6)
                    dur_da.close()
                else:
                    axes[0,col].annotate(var+ " duration files\ndo not exist",
                    xy=(0.1, 0.5), xycoords='axes fraction')
        #plot expsure
        if exp:
            exp.plot_hexbin(gridsize = 40,axis = axes[0,0],linewidths=0.05,extent=extent)
            axes[0,0].set_title('Exposure')
        #determine event ids
        ev_ids = [event_id]
        if isinstance(act_dmgs,Impact):
            ev_ids_observed = [act_dmgs.event_id[act_dmgs.date==date.toordinal()]]
        row = 1
        str_dates = [dt.datetime.strftime(date, "%d%m%Y")]

    elif pl_type == 'temp_ev':
        #plot Hazard at +- 1 day
        #haz.plot_intensity(event=event_id-1,axis=axes[0,0])
        #haz.plot_intensity(event=event_id+1,axis=axes[2,0])

        ev_ids = [imp.event_id[imp.date==(date+dt.timedelta(days=inc_day)).toordinal()]
                  for inc_day in [-1,0,1]]
        str_dates = [dt.datetime.strftime(date+dt.timedelta(days=inc_day), "%d%m%Y")
                     for inc_day in [-1,0,1]]
        if isinstance(act_dmgs,Impact):
            ev_ids_observed = [act_dmgs.event_id[act_dmgs.date==(date+dt.timedelta(days=inc_day)).toordinal()]
                               for inc_day in [-1,0,1]]

        #remove day +1 and day -1 if there is neither modelled nor reported damages
        for i in [2,0]:

            if isinstance(act_dmgs,Impact):
                rep_imp_is_zero = len(ev_ids_observed[i])==0
            elif isinstance(act_dmgs,dict):
                rep_imp_is_zero = str_dates[i] not in act_dmgs.keys()

            mod_imp_is_zero = (len(ev_ids[i])==0) or (imp.at_event[imp.event_id==ev_ids[i]][0] == 0)

            if mod_imp_is_zero and rep_imp_is_zero:
                del ev_ids[i]
                del str_dates[i]
                if isinstance(act_dmgs,Impact):
                    del ev_ids_observed[i]

        fig,axes = plt.subplots(len(ev_ids),3,gridspec_kw={'hspace':0.1,'wspace':0.3},
                                subplot_kw={'projection':projection},
                                figsize=(14,4*len(ev_ids)))
        if len(ev_ids)==1: axes=np.expand_dims(axes,axis=0)
        row = 0

    elif pl_type == 'sel_dates':
        assert(date in sel_dates)
        ev_ids = [haz.event_id[haz.date==date.toordinal()] for date in sel_dates]
        ev_ids_observed = [act_dmgs.event_id[act_dmgs.date==date.toordinal()] for date in sel_dates]
        fig,axes = plt.subplots(len(ev_ids),3,gridspec_kw={'hspace':-0.2,'wspace':0.3},
                                subplot_kw={'projection':projection},figsize=(14,4*len(ev_ids)))
        row = 0

    for i,ev_id in enumerate(ev_ids):

        #plot Hazard
        haz.plot_intensity(event=ev_id,axis=axes[row+i,0],vmin=1e-5,cmap = CMAP_VIR,vmax=vmax_haz)

        #set label
        if haz_var: fig.axes[-1].set_ylabel(haz_var)

        # gl = axes[row+i,0].gridlines(draw_labels=True,transform=projection,
          # #xlocs=lon_grid,ylocs=lat_grid,
          # x_inline=False,y_inline=False,
          # color='k',linestyle='dotted')

        #plot impact (only if it is available)
        if vmax==0:
            vmax_ = None
            n_loop = 2
        else:
            vmax_ = vmax
            n_loop = 1

        for loop in range(n_loop):


            #clear plots if looping through the second time
            if loop==1:
                axes[row+i,1].clear()
                axes[row+i,2].clear()
                if 'cax1' in locals(): cax1.remove()
                if 'cax2' in locals(): cax2.remove()
            try:
                if imp.at_event[imp.event_id==ev_id][0] == 0:
                    raise ValueError("No impact at this date -> continue to except statement")
                if cmap_log:
                    norm = colors.LogNorm(vmin=vmin_log, vmax=vmax_)
                else:
                    norm = colors.Normalize(vmin=1, vmax=vmax_)

                imp.plot_hexbin_impact_exposure(event_id=ev_id,gridsize=50,
                                                axis=axes[row+i,1],linewidths=0.05,
                                                norm=norm,extent=extent,cmap=CMAP_IMPACT_CC,
                                                shapes=shapes,**kwargs)
                axes[row+i,1].set_title(f'Modelled dmg:{imp.at_event[imp.event_id==ev_id][0]:.1e} {imp.unit}')
                cax1 = fig.get_children()[-1]
                vmax1= cax1.dataLim.get_points()[1,1]
            except:
                axes[row+i,1].annotate('Modelled damages are zero\nor unknown',
                xy=(0.1, 0.5), xycoords='axes fraction')
                vmax1=0
                #raise TypeError("imp argument is not CLIMADA Impact class")
                print("Warning: imp argument is not CLIMADA Impact class")
            try:
                ev_id_obs = ev_ids_observed[i]
                if cmap_log:
                    norm1 = colors.LogNorm(vmin=vmin_log, vmax=vmax_)
                else:
                    norm1 = colors.Normalize(vmin=1, vmax=vmax_)
                act_dmgs.plot_hexbin_impact_exposure(event_id=ev_id_obs,gridsize=50,
                                                     axis=axes[row+i,2],linewidths=0.05,
                                                     norm=norm1,extent=extent,cmap=CMAP_IMPACT_CC,
                                                     shapes=shapes,**kwargs)
                axes[row+i,2].set_title(f'Actual dmg:{act_dmgs.at_event[act_dmgs.event_id==ev_id_obs][0]:.1e} {act_dmgs.unit}')
                #TEST EQUAL AXIS
                cax2 = fig.get_children()[-1]
                vmax2= cax2.dataLim.get_points()[1,1]
            except:
                axes[row+i,2].annotate('Actual damages are zero\nor unknown',
                xy=(0.1, 0.5), xycoords='axes fraction')
                vmax2=0

            if  (vmax1!=0 and vmax2!=0 and min(vmax1,vmax2)/max(vmax1,vmax2)<0.1
                 and not cmap_log): #more than 1 OOM difference and not log scale
                vmax_ = max(vmax1,vmax2)/3
            else:
                vmax_ = max(vmax1,vmax2)

    if extent:
        for ax in axes.flatten():
            ax.set_extent(extent)
    if draw_grid:
        for ax in axes.flatten():
            gl=ax.gridlines(draw_labels=True)
            gl.xlabels_top = gl.ylabels_right = False
            gl.xlabel_style = gl.xlabel_style = {'color':'grey'}
    if canton:
        for ax in axes.flatten():
            plot_canton(ax,canton=canton)
    return fig,ax


#define function to convert numpy array data to xr format
def npz_t_xr(date,ncfile,var):
    str_date2 = date.strftime('%Y%m%d')
    npz = 'C:/Users/timo_/Documents/PhD/data/radar_dBZ/CZC_6_6_%s.npy'%str_date2
    arr = np.flip(np.load(npz),axis=[0])
    # arr = np.load(npz)
    ds = xr.Dataset({var: (("chy","chx"),arr)},
                    coords = ncfile.coords).drop('time')
    return ds

#Plotting function for model skill plot
def scatter_from_imp_df(imp_df,unit,xmin,dmg_thresh,eval_dict):
    """plotting function for skill plot from hail_main.py

    Args:
        imp_df (pd.DataFrame): dataframe with modelled and observed damages
        unit (str): unit of impacts
        xmin (float): minimum displayed values
        dmg_thresh (float): threshold used for the skill metrics calculation
        eval_dict (dict): Dictionary of evaluation metrics and hazard,
        exposure and impact info

    Returns:
        plt.Figure: Pyplot figure
    """

    #Set up figure
    fig,axes = plt.subplots(2,2,figsize = (7,7),
                            gridspec_kw={'height_ratios':[5,1],
                                         'width_ratios':[1,5]})
    ax = axes[0,1]
    ax_0model = axes[1,1]
    ax_0model.sharex(ax)
    ax_0rep = axes[0,0]
    ax_0rep.sharey(ax)
    fig.delaxes(axes[1,0])

    #Scatter plot wit label and annotations
    scat=ax.scatter(imp_df.imp_obs,imp_df.imp_modelled,marker='o',
                    facecolors = 'none', edgecolors = 'black')
    neither_zero = (imp_df.imp_modelled>0) & (imp_df.imp_obs>0)
    ax.annotate('n=%d\ndmg=%.1e%s'%(sum(neither_zero),
                                    imp_df.imp_modelled[neither_zero].sum(),unit),
                xy=(0.05, 0.9), xycoords='axes fraction',color='darkblue')

    #Plot events where the modelled damage is zero
    ax_0model.scatter(imp_df.imp_obs[imp_df.imp_modelled==0],
                      imp_df.imp_modelled[imp_df.imp_modelled==0],marker = 'o',
                      facecolors = 'none',edgecolors = 'black')
    tot_0model = sum(imp_df.imp_obs[imp_df.imp_modelled==0])
    str_0model = 'n=%d\ndmg=%.1e%s\n(%.2f%%)'%(sum((imp_df.imp_modelled==0)),
                                               tot_0model,unit,
                                               tot_0model/sum(imp_df.imp_obs))
    ax_0model.annotate(str_0model, xy=(0.68, 0.35), xycoords='axes fraction',
                       color='darkblue')

    #plot events where the reported damage is zero (false alarms)
    ax_0rep.scatter(imp_df.imp_obs[imp_df.imp_obs==0],
                    imp_df.imp_modelled[imp_df.imp_obs==0],marker = 'o',
                    facecolors = 'none',edgecolors = 'black')
    tot_0rep = sum(imp_df.imp_modelled[imp_df.imp_obs==0])
    str_0rep = 'n=%d\ndmg=%.1e%s\n(%.2f%%)'%(sum((imp_df.imp_obs==0)),tot_0rep,
                                             unit,tot_0rep/sum(imp_df.imp_modelled))
    ax_0rep.annotate(str_0rep, xy=(0.01, 0.85), xycoords='axes fraction',color='darkblue')

    ylim = [min(imp_df.imp_modelled[imp_df.imp_modelled!=0])*0.5,
            max(imp_df.imp_modelled)*2]
    ax.set(yscale = 'log',xscale = 'log',ylim = ylim)
    ax_0model.set(xscale = 'log',xlabel = 'reported dmgs [%s]'%unit,
                  xlim=[xmin,imp_df.imp_obs.max()*2])
    ax_0rep.set(yscale = 'log',ylabel = 'modeled dmgs [%s]'%unit,
                ylim=[xmin,imp_df.imp_obs.max()*2])
    for adj,col,lbl in zip([-1,0,1],['salmon','red','salmon'],
                           [None,'1:1 line','+- one order of\nmagnitude']):
        ax.plot([xmin,1e9],[10**(np.log10(xmin)+adj),10**(9+adj)],color = col,label = lbl)

    ax.legend(loc = 'lower right')
    if dmg_thresh:
        ax.add_patch(Rectangle((0, 0), dmg_thresh, dmg_thresh, color="grey", alpha=0.5))
        ax_0model.add_patch(Rectangle((0, ax_0model.get_ylim()[0]), dmg_thresh,
                                      dmg_thresh, color="grey", alpha=0.5))
        ax_0rep.add_patch(Rectangle((ax_0rep.get_xlim()[0], 0), dmg_thresh,
                                    dmg_thresh, color="grey", alpha=0.5))

    fig.suptitle(f'Hazard: {eval_dict["haz_var"]}, Exposure: {eval_dict["exposure"]}, Impf: {eval_dict["impf"]}\n\
                RMSE: {eval_dict["rmse"]:.2e},  RMSF: {eval_dict["rmsf"]:.2f},  weigthed RMSF: {eval_dict["rmsf_weighted"]:.2f}\n\
                POD: {eval_dict["POD"]*100:.1f}%, FAR: {eval_dict["FAR"]*100:.1f}%, within 1 OOM: {eval_dict["p_within_OOM"]*100:.1f}%')
    return fig


def plot_pre_process(exp_dmg,date,poh,haz_poh,ds_PH,event_index,plot_dBZ_dates,
                     poh_level,extent=None):

    #definitions
    poh_cols =['khaki','wheat','darkred','plum','lightsteelblue']
    lines = [mlines.Line2D([], [], color=color) for color in poh_cols]
    labels = ['D-2','D-1','D0','D+1','D+2']
    lines_pt = [mlines.Line2D([], [], color=color, marker=marker, linewidth=0)
                for color,marker in zip(['black','blue','red'],['.','.','*'])]
    labels_pt = ['keep','move','delete']

    #initialize gdf
    if type(exp_dmg) == geopandas.geodataframe.GeoDataFrame:
        exp_gdf = exp_dmg
    else:
        exp_gdf = exp_dmg.gdf
    # if not type(exp_dmg) == climada.entity.exposures.base.Exposures:



    #initialize figure
    fig_crs = ccrs.epsg(2056)
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':fig_crs},figsize = (10,6))

    #plot countours
    for i,day_adj in enumerate([-2,-1,0,1,2]):
        #check that the day exists in the hazard data
        if (0 <= event_index+day_adj < haz_poh.size and
            haz_poh.date[event_index]+ day_adj == haz_poh.date[event_index+day_adj]):

            ax.contour(poh.lon,poh.lat,poh.sel(time=date+pd.Timedelta(day_adj,'d')).fillna(0),
                       transform = ccrs.PlateCarree(),levels = [poh_level],colors = poh_cols[i])
            ax.contour(poh.lon,poh.lat,ds_PH.possible_hail.sel(time=date+pd.Timedelta(day_adj,'d')),
                       transform = ccrs.PlateCarree(),levels = [0.5],
                       colors = poh_cols[i],linestyles = '--',alpha = 0.8)

    shift_str = 'shifted' if 'shifted' in exp_gdf['category'].unique() else 'shif'
    color = exp_gdf['category'].map({'none':'red',shift_str:'blue','day0':'black'})
    size = exp_gdf['category'].map({'none':10, shift_str:1, 'day0':1})
    exp_gdf.to_crs(fig_crs).plot(ax=ax,color = color,#cmap='viridis',cax=cax,legend=True,
                      markersize = size,zorder=10,marker='*')

    if date.strftime('%d%m%Y') in plot_dBZ_dates:
        try:
            #contour_npz(date,poh,ax)
            ds=npz_t_xr(date,poh,"dBZ")
            ax.contourf(ds.lon,ds.lat,ds.dBZ.fillna(0),transform = ccrs.PlateCarree(),
                        levels = [35,40,45,100],colors = ['lightsalmon','red','darkred','none'],alpha=0.2)
        except:
            ValueError("cannot plot dBZ contours, check (hardcoded) path for files")

    if extent: ax.set_extent(extent,crs=fig_crs)
    ax.set_title('%s, #claims: %d, within POH%d: %d%%'%
                  (dt.datetime.strftime(date, "%Y-%m-%d"),exp_gdf.shape[0],poh_level,
                  sum(exp_gdf.POH>=poh_level)/exp_gdf.shape[0]*100))
    plot_canton(ax,canton=['Zürich','Bern','Luzern','Aargau'])
    fig.legend(lines,labels,ncol=2,bbox_to_anchor=(0.5,0.1),title = 'POH %d%%'%poh_level)
    fig.legend(lines_pt,labels_pt,ncol=1,bbox_to_anchor=(0.7,0.1),title ='Damage reports')
    return fig

def plot_haz_hist(haz,bins=30,density=False,**kwargs):

    fig,ax=plt.subplots(1,1)
    ax.hist((haz.intensity>0).sum(axis=1).getA1(),bins=bins,density=density,**kwargs)
    ax.set(yscale='log',xlabel='Area [km$^2$]',ylabel='count')
    if density:
        ax.set_ylabel('density')

    return fig


def plot_skill_plot(ds_KGV_per_event,PAA_lim,MDR_lim,PAA_thresh=1e2,
                    MDR_thresh=1e5,axis0lim=1e-2,case_studies=None,labels=None,
                    color2='tab:blue',gs_fig=None,ret_ax=False,exp_str=''):

    #define cmap
    cmap_post2012 = ListedColormap(["beige", color2])

    ds_KGV_per_event['post2012']=ds_KGV_per_event.date.dt.year>=2013

    #initialize figure
    if gs_fig is None:
        fig,axes = plt.subplots(2,5,figsize=(10,4),
                                gridspec_kw={'height_ratios':[1,0.15],
                                             'width_ratios':[0.15,1,0.3,0.15,1],
                                             'wspace':0.00,'hspace':0.00})
    else:
        gs_in, fig = gs_fig
        gs = gridspec.GridSpecFromSubplotSpec(2,5, subplot_spec=gs_in,wspace=0.0,
                                              hspace=0,width_ratios=[0.15,1,0.3,0.15,1],
                                              height_ratios=[1,0.15])
        axes=np.array([[0,0,0,0,0],[0,0,0,0,0]],dtype=object)
        for x in range(5):
            for y in range(2):
                axes[y,x]= fig.add_subplot(gs[y,x])

    for ax in axes[:,2]:
        ax.set_visible(False)
    axes1 = axes[:,:2]
    axes2 = axes[:,3:]
    ax1 = axes[0,1]
    ax2 = axes[0,4]
    ax3 = axes[1,1]

    for x,y,xlimFull,ylimFull,axes_sel,dmg_thresh in zip(['n_count','imp_observed'],
                                ['n_buildings_MESHS','imp_MESHS'],
                                [PAA_lim,MDR_lim,],
                                [PAA_lim,MDR_lim],
                                (axes1,axes2),
                                [PAA_thresh,MDR_thresh]):
        [[a,b],[c,d]] = axes_sel
        for ax,ylim,yscale,xlim,xscale in zip([a,b,c,d],
                            [ylimFull,ylimFull,(-axis0lim,axis0lim),(-axis0lim,axis0lim)],
                            ['log','log','linear','linear'],
                            [(-axis0lim,axis0lim),xlimFull,(-axis0lim,axis0lim),xlimFull],
                            ['linear','log','linear','log']):
            if ax == c:
                ax.set(xlim=xlim,ylim=ylim) #do not plot the 0,0 axis, as day with neither reported nor observed damage are not included
            else:
                add_legend = None# True if ax == b else False #only works with fillna(0)
                dsNow = ds_KGV_per_event.where((ds_KGV_per_event[x]>xlim[0])&(ds_KGV_per_event[x]<xlim[1]))
                dsNow = dsNow.where((ds_KGV_per_event[y]>ylim[0])&(ds_KGV_per_event[y]<ylim[1]))

                artist=dsNow.plot.scatter(x=x,y=y,hue='post2012',hue_style='discrete',
                                          add_colorbar=False,cmap=cmap_post2012,edgecolor='black',
                                          alpha=1,xscale=xscale,yscale=yscale,ylim=ylim,
                                          xlim=xlim,ax=ax,add_title=False,add_legend=add_legend)
        if case_studies is not None:
            #add case studies
            ds_KGV_per_event.sel(date=case_studies).plot.scatter(x=x,y=y,ax=b,marker='x',
                                                                 color='red',s=50,zorder=10,
                                                                 add_title=False)

        sc.E.make_cut_axis(a,b,c,d)
        #grey out dmg_tresh
        sc.E.grey_out_dmg_thresh(a,b,c,d,dmg_thresh=dmg_thresh)
        #add skill background colors
        sc.E.skill_background_cols(a,b,c,d,min_plot_val=xlimFull[0],
                                   max_plot_val=xlimFull[1],
                                   dmg_thresh=dmg_thresh,alpha=0.4)
        transFigure = fig.transFigure.inverted()
        # for startEnd in [(-0.1,1e9),]
        # for factor,col in zip([0.1,1,10],['grey','black','grey']):
        coord1 = transFigure.transform(b.transData.transform([xlimFull[1],xlimFull[1]]))
        coord2 = transFigure.transform(c.transData.transform([-axis0lim,-axis0lim]))
        line1 = Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                                    transform=fig.transFigure,
                                    ls='--',
                                    color='black')
        fig.lines.extend([line1])

    for factor,col in zip([0.1,10],['grey','grey']):
        ax1.axline(xy1=(1,1*factor),xy2=(10,10*factor), color=col,ls='--')
        ax2.axline(xy1=(1,1*factor),xy2=(10,10*factor), color=col,ls='--')

    if labels is not None:
        assert(len(labels)==2)
        for ax,label in zip([axes1[0,0],axes2[0,0]],labels):
            txt = ax.text(0.15,0.97,label,transform=ax.transAxes,fontsize=11,
                          fontweight='bold',ha='left',va='top')
            txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    #create legend
    if len(np.unique(ds_KGV_per_event['post2012']))>1: #only create if both colors are used
        red_patch = Line2D([],[],color='none',markerfacecolor='wheat',
                           markeredgecolor='black',marker='o', label='2002-2012')
        blue_patch = Line2D([],[],color='none',markerfacecolor='tab:blue',
                            markeredgecolor='black',marker='o', label='2013-2021')
        handles = [blue_patch, red_patch]
        if case_studies is not None:
            handles.append(Line2D([],[],color='none',markerfacecolor='red',
                                  markeredgecolor='red',marker='x', label='Case studies'))
        legend = axes[0,1].legend(handles=handles,loc='upper left')#bbox_to_anchor=(0.5, 0),ncol=1)

    #set labels
    axes[0,0].set_ylabel(f'Modelled No. of {exp_str} damages')
    axes[1,1].set_xlabel(f'Reported No. of {exp_str} damages')
    axes[0,3].set_ylabel(f'Modelled {exp_str} damage [CHF]')
    axes[1,4].set_xlabel(f'Reported {exp_str} damage [CHF]')
    if ret_ax:
        return fig, axes
    else:
        return fig


#impact function plots (from empirical calibration)
def fill_quantile(ds, var, q, col, ax,haz_var, binned=False, dim='b_sample',
                  alpha=0.5,cut_off=None,label=None):
    if cut_off is not None:
        ds=ds.where(ds[haz_var]<=cut_off)
    if label is None:
        label = '%s, Q %.2f-%.2f' % (var, q[0], q[1])

    ax.fill_between(ds[haz_var][1:], ds.quantile(q[0], dim=dim)[var][1:],
                    ds.quantile(q[1], dim=dim)[var][1:],  # step='pre',
                    color=col, alpha=alpha, label=label)

def fill_df_quantile(df,q1,q2,ax,label=None):
    q1_arr=df.quantile(q=q1,axis=1,numeric_only=True)*100
    q2_arr=df.quantile(q=q2,axis=1,numeric_only=True)*100
    if label is None:
        label = 'Q %.2f-%.2f' % (q1, q2)
    ax.fill_between(df.index,q1_arr,q2_arr,color='black',alpha=0.2,
                    label = label)

def impf_plot(df_all,df_roll,df_roll_cut,ds_boot_roll,ds_boot_roll_cut,haz_var,
              impf_emanuel,cut_off,plot_var='PAA',title='',dmg_bin_size=5,
              intensity_label = 'Intensity [?]',color='green'):

    intensity_range = df_roll.index.values
    xlim = [min(intensity_range)-2, max(intensity_range)]

    # plot rolling PAA
    paa_ylim = (0,max(df_roll[plot_var][np.isfinite(df_roll[plot_var])].max(),
                      df_roll.loc[cut_off,plot_var]*1.5)*100)
    fig, ax = plt.subplots()
    ax.plot(df_roll.index[1:], df_roll[plot_var].iloc[1:]*100, color=color, label=plot_var)
    if intensity_range[0] == 0:
        ax.scatter(df_roll.index[0], df_roll[plot_var][0]*100, color=color, alpha=0.5, marker='.')
    fill_quantile(ds_boot_roll*100, plot_var, (0.05, 0.95), f'light{color}', ax,haz_var)

    sc.E.plot_monotone_fit(df_roll_cut*100, plot_var, ax, ls='--')
    sc.E.plot_monotone_fit(ds_boot_roll_cut.quantile(0.95, dim='b_sample').to_dataframe()*100,
                           plot_var, ax, ls='--', color='grey', label='Q 0.05-0.95')
    sc.E.plot_monotone_fit(ds_boot_roll_cut.quantile(0.05, dim='b_sample').to_dataframe()*100,
                           plot_var, ax, ls='--', color='grey', label=None)
    ax.scatter(impf_emanuel.intensity,impf_emanuel.mdd*100,marker='o',
               label='Sigmoidal fit',edgecolor='black',facecolor='none',s=10)
    sc.E.plot_dmg_bin(df_all, ax, color='red', alpha=0.1,bin_size=dmg_bin_size)

    fig.legend(loc='upper left', bbox_to_anchor=ax.get_position())
    ax.set_title(title)
    ax.set(ylabel=f'{plot_var} [%]', xlim=xlim,ylim=paa_ylim, xlabel=intensity_label)
    ax.add_patch(Rectangle((cut_off, ax.get_ylim()[0]), cut_off*3,
                           ax.get_ylim()[1]-ax.get_ylim()[0],
                           color="grey", alpha=0.5))
    return fig,ax

def impf_plot2(df_all,df_roll,df_roll_cut,ds_boot_roll,ds_boot_roll_cut,haz_var,
               impf_emanuel,cut_off,plot_var='PAA',title='',dmg_bin_size=5,
               intensity_label = 'Intensity [?]',color='green',df_boot_emanuel=None,
               quantile = (0.05,0.95),relative_bins=False,figAx=None):

    intensity_range = df_roll.index.values
    xlim = [min(intensity_range)-2, max(intensity_range)]
    # plot rolling PAA
    paa_ylim = (0,max(df_roll[plot_var][np.isfinite(df_roll[plot_var])].max(),
                      ds_boot_roll[plot_var].sel({haz_var:cut_off}).quantile(quantile[1]),#df_roll.loc[cut_off,plot_var]*1.5,
                      impf_emanuel.mdd.max())*100)
    if figAx is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = figAx
    # density plot of bootstrap samples
    iloc_cut_off = np.where(df_roll.index.values==cut_off)[0][0]
    ax.plot(df_roll.index[1:iloc_cut_off+1],
            df_roll[plot_var].iloc[1:iloc_cut_off+1]*100,
            color=color, label=plot_var)
    fill_quantile(ds_boot_roll*100, plot_var, quantile, f'light{color}', ax,
                  haz_var,alpha=0.8,cut_off=cut_off,label='bootstrap CI (5-95%)')

    ax.plot(impf_emanuel.intensity,impf_emanuel.mdd*100,label='Sigmoidal fit',color='black')
    if df_boot_emanuel is not None:
        fill_df_quantile(df_boot_emanuel,0.05,0.95,ax,label='bootstrap CI (5-95%)')
    ax2=sc.E.plot_dmg_bin(df_all, ax, color='tomato', alpha=0.2,
                          bin_size=dmg_bin_size,relative=relative_bins)

    #plot exposed assets
    # sc.E.plot_dmg_bin(df_all, ax, color='blue', alpha=0.1,new_axis=False,bin_size=dmg_bin_size,pl_var='count_all')
    # ax.set(xlim=xlim, xlabel=intensity_label)
    # ax.set(ylabel=f'Proportion of affected assets [%]', xlim=xlim,ylim=paa_ylim, xlabel=intensity_label)
    hdl,lbl = ax.get_legend_handles_labels()
    hdl2,lbl2 = ax2.get_legend_handles_labels()
    fig.legend(hdl+hdl2,lbl+lbl2,loc='upper left', bbox_to_anchor=ax.get_position(),
               framealpha=0.5)
    ax.set_title(title)
    ax.set(ylabel=f'{plot_var} [%]', xlim=xlim,ylim=paa_ylim, xlabel=intensity_label)
    ax.add_patch(Rectangle((cut_off, ax.get_ylim()[0]), cut_off*3,
                           ax.get_ylim()[1]-ax.get_ylim()[0], color="grey",
                           alpha=0.7))
    return fig,ax


def plot_impf_pixel(ds_roll_gc,df_roll,ds_boot_roll,haz_var,cut_off,q1,q2,
                    impf_emanuel=None,plot_var='PAA',title='',color='green',
                    intensity_label='Intensity [?]',pl_type=1):
    ds1 = get_emp_quantile(ds_roll_gc, q1)
    ds2 = get_emp_quantile(ds_roll_gc, q2)
    fig, ax = plt.subplots()

    #Type 1 plot
    if pl_type == 1:
        ax.plot(df_roll.index[1:], df_roll[plot_var][1:]*100, color=color,
                label=plot_var)

        ax.fill_between(ds1[haz_var], ds1[plot_var]*100, ds2[plot_var]* 100,
                        color='orange', alpha=0.3,
                        label='pixel-wise Q%d-Q%d' % (q1*100, q2*100))
        fill_quantile(ds_boot_roll*100, plot_var, (0.05, 0.95), f'light{color}',
                      ax,haz_var,alpha=0.9)

        weighted_PAA = ds_roll_gc[['PAA', 'MDR']].weighted(weights=ds_roll_gc.exp_val.fillna(0)).mean(dim='index')
        ax.plot(df_roll.index[1:], weighted_PAA[plot_var]*100, color='black', label=plot_var)


    #Type 2 plot
    elif pl_type == 2:
        ax.plot(impf_emanuel.intensity,impf_emanuel.mdd*100,
                label='Sigmoidal fit',color='black')

        dsQ0 = get_emp_quantile(ds_roll_gc, 0)
        for quantile in np.flip([0.8,0.85,0.9,0.95]):
            color = (1-(quantile-0.7)*3,0,(quantile-0.7)*3)
            dsQ = get_emp_quantile(ds_roll_gc, quantile)
            ax.fill_between(dsQ0[haz_var], dsQ0[plot_var]*100, dsQ[plot_var]
                            * 100,color=color, alpha=1, label='pixel-wise Q%d' % (quantile*100))

        weighted_PAA = ds_roll_gc[['PAA', 'MDR']].weighted(weights=ds_roll_gc.exp_val.fillna(0)).mean(dim='index')

    #Type 3 plot
    elif pl_type == 3:
        dsMean = get_emp_mean(ds_roll_gc)
        ax.plot(impf_emanuel.intensity,impf_emanuel.mdd*100,
                label='Sigmoidal fit',color='black')
        ax.plot(df_roll.index[1:], df_roll[plot_var][1:]*100, color=color,
                label=plot_var)
        ax.fill_between(ds1[haz_var], dsMean[plot_var]*100, dsMean[plot_var]*100,
                        color=color, alpha=1, label='pixel-wise mean')

        ax.fill_between(ds1[haz_var], dsMean[plot_var]*100/2, dsMean[plot_var]*100*2,
                        color=color, alpha=0.5, label='pixel-wise mean')

        weighted_PAA = ds_roll_gc[['PAA', 'MDR']].weighted(weights=ds_roll_gc.exp_val.fillna(0)).mean(dim='index')
        ax.plot(df_roll.index[1:], weighted_PAA[plot_var]*100, color='black', label=plot_var)


    ax.legend()
    ax.add_patch(Rectangle((cut_off, 0), cut_off*3, ax.get_ylim()[1], color="grey", alpha=0.5))
    ax.set(title=title,ylabel=plot_var,ylim=[0,min(100,ax.get_ylim()[1])],xlabel=intensity_label)
    return fig,ax

def get_emp_quantile(ds_roll_gc, q):
    ds_roll_gc['PAA'] = ds_roll_gc.n_dmgs/ds_roll_gc.n_exp
    ds_roll_gc['MDR'] = ds_roll_gc.dmg_val/ds_roll_gc.exp_val
    ds = ds_roll_gc[['PAA', 'MDR']]
    quant = ds.quantile(q, dim='index')
    return quant

def get_emp_mean(ds_roll_gc):
    ds_roll_gc['PAA'] = ds_roll_gc.n_dmgs/ds_roll_gc.n_exp
    ds_roll_gc['MDR'] = ds_roll_gc.dmg_val/ds_roll_gc.exp_val
    ds = ds_roll_gc[['PAA', 'MDR','exp_val']]
    mean = ds.weighted(weights=ds_roll_gc.exp_val.fillna(0)).mean(dim='index')
    return mean