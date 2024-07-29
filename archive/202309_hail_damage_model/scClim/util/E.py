""" 
Utility function for subproject E (Timo)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import pandas as pd
from climada.engine import Impact
from scipy import sparse
import geopandas as gpd
from shapely.geometry import Point, MultiPoint
import random
import cartopy.crs as ccrs
import sys
from climada import CONFIG
sys.path.append(str(CONFIG.local_data.func_dir))
import scClim as sc
import xarray as xr


def get_hazard_files_TS(haz_var,years,event_def_version,data_dir):
    """get Hazard files (specific to subproject E)

    Args:
        haz_var (str): hazard
        years (np.array): considered years
        event_def_version (int): Event definition version
        data_dir (str,path): data directory

    Raises:
        NotImplementedError: for haz_var= HKE

    Returns:
        paths (paths or xr.Dataset): Data (or path) to use in get_hazard_from_radar()
    """
    #Local paths are not shown in public repository


def calc_and_save_skill_scores(imp_df,dmg_thresh,imp_var,expHazImp_str):
    """
    Calculate skill scores and save them to a csv file

    Args:
        imp_df (pd.DataFrame): modelled impact, observed impact
        dmg_thresh (float): damage threshold
        imp_var (str): name of impact variable (MDR or PAA)
        expHazImp_str (str): String of Exposure Hazard ImpactFunction

    Returns:

    """
    imp_var_str = imp_var.replace('MDR','').replace('PAA','_PAA') #remove MDR from imp_var string (as it is the default)
    #read df_scores from csv
    df_scores = pd.read_csv(str(CONFIG.local_data.out_dir)+'/skill_plots/_scores.csv',index_col=0)

    for y_filter in ['_2012plus','']:
        if y_filter=='_2012plus':
            df_now=imp_df.loc[imp_df.prePost2012,:]
        else:
            df_now=imp_df.copy()

        #caculate skill scores
        rmse,rmsf,rmsf_weighted,FAR,POD,p_within_OOM,considered_events = calc_skill_scores(df_now,dmg_thresh)

        #Save error_measures in csv
        df_scores.loc[f'{expHazImp_str}{imp_var_str}_t{dmg_thresh:.0e}{y_filter}',:] = [
            rmse,rmsf,rmsf_weighted,FAR,POD,p_within_OOM,considered_events]
        if y_filter=='':
            #save df_scores to csv
            df_scores= df_scores.sort_index()
            df_scores.to_csv(str(CONFIG.local_data.out_dir)+'/skill_plots/_scores.csv')
            return rmse,rmsf,rmsf_weighted,FAR,POD,p_within_OOM,considered_events
        
def calc_skill_scores(df_now,dmg_thresh):
        
        #get number of points
        considered_events = sum((df_now.imp_modelled>dmg_thresh) | (df_now.imp_obs>dmg_thresh))

        neither_is_zero = (df_now.imp_modelled!=0) & (df_now.imp_obs!=0) 
        rmse = np.sqrt(np.mean((df_now.imp_modelled[neither_is_zero]-df_now.imp_obs[neither_is_zero])**2))
        rmsf = np.exp(np.sqrt(np.mean(np.log(np.divide(
            df_now.imp_modelled[neither_is_zero],df_now.imp_obs[neither_is_zero]))**2)))
        rmsf_weighted = np.exp(np.sqrt(np.mean(np.log(np.divide(
            df_now.imp_modelled[neither_is_zero],df_now.imp_obs[neither_is_zero]))**2
            *df_now.imp_obs[neither_is_zero])/np.sum(df_now.imp_obs[neither_is_zero])))
        
        above_either_thresh = (df_now.imp_obs>dmg_thresh) | (df_now.imp_modelled>dmg_thresh)
        hits = sum(above_either_thresh & (df_now.imp_obs>=df_now.imp_modelled/10) & 
                   (df_now.imp_obs<df_now.imp_modelled*10))
        misses = sum((df_now.imp_obs>dmg_thresh) & (df_now.imp_modelled<df_now.imp_obs/10))
        false_alarms = sum((df_now.imp_modelled>dmg_thresh) & (df_now.imp_modelled>df_now.imp_obs*10) )

        FAR = (sum((df_now.imp_modelled>dmg_thresh) & (df_now.imp_modelled>df_now.imp_obs*10)))/sum((df_now.imp_modelled>dmg_thresh))
        POD = sum((df_now.imp_obs>dmg_thresh) & (df_now.imp_obs>df_now.imp_modelled/10) &
                  (df_now.imp_obs<df_now.imp_modelled*10))/sum((df_now.imp_obs>dmg_thresh))
        # print(f"POD: {POD}, POD2 = {hits/(hits+misses)}")
        # print(f"FAR: {FAR}, FAR2 = {false_alarms/(hits+false_alarms)}")
        p_within_OOM = sum(above_either_thresh & (df_now.imp_obs>df_now.imp_modelled/10) &
                           (df_now.imp_obs<df_now.imp_modelled*10))/sum(above_either_thresh)
        
        return rmse,rmsf,rmsf_weighted,FAR,POD,p_within_OOM,considered_events


def random_points_in_polygon(number, polygon):
    """Assign random points to polygon

    Args:
        number (int): number of points to assign
        polygon (shapely.polygon): polygon

    Returns:
        Point/MultiPoint: radnomly assingned points
    """
    points = []
    min_x, min_y, max_x, max_y = polygon.bounds
    i= 0
    while i < number:
        point = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y))
        if polygon.contains(point):
            points.append(point)
            i += 1
    if number==1:
        return points[0] 
    if number>1:
        return MultiPoint(points)# returns list of shapely point

def poly_to_pts(gdf,n_pts,epsg,mode='random'):
    """assign random points to geodataframe with polgon geometry

    Args:
        gdf (gpd.GeoDataFrame): gdf
        n_pts (str): 'always_one' to calculate 1 point per polygon
                    otherwise name of column that stores the number 
                    of points to calculate (e.g. n_claims)
        epsg (int): epsg code
        mode (str, optional): Defaults to 'random'.
"""
    gdf_out = gdf.copy()
    #assert that index corresponds to row number
    gdf_out = gdf_out.reset_index().drop(columns='index')

    if not isinstance(gdf_out,gpd.GeoDataFrame):
        gdf_out = gpd.GeoDataFrame(gdf_out,geometry='geometry')
        gdf_out = gdf_out.set_crs(epsg=epsg)


    if n_pts == 'always_one':
        gdf_out['n_points'] = 1
    # elif n_pts == 'n_claims': #n_claims is given at 'c_claims' columns
    #     print('number of points given by n_claims column')
    else:
        gdf_out['n_points'] = gdf_out[n_pts]

    
    if mode=='random':
        rand_points = [random_points_in_polygon(gdf_out.loc[row,'n_points'],gdf_out.loc[row,'geometry'])
                       for row in range(gdf_out.shape[0])]
        gdf_out['random_points'] = rand_points
        gdf_out=gdf_out.set_geometry('random_points')
        gdf_out=gdf_out.drop(columns=['geometry'])
    else:
        NotImplementedError(f'mode {mode} not implemented')
    return(gdf_out)


def filter_imp(imp,sel_ev):
    """Filter impact object by certain events

    Parameters
    ----------
    imp : climada.impact
        imact object to filter
    sel_ev : np.array, dtype=bool
        boolean array of length of event_id with the events to keep

    Returns
    -------
    imp_filtered : climada.impact
    
    """       
    #sel_ev = np.isin(imp.date, dates)
    imp_out = Impact()
    
    for (var_name, var_val) in imp.__dict__.items():
        if var_name == 'eai_exp':
            setattr(imp_out, var_name, var_val) #will be adjusted at end of function
        elif isinstance(var_val, np.ndarray) and var_val.ndim == 1 \
                and var_val.size > 0:
            print(var_name)
            
            setattr(imp_out, var_name, var_val[sel_ev])
        elif isinstance(var_val, sparse.csr_matrix):
            setattr(imp_out, var_name, var_val[sel_ev, :][:, :])
        elif isinstance(var_val, list) and var_val:
            setattr(imp_out, var_name, [var_val[idx] for idx in sel_ev])
        elif var_name == 'centroids':
            setattr(imp_out, var_name, var_val)
        else:
            setattr(imp_out, var_name, var_val)
        
        #define number of years only if frequency is equal, otherwise it will throw an error
        if len(np.unique(imp.frequency))==1: n_years = 1/imp.frequency[0]
        
        imp_out.eai_exp = imp_out.imp_mat.sum(axis=0).getA1()/n_years
        imp_out.aai_agg = imp_out.eai_exp.sum()
    return imp_out


def cut_xarr(da, extent = [5.8, 10.6, 45.7, 47.9]):
    """cut netcdf xarray object to rectangular lat/lon

    Parameters
    ----------
    da : xarray.Dataset or xarray.DataArray
        dataset
    extent : list or str
        lon_min, lon_max, lat_min, lat_max, default = Switzerland

    Returns
    -------
    da : xarray.DataArray, xarray.Dataset
    
    """    
    lon_min, lon_max, lat_min, lat_max = extent
    lon_cond = np.logical_and(da.lon >= lon_min, da.lon <= lon_max)
    lat_cond = np.logical_and(da.lat >= lat_min, da.lat <= lat_max)
    cut_da = da.where(np.logical_and(lat_cond,lon_cond),drop=True)    
    return cut_da
    
    
def plot_dmg_bin(df,ax,bin_size=5,ymin=0,color='green',alpha=0.3,new_axis=True,
                 pl_var = 'count_dmg',ymax=None,relative=False,**kwargs):
    
    #if 0 in df.index:
    bins = pd.cut(df.index,right=False,
                  bins = np.arange(min(df.index),max(df.index)+bin_size,bin_size))
    bin_mids = bins.categories.mid.values
    bin_right = bins.categories.right.values
    bin_l = bins.categories.left.values

    if df.index[0]==0 and df.index[1]>bin_size: #first index is 0 and only contains 0 values (no rolling!)
        bin_mids[0] = 0 #change to zero
        bin_l[0] -= bin_size/2 #center first bin at 0
    df_bins=df.groupby(bins).sum()
    if relative:
        df_bins[pl_var] = df_bins[pl_var]/df_bins[pl_var].sum()
    if new_axis:
        ax2 = ax.twinx()
    else:
        ax2 = ax

    ax2.bar(bin_l,df_bins[pl_var],width=bin_size,color = color,label = '#damages',
            alpha=alpha,align='edge',**kwargs)
    
    if pl_var == 'count_all': #create axis break
        ylabel = "Number of exposed assets"
        ymax = max(df_bins[pl_var][1:])*1.3 if ymax is None else ymax
        # ax2.set_ylim([ymin,ymax])
        ax2.annotate(f'Out of axis: {(df_bins[pl_var][0]):.1e}',xy=(bin_mids[0],ymax),
                     xytext=(bin_right[0]*1.1,ymax*0.9),
                     arrowprops=dict(facecolor=color,shrink=0.05),color=color)

    elif pl_var == 'count_dmg':
        ylabel = "Fraction of damage reports" if relative else "Number of damage reports"
        ymax=max(df_bins[pl_var])*1.1 if ymax is None else ymax
    elif pl_var == 'count_cells':
        ylabel = "Number of grid cells"
        ymax=max(df_bins[pl_var])*1.1 if ymax is None else ymax
    
    ax2.set(ylabel=ylabel,ylim=[ymin,ymax])
    return ax2

def plot_monotone_fit(df,var,ax,color='black',label='monotonic fit',**kwargs):

    if df.index[0]==0:# skip_zero
        assert(np.where(df.index==0) == np.array([0]))
        x=df.index[1:]
        y=df[var][1:]
    else:
        x=df.index
        y=df[var]
    #avoid nan values
    no_nan = ~np.isnan(x) & ~np.isnan(y)
    x = x[no_nan]
    y = y[no_nan]
    
    monotone_fit = sc.smooth_monotonic(x,y,plot=False)
    ax.plot(x, monotone_fit, color=color,label = label,**kwargs)
          
    
####### plotting function for specific plots
#scatter plot of modelled vs observed damages
def make_cut_axis(a,b,c,d):

    #top bottom separation
    a.spines.bottom.set_visible(False)
    b.spines.bottom.set_visible(False)
    c.spines.top.set_visible(False)
    d.spines.top.set_visible(False)
    a.xaxis.tick_top()
    b.xaxis.tick_top()
    a.tick_params(labeltop=False) # don't put tick labels at the top
    b.tick_params(labeltop=False) # don't put tick labels at the top
    c.xaxis.tick_bottom()
    d.xaxis.tick_bottom()

    #left right separation
    a.spines.right.set_visible(False)
    c.spines.right.set_visible(False)
    b.spines.left.set_visible(False)
    d.spines.left.set_visible(False)
    
    a.yaxis.tick_left()
    c.yaxis.tick_left()
    b.yaxis.tick_right()
    d.yaxis.tick_right()
    b.tick_params(labelleft=False) # don't put tick labels left
    d.tick_params(labelleft=False,labelright=False) # don't put tick labels left

    # Now, let's turn towards the cut-out slanted lines.
    # We create line objects in axes coordinates, in which (0,0), (0,1),
    # (1,0), and (1,1) are the four corners of the axes.
    # The slanted lines themselves are markers at those locations, such that the
    # lines keep their angle and position, independent of the axes size or scale
    # Finally, we need to disable clipping.

    v_h = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -v_h), (1, v_h)], markersize=10,
                linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    a.plot([0], [0], transform=a.transAxes, **kwargs)
    b.plot([1], [0], transform=b.transAxes, **kwargs)
    c.plot([0], [1], transform=c.transAxes, **kwargs)
    d.plot([1], [1], transform=d.transAxes, **kwargs)

    kwargs = dict(marker=[(-v_h,-1), (v_h,1)], markersize=10,
                linestyle="none", color='r', mec='k', mew=1, clip_on=False)
    a.plot([1], [1], transform=a.transAxes, **kwargs)
    b.plot([0], [1], transform=b.transAxes, **kwargs)
    c.plot([1], [0], transform=c.transAxes, **kwargs)
    d.plot([0], [0], transform=d.transAxes, **kwargs)

    #additional changes
    a.set(xlabel='',xticks=[0])
    b.set(xlabel='',ylabel='',yticks=[0])
    c.set(xlabel='',ylabel='',yticks=[0],xticks=[0])
    d.set(ylabel='',yticks=[0])

def grey_out_dmg_thresh(a,b,c,d,dmg_thresh=100,alpha=0.4):
    b.add_patch(Rectangle((0, 0), dmg_thresh, dmg_thresh, color="grey", alpha=alpha))
    d.add_patch(Rectangle((0, d.get_ylim()[0]), dmg_thresh, dmg_thresh, color="grey", alpha=alpha))
    a.add_patch(Rectangle((a.get_xlim()[0], 0), dmg_thresh, dmg_thresh, color="grey", alpha=alpha))    
    c.add_patch(Rectangle((c.get_xlim()[0], c.get_ylim()[0]), (c.get_ylim()[1]- c.get_ylim()[0]),
                          c.get_xlim()[1]- c.get_xlim()[0], color="grey", alpha=alpha))

def skill_background_cols(a,b,c,d,min_plot_val,max_plot_val,dmg_thresh=100,alpha=0.3,
                          true_pos='tab:green',false_alarm='tab:red',missed_event='tab:orange'):
    """add background colors for skill score metrics

    Args:
        a-d (axes): axes to add background colors to
        min_plot_val (int): minimum value of "main" plot (axis b)
        dmg_thresh (int, optional): _description_. Defaults to 100.
    """

    #false alarms
    a.add_patch(Rectangle((a.get_xlim()[0], dmg_thresh), 1, 1e18, color=false_alarm, alpha=alpha,zorder=0))
    b.add_patch(Polygon([(min_plot_val,dmg_thresh), (dmg_thresh*0.1, dmg_thresh),
                         (max_plot_val*0.1, max_plot_val),(min_plot_val,max_plot_val)],
                         color=false_alarm, alpha=alpha,zorder=0))
    
    #true positives
    b.add_patch(Polygon([(dmg_thresh, dmg_thresh*0.1), (max_plot_val, max_plot_val*0.1),
                         (max_plot_val, max_plot_val),(max_plot_val*0.1, max_plot_val),
                         (dmg_thresh*0.1, dmg_thresh), (dmg_thresh, dmg_thresh)],
                         color=true_pos, alpha=alpha,zorder=0))

    #missed events
    d.add_patch(Rectangle((dmg_thresh, d.get_ylim()[0]), 1e10, 1, color=missed_event, alpha=alpha,zorder=0))
    b.add_patch(Polygon([(dmg_thresh, min_plot_val), (max_plot_val, min_plot_val),
                         (max_plot_val, max_plot_val*0.1),(dmg_thresh, dmg_thresh*0.1)],
                         color=missed_event, alpha=alpha,zorder=0))