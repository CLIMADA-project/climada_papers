"""
"""
import os
import pickle
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from climada.hazard import TropCyclone, Centroids
from climada.util.save import save
import climada.util.plot as u_plot

MS2KN = 1.9438444924406049
""" m/s to knots factor """

def compute_tc(data_dir, resol):
    """ compute tropical cyclone over all area """
    try:
        abs_path = os.path.join(data_dir, 'sel_hist_syn_1h.p')
        with open(abs_path, 'rb') as f:
            TRACKS = pickle.load(f)
        print('Read', TRACKS.size, 'tracks.')
    except FileNotFoundError:
        print('ERROR: Execute script_exec.py first.')

    # CENTROIDS
    min_lat, max_lat, min_lon, max_lon = 16.99375, 21.95625, -72.48125, -61.66875
    delta = 1.0
    min_lat, max_lat, min_lon, max_lon = min_lat - delta, max_lat + delta, \
                                         min_lon - delta, max_lon + delta

    cent = Centroids()
    cent.coord = (np.mgrid[min_lat : max_lat : complex(0, resol), 
                           min_lon : max_lon : complex(0, resol)]). \
                  reshape(2, resol*resol).transpose()
    cent.id = np.arange(cent.lat.size)
    
    # TROP CYCLONE
    tc = TropCyclone()
    tc.set_from_tracks(TRACKS, cent)
    tc.check()
    save(os.path.join(data_dir, 'tc_all.p'), tc)
    
    return tc

def compute_tc_rp(data_dir, tc):
    """ compute local return periods """
    rp_inten = tc.local_exceedance_inten()
    save(os.path.join(data_dir, 'tc_all_rp_local.p'), rp_inten)
    return rp_inten

def fig05(data_dir, fig_dir=None):
    
    resol = 250

    try:
        abs_path = os.path.join(data_dir, 'tc_all.p')
        with open(abs_path, 'rb') as f:
            tc = pickle.load(f)
        print('Loaded tc_all:', tc.size)
    except FileNotFoundError:
        tc = compute_tc(data_dir, resol)
        
    try:
        abs_path = os.path.join(data_dir, 'tc_all_rp_local.p')
        with open(abs_path, 'rb') as f:
            tc_rp = pickle.load(f)
        print('Loaded tc_all_rp_local:', tc_rp.size)
    except FileNotFoundError:
        tc_rp = compute_tc_rp(data_dir, tc)

    
    f, axes = plt.subplots(2, 2, figsize=(15, 20), 
                           subplot_kw=dict(projection=ccrs.PlateCarree()),
                           squeeze=True, sharex=True, sharey=True)
    for idx, axis in enumerate(axes.flatten()):
        grid = axis.gridlines(draw_labels=True, alpha=0.2)
        grid.xlabels_top = grid.ylabels_right = False
        if idx == 0:
            grid.xlabels_bottom = False
            grid.yformatter = LATITUDE_FORMATTER
        elif idx == 1:
            grid.xlabels_bottom = False
            grid.ylabels_left = False
        elif idx == 2:
            grid.yformatter = LATITUDE_FORMATTER
            grid.xformatter = LONGITUDE_FORMATTER
        elif idx == 3:
            grid.ylabels_left = False
            grid.xformatter = LONGITUDE_FORMATTER    
    
    extent = u_plot._get_borders(tc.centroids.coord)
    tc_rp *= MS2KN
    vmin, vmax = np.min(tc_rp), np.max(tc_rp)
    for array_im, axis, rp in zip(tc_rp, axes.flatten(), [25, 50, 100, 250]):
        axis.set_extent((extent))
        u_plot.add_shapes(axis)

        grid_x, grid_y = np.mgrid[
                extent[0] : extent[1] : complex(0, resol),
                extent[2] : extent[3] : complex(0, resol)]
        grid_im = griddata((tc.centroids.coord[:, 1], tc.centroids.coord[:, 0]),
                           array_im, (grid_x, grid_y)).squeeze()
        im = axis.pcolormesh(grid_x, grid_y, grid_im, vmin=vmin, vmax=vmax, 
                             transform=ccrs.PlateCarree(), cmap='Blues')
        axis.text(-63.1, 22, 'RP = ' + str(rp), fontsize=12)
        if rp == 25:
            axis.text(-70.7, 22.3, 'Turks and',
                      horizontalalignment='left', verticalalignment='center',
                      transform=ccrs.PlateCarree(), fontsize=14)
            axis.text(-70.7, 21.7, 'Caicos islands',
                      horizontalalignment='left', verticalalignment='center',
                      transform=ccrs.PlateCarree(), fontsize=14)
            axis.text(-64.8, 19.2, 'Lesser Antilles',
                      horizontalalignment='left', verticalalignment='center',
                      transform=ccrs.PlateCarree(), fontsize=14)
        
    plt.subplots_adjust(hspace=-0.755, wspace=-0.005)
    f.subplots_adjust(right=0.85)
    cbar_ax = f.add_axes([0.855, 0.353, 0.03, 0.298])
    f.colorbar(im, cax=cbar_ax, orientation='vertical', label='Intensity (kn)')
    
    if fig_dir:
        f.savefig(os.path.join(fig_dir, 'fig05.png'), format='png', bbox_inches='tight')
