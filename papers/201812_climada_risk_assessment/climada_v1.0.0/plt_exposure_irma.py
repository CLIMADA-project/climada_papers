"""
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
from scipy.interpolate import griddata
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from climada.hazard import TCTracks, TropCyclone, Centroids
from climada.entity import BlackMarble, IFTropCyclone, ImpactFuncSet
from climada.engine import Impact
import climada.util.plot as u_plot
from climada.util.scalebar_plot import scale_bar

from constants import YEAR, POLY_VAL, GDP_NLD_ISL, INC_GRP, GDP, CNTRIES

MS2KN = 1.9438444924406049

BUFFER_DEG = -0.051

MAX_VAL = 1.0e7

MIN_VAL = 1.0e2

def exposures():
    """ Generate exposure """
    south_isl = BlackMarble()
    nord_isl = BlackMarble()
    for cntry in CNTRIES:
        if cntry == 'Netherlands':
            ent = BlackMarble()
            ent.set_countries({cntry: ['St. Eustatius', 'Saba']}, YEAR, 
                               res_km=0.1, poly_val=POLY_VAL)
            ent.value = ent.value/ent.value.sum()*GDP_NLD_ISL*(INC_GRP['NLD']+1)
        else:
            ent = BlackMarble()
            ent.set_countries([cntry], YEAR, res_km=0.1, poly_val=POLY_VAL, 
                              **{'gdp': GDP, 'inc_grp': INC_GRP})
    
        if cntry ==  'Turks And Caicos Islands':
            nord_isl.append(ent)
        else:
            south_isl.append(ent)
    
    nord_isl.add_sea((10.0, 0.5))
    south_isl.add_sea((10.0, 0.5))
    
    return south_isl, nord_isl

def add_cntry_names(axis, extent, projection=ccrs.PlateCarree()):
    """Add country names.

    Parameters:
        axis (cartopy.mpl.geoaxes.GeoAxesSubplot): cartopy axis.
        extent (list): geographical limits.
        projection (cartopy.crs projection, optional): geographical projection,
            PlateCarree default.

    """
    shp_file = shapereader.natural_earth(resolution='10m', \
                           category='cultural', name='admin_0_countries')

    shp = shapereader.Reader(shp_file)
    for rec, point in zip(shp.records(), shp.geometries()):
        point_x = point.centroid.xy[0][0]
        point_y = point.centroid.xy[1][0]
        if (point_x <= extent[1]) and (point_x > extent[0]):
            if (point_y <= extent[3]) and (point_y > extent[2]):
                if 'Turks' in rec.attributes['NAME']:
                    axis.text(point_x-0.18, point_y-0.15, rec.attributes['NAME'], \
                    horizontalalignment='left', verticalalignment='top', \
                    transform=projection, fontsize=12)
                elif 'Anguilla' in rec.attributes['NAME']:
                    axis.text(point_x-0.05, point_y+0.05, rec.attributes['NAME'], \
                    horizontalalignment='center', verticalalignment='bottom', \
                    transform=projection, fontsize=12)
                elif 'Barth' in rec.attributes['NAME']:
                    axis.text(point_x-0.05, point_y-0.02, rec.attributes['NAME'], \
                    horizontalalignment='right', verticalalignment='center', \
                    transform=projection, fontsize=12)
                elif 'Briti' in rec.attributes['NAME']:
                    axis.text(point_x, point_y+0.05, rec.attributes['NAME'], \
                    horizontalalignment='center', verticalalignment='bottom', \
                    transform=projection, fontsize=12)
                elif 'U.S.' in rec.attributes['NAME']:
                    axis.text(point_x+0.15, point_y+0.02, rec.attributes['NAME'], \
                    horizontalalignment='center', verticalalignment='bottom', \
                    transform=projection, fontsize=12)
                elif 'Martin' in rec.attributes['NAME']:
                    axis.text(point_x+0.3, point_y, rec.attributes['NAME'], \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)
                elif 'Maarten' in rec.attributes['NAME']:
                    axis.text(point_x-0.4, point_y-0.05, rec.attributes['NAME'], \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)
                elif 'Nev' in rec.attributes['NAME']:
                    axis.text(point_x-0.4, point_y+0.01, 'St. Kitts', \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)
                    axis.text(point_x-0.15, point_y-0.13, 'Nevis', \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)
                elif 'Antigua' in rec.attributes['NAME']:
                    axis.text(point_x-0.1, point_y, 'Antigua', \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)
                    axis.text(point_x-0.34, point_y+0.35, 'Barbuda', \
                    horizontalalignment='center', verticalalignment='center', \
                    transform=projection, fontsize=12)

    if extent[2] < 18.0 < extent[3]:
        axis.text(-63.5, 17.65, 'Saba', \
        horizontalalignment='center', verticalalignment='center', \
        transform=projection, fontsize=12)
        axis.text(-63.4, 17.5, 'St. Eustatius', \
        horizontalalignment='center', verticalalignment='center', \
        transform=projection, fontsize=12)

def irma_tc(exp, data_irma):
    centr = Centroids()
    centr.coord = exp.coord
    centr.id = np.arange(centr.lat.size)

    tc_irma = TropCyclone()
    data_irma.equal_timestep(0.1)
    tc_irma.set_from_tracks(data_irma, centr)

    return tc_irma

def irma_best_track(data_dir):
    file = os.path.join(data_dir, 'ibtracs_global_intp-None_2017242N16333.csv')
    data_irma = TCTracks()
    data_irma.read_ibtracs_csv(file)
    return data_irma

def plot_left(exp, data_irma, tc_irma, ax, scale_pos, cntour_loc, label_loc):
    """ Plot exposed value, irma track and irma wind field. """
    extent = u_plot._get_borders(exp.coord)
    extent = ([extent[0] - BUFFER_DEG, extent[1] + BUFFER_DEG, extent[2] -\
               BUFFER_DEG, extent[3] + BUFFER_DEG])
    ax.set_extent((extent))
    u_plot.add_shapes(ax)

    sel_pos = np.argwhere(exp.value>0)[:,0]
    ax.hexbin(exp.coord[sel_pos, 1], exp.coord[sel_pos, 0], C=exp.value[sel_pos], 
              reduce_C_function=np.average, transform=ccrs.PlateCarree(), 
              gridsize=2000, norm=LogNorm(vmin=MIN_VAL, vmax=MAX_VAL),
              cmap='YlOrRd', vmin=1.0e2, vmax=MAX_VAL)
    ax.set_title('')
    ax.grid(False)
    scale_bar(ax, scale_pos, 10)

    track = data_irma.data[0]
    ax.plot(track.lon.values, track.lat.values, linestyle='solid',
            transform=ccrs.PlateCarree(), lw=2, color='k')
    leg_lines = [Line2D([0], [0], color='k', lw=2), Line2D([0], [0], color='grey', lw=1, ls=':') ]
    leg_names = ['Irma track', 'wind field (kn)']
    if 'bbox' in label_loc:
        ax.legend(leg_lines, leg_names, bbox_to_anchor=label_loc['bbox'], loc=label_loc['loc'])
    else:
        ax.legend(leg_lines, leg_names, loc=label_loc['loc'])

    tc_irma.intensity *= MS2KN
    grid_x, grid_y = np.mgrid[tc_irma.centroids.coord[:, 1].min() : \
                              tc_irma.centroids.coord[:, 1].max() : complex(0, 2000), \
                              tc_irma.centroids.coord[:, 0].min() : \
                              tc_irma.centroids.coord[:, 0].max() : complex(0, 2000)]
    grid_im = griddata((tc_irma.centroids.coord[:, 1], 
                        tc_irma.centroids.coord[:, 0]), 
                        np.array(tc_irma.intensity[0].todense()).squeeze(), 
                        (grid_x, grid_y))
    cs = ax.contour(grid_x, grid_y, grid_im, linewidths=1.0, linestyles=':', \
                    levels=[60, 80, 100, 120], colors=['grey', 'grey', 'grey', 'grey', 'grey'])
    ax.clabel(cs, inline=1, fontsize=10, manual=cntour_loc, rotation=-20, fmt='%1.0f')
    tc_irma.intensity /= MS2KN

def plot_right(irma_tc, exp, ax, scale_pos, plot_line=False):
    """ Plot irma damage in USD. """
    if_exp = ImpactFuncSet()
    if_em = IFTropCyclone()
    if_em.set_emanuel_usa()
    if_exp.add_func(if_em)

    imp_irma = Impact()
    imp_irma.calc(exp, if_exp, irma_tc)
    extent = u_plot._get_borders(exp.coord)
    extent = ([extent[0] - BUFFER_DEG, extent[1] + BUFFER_DEG, 
               extent[2] - BUFFER_DEG, extent[3] + BUFFER_DEG])
    ax.set_extent((extent))
    u_plot.add_shapes(ax)

    sel_pos = np.argwhere(imp_irma.eai_exp>0)[:,0]
    hex_bin = ax.hexbin(imp_irma.coord_exp[sel_pos, 1],  imp_irma.coord_exp[sel_pos, 0], 
                        C=imp_irma.eai_exp[sel_pos], reduce_C_function=np.average,
                        transform=ccrs.PlateCarree(), gridsize=2000, 
                        norm=LogNorm(vmin=MIN_VAL, vmax=MAX_VAL),
                        cmap='YlOrRd', vmin=MIN_VAL, vmax=MAX_VAL)
    ax.set_title('')
    ax.grid(False)
    add_cntry_names(ax, extent)
    scale_bar(ax, scale_pos, 10)

    if plot_line:
        x1, y1 = [-64.57, -64.82], [18.28, 18.47]
        ax.plot(x1, y1, linewidth=1.0, color='grey', linestyle='--')

    return hex_bin

def plot_percen(irma_tc, exp, if_exp, axs):
    """ Plot irma damage in %. """
    # south
    extent = u_plot._get_borders(exp.coord)
    extent = ([extent[0] - BUFFER_DEG, extent[1] + BUFFER_DEG, 
               extent[2] - BUFFER_DEG, extent[3] + BUFFER_DEG])
    axs.set_extent((extent))
    u_plot.add_shapes(axs)

    imp_irma = Impact()
    imp_irma.calc(exp, if_exp, irma_tc)
    imp_irma.eai_exp[exp.value > 0] = \
        imp_irma.eai_exp[exp.value > 0]/exp.value[exp.value > 0]*100
    imp_irma.eai_exp[exp.value == 0] = 0.
    sel_exp = imp_irma.eai_exp > 0
    im = axs.hexbin(exp.coord[sel_exp, 1], exp.coord[sel_exp, 0], 
                    C=imp_irma.eai_exp[sel_exp], reduce_C_function=np.average,
                    transform=ccrs.PlateCarree(), gridsize=2000, cmap='YlOrRd', 
                    vmin=0, vmax=50)
    axs.set_title('')
    axs.grid(False)
    scale_bar(axs, (0.90, 0.90), 10)
    
    return im

def irma_percen(irma_tc_s, irma_tc_n, exp_s, exp_n):
    """ Plot irma damage in % in lesser antilles and TCA. """
    if_exp = ImpactFuncSet()
    if_em = IFTropCyclone()
    if_em.set_emanuel_usa()
    if_exp.add_func(if_em)

    fig, axs = plt.subplots(1, 2, figsize=(16, 20), 
                            subplot_kw=dict(projection=ccrs.PlateCarree()),
                            squeeze=True, sharex=False, sharey=False)

    for idx, axis in enumerate(axs.flatten()):
        grid = axis.gridlines(draw_labels=True, alpha=0.2)
        grid.xlabels_top = grid.ylabels_right = False
        grid.xformatter = LONGITUDE_FORMATTER
        grid.yformatter = LATITUDE_FORMATTER

    # south
    plot_percen(irma_tc_s, exp_s, if_exp, axs[0])

    # nord
    im = plot_percen(irma_tc_n, exp_n, if_exp, axs[1])

    plt.subplots_adjust(wspace=0.14)
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.885, 0.4285, 0.03, 0.148])
    fig.colorbar(im, cax=cbar_ax, orientation='vertical', label='% Damage')

    return fig

def irma_and_exposures(exp_s, exp_t, tc_irma_s, tc_irma_t, data_irma):
    """ Plot exposure, irma track and contour lines, and irma damage. """
    # Lesser Antilles
    f1, axes1 = plt.subplots(1, 2, figsize=(16, 20), 
                             subplot_kw=dict(projection=ccrs.PlateCarree()),
                             squeeze=True, sharex=False, sharey=True)

    for idx, axis in enumerate(axes1.flatten()):
        grid = axis.gridlines(draw_labels=True, alpha=0.2)
        grid.xlabels_top = grid.ylabels_right = False
        if idx == 0:
            grid.xformatter = LONGITUDE_FORMATTER
            grid.yformatter = LATITUDE_FORMATTER
        elif idx == 1:
            grid.xformatter = LONGITUDE_FORMATTER
            grid.ylabels_left = False

    cntour_loc = [(-64.6, 17.4), (-63.6, 17.4), (-64.3, 18), (-63.9, 18), \
                  (-63.8, 18.4), (-62.4, 18.4)]
    label_loc = dict()
    label_loc['loc'] = 3
    plot_left(exp_s, data_irma, tc_irma_s, axes1[0], (0.90, 0.90), cntour_loc, 
              label_loc)
    im = plot_right(tc_irma_s, exp_s, axes1[1], (0.90, 0.90), plot_line=True)

    plt.subplots_adjust(wspace=0)
    f1.subplots_adjust(right=0.88)
    cbar_ax = f1.add_axes([0.885, 0.4221, 0.03, 0.161])
    f1.colorbar(im, cax=cbar_ax, orientation='vertical', label='Value (USD)')

    # Turks and Caicos Isl.
    f2, axes2 = plt.subplots(1, 2, figsize=(16, 20), 
                             subplot_kw=dict(projection=ccrs.PlateCarree()),
                             squeeze=True, sharex=False, sharey=True)

    for idx, axis in enumerate(axes2.flatten()):
        grid = axis.gridlines(draw_labels=True, alpha=0.2)
        grid.xlabels_top = grid.ylabels_right = False
        if idx == 0:
            grid.yformatter = LATITUDE_FORMATTER
            grid.xformatter = LONGITUDE_FORMATTER
        elif idx == 1:
            grid.ylabels_left = False
            grid.xformatter = LONGITUDE_FORMATTER

    cntour_loc = [(-72, 21.65), (-71.6, 21.67), (-71.64, 21.9)]
    label_loc = dict()
    label_loc['loc'] = 1
    label_loc['bbox'] = [0.985, 0.885]
    plot_left(exp_t, data_irma, tc_irma_t, axes2[0], (0.90, 0.90), cntour_loc, 
              label_loc)
    im = plot_right(tc_irma_t, exp_t, axes2[1], (0.90, 0.90))

    plt.subplots_adjust(wspace=0)
    f2.subplots_adjust(right=0.88)
    cbar_ax = f2.add_axes([0.885, 0.4235, 0.03, 0.158])
    f2.colorbar(im, cax=cbar_ax, orientation='vertical', label='Value (USD)')
    return f1, f2

def fig03_fig04(ibtracs_dir, data_dir, fig_dir=None):
    """ Generate figures 03 and 04 """
    south_isl, nord_isl = exposures()
    
    data_irma = irma_best_track(ibtracs_dir)
    tc_irma_s = irma_tc(south_isl, data_irma)
    
    data_irma = irma_best_track(ibtracs_dir)
    tc_irma_t = irma_tc(nord_isl, data_irma)
    
    data_irma = irma_best_track(ibtracs_dir)

    fig03a, fig03b = irma_and_exposures(south_isl, nord_isl, tc_irma_s, 
                                        tc_irma_t, data_irma)
    if fig_dir:
        fig03a.savefig(os.path.join(fig_dir, 'fig03a.png'), format='png', bbox_inches='tight')
        fig03b.savefig(os.path.join(fig_dir, 'fig03b.png'), format='png', bbox_inches='tight')

    fig04 = irma_percen(tc_irma_s, tc_irma_t, south_isl, nord_isl)
    if fig_dir:
        fig04.savefig(os.path.join(fig_dir, 'fig04.png'), format='png', bbox_inches='tight')
