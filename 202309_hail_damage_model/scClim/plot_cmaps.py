# -*- coding: utf-8 -*-
"""
Customized colormaps for plotting

"""

# modules
import matplotlib.colors
import numpy as np
import colorcet as cc
import copy

# colormab perc2_9lev
# source: https://www.ncl.ucar.edu/Document/Graphics/ColorTables/perc2_9lev.shtml
rgbtable = np.array([[215., 227., 238.],
                     [181., 202., 255.],
                     [143., 179., 255.],
                     [127., 151., 255.],
                     [171., 207.,  99.],
                     [232., 245., 158.],
                     [255., 209.,  33.],
                     [255., 250.,  20.],
                     [255., 163.,  10.],
                     [255.,  76.,   0.]])
N = len(rgbtable)
cdict={
    'red':    [(i/(N-1),rgbtable[i,0]/255,rgbtable[i,0]/255) for i in range(0,N)],
    'green':  [(i/(N-1),rgbtable[i,1]/255,rgbtable[i,1]/255) for i in range(0,N)],
    'blue':   [(i/(N-1),rgbtable[i,2]/255,rgbtable[i,2]/255) for i in range(0,N)]
}
CMAP_perc2_9lev = matplotlib.colors.LinearSegmentedColormap('camp',cdict)

# colormab precip_11lev
# source: https://www.ncl.ucar.edu/Document/Graphics/ColorTables/precip_11lev.shtml
rgbtable = np.array([#[255., 255., 255.], #remove white
                     [237., 250., 194.],
                     [205., 255., 205.],
                     [153., 240., 178.],
                     [ 83., 189., 159.],
                     [ 50., 166., 150.],
                     [ 50., 150., 180.],
                     [  5., 112., 176.],
                     [  5.,  80., 140.],
                     [ 10.,  31., 150.],
                     [ 44.,   2.,  70.],
                     [106.,  44.,  90.]])
N = len(rgbtable)
cdict={
    'red':    [(i/(N-1),rgbtable[i,0]/255,rgbtable[i,0]/255) for i in range(0,N)],
    'green':  [(i/(N-1),rgbtable[i,1]/255,rgbtable[i,1]/255) for i in range(0,N)],
    'blue':   [(i/(N-1),rgbtable[i,2]/255,rgbtable[i,2]/255) for i in range(0,N)]
}
CMAP_precip_11lev = matplotlib.colors.LinearSegmentedColormap('camp',cdict)

# colormab precip3_16lev
# source: https://www.ncl.ucar.edu/Document/Graphics/ColorTables/precip3_16lev.shtml
rgbtable = np.array([[255., 255., 255.],
                     [214., 226., 255.],
                     [181., 201., 255.],
                     [142., 178., 255.],
                     [127., 150., 255.],
                     [ 99., 112., 247.],
                     [  0.,  99., 255.],
                     [  0., 150., 150.],
                     [  0., 198.,  51.],
                     [ 99., 255.,   0.],
                     [150., 255.,   0.],
                     [198., 255.,  51.],
                     [255., 255.,   0.],
                     [255., 198.,   0.],
                     [255., 160.,   0.],
                     [255., 124.,   0.],
                     [255.,  25.,   0.]])
N = len(rgbtable)
cdict={
    'red':    [(i/(N-1),rgbtable[i,0]/255,rgbtable[i,0]/255) for i in range(0,N)],
    'green':  [(i/(N-1),rgbtable[i,1]/255,rgbtable[i,1]/255) for i in range(0,N)],
    'blue':   [(i/(N-1),rgbtable[i,2]/255,rgbtable[i,2]/255) for i in range(0,N)]
}
CMAP_precip3_16lev = matplotlib.colors.LinearSegmentedColormap('camp',cdict)

# COLORMAP DAMAGES TRUNCATED / PERCEPTUALLY UNIFORM
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

CMAP_IMPACT_CC = truncate_colormap(cc.cm.fire_r, 0.1, 1)
CMAP_IMPACT_CC.set_under('white',alpha=0)
CMAP_IMPACT_CHF = truncate_colormap(cc.cm.rainbow4, 0.5, 1)

CMAP_HAZARD_HAILCAST = copy.copy(cc.cm.bmy)