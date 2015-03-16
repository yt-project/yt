"""


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np
from yt.extern.six.moves import zip as izip

import matplotlib
import matplotlib.colors as cc
import matplotlib.cm as mcm
from . import _colormap_data as _cm

def is_colormap(cmap):
    return isinstance(cmap,cc.Colormap)

def check_color(name):
    try:
        ss = cc.colorConverter.to_rgb(name)
        return True
    except ValueError:
        return False

yt_colormaps = {}

def add_cmap(name, cdict):
    yt_colormaps[name] = \
        cc.LinearSegmentedColormap(name,cdict,256)
    mcm.datad[name] = cdict
    mcm.__dict__[name] = cdict
    try: # API compatibility
        mcm.register_cmap(name, yt_colormaps[name])
    except AttributeError:
        pass
    

# The format is as follows:
#   First number is the number at which we are defining a color breakpoint
#   Second number is the (0..1) number to interpolate to when coming *from below*
#   Third number is the (0..1) number to interpolate to when coming *from above*

# Next up is boilerplate -- the name, the colormap dict we just made, and the
# number of segments we want.  This is probably fine as is.

cdict = {'red':   ((0.0, 80/256., 80/256.),
                   (0.2, 0.0, 0.0),
                   (0.4, 0.0, 0.0),
                   (0.6, 256/256., 256/256.),
                   (0.95, 256/256., 256/256.),
                   (1.0, 150/256., 150/256.)),
         'green': ((0.0, 0/256., 0/256.),
                   (0.2, 0/256., 0/256.),
                   (0.4, 130/256., 130/256.),
                   (0.6, 256/256., 256/256.),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 80/256., 80/256.),
                   (0.2, 220/256., 220/256.),
                   (0.4, 0.0, 0.0),
                   (0.6, 20/256., 20/256.),
                   (1.0, 0.0, 0.0))}

add_cmap('bds_highcontrast', cdict)
add_cmap('algae', cdict)

# Set the default colormap to be algae.
matplotlib.rc('image', cmap="algae")

# This next colormap was designed by Tune Kamae and converted here by Matt
_vs = np.linspace(0,1,255)
_kamae_red = np.minimum(255,
                113.9*np.sin(7.64*(_vs**1.705)+0.701)-916.1*(_vs+1.755)**1.862 \
              + 3587.9*_vs+2563.4)/255.0
_kamae_grn = np.minimum(255,
                70.0*np.sin(8.7*(_vs**1.26)-2.418)+151.7*_vs**0.5+70.0)/255.0
_kamae_blu = np.minimum(255,
                194.5*_vs**2.88+99.72*np.exp(-77.24*(_vs-0.742)**2.0)
              + 45.40*_vs**0.089+10.0)/255.0

cdict = {'red':zip(_vs,_kamae_red,_kamae_red),
         'green':zip(_vs,_kamae_grn,_kamae_grn),
         'blue':zip(_vs,_kamae_blu,_kamae_blu)}
add_cmap('kamae', cdict)

# This one is a simple black & green map

cdict = {'red':   ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),
         'blue':  ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0))}

add_cmap('black_green', cdict)

# This one is a variant of a colormap commonly
# used for X-ray observations by Maxim Markevitch

cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.3, 0.0, 0.0),
                 (0.352, 0.245, 0.245),
                 (0.42, 0.5, 0.5),
                 (0.51, 0.706, 0.706),
                 (0.613, 0.882, 0.882),
                 (0.742, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),
                   (0.585, 0.0, 0.0),
                   (0.613, 0.196, 0.196),
                   (0.693, 0.48, 0.48),
                   (0.785, 0.696, 0.696),
                   (0.885, 0.882, 0.882),
                   (1.0, 1.0, 1.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.136, 0.0, 0.0),
                  (0.136, 0.373, 0.373),
                  (0.391, 1.0, 1.0),
                  (1.0, 1.0, 1.0))}

add_cmap("purple_mm", cdict)

# This one comes from
# http://permalink.gmane.org/gmane.comp.python.matplotlib.devel/10518
# and is an implementation of http://arxiv.org/abs/1108.5083
#

# cubehelix parameters
_gamma_cubehelix = 1.0
_s_cubehelix = 0.5
_r_cubehelix = -1.5
_h_cubehelix = 1.0

_cubehelix_data = {
        'red': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.14861 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) + 1.78277 * np.sin(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'green': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.29227 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) - 0.90649 * np.sin(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'blue': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (1.97294 * np.cos(2 * np.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
}

add_cmap("cubehelix", _cubehelix_data)

# Add colormaps in _colormap_data.py that weren't defined here
_vs = np.linspace(0,1,256)
for k,v in list(_cm.color_map_luts.items()):
    if k not in yt_colormaps and k not in mcm.cmap_d:
        cdict = { 'red': zip(_vs,v[0],v[0]),
                  'green': zip(_vs,v[1],v[1]),
                  'blue': zip(_vs,v[2],v[2]) }
        add_cmap(k, cdict)

def _extract_lookup_table(cmap_name):
    cmap = mcm.get_cmap(cmap_name)
    if not cmap._isinit: cmap._init()
    r = cmap._lut[:-3, 0]
    g = cmap._lut[:-3, 1]
    b = cmap._lut[:-3, 2]
    a = np.ones(b.shape)
    return [r, g, b, a]

def show_colormaps(subset = "all", filename=None):
    """
    Displays the colormaps available to yt.  Note, most functions can use
    both the matplotlib and the native yt colormaps; however, there are 
    some special functions existing within image_writer.py (e.g. write_image()
    write_bitmap(), etc.), which cannot access the matplotlib
    colormaps.

    In addition to the colormaps listed, one can access the reverse of each 
    colormap by appending a "_r" to any map.
    
    Parameters
    ----------

    subset : string, opt

        valid values : "all", "yt_native"
        default : "all"

        As mentioned above, a few functions can only access yt_native 
        colormaps.  To display only the yt_native colormaps, set this
        to "yt_native".  

    filename : string, opt

        default: None

        If filename is set, then it will save the colormaps to an output
        file.  If it is not set, it will "show" the result interactively.
    """
    import pylab as pl

    a=np.outer(np.arange(0,1,0.01), np.ones(10))
    if (subset == "all"):
        maps = [ m for m in pl.cm.datad if (not m.startswith("idl")) & (not m.endswith("_r"))]
    if (subset == "yt_native"):
        maps = [ m for m in _cm.color_map_luts if (not m.startswith("idl")) & (not m.endswith("_r"))]
    maps = list(set(maps))
    maps.sort()
    # scale the image size by the number of cmaps
    pl.figure(figsize=(2.*len(maps)/10.,6))
    pl.subplots_adjust(top=0.7,bottom=0.05,left=0.01,right=0.99)
    l = len(maps)+1
    for i,m in enumerate(maps):
        pl.subplot(1,l,i+1)
        pl.axis("off")
        pl.imshow(a, aspect='auto',cmap=pl.get_cmap(m),origin="lower")      
        pl.title(m,rotation=90, fontsize=10, verticalalignment='bottom')
    if filename is not None:
        pl.savefig(filename, dpi=100, facecolor='gray') 
    else:  
        pl.show()
