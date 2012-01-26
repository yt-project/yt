"""
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation:  UC Boulder, KIPAC
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
License:
  Copyright (C) 2008-2011 Britton Smith, Matthew Turk  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as na

import matplotlib
import matplotlib.colors as cc
import matplotlib.cm as mcm
import _colormap_data as _cm

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
_vs = na.linspace(0,1,255)
_kamae_red = na.minimum(255,
                113.9*na.sin(7.64*(_vs**1.705)+0.701)-916.1*(_vs+1.755)**1.862 \
              + 3587.9*_vs+2563.4)/255.0
_kamae_grn = na.minimum(255,
                70.0*na.sin(8.7*(_vs**1.26)-2.418)+151.7*_vs**0.5+70.0)/255.0
_kamae_blu = na.minimum(255,
                194.5*_vs**2.88+99.72*na.exp(-77.24*(_vs-0.742)**2.0)
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
        'red': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.14861 * na.cos(2 * na.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) + 1.78277 * na.sin(2 * na.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'green': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (-0.29227 * na.cos(2 * na.pi * (_s_cubehelix / 3 + _r_cubehelix * x)) - 0.90649 * na.sin(2 * na.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
        'blue': lambda x: x**_gamma_cubehelix + (_h_cubehelix * x**_gamma_cubehelix * (1 - x**_gamma_cubehelix) / 2) * (1.97294 * na.cos(2 * na.pi * (_s_cubehelix / 3 + _r_cubehelix * x))),
}

add_cmap("cubehelix", _cubehelix_data)

# Add colormaps in _colormap_data.py that weren't defined here
_vs = na.linspace(0,1,255)
for k,v in _cm.color_map_luts.iteritems():
    if k not in yt_colormaps:
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
    a = na.ones(b.shape)
    return [r, g, b, a]
