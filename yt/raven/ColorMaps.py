"""
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation:  UC Boulder, KIPAC
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
License:
  Copyright (C) 2008-2009 Britton Smith, Matthew Turk  All Rights Reserved.

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

from yt.raven import *
import matplotlib.colors as cc
import matplotlib.cm as mcm

def check_color(name):
    try:
        ss = cc.colorConverter.to_rgb(name)
        return True
    except ValueError:
        return False

raven_colormaps = {}

def add_cmap(name, cdict):
    raven_colormaps[name] = \
        cc.LinearSegmentedColormap(name,cdict,256)
    mcm.register_cmap(name, raven_colormaps[name])

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
