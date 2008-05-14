"""
@author: Britton Smith, Matthew Turk
@organization: UC Boulder, KIPAC
@contact: brittonsmith@gmail.com, mturk@slac.stanford.edu
@license:
  Copyright (C) 2008 Britton Smith, Matthew Turk  All Rights Reserved.

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

import matplotlib.colors as cc

raven_colormaps = {}

def add_cmap(name, cdict):
    raven_colormaps[name] = \
        cc.LinearSegmentedColormap(name,cdict,256)

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
                   (0.8, 200/256., 200/256.),
                   (1.0, 230/256., 230/256.)),
         'green': ((0.0, 0/256., 0/256.),
                   (0.2, 0/256., 0/256.),
                   (0.4, 130/256., 130/256.),
                   (0.6, 256/256., 256/256.),
                   (0.8, 100/256., 100/256.),
                   (1.0, 0.0, 0.0)),
         'blue':  ((0.0, 80/256., 80/256.),
                   (0.2, 220/256., 220/256.),
                   (0.4, 0.0, 0.0),
                   (0.6, 20/256., 20/256.),
                   (0.8, 0/256., 0/256.),
                   (1.0, 0.0, 0.0))}

add_cmap('bds_highcontrast', cdict)
