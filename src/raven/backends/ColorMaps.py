"""
ColorMaps
---------

Here we grab all of the color maps we can find in the distribution.

Uses a method inspired by
U{http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps}.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

from matplotlib.colors import LinearSegmentedColormap

ravenColorMaps = {}

files = [] # Temporarily, we will not do anything here

for file in files:
    cdict = {'red':[],'green':[],'blue':[]}
    # Put in a checker here
    lines = open(file).readlines()
    for lineI in range(len(lines)):
        line = lines[lineI]
        if line.startswith("#"): continue
        r,g,b = line.split()[:3]
        cdict['red'].append(float(r))
        cdict['green'].append(float(g))
        cdict['blue'].append(float(b))
    
    ravenColorMaps[os.basename(file)] = LinearSegmentedColormap(cdict, len(cdict['r']))
