"""
ColorMaps
---------

Here we grab all of the color maps we can find in the distribution.

Uses a method inspired by
U{http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps}.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
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
