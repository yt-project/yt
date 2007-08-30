# Simple example of how to use matplotlib for plotting.
# A lot of this will become automated, but here's some lower level stuff
# to give an idea of what it's all about.

# Note that this just plots three linked AMR plots.  (They can be projections,
# too!)  We could just as easily make a fourth plot that is of any type.
# Color bar support will be added soon.  For now it is not clear to me the best
# way to proceed.

from yt import ytcfg
ytcfg["raven","backend"] = "MPL"

# Set the parameters of our movie
maxwidth = (1.0,'1')
minwidth = (10.0,"rsun")
numframes = 400

filename_template = "frames/frame%04i.png"
fieldName = "NumberDensity"
hierarchy_filename = "DataDump0043.dir/DataDump0043"

# Import the modules we need
import yt.lagos as lagos
from yt.arraytypes import *
import yt.raven as raven
import matplotlib
import matplotlib.cm as cm
import matplotlib.colorbar as colorbar
import matplotlib.contour as contour
import pylab
from math import log10

# Standard data loading
a=lagos.EnzoStaticOutput(hierarchy_filename)
v,c=a.hierarchy.findMax("NumberDensity")
centers = [(c[1],c[2]),(c[0],c[2]),(c[0],c[1])]

# We want a square figure
fig = pylab.figure(figsize=(8,8))

# This makes them all adjacent
fig.subplots_adjust(hspace=0,wspace=0,bottom=0.0, top=1.0, left=0.0, right=1.0)

# This sets the color map
pylab.prism()

# These are the inital bounds
absmin = 1e+30
absmax = 1e-30

# Set up lists...
slices = []
axes=[]
ims = []

# Now, for each axis, do the thingie
for i in range(3):
    slices.append(a.h.slice(i, c[i], [fieldName], center=c))
    # Subplots are 1-indexed, so we do i+1
    axes.append(fig.add_subplot(2,2,i+1, aspect='equal'))
    ims.append(raven.be.SlicePlot(slices[-1], fieldName, fig, axes[-1], False))
    absmin = min((slices[-1][fieldName]).min(), absmin)
    absmax = max((slices[-1][fieldName]).max(), absmax)
    #axes[-1].set_xticks(()) # we don't want any of these things
    #axes[-1].set_yticks(())
    #axes[-1].set_xlabel("")
    #axes[-1].set_ylabel("")

axes.append(fig.add_subplot(2,2,4, aspect='equal'))
axes[-1].axesFrame.set_visible(False)
axes[-1].set_axis_off()

# For easier access
ac = zip(axes, centers, ims)
i=0

#ss = "Your Text Here"
#text=fig.text(0.55, 0.25, ss, size='large')

# We want to do this log-spaced.
# Note that we make numframes complex, so that it is regarded by mgrid
# as a number of steps to take, rather than a step-size interval.
for w in na.mgrid[log10(maxwidth[0]/a[maxwidth[1]])\
                 :log10(minwidth[0]/a[minwidth[1]])\
                 :numframes*1j]:
    width = 10**w
    # Zoom in...
    for ax, cc, im in ac:
        im.set_width(width, '1')
        im.set_zlim(absmin, absmax)
    #for im in ims: im.set_zlim(absmin, absmax)
    print "About to save with limits", absmin, absmax
    s = filename_template % (i)
    fig.savefig(s)
    print "Saved %s" % (s)
    i+=1
