#!python2.5

import optparse
import math
import os, os.path
import sys

my_opts = optparse.OptionParser()
my_opts.add_option("-w", "--width",
                   action="store", type="float",
                   dest="width", default=1.0,
                   help="Width in specified units")
my_opts.add_option("-u", "--unit",
                   action="store", type="string",
                   dest="unit", default='1',
                   help="Desired units")
my_opts.add_option("-k", "--slab-width",
                   action="store", type="float",
                   dest="slab_width", default=1.0,
                   help="Slab width in specified units")
my_opts.add_option("-g", "--slab-unit",
                   action="store", type="string",
                   dest="slab_unit", default='1',
                   help="Desired units for the slab")
my_opts.add_option("-b", "--basename",
                   action="store", type="string",
                   dest="basename", default="galaxy",
                   help="Basename of parameter files")
my_opts.add_option("-c", "--center",
                   action="store", type="float",
                   dest="center", default=None,
                   nargs=3,
                   help="Center (-1,-1,-1 for max)")
my_opts.add_option("-z", "--zlim",
                   action="store", type="float",
                   dest="zlim", default=None,
                   nargs=2,
                   help="Color limits (min, max)")
my_opts.add_option("-a", "--axis",
                   action="store", type="int",
                   dest="axis", default=4,
                   help="Axis (4 for all three)")
my_opts.add_option("-l", "--log",
                   action="store_true",
                   dest="takelog", default=False,
                   help="Take the log of the field?")
my_opts.add_option("-f", "--field",
                   action="store", type="string",
                   dest="field", default=None,
                   help="Field to color by")
my_opts.add_option("-s", "--skip",
                   action="store", type="int",
                   dest="skip", default=1,
                   help="Skip factor for outputs")
my_opts.add_option("-m", "--colormap",
                   action="store", type="string",
                   dest="cmap", default="jet",
                   help="Colormap name")
my_opts.add_option("-o", "--output",
                   action="store", type="string",
                   dest="output", default="particle_frames/",
                   help="Folder in which to place output frames")
my_opts.add_option("-t", "--particle-type",
                   action="store", type="int",
                   dest="ptype", default=2,
                   help="Particle type to select")
my_opts.add_option("-e", "--age-cut",
                   action="store", type="float",
                   dest="age_filter", default=None,
                   nargs=2,
                   help="Bounds for the field to select")


opts, args = my_opts.parse_args()

import yt.lagos as lagos
import yt.raven as raven
import yt.fido as fido
import numpy as na
mylog = fido.mylog



if not os.path.isdir(opts.output):
    os.mkdir(opts.output)

try:
    first = int(args[0])
    last = int(args[1])
except:
    mylog.error("Hey, sorry, but you need to specify the first and last outputs you want to look at.")
    sys.exit()

print first, last, range(first, last+1, opts.skip)

import pylab
from matplotlib.colors import LogNorm, Normalize

if opts.takelog:
    myNorm = LogNorm
else:
    myNorm = Normalize

for n in range(first,last+1,opts.skip): # This is the array of galaxy outputs we want
    # Now we figure out where this file is
    bn_try = opts.basename + "%04i" % n
    if os.path.isfile(bn_try):
        fn = bn_try
    else:
        fn = os.path.join(bn_try + ".dir", bn_try)
    mylog.info("Now attempting to zoom in on %s", fn)
    try:
        a = lagos.EnzoStaticOutput(fn)
        min_dx = a.h.get_smallest_dx()
    except:
        mylog.warning("Something messed up!  Are you sure you gave good info?")
        continue

    if opts.center == None:
        mylog.info("No center fed in; seeking.")
        v, center = a.h.find_max("Density")
    else:
        mylog.info("Center fed in; not seeking.")
        center = opts.center
    mylog.info("Setting center to %0.3f %0.3f %0.3f", 
                center[0], center[1], center[2])

    # Now we deal with all the plotting.

    if opts.axis == 4:
        axes = range(3)
    else:
        axes = [opts.axis]

    gI = na.where(a.h.gridNumberOfParticles.ravel() > 0)

    npart = na.sum(a.h.gridNumberOfParticles)
    pos = na.zeros((npart,4)) # three plus fielD
    ni = 0
    mylog.info("Getting particles from %s grids", len(gI[0]))
    for grid in a.h.grids[gI]:
        pi = na.where(grid["particle_type"] == opts.ptype)
        x = grid["particle_position_x"][pi]
        y = grid["particle_position_y"][pi]
        z = grid["particle_position_z"][pi]
        pos[ni:ni+x.shape[0],0] = x
        pos[ni:ni+x.shape[0],1] = y
        pos[ni:ni+x.shape[0],2] = z
        if opts.field:
            if opts.field == "age":
                pos[ni:ni+x.shape[0],3] = \
                  (a["InitialTime"] - grid["creation_time"][pi]) \
                  * a["years"]
            else:
                pos[ni:ni+x.shape[0],3] = grid[opts.field][pi]
        ni += x.shape[0]

    for ax in axes:
        mylog.info("Adding plot for axis %i", ax)
        pylab.clf()
        xa = lagos.x_dict[ax]
        ya = lagos.y_dict[ax]
        
        z1 = center[ax] - opts.slab_width/a[opts.slab_unit] * 0.5
        z2 = center[ax] + opts.slab_width/a[opts.slab_unit] * 0.5

        newpos = pos[na.where((pos[:,ax] > z1)
                            & (pos[:,ax] < z2))]
        if newpos.shape[0] == 0:
            mylog.info("Hey, no particles match the criteria.  Moving on.")
            continue
        if opts.age_filter and opts.field == "age":
            newpos = newpos[na.where((newpos[:,3] > opts.age_filter[0])
                                   & (newpos[:,3] < opts.age_filter[1]))]
        mylog.info("Constraining to %s - %s, which gives us %s particles",
                    z1,z2,newpos.shape[0])
        if opts.field:
            #print "MM", \
                  #newpos[:,xa].min(), newpos[:,xa].max(), \
                  #newpos[:,3].min(), newpos[:,3].max()
            pylab.scatter(newpos[:,xa].ravel(),
                          newpos[:,ya].ravel(),
                          c=newpos[:,3].ravel(),
                          faceted=False,s=1.0,
                          alpha=1.0,
                          norm=myNorm())
            if opts.zlim:
                pylab.clim(opts.zlim[0], opts.zlim[1])
            pylab.colorbar()
        else:
            pylab.scatter(newpos[:,xa].ravel(),
                          newpos[:,ya].ravel(),
                          c='k',faceted=False,s=1.0,
                          alpha=1.0)
            
        x1 = center[xa] - \
             opts.width / (2.0*a[opts.unit])
        x2 = center[xa] + \
             opts.width / (2.0*a[opts.unit])
        y1 = center[ya] - \
             opts.width / (2.0*a[opts.unit])
        y2 = center[ya] + \
             opts.width / (2.0*a[opts.unit])
        pylab.xlim(x1, x2)
        pylab.ylim(y1, y2)
        pylab.savefig(os.path.join(opts.output,"particle_frame%06i_%s.png"
                      % (n,lagos.axis_names[ax])))
