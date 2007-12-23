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
my_opts.add_option("-b", "--basename",
                   action="store", type="string",
                   dest="basename", default="galaxy",
                   help="Basename of parameter files")
my_opts.add_option("-p", "--projection",
                   action="store_true", 
                   dest="projection", default=False,
                   help="Make a zoomin of a projection")
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
my_opts.add_option("-f", "--field",
                   action="store", type="string",
                   dest="field", default="Density",
                   help="Field to display")
my_opts.add_option("-g", "--weight",
                   action="store", type="string",
                   dest="weight", default=None,
                   help="Field to weight projections with")
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
                   dest="output", default="frames/",
                   help="Folder in which to place output frames")

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

for n in range(first,last+1,opts.skip): # This is the array of galaxy outputs we want
    # Now we figure out where this file is
    bn_try = opts.basename + "%04i" % n
    if os.path.isfile(bn_try):
        fn = bn_try
    else:
        fn = os.path.join(bn_try + ".dir", bn_try)
    mylog.info("Now attempting to make frame from %s", fn)
    try:
        a = lagos.EnzoStaticOutput(fn)
        min_dx = a.h.get_smallest_dx()
    except:
        mylog.warning("Something messed up!  Are you sure you gave good info?")
        continue

    pc=raven.PlotCollection(a)
    if opts.center == None:
        mylog.info("No center fed in; seeking.")
        v, center = a.h.find_max("Density")
    else:
        mylog.info("Center fed in; not seeking.")
        center = opts.center
        a.h.center = na.array(center)
    mylog.info("Setting center to %0.3f %0.3f %0.3f", 
                center[0], center[1], center[2])
    if opts.axis == 4:
        axes = range(3)
    else:
        axes = [opts.axis]
    for ax in axes:
        mylog.info("Adding plot for axis %i", ax)
        if opts.projection: pc.add_projection(opts.field, ax,
                                weight_field=opts.weight, center=center)
        else: pc.add_slice(opts.field, ax, center=center)
    pc.set_width(opts.width, opts.unit)
    pc.set_cmap(opts.cmap)
    if opts.zlim: pc.set_zlim(*opts.zlim)
    pc.save(os.path.join(opts.output,"frame%06i" % (n)))
