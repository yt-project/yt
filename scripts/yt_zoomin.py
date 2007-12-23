#!python2.5

import yt.lagos as lagos
import yt.raven as raven
import yt.fido as fido
import numpy as na
import optparse
import math
import os, os.path

mylog = fido.mylog

my_opts = optparse.OptionParser()
my_opts.add_option("-w", "--max-width",
                   action="store", type="float",
                   dest="max_width", default=1.0,
                   help="Maximum width in code units")
my_opts.add_option("-z", "--min-width",
                   action="store", type="float",
                   dest="min_width", default=50,
                   help="Minimum width in units of smallest dx, defaulting to 50")
my_opts.add_option("-p", "--projection",
                   action="store_true", 
                   dest="projection", default=False,
                   help="Make a zoomin of a projection")
my_opts.add_option("-a", "--axis",
                   action="store", type="int",
                   dest="axis", default=4,
                   help="Axis (4 for all three)")
my_opts.add_option("-f", "--field",
                   action="store", type="string",
                   dest="field", default="Density",
                   help="Field to zoom in on")
my_opts.add_option("-g", "--weight",
                   action="store", type="string",
                   dest="weight", default=None,
                   help="Field to weight projections with")
my_opts.add_option("-n", "--nframes",
                   action="store", type="int",
                   dest="nframes", default=100,
                   help="Number of frames to generate")
my_opts.add_option("-o", "--output",
                   action="store", type="string",
                   dest="output", default="frames/",
                   help="Folder in which to place output frames")

opts, args = my_opts.parse_args()

for arg in args:
    mylog.info("Now attempting to zoom in on %s", arg)
    try:
        a = lagos.EnzoStaticOutput(arg)
        min_dx = a.h.get_smallest_dx()
        min_width = min_dx * opts.min_width
    except:
        mylog.warning("Something messed up!  Are you sure you specified the file correctly?")
        continue
    if opts.axis == 4:
        axes = range(3)
    else:
        axes = [opts.axis]
    pc = raven.PlotCollection(a)
    for ax in axes:
        mylog.info("Adding plot for axis %i", ax)
        if opts.projection: pc.add_projection(opts.field, ax,
                                weight_field=opts.weight)
        else: pc.add_slice(opts.field, ax)
    pc.set_width(opts.max_width,'1')
    # Check the output directory
    if not os.path.isdir(opts.output):
        os.mkdir(opts.output)
    # Figure out our zoom factor
    # Recall that factor^nframes = min_width / max_width
    # so factor = (log(min/max)/log(nframes))
    mylog.info("min_width: %0.3e max_width: %0.3e nframes: %0.3e",
               min_width, opts.max_width, opts.nframes)
    factor=10**(math.log10(min_width/opts.max_width)/opts.nframes)
    mylog.info("Zoom factor: %0.3e", factor)
    w = 1.0
    for i in range(opts.nframes):
        mylog.info("Setting width to %0.3e", w)
        mylog.info("Saving frame %06i",i)
        pc.set_width(w,"1")
        pc.save(os.path.join(opts.output,"%s_frame%06i" % (a,i)))
        w *= factor
