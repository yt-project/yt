"""
A means of running standalone commands with a shared set of options.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

from yt.mods import *
from yt.recipes import _fix_pf
import optparse, os, os.path, math

_common_options = dict(
    axis    = dict(short="-a", long="--axis",
                   action="store", type="int",
                   dest="axis", default=4,
                   help="Axis (4 for all three)"),
    log     = dict(short="-l", long="--log",
                   action="store_true",
                   dest="takelog", default=True,
                   help="Take the log of the field?"),
    field   = dict(short="-f", long="--field",
                   action="store", type="string",
                   dest="field", default="Density",
                   help="Field to color by"),
    weight  = dict(short="-g", long="--weight",
                   action="store", type="string",
                   dest="weight", default=None,
                   help="Field to weight projections with"),
    cmap    = dict(short="", long="--colormap",
                   action="store", type="string",
                   dest="cmap", default="jet",
                   help="Colormap name"),
    zlim    = dict(short="-z", long="--zlim",
                   action="store", type="float",
                   dest="zlim", default=None,
                   nargs=2,
                   help="Color limits (min, max)"),
    width   = dict(short="-w", long="--width",
                   action="store", type="float",
                   dest="width", default=1.0,
                   help="Width in specified units"),
    unit    = dict(short="-u", long="--unit",
                   action="store", type="string",
                   dest="unit", default='1',
                   help="Desired units"),
    center  = dict(short="-c", long="--center",
                   action="store", type="float",
                   dest="center", default=None,
                   nargs=3,
                   help="Center (-1,-1,-1 for max)"),
    bn      = dict(short="-b", long="--basename",
                   action="store", type="string",
                   dest="basename", default="galaxy",
                   help="Basename of parameter files"),
    output  = dict(short="-o", long="--output",
                   action="store", type="string",
                   dest="output", default="frames/",
                   help="Folder in which to place output images"),
    skip    = dict(short="-s", long="--skip",
                   action="store", type="int",
                   dest="skip", default=1,
                   help="Skip factor for outputs"),
    proj    = dict(short="-p", long="--projection",
                   action="store_true", 
                   dest="projection", default=False,
                   help="Use a projection rather than a slice"),
    maxw    = dict(short="", long="--max-width",
                   action="store", type="float",
                   dest="max_width", default=1.0,
                   help="Maximum width in code units"),
    minw    = dict(short="", long="--min-width",
                   action="store", type="float",
                   dest="min_width", default=50,
                   help="Minimum width in units of smallest dx (default: 50)"),
    nframes = dict(short="-n", long="--nframes",
                   action="store", type="int",
                   dest="nframes", default=100,
                   help="Number of frames to generate"),
    slabw   = dict(short="", long="--slab-width",
                   action="store", type="float",
                   dest="slab_width", default=1.0,
                   help="Slab width in specified units"),
    slabu   = dict(short="-g", long="--slab-unit",
                   action="store", type="string",
                   dest="slab_unit", default='1',
                   help="Desired units for the slab"),
    ptype   = dict(short="", long="--particle-type",
                   action="store", type="int",
                   dest="ptype", default=2,
                   help="Particle type to select"),
    agecut  = dict(short="", long="--age-cut",
                   action="store", type="float",
                   dest="age_filter", default=None,
                   nargs=2,
                   help="Bounds for the field to select"),
    )

def _add_options(parser, *options):
    for opt in options:
        oo = _common_options[opt].copy()
        parser.add_option(oo.pop("short"), oo.pop("long"), **oo)

def _get_parser(*options):
    parser = optparse.OptionParser()
    _add_options(parser, *options)
    return parser

# Now we define our functions, each of which will be an 'entry_point' when
# passed to setuptools.

def zoomin():
    parser = _get_parser("maxw", "minw", "proj", "axis", "field", "weight",
                             "zlim", "nframes", "output", "cmap")
    opts, args = parser.parse_args()

    for arg in args:
        pf = _fix_pf(arg)
        min_width = opts.min_width * pf.h.get_smallest_dx()
        if opts.axis == 4:
            axes = range(3)
        else:
            axes = [opts.axis]
        pc = PlotCollection(pf)
        for ax in axes: 
            if opts.projection: pc.add_projection(opts.field, ax,
                                    weight_field=opts.weight)
            else: pc.add_slice(opts.field, ax)
        pc.set_width(opts.max_width,'1')
        pc.set_cmap(opts.cmap)
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
            if opts.zlim: pc.set_zlim(*opts.zlim)
            pc.save(os.path.join(opts.output,"%s_frame%06i" % (pf,i)))
            w *= factor

def timeseries():
    parser = _get_parser("width", "unit", "bn", "proj", "center",
                         "zlim", "axis", "field", "weight", "skip",
                         "cmap", "output")
    opts, args = parser.parse_args()

    try:
        first = int(args[0])
        last = int(args[1])
    except:
        mylog.error("Hey, sorry, but you need to specify the indices of the first and last outputs you want to look at.")
        sys.exit()

    for n in range(first,last+1,opts.skip):
        # Now we figure out where this file is
        bn_try = "%s%04i" % (opts.basename, n)
        try:
            pf = _get_pf(bn_try)
        except IOError:
            pf = _get_pf("%s.dir/%s" % (bn_try, bn_try))
        pc=raven.PlotCollection(pf)
        center = opts.center
        if center is None or opts.center == (-1,-1,-1):
            mylog.info("No center fed in; seeking.")
            v, center = a.h.find_max("Density")
        center = na.array(center)
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
        pc.save(os.path.join(opts.output,"%s" % (pf)))
