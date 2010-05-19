"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation:  UCSD
License:
  Copyright (C) 2010 Matthew Turk  All Rights Reserved.

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
from yt.funcs import *
import _colormap_data as cmd
import yt.amr_utils as au

def write_image(image, filename, color_bounds = None, cmap_name = "algae"):
    if color_bounds is None:
        mi = na.nanmin(image[~na.isinf(image)])
        ma = na.nanmax(image[~na.isinf(image)])
        color_bounds = mi, ma
    image = (image - color_bounds[0])/(color_bounds[1] - color_bounds[0])
    to_plot = map_to_colors(image, cmap_name)
    au.write_png(to_plot, filename)
    return to_plot

def map_to_colors(buff, cmap_name):
    if cmap_name not in cmd.color_map_luts:
        print "Your color map was not found in the extracted colormap file."
        raise KeyError(cmap_name)
    lut = cmd.color_map_luts[cmap_name]
    x = na.mgrid[0.0:1.0:lut[0].shape[0]*1j]
    shape = buff.shape
    mapped = na.dstack(
            [(na.interp(buff, x, v)*255) for v in lut ]).astype("uint8")
    return mapped.copy("C")

def strip_colormap_data(fn = "color_map_data.py",
            cmaps = ("jet", "algae", "hot", "gist_stern")):
    import yt.raven, pprint
    import yt.raven.ColorMaps as rcm
    f = open(fn, "w")
    f.write("### Auto-generated colormap tables, taken from Matplotlib ###\n\n")
    f.write("from numpy import array\n")
    f.write("color_map_luts = {}\n\n\n")
    if cmaps is None: cmaps = yt.raven.ColorMaps
    for cmap_name in sorted(cmaps):
        print "Stripping", cmap_name
        vals = rcm._extract_lookup_table(cmap_name)
        f.write("### %s ###\n\n" % (cmap_name))
        f.write("color_map_luts['%s'] = \\\n" % (cmap_name))
        f.write("   (\n")
        for v in vals:
            f.write(pprint.pformat(v, indent=3))
            f.write(",\n")
        f.write("   )\n\n")
    f.close()

def splat_points(image, points_x, points_y,
                 contribution = None, transposed = False):
    if contribution is None:
        contribution = 100.0
    val = contribution * 1.0/points_x.size
    if transposed:
        points_y = 1.0 - points_y
        points_x = 1.0 - points_x
    im = image.copy()
    au.add_points_to_image(im, points_x, points_y, val)
    return im
