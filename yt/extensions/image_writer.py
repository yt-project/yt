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
import types

def _scale_image(image):
    if isinstance(image, na.ndarray) and image.dtype == na.uint8:
        return image
    if isinstance(image, (types.TupleType, types.ListType)):
        image, mi, ma = image
    else:
        mi, ma = image.min(), image.max()
    image = (na.clip((image-mi)/(ma-mi) * 255, 0, 255)).astype('uint8')
    return image

def multi_image_composite(fn, red_channel, blue_channel,
                          green_channel = None, alpha_channel = None):
    """
    Accepts at least a red and a blue array, of shape (N,N) each, that are
    optionally scaled and composited into a final image, written into *fn*.
    Can also accept green and alpha.  For instance:

        >>> red_channel = na.log10(frb["Temperature"])
        >>> blue_channel = na.log10(frb["Density"])
        >>> iw.multi_image_composite("multi_channel1.png", red_channel, blue_channel)

    Optionally, you can also specify a tuple that includes scaling information,
    in the form of (array_to_plot, min_value_to_scale, max_value_to_scale).
    """
    red_channel = _scale_image(red_channel)
    blue_channel = _scale_image(blue_channel)
    if green_channel is None:
        green_channel = na.zeros(red_channel.shape, dtype='uint8')
    else:
        green_channel = _scale_image(green_channel)
    if alpha_channel is None:
        alpha_channel = na.zeros(red_channel.shape, dtype='uint8') + 255
    else:
        alpha_channel = _scale_image(alpha_channel) 
    image = na.array([red_channel, green_channel, blue_channel, alpha_channel])
    image = image.transpose().copy() # Have to make sure it's contiguous 
    au.write_png(image, fn)

def write_bitmap(bitmap_array, filename, max_val = None):
    """
    This accepts a three- or four-channel *bitmap_array*.  If the image is not
    already uint8, it will be scaled and converted.  If it is not four channel,
    a fourth alpha channel will be added and set to fully opaque.  The
    resultant image will be directly written to *filename* as a PNG with no
    colormap applied.  *max_val* is a value used if the array is passed in as
    anything other than uint8; it will be the value used for scaling and
    clipping when the array is converted.

    Additionally, the minimum is assumed to be zero; this makes it primarily
    suited for the results of volume rendered images, rather than misaligned
    projections.
    """
    if bitmap_array.dtype != na.uint8:
        if max_val is None: max_val = bitmap_array.max()
        bitmap_array = na.clip(bitmap_array / max_val, 0.0, 1.0) * 255
        bitmap_array = bitmap_array.astype("uint8")
    if len(bitmap_array.shape) != 3 or bitmap_array.shape[-1] not in (3,4):
        raise RuntimeError
    if bitmap_array.shape[-1] == 3:
        s1, s2 = bitmap_array.shape[:2]
        alpha_channel = 255*na.ones((s1,s2,1), dtype='uint8')
        bitmap_array = na.concatenate([bitmap_array, alpha_channel], axis=-1)
    au.write_png(bitmap_array.copy(), filename)

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
