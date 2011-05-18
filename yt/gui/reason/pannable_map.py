"""
A simple leaflet-based pannable map server

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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
import os
import numpy as na

from yt.visualization.image_writer import apply_colormap
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.utilities.amr_utils import write_png_to_string

import yt.utilities.bottle as bottle

from yt.funcs import *
local_dir = os.path.dirname(__file__)

def get_realistic_bounds(field, pdx, width):
    # We don't want to scale to anything that's as big or bigger than the
    # possible viewing area.  This should be unique across tiles on this level.
    min_dx = width / 256
    ff = na.log10(field[(pdx > min_dx / 64.0) & (pdx < width * 2.0)])
    return ff.min(), ff.max()

class PannableMapServer(object):
    def __init__(self, data):
        self.data = data
        self.pf = data.pf
        bottle.route("/map/:L/:x/:y.png")(self.map)
        bottle.route("/")(self.index)
        bottle.route("/index.html")(self.index)
        bottle.route("/static/:filename#.+#")(self.static)

    def map(self, L, x, y):
        dd = 1.0 / (2.0**(int(L)-1))
        relx = int(x) * dd
        rely = int(y) * dd
        DW = (self.pf.domain_left_edge + self.pf.domain_right_edge)/2.0
        xp = self.pf.domain_left_edge[0] + relx * DW[0]
        yp = self.pf.domain_left_edge[1] + rely * DW[1]
        frb = FixedResolutionBuffer(self.data, (xp, xp+dd*DW[0], yp, yp+dd*DW[1]), (256, 256))
        cmi, cma = get_realistic_bounds(self.data["Density"], self.data["pdx"], dd*DW[0])
        to_plot = apply_colormap(na.log10(frb["Density"]), color_bounds = (cmi, cma))
        rv = write_png_to_string(to_plot)
        return rv

    def index(self):
        return bottle.static_file("map_index.html", root=self._path("html"))

    def static(self, filename):
        return bottle.static_file(filename, root=self._path("html/leaflet"))

    def _path(self, fn):
        return os.path.join(local_dir, fn) 
