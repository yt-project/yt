"""
A simple leaflet-based pannable map server

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
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
from yt.utilities.amr_utils import write_png_to_string, get_color_bounds

import yt.utilities.bottle as bottle

from yt.funcs import *
local_dir = os.path.dirname(__file__)

def exc_writeout(f):
    import traceback
    @wraps(f)
    def func(*args, **kwargs):
        try:
            rv = f(*args, **kwargs)
            return rv
        except Exception as e:
            traceback.print_exc(None, open("temp.exc", "w"))
            raise
    return func

class PannableMapServer(object):
    _widget_name = "pannable_map"
    def __init__(self, data, field, route_prefix = ""):
        self.data = data
        self.pf = data.pf
        self.field = field
        
        bottle.route("%s/map/:L/:x/:y.png" % route_prefix)(self.map)
        bottle.route("%s/" % route_prefix)(self.index)
        bottle.route("%s/index.html" % route_prefix)(self.index)
        bottle.route("%s/static/:filename#.+#" % route_prefix)(self.static)
        # This is a double-check, since we do not always mandate this for
        # slices:
        self.data[self.field] = self.data[self.field].astype("float64")

    def map(self, L, x, y):
        mylog.debug("Hello!")
        dd = 1.0 / (2.0**(int(L)))
        relx = int(x) * dd
        rely = int(y) * dd
        DW = (self.pf.domain_right_edge - self.pf.domain_left_edge)
        xl = self.pf.domain_left_edge[0] + relx * DW[0]
        yl = self.pf.domain_left_edge[1] + rely * DW[1]
        xr = xl + dd*DW[0]
        yr = yl + dd*DW[1]
        frb = FixedResolutionBuffer(self.data, (xl, xr, yl, yr), (256, 256))
        cmi, cma = get_color_bounds(self.data['px'], self.data['py'],
                                    self.data['pdx'], self.data['pdy'],
                                    self.data[self.field],
                                    self.pf.domain_left_edge[0],
                                    self.pf.domain_right_edge[0],
                                    self.pf.domain_left_edge[1],
                                    self.pf.domain_right_edge[1],
                                    dd*DW[0] / (64*256),
                                    dd*DW[0])
        if self.pf.field_info[self.field].take_log:
            cmi = na.log10(cmi)
            cma = na.log10(cma)
            to_plot = apply_colormap(na.log10(frb[self.field]), color_bounds = (cmi, cma))
        else:
            to_plot = apply_colormap(frb[self.field], color_bounds = (cmi, cma))
        rv = write_png_to_string(to_plot)
        return rv

    def index(self):
        return bottle.static_file("map_index.html", root=self._path("html"))

    def static(self, filename):
        return bottle.static_file(filename, root=self._path("html/leaflet"))

    def _path(self, fn):
        return os.path.join(local_dir, fn) 
