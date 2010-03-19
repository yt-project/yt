"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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
from yt.raven import FixedResolutionBuffer, ObliqueFixedResolutionBuffer
from yt.lagos import data_object_registry, AMRProjBase, AMRSliceBase, \
                     x_dict, y_dict

class VariableMeshPanner(object):
    _buffer = None

    def __init__(self, source, size, field, callback = None, pf = None):
        if not isinstance(source, (AMRProjBase, AMRSliceBase)):
            raise RuntimeError
        if callback is None:
            callback = lambda a: None
        self.callback = callback
        self.size = size
        self.source = source
        self.field = field
        self.xlim, self.ylim = self.bounds

    @property
    def bounds(self):
        DLE, DRE = self.pf["DomainLeftEdge"], self.pf["DomainRightEdge"]
        ax = self.source.axis
        xax, yax = x_dict[ax], y_dict[ax]
        xbounds = DLE[xax], DRE[xax]
        ybounds = DLE[yax], DRE[yax]
        return (xbounds, ybounds)

    @property
    def width(self):
        Wx = self.xlim[1] - self.xlim[0]
        Wy = self.ylim[1] - self.ylim[0]
        return (Wx, Wy)

    def zoom(self, factor):
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        nWx, nWy = Wx/factor, Wy/factor
        self.xlim = (centerx - nWx*0.5, centerx + nWx*0.5)
        self.ylim = (centery - nWy*0.5, centery + nWy*0.5)
        self.callback(self.buffer)

    def pan(self, deltas):
        self.xlim = (self.xlim[0] + deltas[0], self.xlim[1] + deltas[0])
        self.ylim = (self.ylim[0] + deltas[1], self.ylim[1] + deltas[1])
        self.callback(self.buffer)

    def pan_x(self, delta):
        self.pan( (delta, 0.0) )

    def pan_y(self, delta):
        self.pan( (0.0, delta) )

    def pan_rel(self, deltas):
        Wx, Wy = self.width
        px = deltas[0] * Wx
        py = deltas[1] * Wy
        self.pan( (px, py) )

    def pan_rel_x(self, delta):
        self.pan_rel( (delta, 0.0) )

    def pan_rel_y(self, delta):
        self.pan_rel( (0.0, delta) )

    @property
    def buffer(self):
        my_bounds = self.xlim + self.ylim # Tuples concatenate
        # If the buffer is None, we just kill it regardless.
        if self._buffer is None:
            self._regenerate_buffer()
        elif self._buffer.bounds != my_bounds:
            self._regenerate_buffer()
        return self._buffer[self.field]

    def _regenerate_buffer(self):
        new_buffer = FixedResolutionBuffer(
            self.source, self.xlim + self.ylim,
            self.size)
        self._buffer = new_buffer

    def set_low_high(self, low, high):
        print "Setting low, high", low, high
        self.xlim = (low[0], high[0])
        self.ylim = (low[1], high[1])
        b = na.log10(self.buffer)
        mi, ma = b.min(), b.max()
        print "Returning buffer with extrema",
        print mi, ma
        b = (b - mi)/(ma - mi)
        return b

data_object_registry["image_panner"] = VariableMeshPanner
