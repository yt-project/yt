"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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
import types, os
from yt.visualization.fixed_resolution import \
    FixedResolutionBuffer, ObliqueFixedResolutionBuffer
from yt.data_objects.data_containers import \
    data_object_registry, AMRProjBase, AMRSliceBase
from yt.utilities.definitions import \
    x_dict, y_dict
from yt.funcs import *

class VariableMeshPanner(object):
    _buffer = None
    _hold = False

    def __init__(self, source, size, field, callback = None,
                 viewport_callback = None):
        """
        This class describes a meta-view onto a 2D data interface.  You provide
        it with a slice or a projection (*source*), a tuple of the x and y
        resolution of the final image you want (*size) and a *field* to
        display, and then it will generate a new buffer each time you call its
        methods for moving the window around and changing its size.  *callback*
        is a routine called with the new buffer each time the buffer changes
        and *viewport_callback* is called with the new *xlim* and *ylim* values
        each time the viewport changes.
        """
        #if not isinstance(source, (AMRProjBase, AMRSliceBase)):
        #    raise RuntimeError
        if callback is None:
            callback = lambda a: None
        self.callback = callback
        if viewport_callback is None:
            viewport_callback = lambda a, b: None
        self.viewport_callback = viewport_callback
        self.size = size
        self.source = source
        self.field = field
        self.xlim, self.ylim = self.bounds

    def _run_callbacks(self):
        if self._hold: return
        self.callback(self.buffer)
        self.viewport_callback(self.xlim, self.ylim)

    @property
    def bounds(self):
        if not hasattr(self, 'pf'): self.pf = self.source.pf
        DLE, DRE = self.pf.domain_left_edge, self.pf.domain_right_edge
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
        """
        This zooms the window by *factor*.
        """
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        nWx, nWy = Wx/factor, Wy/factor
        self.xlim = (centerx - nWx*0.5, centerx + nWx*0.5)
        self.ylim = (centery - nWy*0.5, centery + nWy*0.5)
        self._run_callbacks()

    def pan(self, deltas):
        """
        This accepts a tuple of *deltas*, composed of (delta_x, delta_y) that
        will pan the window by those values in absolute coordinates.
        """
        self.xlim = (self.xlim[0] + deltas[0], self.xlim[1] + deltas[0])
        self.ylim = (self.ylim[0] + deltas[1], self.ylim[1] + deltas[1])
        self._run_callbacks()

    def pan_x(self, delta):
        """
        This pans the window by *delta* in the x direction.
        """
        self.pan( (delta, 0.0) )

    def pan_y(self, delta):
        """
        This pans the window by *delta* in the y direction.
        """
        self.pan( (0.0, delta) )

    def pan_rel(self, deltas):
        """
        This accepts a tuple of *deltas*, composed of (delta_x, delta_y) that
        will pan the window by those values in relative values, relative to the
        current window view.
        """
        Wx, Wy = self.width
        px = deltas[0] * Wx
        py = deltas[1] * Wy
        self.pan( (px, py) )

    def pan_rel_x(self, delta):
        """
        This pans the window by *delta* in the x direction, where *delta* is
        relative to the current window view.
        """
        self.pan_rel( (delta, 0.0) )

    def pan_rel_y(self, delta):
        """
        This pans the window by *delta* in the y direction, where *delta* is
        relative to the current window view.
        """
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
        """
        This accepts *low* in the format (x_min, y_min) and *high* in the
        format (x_max, y_max) and changes the viewport appropriately.  In
        breaking from tradition, it also returns the buffer.
        """
        self.xlim = (low[0], high[0])
        self.ylim = (low[1], high[1])
        return na.log10(self.buffer)

    def set_width(self, width):
        """
        This sets the width based on the current center.
        """
        if not iterable(width): width = (width, width)
        Wx, Wy = self.width
        centerx = self.xlim[0] + Wx*0.5
        centery = self.ylim[0] + Wy*0.5
        self.xlim = (centerx - width[0]/2.0,
                     centerx + width[0]/2.0)
        self.ylim = (centery - width[1]/2.0,
                     centery + width[1]/2.0)
        self._run_callbacks()

    def set_limits(self, xlim, ylim):
        """
        This accepts a new *xlim* and *ylim*.
        """
        self.xlim = xlim
        self.ylim = ylim
        self._run_callbacks()

    def set_center(self, center):
        if len(center) == 2:
            centerx, centery = center
        elif len(center) == 3:
            centerx = center[x_dict[self.source.axis]]
            centery = center[y_dict[self.source.axis]]
        else:
            raise RuntimeError
        Wx, Wy = self.width
        self.xlim = (centerx - Wx*0.5, centerx + Wx*0.5)
        self.ylim = (centery - Wy*0.5, centery + Wy*0.5)
        self._run_callbacks()

data_object_registry["image_panner"] = VariableMeshPanner

class WindowedVariableMeshPanner(VariableMeshPanner):

    def __init__(self, source, full_size, my_size, start_indices,
                 field, callback = None):
        """
        This image panner accepts a *full_size*, which describes the full size
        of the image which it will be displaying a portion of.  *my_size* is
        the size that this window will be responsible for, and *start_indices*
        is a tuple of indices that this window begins at.  *field* and
        *callback* function as in the vanilla image panner.
        """
        self.my_size = my_size
        self.start_indices = start_indices
        VariableMeshPanner.__init__(self, source, full_size, field, callback)

    def _regenerate_buffer(self):
        dx = (self.xlim[1] - self.xlim[0])/self.size[0]
        dy = (self.ylim[1] - self.ylim[0])/self.size[1]
        my_lim = (self.xlim[0] + dx*self.start_indices[0],
                  self.xlim[0] + dx*(self.start_indices[0] + self.my_size[0]),
                  self.ylim[0] + dx*self.start_indices[1],
                  self.ylim[0] + dx*(self.start_indices[1] + self.my_size[1]))
        new_buffer = FixedResolutionBuffer(self.source, my_lim, self.my_size)
        self._buffer = new_buffer

data_object_registry["windowed_image_panner"] = WindowedVariableMeshPanner

class MultipleWindowVariableMeshPanner(object):
    def __init__(self, windows):
        """
        This panner is exclusively a controller.  It accepts a list of
        *windows*, to which it will issue commands.  It knows nothing other
        than their interfaces.
        """
        self.windows = windows

    def zoom(self, factor):
        for w in self.windows: w.zoom(factor)

    def pan(self, deltas):
        for w in self.windows: w.pan(factor)

    def pan_x(self, delta):
        for w in self.windows: w.pan_x(delta)

    def pan_y(self, delta):
        for w in self.windows: w.pan_y(delta)

    def pan_rel(self, deltas):
        for w in self.windows: w.pan_rel(deltas)

    def pan_rel_x(self, delta):
        for w in self.windows: w.pan_rel_x(delta)

    def pan_rel_y(self, delta):
        for w in self.windows: w.pan_rel_y(delta)

    def set_limits(self, xlim, ylim):
        for w in self.windows: w.set_limits(xlim, ylim)

class ImageSaver(object):
    def __init__(self, tile_id):
        """
        This is a sample callback for an image panner that will save to a file
        named ``wimage_%03i.png`` where the final variable is specified here as
        *tile_id*.
        """
        self.tile_id = tile_id

        import matplotlib;matplotlib.use("Agg");import pylab
        self.pylab = pylab
        self.pylab.clf()
        fig = pylab.gcf()
        fig.subplots_adjust(left=0.0, bottom=0.0,
                            right=1.0, top=1.0,
                            wspace=0.0, hspace=0.0)
        fig.set_dpi(100.0)
        fig.set_size_inches((5.12, 5.12))

    def __call__(self, val):
        self.pylab.clf()
        self.pylab.imshow(na.log10(val), interpolation='nearest')
        self.pylab.savefig("wimage_%03i.png" % self.tile_id)

class TransportAppender(object):
    def __init__(self, transport):
        """
        This is a simple callback that appends the latest buffer to the
        supplied *transport* object.
        """
        self.transport = transport

    def __call__(self, val):
        from yt.utilities.lib import write_png_to_string
        from yt.visualization.image_writer import map_to_colors
        image = na.log10(val)
        mi = na.nanmin(image[~na.isinf(image)])
        ma = na.nanmax(image[~na.isinf(image)])
        color_bounds = mi, ma
        image = (image - color_bounds[0])/(color_bounds[1] - color_bounds[0])
        to_plot = map_to_colors(image, "algae")
        to_plot = na.clip(to_plot, 0, 255)
        s = write_png_to_string(to_plot)
        response_body = "data:image/png;base64," + base64.encodestring(s)
        tf.close()
        self.transport.append(response_body)

