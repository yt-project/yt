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
import types, os
from yt.raven import FixedResolutionBuffer, ObliqueFixedResolutionBuffer
from yt.lagos import data_object_registry, AMRProjBase, AMRSliceBase, \
                     x_dict, y_dict

class SourceTypeException(Exception):
    def __init__(self, source):
        self.source = source

    def __repr__(self):
        return "Wrong type: %s, %s" % (self.source, type(self.source))

class VariableMeshPanner(object):
    _buffer = None

    def __init__(self, source, size, field, callback = None):
        if not isinstance(source, (AMRProjBase, AMRSliceBase)):
            print "PROBLEM WITH SOURCE", source, type(source)
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

class WindowedVariableMeshPanner(VariableMeshPanner):

    def __init__(self, source, full_size, my_size, start_indices,
                 field, callback = None):
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

class RemoteWindowedVariableMeshController(MultipleWindowVariableMeshPanner):
    def __init__(self, source, mec = None):
        if mec is None:
            from IPython.kernel.client import get_multiengine_client
            mec = get_multiengine_client()
        self.mec = mec
        self.mec.execute("import yt.extensions.image_panner")
        self._var_name = "_image_panner_%s" % (id(self))
        self._pf_name = "_pf_%s" % (id(self))
        self._source_name = "_source_%s" % (id(self))
        self.source = source
        self.mec.execute("from yt.mods import *")
        self.mec.execute("from yt.funcs import iterable")
        self.mec.push({self._pf_name: self.pf})
        self.mec.execute("%s.h" % self._pf_name)
        self.mec.push({self._source_name: self.source})
        # Now, because the double pickling tosses a PF hash reference inside
        # the unpickled object, we work around it a little
        self.mec.execute("while iterable(%s): %s = %s[1]" % (
            self._source_name, self._source_name, self._source_name))
        self.windows = []

    def add_window(self, *args, **kwargs):
        engine_id = len(self.windows)
        an = "_args_%s" % id(self)
        kn = "_kwargs_%s" % id(self)
        if 'callback' not in kwargs:
            kwargs['callback'] = ImageSaver(engine_id)
        self.mec.push({an: args, kn: kwargs}, engine_id)
        exec_string = "%s = %s.h.windowed_image_panner(%s, *%s, **%s)" % (
            self._var_name, self._pf_name, self._source_name, an, kn)
        self.mec.execute(exec_string, engine_id)
        self.windows.append(WindowedVariableMeshPannerProxy(
            self.mec, engine_id, self._var_name, id(self)))

data_object_registry["remote_image_panner"] = RemoteWindowedVariableMeshController

_wrapped_methods = ["zoom", "pan", "pan_x", "pan_y", "pan_rel",
                     "pan_rel_x", "pan_rel_y", "set_low_high"]

class WindowedVariableMeshPannerProxy(object):
    class __metaclass__(type):
        def __new__(cls, name, b, d):
            # We add on a bunch of proxy functions
            def return_proxy(fname):
                def func(self, *args, **kwargs):
                    vn = "_ret_%s" % self._cid
                    an = "_args_%s" % self._cid
                    kn = "_kwargs_%s" % self._cid
                    self.mec.push({an: args, kn: kwargs}, self.engine_id)
                    exec_string = "%s = %s.%s(*%s, **%s)" % (
                        vn, self._var_name, fname, an, kn)
                    print "Executing %s on %s" % (exec_string, self.engine_id)
                    self.mec.execute(exec_string, self.engine_id)
                    return self.mec.pull(vn, self.engine_id)
                return func
            new_dict = {}
            new_dict.update(d)
            for f in _wrapped_methods:
                print "Constructing proxy for", f
                new_dict[f] = return_proxy(f)
            return type.__new__(cls, name, b, new_dict)

    def __init__(self, mec, engine_id, var_name, cid):
        # mec here is, optionally, an instance of MultiEngineClient
        self._var_name = var_name
        self._cid = cid
        self.engine_id = engine_id
        self.mec = mec

    @property
    def bounds(self):
        vn = "_ret_%s" % self._cid
        self.mec.execute("%s = %s.bounds" % (vn, self._var_name),
                         self.engine_id)
        return self.mec.pull(vn, self.engine_id)

    @property
    def width(self):
        vn = "_ret_%s" % self._cid
        self.mec.execute("%s = %s.width" % (vn, self._var_name),
                         self.engine_id)
        return self.mec.pull(vn, self.engine_id)

    @property
    def buffer(self):
        vn = "_ret_%s" % self._cid
        self.mec.execute("%s = %s.buffer" % (vn, self._var_name),
                         self.engine_id)
        return self.mec.pull(vn, self.engine_id)

    def _regenerate_buffer(self):
        return

    def _run_callback(self):
        self.mec.execute("%s._regenerate_buffer()" % self._var_name,
                         self.engine_id)
        self.mec.execute("%s.callback(%s.buffer)" % (
            self._var_name, self._var_name), self.engine_id)

class ImageSaver(object):
    def __init__(self, tile_id):
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

class PanningCeleritasStreamer(object):
    _initialized = False
    def __init__(self, tile_id, cmap = "algae", port = 9988,
                 zlim = (0.0, 1.0)):
        self.tile_id = tile_id
        self._port = port
        self.cmap = cmap
        self.zlim = zlim

    def initialize(self, shape):
        if isinstance(self.cmap, types.StringTypes):
            import matplotlib.cm
            self.cmap = matplotlib.cm.get_cmap(self.cmap)

        import celeritas_streamer
        self.cs = celeritas_streamer.CeleritasStream()
        #print "Setting shape: %s and port: %s in %s" % (
        #    shape, self._port, os.getpid())
        self.cs.setSize(*shape)
        self.cs.setLocalPort(self._port)
        self.cs.initialize()
        self._initialized = True

    def __call__(self, val):
        if not self._initialized: self.initialize(val.shape)
        vv = val.copy()
        na.subtract(vv, self.zlim[0], vv)
        na.divide(vv, (self.zlim[1]-self.zlim[0]), vv)
        new_buf = self.cmap(vv)
        na.multiply(new_buf, 255.0, new_buf)
        new_buf = new_buf[:,:,:3].astype('uint8')
        self.cs.readFromRGBMemAndSend(new_buf)
