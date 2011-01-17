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

class RemoteWindowedVariableMeshController(MultipleWindowVariableMeshPanner):
    def __init__(self, source, mec = None):
        """
        This panner controls remote windowed panners.  It requires a *source*,
        which will be pickled and sent to the remote panners, which it will
        create as requested.  If not supplied with a *mec* (an IPython
        MultiEngineClient) it will create one itself.
        """
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
        """
        This will create a new remote WindowedVariableMeshImagePanner.  The
        *args* and *kwargs* supplied here will be passed over, but the *source*
        argument is implicitly handled by this routine.
        """
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
                     "pan_rel_x", "pan_rel_y", "set_limits"]

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

class ProxySource(object):
    # This proxies only the things we know we need
    # Note that we assume we will only have a single engine.
    def __init__(self, mec, idnum, source_varname):
        self.mec = mec
        self.idnum = idnum
        self.source_varname = source_varname
        self.mec.execute("_tmp_%s = %s.axis" % (
            self.idnum, self.source_varname))
        self.axis = self.mec.pull("_tmp_%s" % self.idnum)[0]

    def keys(self):
        self.mec.execute("_tmp_%s = %s.keys()" % (
            self.idnum, self.source_varname))
        keys = self.mec.pull("_tmp_%s" % self.idnum)[0]
        dd = dict( (k, None) for k in keys )
        return dd

    @property
    def pf(self):
        self.mec.execute("_tmp_%s = %s.pf.domain_left_edge" % (
            self.idnum, self.source_varname))
        DLE = self.mec.pull("_tmp_%s" % self.idnum)[0]
        self.mec.execute("_tmp_%s = %s.pf.domain_right_edge" % (
            self.idnum, self.source_varname))
        DRE = self.mec.pull("_tmp_%s" % self.idnum)[0]
        return dict(DomainLeftEdge = DLE, DomainRightEdge = DRE)

class ProxyFixedResolutionBuffer(dict):
    pass

class NonLocalDataImagePanner(VariableMeshPanner):
    def __init__(self, mec, source_varname, size, field,
                 callback = None, viewport_callback = None):
        self.source_varname = source_varname
        self._var_name = "_image_panner_%s" % (id(self))
        self.mec = mec
        self.mec.execute("import yt.extensions.image_panner")
        self.mec.execute("%s = yt.extensions.image_panner.VariableMeshPanner(" % (
                        self._var_name) +
                          "%s, (%s, %s), '%s')" % (
                        source_varname, size[0], size[1], field))

        ps = ProxySource(mec, id(self), source_varname)
        self._prfb = ProxyFixedResolutionBuffer()

        VariableMeshPanner.__init__(self, ps, size, field,
                        callback, viewport_callback)

    def _regenerate_buffer(self):
        args = (self.xlim, self.ylim)
        self.mec.push({'_tmp_%s' % id(self) : args}, block=False)
        self.mec.execute("%s.set_limits(*_tmp_%s)" % (self._var_name, id(self)),
                         block=False)
        self.mec.execute("_tmp_%s = %s.buffer" % (id(self), self._var_name),
                         block=False)
        self._prfb[self.field] = self.mec.pull("_tmp_%s" % (id(self)))[0]
        self._prfb.bounds = self.xlim + self.ylim
        self._buffer = self._prfb

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
        from yt.utilities.amr_utils import write_png_to_file
        from yt.visualization.image_writer import map_to_colors
        image = na.log10(val)
        mi = na.nanmin(image[~na.isinf(image)])
        ma = na.nanmax(image[~na.isinf(image)])
        color_bounds = mi, ma
        image = (image - color_bounds[0])/(color_bounds[1] - color_bounds[0])
        to_plot = map_to_colors(image, "algae")
        to_plot = na.clip(to_plot, 0, 255)
        tf = tempfile.TemporaryFile()
        write_png_to_file(to_plot, tf)
        tf.seek(0)
        s = tf.read()
        response_body = "data:image/png;base64," + base64.encodestring(s)
        tf.close()
        self.transport.append(response_body)

class PanningCeleritasStreamer(object):
    _initialized = False
    def __init__(self, tile_id, cmap = "algae", port = 9988,
                 zlim = (0.0, 1.0), take_log = True):
        """
        This is an in-development mechanism for supplying buffers to a
        Celeritas server.
        """
        self.tile_id = tile_id
        self._port = port
        self.cmap = cmap
        self.zlim = zlim
        self.take_log = True

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
        if self.take_log:
            vv = na.log10(val)
        else:
            vv = val.copy()
        na.subtract(vv, self.zlim[0], vv)
        na.divide(vv, (self.zlim[1]-self.zlim[0]), vv)
        new_buf = self.cmap(vv)[:,:,:3]
        na.multiply(new_buf, 255.0, new_buf)
        new_buf = new_buf.astype('uint8')
        self.cs.readFromRGBMemAndSend(new_buf)
