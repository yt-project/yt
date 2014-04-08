"""
FITS-specific miscellaneous functions
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import aplpy
from yt.utilities.fits_image import FITSImageBuffer
from yt.funcs import fix_axis
import astropy.wcs as pywcs

axis_wcs = [[1,2],[0,2],[0,1]]

plot_method_list = ["recenter","show_colorscale","show_grayscale",
                    "refresh","add_colorbar","remove_colorbar"]

def plot_method(method, plots):
    def _method(*args, **kwargs):
        for plot in plots.values():
            getattr(plot, method)(*args, **kwargs)
        return
    return _method

class FITSPlot(object):
    def __init__(self, ds, data, axis, fields, **kwargs):
        self.ds = ds
        self.fields = fields
        self.plots = {}
        w = pywcs.WCS(naxis=2)
        w.wcs.crpix = self.ds.wcs.wcs.crpix[axis_wcs[axis]]
        w.wcs.cdelt = self.ds.wcs.wcs.cdelt[axis_wcs[axis]]
        w.wcs.crval = self.ds.wcs.wcs.crval[axis_wcs[axis]]
        w.wcs.cunit = [str(self.ds.wcs.wcs.cunit[idx]) for idx in axis_wcs[axis]]
        w.wcs.ctype = [self.ds.wcs.wcs.ctype[idx] for idx in axis_wcs[axis]]
        self.buffer = FITSImageBuffer(data, fields=fields, wcs=w)
        for field in self.fields:
            self.plots[field] = aplpy.FITSFigure(self.buffer[field],
                                                 **kwargs)
            self.plots[field].set_auto_refresh(False)
        self._setup_plot_methods()
        self.set_font(family="serif", size=15)

    def _setup_plot_methods(self):
        for method in plot_method_list:
            self.__dict__[method] = plot_method(method, self.plots)

    def __getitem__(self, key):
        return self.plots[key]

    def keys(self):
        return self.plots.keys()

    def values(self):
        return self.plots.values()

    def items(self):
        return self.plots.items()

    def set_font(self, **kwargs):
        for plot in self.keys():
            self[plot].axis_labels.set_font(**kwargs)
            self[plot].tick_labels.set_font(**kwargs)

    def set_stretch(self, name, stretch):
        self[name].show_colorscale(stretch=stretch)

    def set_zlim(self, name, zmin, zmax):
        self[name].show_colorscale(vmin=zmin, vmax=zmax)

class FITSSlicePlot(FITSPlot):
    def __init__(self, ds, axis, fields, coord=None, **kwargs):
        axis = fix_axis(axis)
        if coord is None:
            coord = ds.domain_center.ndarray_view()[axis]
        slc = ds.slice(axis, coord)
        data = {}
        for field in fields:
            data[field] = slc[field].reshape(ds.domain_dimensions[axis_wcs[axis]]).transpose()
        super(FITSSlicePlot, self).__init__(ds, data, axis, fields, **kwargs)

class FITSProjectionPlot(FITSPlot):
    def __init__(self, ds, axis, fields, weight_field=None, data_source=None, **kwargs):
        axis = fix_axis(axis)
        prj = ds.proj(fields[0], axis, weight_field=weight_field, data_source=data_source)
        data = {}
        for field in fields:
            data[field] = prj[field].reshape(ds.domain_dimensions[axis_wcs[axis]]).transpose()
        super(FITSProjectionPlot, self).__init__(ds, data, axis, fields, **kwargs)

