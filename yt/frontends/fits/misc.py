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

import __builtin__
import aplpy
from yt.utilities.fits_image import FITSImageBuffer
from yt.funcs import fix_axis, ensure_list
import astropy.wcs as pywcs
from yt.utilities.exceptions import \
    YTNotInsideNotebook
import matplotlib.pyplot as plt

axis_wcs = [[1,2],[0,2],[0,1]]

plot_method_list = ["recenter","refresh","add_colorbar",
                    "remove_colorbar"]

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
            self.plots[field] = aplpy.FITSFigure(self.buffer[field], **kwargs)
            self.plots[field].set_auto_refresh(False)
        self._setup_plot_methods()
        self.set_font(family="serif", size=15)
        for v in self.values():
            v.show_colorscale()
        plt.close("all")

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

    def show(self):
        r"""This will send any existing plots to the IPython notebook.
        function name.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will send any existing plots to the
        notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Examples
        --------

        >>> from yt.mods import SlicePlot
        >>> slc = SlicePlot(pf, "x", ["Density", "VelocityMagnitude"])
        >>> slc.show()

        """
        if "__IPYTHON__" in dir(__builtin__):
            from IPython.display import display
            for k, v in sorted(self.plots.iteritems()):
                display(v._figure)
        else:
            raise YTNotInsideNotebook

class FITSSlicePlot(FITSPlot):
    def __init__(self, ds, axis, fields, coord=None, field_parameters=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis)
        if coord is None:
            coord = ds.domain_center.ndarray_view()[axis]
        slc = ds.slice(axis, coord, field_parameters=field_parameters)
        data = {}
        for field in fields:
            data[field] = slc.to_frb((1.0,"unitary"), ds.domain_dimensions[axis_wcs[axis]])[field]
        super(FITSSlicePlot, self).__init__(ds, data, axis, fields, **kwargs)

class FITSProjectionPlot(FITSPlot):
    def __init__(self, ds, axis, fields, weight_field=None, data_source=None,
                 field_parameters=None, **kwargs):
        fields = ensure_list(fields)
        axis = fix_axis(axis)
        prj = ds.proj(fields[0], axis, weight_field=weight_field, data_source=data_source)
        data = {}
        for field in fields:
            data[field] = prj.to_frb((1.0,"unitary"), ds.domain_dimensions[axis_wcs[axis]])[field]
        super(FITSProjectionPlot, self).__init__(ds, data, axis, fields, **kwargs)

