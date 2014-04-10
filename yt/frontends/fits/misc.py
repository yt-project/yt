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

import numpy as np
from yt.funcs import fix_axis, ensure_list, iterable
from yt.visualization.plot_window import AxisAlignedSlicePlot, \
    OffAxisSlicePlot, ProjectionPlot, OffAxisProjectionPlot

def force_aspect(ax,aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def set_onaxis_wcs(pw):
    return
    if pw.axis == pw.ds.ra_axis:
        xname = "Dec"
        yname = pw.ds.vel_name
        xunit = str(pw.ds.wcs_2d.wcs.cunit[1])
        yunit = str(pw.ds.wcs_1d.wcs.cunit[0])
    elif pw.axis == pw.ds.dec_axis:
        xname = "RA"
        yname = pw.ds.vel_name
        xunit = str(pw.ds.wcs_2d.wcs.cunit[0])
        yunit = str(pw.ds.wcs_1d.wcs.cunit[0])
    elif pw.axis == pw.ds.vel_axis:
        xname = "RA"
        yname = "Dec"
        xunit = str(pw.ds.wcs_2d.wcs.cunit[0])
        yunit = str(pw.ds.wcs_2d.wcs.cunit[1])

    for k,v in pw.plots.iteritems():
        v.axes.set_xlabel(r"%s (%s)" % (xname, xunit))
        v.axes.set_ylabel(r"%s (%s)" % (yname, yunit))
        v.axes.set_aspect('auto')

class FITSSlicePlot(AxisAlignedSlicePlot):

    def __init__(self, ds, axis, fields, set_wcs=True, **kwargs):

        if isinstance(axis, basestring):
            if axis in ds.axis_names:
                axis = ds.axis_names[axis]
        self.axis = fix_axis(axis)
        self.ds = ds
        self.set_wcs = set_wcs
        super(FITSSlicePlot, self).__init__(ds, axis, fields, origin="native", **kwargs)
        self.set_axes_unit("pixel")

    def _set_wcs(self):
        if self.set_wcs:
            set_onaxis_wcs(self)

    def show(self):
        self._set_wcs()
        super(FITSSlicePlot, self).show()

    def save(self, *args, **kwargs):
        self._set_wcs()
        super(FITSSlicePlot, self).save(*args, **kwargs)

class FITSOffAxisSlicePlot(OffAxisSlicePlot):

    def __init__(self, ds, normal, fields, set_wcs=True, **kwargs):

        self.ds = ds
        my_normal = normal
        if ds.xyv_data:
            if len(normal) > 2:
                raise NotImplementedError("Normal vector must be in two dimensions for this dataset!")
            my_normal = np.zeros((3))
            my_normal[ds.ra_axis] = normal[0]
            my_normal[ds.dec_axis] = normal[1]

        super(FITSOffAxisSlicePlot, self).__init__(ds, my_normal, fields, **kwargs)
        self.set_axes_unit("pixel")

class FITSProjectionPlot(ProjectionPlot):

    def __init__(self, ds, axis, fields, set_wcs=True, **kwargs):

        self.ds = ds
        if isinstance(axis, basestring):
            if axis in ds.axis_names:
                axis = ds.axis_names[axis]
        self.axis = fix_axis(axis)
        self.set_wcs = set_wcs

        super(FITSProjectionPlot, self).__init__(ds, axis, fields, origin="native", **kwargs)
        self.set_axes_unit("pixel")

    def _set_wcs(self):
        if self.set_wcs:
            set_onaxis_wcs(self)

    def show(self):
        self._set_wcs()
        super(FITSProjectionPlot, self).show()

    def save(self, *args, **kwargs):
        self._set_wcs()
        super(FITSProjectionPlot, self).save(*args, **kwargs)

class FITSOffAxisProjectionPlot(OffAxisProjectionPlot):

    def __init__(self, ds, normal, fields, set_wcs=True, **kwargs):

        self.ds = ds
        my_normal = normal
        if ds.xyv_data:
            if len(normal) > 2:
                raise ValueError("Normal vector must be in two dimensions for this dataset!")
            my_normal = np.zeros((3))
            my_normal[ds.ra_axis] = normal[0]
            my_normal[ds.dec_axis] = normal[1]

        super(FITSOffAxisProjectionPlot, self).__init__(ds, my_normal, fields, axes_unit="pixel", **kwargs)


