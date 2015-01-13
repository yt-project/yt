"""
Miscellaneous FITS routines
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
import base64
from yt.extern.six.moves import StringIO
from yt.fields.derived_field import ValidateSpatial
from yt.utilities.on_demand_imports import _astropy
from yt.funcs import mylog, get_image_suffix
from yt.visualization._mpl_imports import FigureCanvasAgg

import os

def _make_counts(emin, emax):
    def _counts(field, data):
        e = data["event_energy"].in_units("keV")
        mask = np.logical_and(e >= emin, e < emax)
        x = data["event_x"][mask]
        y = data["event_y"][mask]
        z = np.ones(x.shape)
        pos = np.array([x,y,z]).transpose()
        img = data.deposit(pos, method="count")
        if data.has_field_parameter("sigma"):
            sigma = data.get_field_parameter("sigma")
        else:
            sigma = None
        if sigma is not None and sigma > 0.0:
            kern = _astropy.conv.Gaussian2DKernel(stddev=sigma)
            img[:,:,0] = _astropy.conv.convolve(img[:,:,0], kern)
        return data.ds.arr(img, "counts/pixel")
    return _counts

def setup_counts_fields(ds, ebounds, ftype="gas"):
    r"""
    Create deposited image fields from X-ray count data in energy bands.

    Parameters
    ----------
    ds : Dataset
        The FITS events file dataset to add the counts fields to.
    ebounds : list of tuples
        A list of tuples, one for each field, with (emin, emax) as the
        energy bounds for the image.
    ftype : string, optional
        The field type of the resulting field. Defaults to "gas".

    Examples
    --------
    >>> ds = yt.load("evt.fits")
    >>> ebounds = [(0.1,2.0),(2.0,3.0)]
    >>> setup_counts_fields(ds, ebounds)
    """
    for (emin, emax) in ebounds:
        cfunc = _make_counts(emin, emax)
        fname = "counts_%s-%s" % (emin, emax)
        mylog.info("Creating counts field %s." % fname)
        ds.add_field((ftype,fname), function=cfunc,
                     units="counts/pixel",
                     validators = [ValidateSpatial()],
                     display_name="Counts (%s-%s keV)" % (emin, emax))

def ds9_region(ds, reg, obj=None, field_parameters=None):
    r"""
    Create a data container from a ds9 region file. Requires the pyregion
    package (http://leejjoon.github.io/pyregion/) to be installed.

    Parameters
    ----------
    ds : FITSDataset
        The Dataset to create the region from.
    reg : string
        The filename of the ds9 region, or a region string to be parsed.
    obj : data container, optional
        The data container that will be used to create the new region.
        Defaults to ds.all_data.
    field_parameters : dictionary, optional
        A set of field parameters to apply to the region.

    Examples
    --------

    >>> ds = yt.load("m33_hi.fits")
    >>> circle_region = ds9_region(ds, "circle.reg")
    >>> print circle_region.quantities.extrema("flux")
    """
    import pyregion
    if os.path.exists(reg):
        r = pyregion.open(reg)
    else:
        r = pyregion.parse(reg)
    reg_name = reg
    filter = r.get_filter(header=ds.wcs_2d.to_header())
    nx = ds.domain_dimensions[ds.lon_axis]
    ny = ds.domain_dimensions[ds.lat_axis]
    mask = filter.mask((ny,nx)).transpose()
    def _reg_field(field, data):
        i = data["xyz"[ds.lon_axis]].ndarray_view().astype("int")-1
        j = data["xyz"[ds.lat_axis]].ndarray_view().astype("int")-1
        new_mask = mask[i,j]
        ret = data["zeros"].copy()
        ret[new_mask] = 1.
        return ret
    ds.add_field(("gas",reg_name), function=_reg_field)
    if obj is None:
        obj = ds.all_data()
    if field_parameters is not None:
        for k,v in field_parameters.items():
            obj.set_field_parameter(k,v)
    return obj.cut_region(["obj['%s'] > 0" % (reg_name)])

class PlotWindowWCS(object):
    r"""
    Use the wcsaxes library to plot celestial coordinates on the axes of a
    on-axis PlotWindow plot. See http://wcsaxes.readthedocs.org for details.

    Parameters
    ----------
    pw : on-axis PlotWindow instance
        The PlotWindow instance to add celestial axes to.
    """
    def __init__(self, pw):
        from wcsaxes import WCSAxes
        if pw.oblique:
            raise NotImplementedError("WCS axes are not implemented for oblique plots.")
        if not hasattr(pw.ds, "wcs_2d"):
            raise NotImplementedError("WCS axes are not implemented for this dataset.")
        if pw.data_source.axis != pw.ds.spec_axis:
            raise NotImplementedError("WCS axes are not implemented for this axis.")
        self.plots = {}
        self.pw = pw
        for f in pw.plots:
            rect = pw.plots[f]._get_best_layout()[1]
            fig = pw.plots[f].figure
            ax = fig.axes[0]
            wcs_ax = WCSAxes(fig, rect, wcs=pw.ds.wcs_2d, frameon=False)
            fig.add_axes(wcs_ax)
            wcs = pw.ds.wcs_2d.wcs
            xax = pw.ds.coordinates.x_axis[pw.data_source.axis]
            yax = pw.ds.coordinates.y_axis[pw.data_source.axis]
            xlabel = "%s (%s)" % (wcs.ctype[xax].split("-")[0],
                                  wcs.cunit[xax])
            ylabel = "%s (%s)" % (wcs.ctype[yax].split("-")[0],
                                  wcs.cunit[yax])
            fp = pw._font_properties
            wcs_ax.coords[0].set_axislabel(xlabel, fontproperties=fp)
            wcs_ax.coords[1].set_axislabel(ylabel, fontproperties=fp)
            wcs_ax.coords[0].ticklabels.set_fontproperties(fp)
            wcs_ax.coords[1].ticklabels.set_fontproperties(fp)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            wcs_ax.set_xlim(pw.xlim[0].value, pw.xlim[1].value)
            wcs_ax.set_ylim(pw.ylim[0].value, pw.ylim[1].value)
            wcs_ax.coords.frame._update_cache = []
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            self.plots[f] = fig

    def keys(self):
        return self.plots.keys()

    def values(self):
        return self.plots.values()

    def items(self):
        return self.plots.items()

    def __getitem__(self, key):
        for k in self.keys():
            if k[1] == key:
                return self.plots[k]

    def show(self):
        return self

    def save(self, name=None, mpl_kwargs=None):
        if mpl_kwargs is None:
            mpl_kwargs = {}
        mpl_kwargs["bbox_inches"] = "tight"
        self.pw.save(name=name, mpl_kwargs=mpl_kwargs)

    def _repr_html_(self):
        ret = ''
        for k, v in self.plots.iteritems():
            canvas = FigureCanvasAgg(v)
            f = StringIO()
            canvas.print_figure(f)
            f.seek(0)
            img = base64.b64encode(f.read())
            ret += '<img src="data:image/png;base64,%s"><br>' % img
        return ret
