"""
This is an interface to MatPlotLib <http://matplotlib.sf.net> to plot
irregularly shaped grids, with the presumption that at any point we could have
data that is "hidden" in deeper levels of refinement.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

from yt.raven import *
from yt.funcs import *

import _MPL

# We only get imported if matplotlib was imported successfully

def ClusterFilePlot(cls, x, y, xlog=None, ylog=None, fig=None, filename=None,
                    format="png", xbounds = None, ybounds = None):
    """

    """
    if not fig:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = matplotlib.figure.Figure(figsize=(8,8))
        canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    if not iterable(cls):
        cls = [cls]
    if xlog == None:
        if lagos.CFfieldInfo.has_key(x):
            xlog = lagos.CFfieldInfo[x][2]
    if ylog == None:
        if lagos.CFfieldInfo.has_key(y):
            ylog = lagos.CFfieldInfo[y][2]
    if xlog and ylog:
        pp=ax.loglog
    elif xlog and not ylog:
        pp=ax.semilogx
    elif ylog and not xlog:
        pp=ax.semilogy
    else:
        pp=ax.plot

    fig.hold(True)
    colors = 'krbgm' * 10
    for cl, cc in zip(cls, colors):
        #pp(cl[x],cl[y], lw=2.5)
        pp(cl[x], cl[y], lw=2.5, color=cc)
    if lagos.CFfieldInfo.has_key(x):
        ax.set_xlabel(lagos.CFfieldInfo[x][1], fontsize=18)
        print lagos.CFfieldInfo[x][1]
    if lagos.CFfieldInfo.has_key(y):
        ax.set_ylabel(lagos.CFfieldInfo[y][1], fontsize=18)
        print lagos.CFfieldInfo[y][1]
    if xbounds:
        ax.set_xlim(xbounds)
    if ybounds:
        ax.set_ylim(ybounds)
    ax.axesFrame.set_linewidth(2)
    for tickLabel in ax.get_xticklabels() + ax.get_yticklabels():
        tickLabel.set_fontsize(14)
    if filename:
        canvas.print_figure(filename, format=format)
    return fig

engineVals = {}

def Initialize(*args, **kwargs):
    engineVals["initialized"] = True
    if not kwargs.has_key("canvas"):
        from matplotlib.backends.backend_agg \
                import FigureCanvasAgg as FigureCanvas
    else:
        FigureCanvas = kwargs["canvas"]
    engineVals["canvas"] = FigureCanvas
    return

def CleanUp(*args, **kwargs):
    pass

class RavenPlot:

    datalabel = None
    colorbar = None
    def __init__(self, data, fields, figure = None, axes=None, size=(10,8)):
        self.data = data
        self.fields = fields
        self.size = size
        self.set_autoscale(True)
        self.im = defaultdict(lambda: "")
        self["ParameterFile"] = "%s" % self.data.pf
        self.axis_names = {}
        self._ax_max = self.data.pf["DomainRightEdge"]
        if not figure:
            self._figure = matplotlib.figure.Figure(size)
        else:
            self._figure = figure
        if not figure:
            self._axes = self._figure.add_subplot(1,1,1)
        else:
            self._axes = axes
        self._callbacks = []

    def set_autoscale(self, val):
        self.do_autoscale = val

    def __getitem__(self, item):
        return self.data[item] # Should be returned in CGS

    def save_image(self, prefix, format, submit=None, override=False):
        """
        Save this plot image.  Will generate a filename based on the *prefix*,
        *format*.  *submit* will govern the submission to the Deliverator and
        *override* will force no filename generation beyond the prefix.
        """
        self._redraw_image()
        if not override:
            self._generate_prefix(prefix)
            my_prefix = self.prefix
        else:
            my_prefix = prefix
        fn = ".".join([my_prefix, format])
        canvas = engineVals["canvas"](self._figure)
        #self._figure.savefig(fn, format)
        canvas.print_figure(fn)
        self["Type"] = self._type_name
        self["GeneratedAt"] = self.data.hierarchy["CurrentTimeIdentifier"]
        return fn

    def _redraw_image(self):
        pass

    def _generate_prefix(self):
        pass

    def set_xlim(self, xmin, xmax):
        """
        Set the x boundaries of this plot.
        """
        self._axes.set_xlim(xmin, xmax)

    def set_ylim(self, ymin, ymax):
        """
        Set the y boundaries of this plot.
        """
        self._axes.set_ylim(ymin, ymax)

    def set_zlim(self, zmin, zmax):
        """
        Set the z boundaries of this plot.
        """
        self.norm.autoscale(na.array([zmin,zmax]))
        self.image.changed()
        if self.colorbar is not None:
            _notify(self.image, self.colorbar)

    def set_cmap(self, cmap):
        """
        Change the colormap of this plot to *cmap*.
        """
        if isinstance(cmap, types.StringTypes):
            if str(cmap) in raven_colormaps:
                cmap = raven_colormaps[str(cmap)]
            elif hasattr(matplotlib.cm, cmap):
                cmap = getattr(matplotlib.cm, cmap)
        self.cmap = cmap

    def __setitem__(self, item, val):
        self.im[item] = val

    def add_callback(self, func):
        """
        Add *func* as a callback to this plot.  *func* will be called with this
        plot as its first argument every time the plot is redrawn.  Returns the
        id of the callback (for use with :meth:`remove_callback`.)
        """
        self._callbacks.append(func)
        return len(self._callbacks)-1

    def remove_callback(self, id):
        """
        Given an *id*, remove that index in the callbacks list.
        """
        self._callbacks[id] = lambda a: None

    def _run_callbacks(self):
        for cb in self._callbacks:
            cb(self)

    def set_label(self, label):
        """
        Set the datalabel to *label*.  (This has different meanings based on
        the plot.)
        """
        self.datalabel = label
        if self.colorbar != None: self.colorbar.set_label(str(label))

class VMPlot(RavenPlot):
    _antialias = True
    def __init__(self, data, field, figure = None, axes = None,
                 use_colorbar = True, size=None):
        fields = ['X', 'Y', field, 'X width', 'Y width']
        if not size:
            size = (10,8)
            if not use_colorbar: size=(8,8)
        RavenPlot.__init__(self, data, fields, figure, axes, size=size)
        self._figure.subplots_adjust(hspace=0, wspace=0, bottom=0.0,
                                    top=1.0, left=0.0, right=1.0)
        self.xmin = 0.0
        self.ymin = 0.0
        self.xmax = 1.0
        self.ymax = 1.0
        self.cmap = None
        if self.data.axis < 3:
            self._x_max = self._ax_max[lagos.x_dict[self.data.axis]]
            self._y_max = self._ax_max[lagos.y_dict[self.data.axis]]
        self.__setup_from_field(field)
        self.__init_temp_image(use_colorbar)

    def __setup_from_field(self, field):
        self.set_log_field(field in lagos.log_fields
                           or lagos.fieldInfo[field].take_log)
        self.axis_names["Z"] = field

    def set_log_field(self, val):
        if val:
            self.log_field = True
            self.norm = matplotlib.colors.LogNorm()
            ttype = matplotlib.ticker.LogFormatter
        else:
            self.log_field = False
            self.norm = matplotlib.colors.Normalize()
            ttype = matplotlib.ticker.ScalarFormatter
        if self.colorbar:
            self.colorbar.set_norm(self.norm)
            self.colorbar.formatter = ttype()

    def __init_temp_image(self, setup_colorbar):
        temparray = na.ones(self.size)
        self.image = \
            self._axes.imshow(temparray, interpolation='nearest',
                             norm = self.norm, aspect=1.0, picker=True,
                             origin='lower')
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if setup_colorbar:
            self.colorbar = self._figure.colorbar(self._axes.images[-1], \
                                                  extend='neither', \
                                                  shrink=0.95)
        else:
            self.colorbar = None
        self.set_width(1,'1')

    def _get_buff(self, width=None):
        x0, x1 = self.xlim
        y0, y1 = self.ylim
        if width is None:
            l, b, width, height = _get_bounds(self._axes.bbox)
        else:
            height = width
        self.pix = (width,height)
        # 'px' == pixel x, or x in the plane of the slice
        # 'x' == actual x
        aa = int(self._antialias)
        buff = _MPL.Pixelize(self.data['px'],
                            self.data['py'],
                            self.data['pdx'],
                            self.data['pdy'],
                            self[self.axis_names["Z"]],
                            int(width), int(width),
                            (x0, x1, y0, y1),aa).transpose()
        return buff

    def _redraw_image(self, *args):
        self._axes.clear() # To help out the colorbar
        buff = self._get_buff()
        mylog.debug("Received buffer of min %s and max %s (%s %s)",
                    buff.min(), buff.max(),
                    self[self.axis_names["Z"]].min(),
                    self[self.axis_names["Z"]].max())
        if self.log_field:
            bI = na.where(buff > 0)
            newmin = buff[bI].min()
            newmax = buff[bI].max()
        else:
            newmin = buff.min()
            newmax = buff.max()
        if self.do_autoscale:
            self.norm.autoscale(na.array((newmin,newmax)))
        self.image = \
            self._axes.imshow(buff, interpolation='nearest', norm = self.norm,
                            aspect=1.0, picker=True, origin='lower')
        self._reset_image_parameters()
        self._run_callbacks()

    def _reset_image_parameters(self):
        self._axes.set_xticks(())
        self._axes.set_yticks(())
        self._axes.set_ylabel("")
        self._axes.set_xlabel("")
        if self.cmap:
            self.image.set_cmap(self.cmap)
        if self.colorbar != None:
            self.image.set_norm(self.norm)
            self.colorbar.set_norm(self.norm)
            if self.do_autoscale: _notify(self.image, self.colorbar)
        self.autoset_label()

    def set_xlim(self, xmin, xmax):
        self.xlim = (xmin,xmax)

    def set_ylim(self, ymin, ymax):
        self.ylim = (ymin,ymax)

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name, \
            lagos.axis_names[self.data.axis], self.axis_names['Z']])
        self["Field1"] = self.axis_names["Z"]
        self["Field2"] = None
        self["Field3"] = None

    def set_width(self, width, unit):
        self["Unit"] = str(unit)
        self["Width"] = float(width)
        if isinstance(unit, types.StringType):
            unit = self.data.hierarchy[unit]
        self.width = width / unit
        self._refresh_display_width()

    def _refresh_display_width(self, width=None):
        if width:
            self.width = width
        else:
            width = self.width
        if iterable(width):
            width_x, width_y = width
        else:
            width_x = width
            width_y = width
        l_edge_x = self.data.center[lagos.x_dict[self.data.axis]] - width_x/2.0
        r_edge_x = self.data.center[lagos.x_dict[self.data.axis]] + width_x/2.0
        l_edge_y = self.data.center[lagos.y_dict[self.data.axis]] - width_y/2.0
        r_edge_y = self.data.center[lagos.y_dict[self.data.axis]] + width_y/2.0
        self.set_xlim(max(l_edge_x,0.0), min(r_edge_x,self._x_max))
        self.set_ylim(max(l_edge_y,0.0), min(r_edge_y,self._y_max))
        self._redraw_image()

    def autoscale(self):
        zmin = self._axes.images[-1]._A.min()
        zmax = self._axes.images[-1]._A.max()
        self.set_zlim(zmin, zmax)

    def switch_y(self, *args, **kwargs):
        pass

    def switch_x(self, *args, **kwargs):
        pass

    def switch_z(self, field):
        self.set_log_field(field in lagos.log_fields
                           or lagos.fieldInfo[field].take_log)
        self.axis_names["Z"] = field
        self._redraw_image()

    def selfSetup(self):
        pass

class SlicePlot(VMPlot):
    _type_name = "Slice"

    def autoset_label(self):
        if self.datalabel != None:
            self.colorbar.set_label(str(self.datalabel))
            return
        field_name = self.axis_names["Z"]
        data_label = r"$\rm{%s}" % field_name.replace("_","\hspace{0.5}")
        if lagos.fieldInfo.has_key(field_name):
            data_label += r"\/\/ (%s)" % (lagos.fieldInfo[field_name].get_units())
        data_label += r"$"
        if self.colorbar != None: self.colorbar.set_label(str(data_label))

class ProjectionPlot(VMPlot):

    _type_name = "Projection"
    def autoset_label(self):
        if self.datalabel != None:
            self.colorbar.set_label(str(self.datalabel))
            return
        field_name = self.axis_names["Z"]
        data_label = r"$\rm{%s}" % field_name.replace("_","\hspace{0.5}")
        if lagos.fieldInfo.has_key(field_name):
            data_label += r"\/\/ (%s)" % (lagos.fieldInfo[field_name].get_projected_units())
        data_label += r"$"
        if self.colorbar != None: self.colorbar.set_label(str(data_label))

    def switch_z(self, field):
        mylog.warning("Choosing not to change the field of a projection instance")


class CuttingPlanePlot(SlicePlot):

    _type_name = "CuttingPlane"
    def _get_buff(self, width=None):
        px_min, px_max = self.xlim
        py_min, py_max = self.ylim
        if width is None:
            l, b, width, height = _get_bounds(self._axes.bbox)
        else:
            height = width
        self.pix = (width,height)
        indices = na.argsort(self.data['dx'])[::-1]
        buff = _MPL.CPixelize( self.data['x'], self.data['y'], self.data['z'],
                               self.data['px'], self.data['py'],
                               self.data['pdx'], self.data['pdy'], self.data['pdz'],
                               self.data.center, self.data._inv_mat, indices,
                               self.data[self.axis_names['Z']],
                               int(width), int(width),
                               (px_min, px_max, py_min, py_max)).transpose()
        return buff

    def _refresh_display_width(self, width=None):
        
        if width:
            self.width = width
        else:
            width = self.width
        if iterable(width):
            width_x, width_y = width
        else:
            width_x = width
            width_y = width
        l_edge_x = -width_x/2.0
        r_edge_x = +width_x/2.0
        l_edge_y = -width_y/2.0
        r_edge_y = +width_y/2.0
        self.set_xlim(l_edge_x, r_edge_x) # We have no real limits
        self.set_ylim(l_edge_y, r_edge_y) # At some point, perhaps calculate them?
        self._redraw_image()

class ProfilePlot(RavenPlot):
    def setup_bins(self, field, func=None):
        if field in lagos.fieldInfo and lagos.fieldInfo[field].take_log:
            log_field = True
            if func: func('log')
        else:
            log_field = False
            if func: func('linear')
        mylog.debug("Field: %s, log_field: %s", field, log_field)
        return log_field

    def autoset_label(self, field, func):
        dataLabel = r"$\rm{%s}" % (field.replace("_","\hspace{0.5}"))
        if field in lagos.fieldInfo:
            dataLabel += r" (%s)" % (lagos.fieldInfo[field].get_units())
        dataLabel += r"$"
        func(str(dataLabel))


class Profile1DPlot(ProfilePlot):

    _type_name = "Profile1D"
    def __init__(self, data, fields, id, ticker=None, 
                 figure=None, axes=None, plot_options=None):
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]

        if plot_options is None: plot_options = {}
        self.plot_options = plot_options

        self._log_x = self.data._x_log
        self._log_y = self.setup_bins(self.fields[1])

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name,
                       str(self._semi_unique_id),
                       self.axis_names['X'], self.axis_names['Y']])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]

    def _redraw_image(self):
        self._axes.clear()
        if not self._log_x and not self._log_y:
            func = self._axes.plot
        elif self._log_x and not self._log_y:
            func = self._axes.semilogx
        elif not self._log_x and self._log_y:
            func = self._axes.semilogy
        elif self._log_x and self._log_y:
            func = self._axes.loglog
        indices = na.argsort(self.data[self.fields[0]])
        func(self.data[self.fields[0]][indices],
             self.data[self.fields[1]][indices],
             **self.plot_options)
        self.autoset_label(self.fields[0], self._axes.set_xlabel)
        self.autoset_label(self.fields[1], self._axes.set_ylabel)
        self._run_callbacks()

    def set_log_field(self, val):
        if val:
            self._log_y = True
        else:
            self._log_y = False

    def switch_x(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[0] = field
        self.axis_names["X"] = field
        if field not in self.data.keys():
            self.data.add_fields(field, weight, accumulation)
        self._log_x = self.setup_bins(self.fields[0])
    
    def switch_z(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[1] = field
        self.axis_names["Y"] = field
        if field not in self.data.keys():
            self.data.add_fields(field, weight, accumulation)
        self._log_y = self.setup_bins(self.fields[1])
    switch_y = switch_z # Compatibility...

class PhasePlot(ProfilePlot):

    _type_name = "Profile2D"
    def __init__(self, data, fields, id, ticker=None, cmap=None,
                 figure=None, axes=None):
        self._semi_unique_id = id
        RavenPlot.__init__(self, data, fields, figure, axes)
        self.ticker = ticker
        self.image = None
        self.set_cmap(cmap)

        self.axis_names["X"] = fields[0]
        self.axis_names["Y"] = fields[1]
        self.axis_names["Z"] = fields[2]

        self.x_bins = self.data[self.fields[0]]
        self.y_bins = self.data[self.fields[1]]
        self._log_x = self.data._x_log
        self._log_y = self.data._y_log
        self._log_z = self.setup_bins(self.fields[2])
        self.__init_colorbar()

    def __init_colorbar(self):
        temparray = na.ones((self.x_bins.size, self.y_bins.size))
        self.norm = matplotlib.colors.Normalize()
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins,
                                      temparray, shading='flat',
                                      norm=self.norm)
        self.colorbar = self._figure.colorbar(self.image,
                                    extend='neither', shrink=0.95,
                                    format="%0.2e" )

    def set_cmap(self, cmap):
        RavenPlot.set_cmap(self, cmap)
        if self.image != None and self.cmap != None:
            self.image.set_cmap(self.cmap)

    def switch_z(self, field, weight="CellMassMsun", accumulation=False):
        self.fields[2] = field
        self.axis_names["Z"] = field
        if field not in self.data.keys(): self.data.add_fields(field, weight, accumulation)
        self._log_z = self.setup_bins(self.fields[2])

    def set_log_field(self, val):
        if val:
            self._log_z = True
            self.norm = matplotlib.colors.LogNorm()
            ttype = matplotlib.ticker.LogFormatter
        else:
            self._log_z = False
            self.norm = matplotlib.colors.Normalize()
            ttype = matplotlib.ticker.ScalarFormatter
        if self.colorbar:
            self.colorbar.set_norm(self.norm)
            self.colorbar.formatter = ttype()

    def _redraw_image(self):
        vals = self.data[self.fields[2]].transpose()
        used_bin = self.data["UsedBins"].transpose()
        vmin = na.nanmin(vals[used_bin])
        vmax = na.nanmax(vals[used_bin])
        if self._log_z:
            # We want smallest non-zero vmin
            self.norm=matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax,
                                                clip=False)
            self.ticker = matplotlib.ticker.LogLocator()
            vI = na.where(vals > 0)
            newmin = vals[vI].min()
            newmax = vals[vI].max()
            self.norm.autoscale(na.array((newmin,newmax)))
        else:
            self.norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax,
                                                  clip=False)
            self.ticker = matplotlib.ticker.MaxNLocator()
        self.colorbar.set_norm(self.norm)
        if self.cmap == None:
            self.cmap = matplotlib.cm.get_cmap()
        self.cmap.set_bad("w")
        self.cmap.set_under("w")
        self.cmap.set_over("w")
        self._axes.clear()
        self._axes.set_xscale("linear")
        self._axes.set_yscale("linear")
        self.image = self._axes.pcolormesh(self.x_bins, self.y_bins, \
                                      vals, shading='flat', \
                                      norm=self.norm, cmap=self.cmap)
        self._axes.set_xscale({0:"linear",1:"log"}[int(self._log_x)])
        self._axes.set_yscale({0:"linear",1:"log"}[int(self._log_y)])
        self.vals = vals

        _notify(self.image, self.colorbar)
        self.autoset_label(self.fields[0], self._axes.set_xlabel)
        self.autoset_label(self.fields[1], self._axes.set_ylabel)
        self.autoset_label(self.fields[2], self.colorbar.set_label)
        self._run_callbacks()

    def _generate_prefix(self, prefix):
        self.prefix = "_".join([prefix, self._type_name,
            str(self._semi_unique_id),
            self.axis_names['X'], self.axis_names['Y'],
            self.axis_names['Z'], ])
        self["Field1"] = self.axis_names["X"]
        self["Field2"] = self.axis_names["Y"]
        self["Field3"] = self.axis_names["Z"]

    def autoset_label(self, field_name, func):
        if self.datalabel != None:
            func(str(self.datalabel))
            return
        data_label = r"$\rm{%s}" % field_name.replace("_"," ")
        if field_name in lagos.fieldInfo:
            data_label += r"\/\/ (%s)" % (lagos.fieldInfo[field_name].get_units())
        data_label += r"$"
        func(str(data_label))

    def set_width(self, width, unit):
        mylog.warning("Choosing not to change the width of a phase plot instance")


# Now we provide some convenience functions to get information about plots.
# With Matplotlib 0.98.x, the 'transforms' branch broke backwards
# compatibility.  Despite that, the various packagers are plowing ahead with
# packaging 0.98.x with new distributions of python software.  So I guess
# we have to support it.

_compatibility_functions = ["_get_bounds","_notify"]

_mpl98_get_bounds = lambda bbox: bbox.bounds
_mpl9x_get_bounds = lambda bbox: bbox.get_bounds()
_mpl98_notify = lambda im,cb: cb.update_bruteforce(im)
_mpl9x_notify = lambda im,cb: cb.notify(im)

# This next function hurts, because it relies on the fact that
# we're only differentiating between 0.9[01] and 0.98.

_mpl_version = float(matplotlib.__version__[:4])

if _mpl_version < 0.98:
    _prefix = '_mpl9x'
    mylog.debug("Turning on matplotlib 0.9X compat (%s)",
                matplotlib.__version__)
else:
    _prefix = '_mpl98'
    mylog.debug("Turning on matplotlib 0.98 compat (%s)",
                matplotlib.__version__)

for fn in _compatibility_functions:
    exec("%s = %s%s" % (fn, _prefix, fn))
