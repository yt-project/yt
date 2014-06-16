"""
Easy plotting.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from ._mpl_imports import *
from yt.data_objects.profiles import BinnedProfile1D

plot_type_registry = {}

class EasyPlot(object):
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            if hasattr(cls, "_type_name"):
                plot_type_registry[cls._type_name] = cls

    def __init__(self):
        pass

    def save(self, fn):
        canvas = FigureCanvasAgg(self.figure)
        canvas.print_figure(fn)

class EasyPDFPlot(EasyPlot):
    _type_name = "PDF"

    def __init__(self, data_source, x_field, x_log = True,
                 n_bins = 128, y_log = True,
                 plot_args = None,
                 figure_args = None):
        if plot_args is None: plot_args = {}
        if figure_args is None: figure_args = {}
        self.x_field = x_field
        self.data_source = data_source
        # Now we just make the plot
        x_min, x_max = self.data_source.quantities["Extrema"](
                x_field, non_zero = x_log)[0]
        self.profile = BinnedProfile1D(self.data_source,
            n_bins, self.x_field, x_min, x_max, x_log)
        self.profile.add_fields(["CellMassMsun"], weight=None)
        self.profile["CellMassMsun"] /= self.profile["CellMassMsun"].sum()
        self.figure = matplotlib.figure.Figure(**figure_args)
        self.axes = self.figure.add_subplot(1,1,1)
        if y_log and x_log: f = self.axes.loglog
        elif y_log: f = self.axes.semilogy
        elif x_log: f = self.axes.semilogx
        else: f = self.axes.plot
        self.plot = f(self.profile[x_field], self.profile["CellMassMsun"],
                      **plot_args)
        self.axes.set_xlabel(data_source.ds.field_info[x_field].get_label())
