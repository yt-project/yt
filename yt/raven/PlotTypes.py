"""
All of the base-level stuff for plotting.

Think of this as a way of getting rid of EnzoHippo.  We should have access to
all of the engine-appropriate methods.

@todo: Implement regions

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

class PlotCollection:
    def __init__(self, pf, deliverator_id=-1, center=None):
        be.Initialize()
        self.plots = []
        self._run_id = deliverator_id
        self.pf = pf
        if not center:
            v,self.c = pf.h.findMax("Density") # @todo: ensure no caching
        else:
            self.c = center
        if deliverator_id > 0:
            self.submit = True
            self._run_id = deliverator_id
            r=deliveration.SubmitParameterFile(\
                deliverator_id, self.pf)
            mylog.debug("Received response '%s'", r)
            self._http_prefix = ytcfg["raven","httpPrefix"] % self.pf
        else:
            self.submit = False
    def save(self, basename, format="png", override=False):
        fn = []
        for plot in self.plots:
            fn.append(plot.save_image(basename, \
                      format=format, submit=self._run_id,
                      override=override))
            if self.submit:
                im = plot.im.copy()
                im["Filename"] = self._http_prefix + "/" \
                                + os.path.basename(fn[-1])
                im["RunID"] = self._run_id
                deliveration.SubmitImage(\
                      self.pf.hierarchy, im)
            mylog.info("Saved %s", fn[-1])
        return fn
    def set_xlim(self, xmin, xmax):
        for plot in self.plots:
            plot.set_xlim(xmin, xmax)
    def set_ylim(self, ymin, ymax):
        for plot in self.plots:
            plot.set_ylim(ymin, ymax)
    def autoscale(self):
        for plot in self.plots:
            plot.set_autoscale(True)
    def set_zlim(self, zmin, zmax):
        for plot in self.plots:
            plot.set_autoscale(False)
            plot.set_zlim(zmin, zmax)
    def set_lim(self, lim):
        for plot in self.plots:
            plot.set_xlim(lim[:2])
            plot.set_ylim(lim[2:])
    def set_width(self, width, unit):
        for plot in self.plots:
            plot.set_width(width, unit)
    def set_cmap(self, cmap):
        for plot in self.plots:
            plot.set_cmap(cmap)
    def switch_field(self, field):
        for plot in self.plots:
            plot.switch_z(field)
    switch_z = switch_field
    # Now we get around to adding the plots we want
    def _add_plot(self, plot):
        # Accept some instance that's already been created
        # Handy for the subplot stuff that matplotlib is good at
        # And, as long as the 'plot' object is duck-typed, we should be fine
        # with it, right?
        self.plots.append(plot)
        return plot
    def add_slice(self, field, axis, coord=None, center=None,
                 use_colorbar=True,
                 figure = None, axes = None):
        # Make the slice object here
        # Pass it in to the engine to get a SlicePlot, which we then append
        # onto self.plots
        if center == None:
            center = self.c
        if coord == None:
            coord = center[axis]
        slice = self.pf.hierarchy.slice(axis, coord, field, center)
        p = self._add_plot(be.SlicePlot(slice, field, use_colorbar=use_colorbar,
                                      axes=axes, figure=figure))
        p["Axis"] = lagos.axis_names[axis]
        return p
    def add_projection(self, field, axis, weight_field=None,
                      center=None, use_colorbar=True,
                      figure = None, axes = None, fig_size=None):
        # Make the projection object here
        # Pass it in to the engine to get a SlicePlot, which we then append
        # onto self.plots
        if center == None:
            center = self.c
        proj = self.pf.hierarchy.proj(axis, field, weight_field, center=center)
        p = self._add_plot(be.ProjectionPlot(proj, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        p["Axis"] = lagos.axis_names[axis]
        return p
    def add_twophase_sphere(self, radius, unit, fields, center=None, cmap=None):
        if center == None:
            center = self.c
        fields = fields[:2] + ["CellsPerBin"] + fields[2:]
        r = radius/self.pf[unit]
        sphere = self.pf.hierarchy.sphere(center, r, fields)
        p = self._add_plot(be.PhasePlot(sphere, fields, width=radius,
                                      unit=unit, cmap=cmap))
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p
    def add_threephase_sphere(self, radius, unit, fields, center=None, cmap=None,
                            weight="CellMass"):
        if center == None:
            center = self.c
        r = radius/self.pf[unit]
        sphere = self.pf.hierarchy.sphere(center, r, fields)
        p = self._add_plot(be.PhasePlot(sphere, fields, width=radius,
                                      unit=unit, cmap=cmap, weight=weight))
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p
