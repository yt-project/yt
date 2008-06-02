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
    __id_counter = 0
    def __init__(self, pf, deliverator_id=-1, center=None):
        PlotTypes.Initialize()
        self.plots = []
        self._run_id = deliverator_id
        self.pf = pf
        if center == None:
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
        p = self._add_plot(PlotTypes.SlicePlot(slice, field, use_colorbar=use_colorbar,
                                      axes=axes, figure=figure))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_cutting_plane(self, field, normal,
                          center=None, use_colorbar=True,
                          figure = None, axes = None, fig_size=None, obj=None):
        if center == None:
            center = self.c
        if not obj:
            cp = self.pf.hierarchy.cutting(normal, center, field)
        else:
            cp = obj
        p = self._add_plot(PlotTypes.CuttingPlanePlot(cp, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        p["Axis"] = "CuttingPlane"
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
        p = self._add_plot(PlotTypes.ProjectionPlot(proj, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_profile_object(self, object, fields, cmap=None,
                           weight="CellMassMsun", accumulation=False,
                           x_bins=64, x_log=True, x_bounds=None,
                           lazy_reader=False, id=None):
        if x_bounds is None:
            x_min, x_max = object[fields[0]].min(), object[fields[0]].max()
        else:
            x_min, x_max = x_bounds
        profile = lagos.BinnedProfile1D(object,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     lazy_reader)
        if len(fields) > 1:
            profile.add_fields(fields[1], weight=weight, accumulation=accumulation)
        # These next two lines are painful.
        profile.pf = self.pf
        profile.hierarchy = self.pf.hierarchy
        if id is None: id = self._get_new_id()
        p = self._add_plot(PlotTypes.Profile1DPlot(profile, fields, 
                                                   id, cmap=cmap))
        return p

    def add_profile_sphere(self, radius, unit, fields, **kwargs):
        center = kwargs.pop("center",self.c)
        r = radius/self.pf[unit]
        if 'sphere' in kwargs:
            sphere = kwargs.pop('sphere')
        else:
            ftg = fields[:]
            if kwargs.get("lazy_reader",False): ftg = []
            sphere = self.pf.hierarchy.sphere(center, r, ftg)
        p = self.add_profile_object(sphere, fields, **kwargs)
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p

    def add_phase_object(self, object, fields, cmap=None,
                               weight="CellMassMsun", accumulation=False,
                               x_bins=64, x_log=True, x_bounds=None,
                               y_bins=64, y_log=True, y_bounds=None,
                               lazy_reader=False, id=None):
        if x_bounds is None:
            x_min, x_max = object[fields[0]].min(), object[fields[0]].max()
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = object[fields[1]].min(), object[fields[1]].max()
        else:
            y_min, y_max = y_bounds
        profile = lagos.BinnedProfile2D(object,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     y_bins, fields[1], y_min, y_max, y_log,
                                     lazy_reader)
        if len(fields) > 2:
            profile.add_fields(fields[2], weight=weight, accumulation=accumulation)
        # These next two lines are painful.
        profile.pf = self.pf
        profile.hierarchy = self.pf.hierarchy
        if id is None: id = self._get_new_id()
        p = self._add_plot(PlotTypes.PhasePlot(profile, fields, 
                                               id, cmap=cmap))
        return p

    def add_phase_sphere(self, radius, unit, fields, **kwargs):
        center = kwargs.pop("center",self.c)
        r = radius/self.pf[unit]
        if 'sphere' in kwargs:
            sphere = kwargs.pop('sphere')
        else:
            ftg = fields[:]
            if kwargs.get("lazy_reader",False): ftg = []
            sphere = self.pf.hierarchy.sphere(center, r, ftg)
        p = self.add_phase_object(sphere, fields, **kwargs)
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p

    def _get_new_id(self):
        self.__id_counter += 1
        return self.__id_counter-1


    def clear_plots(self):
        for i in range(len(self.plots)):
            del self.plots[i].data
            del self.plots[i]
