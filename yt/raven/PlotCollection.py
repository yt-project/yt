"""
All of the base-level stuff for plotting.

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

class PlotCollection:
    __id_counter = 0
    def __init__(self, pf, deliverator_id=-1, center=None):
        """
        Generate a collection of linked plots using *pf* as a source,
        optionally submitting to the deliverator with *deliverator_id*
        and with *center*, which will otherwise be taken to be the point of
        maximum density.
        """
        PlotTypes.Initialize()
        self.plots = []
        self._run_id = deliverator_id
        self.pf = pf
        if center == None:
            v,self.c = pf.h.find_max("Density") # @todo: ensure no caching
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
        mylog.info("Created plot collection with default plot-center = %s",
                    list(self.c))

    def save(self, basename, format="png", override=False):
        """
        Same plots with automatically generated names, prefixed with *basename*
        (including directory path) unless *override* is specified, and in
        *format*.
        """
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
        """
        Set the x boundaries of all plots.
        """
        for plot in self.plots:
            plot.set_xlim(xmin, xmax)

    def set_ylim(self, ymin, ymax):
        """
        Set the y boundaries of all plots.
        """
        for plot in self.plots:
            plot.set_ylim(ymin, ymax)

    def set_zlim(self, zmin, zmax):
        """
        Set the limits of the colorbar.
        """
        for plot in self.plots:
            plot.set_autoscale(False)
            plot.set_zlim(zmin, zmax)

    def set_lim(self, lim):
        """
        Shorthand for setting x,y at same time.
        *lim* should be formatted as (xmin,xmax,ymin,ymax)
        """
        for plot in self.plots:
            plot.set_xlim(*lim[:2])
            plot.set_ylim(*lim[2:])

    def autoscale(self):
        """
        Turn back on autoscaling.
        """
        for plot in self.plots:
            plot.set_autoscale(True)

    def set_width(self, width, unit):
        """
        Set the witdh of the slices, cutting planes and projections to be
        *width* *units*
        """
        for plot in self.plots:
            plot.set_width(width, unit)

    def set_cmap(self, cmap):
        """
        Change the colormap of all plots to *cmap*.
        """
        for plot in self.plots:
            plot.set_cmap(cmap)

    def switch_field(self, field):
        """
        Change all the fields displayed to be *field*
        """
        for plot in self.plots:
            plot.switch_z(field)
    switch_z = switch_field

    def _add_plot(self, plot):
        """
        Accept some *plot* instance that's already been created.
        Handy for the subplot stuff that matplotlib is good at.
        And, as long as the 'plot' object is duck-typed, we should be fine
        with it, right?
        """
        self.plots.append(plot)
        return plot

    def add_slice(self, field, axis, coord=None, center=None,
                 use_colorbar=True, figure = None, axes = None, fig_size=None,
                 periodic = False):
        """
        Generate a slice through *field* along *axis*, optionally at
        [axis]=*coord*, with the *center* attribute given (some 
        degeneracy with *coord*, but not complete), with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)
        """
        if center == None:
            center = self.c
        if coord == None:
            coord = center[axis]
        slice = self.pf.hierarchy.slice(axis, coord, field, center)
        p = self._add_plot(PlotTypes.SlicePlot(slice, field, use_colorbar=use_colorbar,
                         axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        mylog.info("Added slice of %s at %s = %s with 'center' = %s", field,
                    axis_names[axis], coord, list(center))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_cutting_plane(self, field, normal,
                          center=None, use_colorbar=True,
                          figure = None, axes = None, fig_size=None, obj=None):
        """
        Generate a cutting plane of *field* with *normal*, centered at *center*
        (defaults to PlotCollection center) with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches).  If so desired,
        *obj* is a pre-existing cutting plane object.
        """
        if center == None:
            center = self.c
        if not obj:
            cp = self.pf.hierarchy.cutting(normal, center, field)
        else:
            cp = obj
        p = self._add_plot(PlotTypes.CuttingPlanePlot(cp, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        mylog.info("Added plane of %s with 'center' = %s and normal = %s", field,
                    list(center), list(normal))
        p["Axis"] = "CuttingPlane"
        return p

    def add_projection(self, field, axis, weight_field=None,
                      center=None, use_colorbar=True,
                      figure = None, axes = None, fig_size=None,
                      periodic = False):
        """
        Generate a projection of *field* along *axis*, optionally giving
        a *weight_field*-weighted average with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)
        """
        if center == None:
            center = self.c
        proj = self.pf.hierarchy.proj(axis, field, weight_field, center=center)
        p = self._add_plot(PlotTypes.ProjectionPlot(proj, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_profile_object(self, data_source, fields, 
                           weight="CellMassMsun", accumulation=False,
                           x_bins=64, x_log=True, x_bounds=None,
                           lazy_reader=False, id=None,
                           axes=None, figure=None):
        """
        Use an existing data object, *data_source*, to be the source of a
        one-dimensional profile.  *fields* will define the x and y bin-by
        fields, *weight* is used to weight the y value, *accumulation* determines
        if y is summed along x, *x_bins*, *x_log* and *x_bounds* define the means of
        choosing the bins.  *id* is used internally to differentiate between
        multiple plots in a single collection.  *lazy_reader* determines the
        memory-conservative status.
        """
        if x_bounds is None:
            x_min, x_max = data_source[fields[0]].min(), data_source[fields[0]].max()
        else:
            x_min, x_max = x_bounds
        profile = lagos.BinnedProfile1D(data_source,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     lazy_reader)
        if len(fields) > 1:
            profile.add_fields(fields[1], weight=weight, accumulation=accumulation)
        # These next two lines are painful.
        profile.pf = self.pf
        profile.hierarchy = self.pf.hierarchy
        if id is None: id = self._get_new_id()
        p = self._add_plot(PlotTypes.Profile1DPlot(profile, fields, id,
                                                   axes=axes, figure=figure))
        return p

    def add_profile_sphere(self, radius, unit, fields, **kwargs):
        """
        Generate a spherical 1D profile, given only a *radius*, a *unit*,
        and at least two *fields*.  Any remaining *kwargs* will be passed onto
        :meth:`add_profile_object`.
        """
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

    def add_phase_object(self, data_source, fields, cmap=None,
                               weight="CellMassMsun", accumulation=False,
                               x_bins=64, x_log=True, x_bounds=None,
                               y_bins=64, y_log=True, y_bounds=None,
                               lazy_reader=False, id=None,
                               axes = None, figure = None):
        """
        Given a *data_source*, and *fields*, automatically generate a 2D
        profile and plot it.  *id* is used internally to add onto the prefix,
        and will be automatically generated if not given. Remainder of
        arguments are identical to :meth:`add_profile_object`.
        """
        if x_bounds is None:
            x_min, x_max = data_source[fields[0]].min(), data_source[fields[0]].max()
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = data_source[fields[1]].min(), data_source[fields[1]].max()
        else:
            y_min, y_max = y_bounds
        profile = lagos.BinnedProfile2D(data_source,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     y_bins, fields[1], y_min, y_max, y_log,
                                     lazy_reader)
        # These next two lines are painful.
        profile.pf = self.pf
        profile.hierarchy = self.pf.hierarchy
        if id is None: id = self._get_new_id()
        p = self._add_plot(PlotTypes.PhasePlot(profile, fields, 
                                               id, cmap=cmap,
                                               figure=figure, axes=axes))
        if len(fields) > 2:
            # This will add it to the profile object
            p.switch_z(fields[2], weight=weight, accumulation=accumulation)
        return p

    def add_phase_sphere(self, radius, unit, fields, **kwargs):
        """
        Given a *radius* and *unit*, generate a 2D profile from a sphere, with
        *fields* as the x,y,z.  Automatically weights z by CellMassMsun.
        *kwargs* get passed onto :meth:`add_phase_object`.
        """
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
        """
        Delete all plots and their attendant data.
        """
        for i in range(len(self.plots)):
            del self.plots[-1].data
            del self.plots[-1]

def wrap_pylab_newplot(func):
    @wraps(func)
    def pylabify(self, *args, **kwargs):
        # Let's assume that axes and figure are not in the positional
        # arguments -- probably safe!
        new_fig = self.pylab.figure()
        kwargs['axes'] = self.pylab.gca()
        kwargs['figure'] = self.pylab.gcf()
        retval = func(self, *args, **kwargs)
        retval._redraw_image()
        retval._fig_num = new_fig.number
        self.pylab.show()
        return retval
    return pylabify

def wrap_pylab_show(func):
    @wraps(func)
    def pylabify(self, *args, **kwargs):
        retval = func(self, *args, **kwargs)
        self.pylab.show()
        return retval
    return pylabify

class PlotCollectionInteractive(PlotCollection):
    add_slice = wrap_pylab_newplot(PlotCollection.add_slice)
    add_cutting_plane = wrap_pylab_newplot(PlotCollection.add_cutting_plane)
    add_projection = wrap_pylab_newplot(PlotCollection.add_projection)
    add_profile_object = wrap_pylab_newplot(PlotCollection.add_profile_object)
    add_phase_object = wrap_pylab_newplot(PlotCollection.add_phase_object)
    
    set_xlim = wrap_pylab_show(PlotCollection.set_xlim)
    set_ylim = wrap_pylab_show(PlotCollection.set_ylim)
    set_zlim = wrap_pylab_show(PlotCollection.set_zlim)
    set_lim = wrap_pylab_show(PlotCollection.set_lim)
    autoscale = wrap_pylab_show(PlotCollection.autoscale)
    set_width = wrap_pylab_show(PlotCollection.set_width)
    set_cmap = wrap_pylab_show(PlotCollection.set_cmap)
    switch_field = wrap_pylab_show(PlotCollection.switch_field)

    def __init__(self, *args, **kwargs):
        import pylab
        self.pylab = pylab
        PlotCollection.__init__(self, *args, **kwargs)

    def redraw(self):
        for plot in self.plots:
            plot._redraw_image()
        self.pylab.show()

    def clear_plots(self):
        for plot in self.plots:
            self.pylab.figure(pylab._fig_num)
            self.pylab.clf()
        PlotCollection.clear_data(self)
