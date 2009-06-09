"""
All of the base-level stuff for plotting.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

# No better place to put this
def concatenate_pdfs(output_fn, input_fns):
    from pyPdf import PdfFileWriter, PdfFileReader
    outfile = PdfFileWriter()
    for fn in input_fns:
        infile = PdfFileReader(open(fn, 'rb'))
        outfile.addPage(infile.getPage(0))
    outfile.write(open(output_fn, "wb"))

class PlotCollection(object):
    __id_counter = 0
    def __init__(self, pf, center=None, deliverator_id=-1):
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
            self.c = na.array(center, dtype='float64')
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

    def __iter__(self):
        for p in self.plots:
            yield p

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

    def add_slice(self, *args, **kwargs):
        """
        Generate a slice through *field* along *axis*, optionally at
        [axis]=*coord*, with the *center* attribute given (some 
        degeneracy with *coord*, but not complete), with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)
        """
        return self.__add_slice(PlotTypes.SlicePlot, *args, **kwargs)

    def add_slice_interpolated(self, *args, **kwargs):
        """
        Generate a slice through *field* along *axis*, optionally at
        [axis]=*coord*, with the *center* attribute given (some 
        degeneracy with *coord*, but not complete), with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)

        The slice will be interpolated using the delaunay module, with natural
        neighbor interpolation.
        """
        return self.__add_slice(PlotTypes.SlicePlotNaturalNeighbor, *args, **kwargs)

    def __add_slice(self, ptype, field, axis, coord=None, center=None,
                 use_colorbar=True, figure = None, axes = None, fig_size=None,
                 periodic = True, data_source = None, **kwargs):
        if center == None:
            center = self.c
        if coord == None:
            coord = center[axis]
        if data_source is None:
            data_source = self.pf.hierarchy.slice(axis, coord, field, center, **kwargs)
        p = self._add_plot(ptype(data_source, field, use_colorbar=use_colorbar,
                         axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        mylog.info("Added slice of %s at %s = %s with 'center' = %s", field,
                    axis_names[axis], coord, list(center))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_particles(self, axis, width, p_size=1.0, col='k', stride=1.0,
                      data_source=None):
        LE = self.pf["DomainLeftEdge"].copy()
        RE = self.pf["DomainRightEdge"].copy()
        LE[axis] = self.c[axis] - width/2.0
        RE[axis] = self.c[axis] + width/2.0
        if data_source is None: data_source = self.pf.h.region(self.c, LE, RE)
        p = self._add_plot(PlotTypes.ParticlePlot(data_source, axis,
                                        width, p_size, col, stride))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_cutting_plane(self, field, normal,
                          center=None, use_colorbar=True,
                          figure = None, axes = None, fig_size=None, obj=None,
                           **kwargs):
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
            cp = self.pf.hierarchy.cutting(normal, center, field, **kwargs)
        else:
            cp = obj
        p = self._add_plot(PlotTypes.CuttingPlanePlot(cp, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        mylog.info("Added plane of %s with 'center' = %s and normal = %s", field,
                    list(center), list(normal))
        p["Axis"] = "CuttingPlane"
        return p

    def add_projection(self, *args, **kwargs):
        """
        Generate a projection of *field* along *axis*, optionally giving
        a *weight_field*-weighted average with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)
        """
        return self._add_projection(PlotTypes.ProjectionPlot, *args, **kwargs)

    def add_projection_interpolated(self, *args, **kwargs):
        """
        Generate a projection of *field* along *axis*, optionally giving
        a *weight_field*-weighted average with *use_colorbar*
        specifying whether the plot is naked or not and optionally
        providing pre-existing Matplotlib *figure* and *axes* objects.
        *fig_size* in (height_inches, width_inches)

        The projection will be interpolated using the delaunay module, with
        natural neighbor interpolation.
        """
        return self._add_projection(PlotTypes.ProjectionPlotNaturalNeighbor, *args, **kwargs)

    def _add_projection(self, ptype, field, axis, weight_field=None,
                      center=None, use_colorbar=True,
                      figure = None, axes = None, fig_size=None,
                      periodic = True, data_source = None, **kwargs):
        if center == None:
            center = self.c
        if data_source is None:
            data_source = self.pf.hierarchy.proj(axis, field, weight_field,
                                center=center, **kwargs)
        p = self._add_plot(ptype(data_source, field,
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
            x_min, x_max = data_source.quantities["Extrema"](
                            fields[0], lazy_reader=lazy_reader)[0]
        else:
            x_min, x_max = x_bounds
        profile = lagos.BinnedProfile1D(data_source,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     lazy_reader)
        if len(fields) > 1:
            profile.add_fields(fields[1], weight=weight, accumulation=accumulation)
        # These next two lines are painful.
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
            x_min, x_max = data_source.quantities["Extrema"](
                                    fields[0], lazy_reader=lazy_reader)[0]
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = data_source.quantities["Extrema"](
                                    fields[1], lazy_reader=lazy_reader)[0]
        else:
            y_min, y_max = y_bounds
        profile = lagos.BinnedProfile2D(data_source,
                                     x_bins, fields[0], x_min, x_max, x_log,
                                     y_bins, fields[1], y_min, y_max, y_log,
                                     lazy_reader)
        # These next two lines are painful.
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

    def add_fixed_resolution_plot(self, frb, field, center=None, use_colorbar=True,
                      figure = None, axes = None, fig_size=None, **kwargs):
        p = self._add_plot(PlotTypes.FixedResolutionPlot(frb, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        p["Axis"] = "na"
        return p

    def add_ortho_ray(self, axis, coords, field, axes = None,
                      figure = None, **kwargs):
        data_source = self.pf.h.ortho_ray(axis, coords, field)
        p = self._add_plot(PlotTypes.LineQueryPlot(data_source,
                [axis_names[axis], field], self._get_new_id(),
                figure, axes, plot_options=kwargs))
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

    @rootonly
    def save_book(self, filename):
        from pyPdf import PdfFileWriter, PdfFileReader
        outfile = PdfFileWriter()
        fns = self.save("__temp", format="pdf")
        concatenate_pdfs(filename, fns)
        for fn in fns: os.unlink(fn)

def wrap_pylab_newplot(func):
    @wraps(func)
    def pylabify(self, *args, **kwargs):
        # Let's assume that axes and figure are not in the positional
        # arguments -- probably safe!
        new_fig = self.pylab.figure()
        try:
            new_fig.canvas.set_window_title("%s" % (self.pf))
        except AttributeError:
            pass
        kwargs['axes'] = self.pylab.gca()
        kwargs['figure'] = self.pylab.gcf()
        retval = func(self, *args, **kwargs)
        retval._redraw_image()
        retval._fig_num = new_fig.number
        self.pylab.show()
        self.pylab.draw()
        return retval
    return pylabify

def wrap_pylab_show(func):
    @wraps(func)
    def pylabify(self, *args, **kwargs):
        retval = func(self, *args, **kwargs)
        fig_num = self.pylab.gcf().number
        for p in self.plots:
            self.pylab.figure(p._fig_num)
            self.pylab.draw()
        self.pylab.figure(fig_num)
        return retval
    return pylabify

class _Interactify(type):
    # All inherited methods get wrapped if they start with add_ or set_
    # So anything inheriting this automatically gets set up; additional
    # wrappings can be done manually.  Note that this does NOT modify
    # methods that are only in the subclasses.
    def __init__(cls, name, bases, d):
        super(_Interactify, cls).__init__(name, bases, d)
        for base in bases:
            for attrname in dir(base):
                if attrname in d: continue # If overridden, don't reset
                attr = getattr(cls, attrname)
                if type(attr) == types.MethodType:
                    if attrname.startswith("add_"):
                        setattr(cls, attrname, wrap_pylab_newplot(attr))
                    elif attrname.startswith("set_"):
                        setattr(cls, attrname, wrap_pylab_show(attr))

class PlotCollectionInteractive(PlotCollection):
    __metaclass__ = _Interactify

    autoscale = wrap_pylab_show(PlotCollection.autoscale)
    switch_field = wrap_pylab_show(PlotCollection.switch_field)

    def __init__(self, *args, **kwargs):
        import pylab
        self.pylab = pylab
        super(PlotCollectionInteractive, self).__init__(*args, **kwargs)

    def redraw(self):
        for plot in self.plots:
            plot._redraw_image()
        self.pylab.show()

    def clear_plots(self):
        for plot in self.plots:
            self.pylab.figure(plot._fig_num)
            self.pylab.clf()
        PlotCollection.clear_plots(self)

def get_multi_plot(nx, ny, colorbar = 'vertical', bw = 4, dpi=300):
    """
    This returns *nx* and *ny* axes on a single figure, set up so that the
    *colorbar* can be placed either vertically or horizontally in a bonus
    column or row, respectively.  The axes all have base width of *bw* inches.
    """
    PlotTypes.Initialize()
    hf, wf = 1.0/ny, 1.0/nx
    fudge_x = fudge_y = 1.0
    if colorbar.lower() == 'vertical':
        fudge_x = nx/(0.25+nx)
        fudge_y = 1.0
    elif colorbar.lower() == 'horizontal':
        fudge_x = 1.0
        fudge_y = ny/(0.40+ny)
    fig = matplotlib.figure.Figure((bw*nx/fudge_x, bw*ny/fudge_y), dpi=dpi)
    fig.set_canvas(be.engineVals["canvas"](fig))
    fig.subplots_adjust(wspace=0.0, hspace=0.0,
                        top=1.0, bottom=0.0,
                        left=0.0, right=1.0)
    tr = []
    print fudge_x, fudge_y
    for j in range(ny):
        tr.append([])
        for i in range(nx):
            left = i*wf*fudge_x
            bottom = fudge_y*(1.0-(j+1)*hf) + (1.0-fudge_y)
            ax = fig.add_axes([left, bottom, wf*fudge_x, hf*fudge_y])
            tr[-1].append(ax)
    cbars = []
    if colorbar.lower() == 'horizontal':
        for i in range(nx):
            # left, bottom, width, height
            # Here we want 0.10 on each side of the colorbar
            # We want it to be 0.05 tall
            # And we want a buffer of 0.15
            ax = fig.add_axes([wf*(i+0.10)*fudge_x, hf*fudge_y*0.20,
                               wf*(1-0.20)*fudge_x, hf*fudge_y*0.05])
            cbars.append(ax)
    elif colorbar.lower() == 'vertical':
        for j in range(nx):
            ax = fig.add_axes([wf*(nx+0.05)*fudge_x, hf*fudge_y*(ny-(j+0.95)),
                               wf*fudge_x*0.05, hf*fudge_y*0.90])
            ax.clear()
            cbars.append(ax)
    return fig, tr, cbars
