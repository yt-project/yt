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
    def __init__(self, pf, center=None):
        r"""The primary interface for creating plots.

        The PlotCollection object was created to ease the creation of multiple
        slices, projections and so forth made from a single parameter file.
        The concept is that when the width on one image changes, it should
        change on all the others.  The PlotCollection can create all plot types
        available in yt.

        Parameters
        ----------
        pf : `StaticOutput`
            The parameter file from which all the plots will be created.
        center : array_like, optional
            The 'center' supplied to plots like sphere plots, slices, and so
            on.  Should be 3 elements.  Defaults to the point of maximum
            density.
        Long_variable_name : {'hi', 'ho'}, optional
            Choices in brackets, default first when optional.

        Notes
        -----
        This class is the primary entry point to creating plots, but it is not
        the only entry point.  Additionally, creating a PlotCollection should
        be a "cheap" operation.

        You may iterate over the plots in the PlotCollection, via something
        like:

        >>> pc = PlotCollection(pf)
        >>> for p in pc: print p

        Examples
        --------

        >>> pc = PlotCollection(pf, center=[0.5, 0.5, 0.5])
        >>> pc.add_slice("Density", 0)
        >>> pc.save()

        """
        PlotTypes.Initialize()
        self.plots = []
        self.pf = pf
        if center == None:
            v,self.c = pf.h.find_max("Density") # @todo: ensure no caching
        elif center == "center" or center == "c":
            self.c = (pf["DomainRightEdge"] + pf["DomainLeftEdge"])/2.0
        else:
            self.c = na.array(center, dtype='float64')
        mylog.info("Created plot collection with default plot-center = %s",
                    list(self.c))

    def __iter__(self):
        for p in self.plots:
            yield p

    def save(self, basename=None, format="png", override=False, force_save=False):
        r"""Save out all the plots hanging off this plot collection, using
        generated names.

        This function will create names for every plot that belongs to the
        PlotCollection and save them out.  The names will, by default, be
        prefixed with the name of the affiliate parameter file, and the file
        names should indicate clearly what each plot represents.

        Parameters
        ----------
        basename : string, optional
            The prefix for all of the plot filenames.
        format : string, optional
            The plot file format.  Can be 'png', 'pdf', 'eps', 'jpg', or
            anything else that matplotlib understands.
        override : boolean
            If this is true, then no generated filenames will be appended to
            the base name.  You probably don't want this.
        force_save : boolean
            In parallel, only the root task (proc 0) saves an image, unless
            this is set to True.

        Returns
        -------
        items : string
            This function returns a list of the filenames created.

        Examples
        --------

        >>> fns = pc.save()
        >>> for fn in fns: print "Saved", fn
        """
        if basename is None: basename = str(self.pf)
        fn = []
        for plot in self.plots:
            fn.append(plot.save_image(basename, format=format, 
                      override=override, force_save=force_save))
            mylog.info("Saved %s", fn[-1])
        return fn

    def set_xlim(self, xmin, xmax):
        r"""Set the x-limits of all plots.

        set_xlim on all plots is called with the parameters passed to this
        function.

        Parameters
        ----------
        xmin : float
            The left boundary for the x axis.
        xmax : float
            The right boundary for the x axis.
        """
        for plot in self.plots:
            plot.set_xlim(xmin, xmax)

    def set_ylim(self, ymin, ymax):
        r"""Set the y-limits of all plots.

        set_ylim on all plots is called with the parameters passed to this
        function.

        Parameters
        ----------
        ymin : float
            The left boundary for the x axis.
        ymax : float
            The right boundary for the x axis.
        """
        for plot in self.plots:
            plot.set_ylim(ymin, ymax)

    def set_zlim(self, zmin, zmax, *args, **kwargs):
        """
        Set the limits of the colorbar. 'min' or 'max' are possible inputs 
        when combined with dex=value, where value gives the maximum number of 
        dex to go above/below the min/max.  If value is larger than the true
        range of values, min/max are limited to true range.

        Only ONE of the following options can be specified. If all 3 are
        specified, they will be used in the following precedence order:
            ticks - a list of floating point numbers at which to put ticks
            minmaxtick - display DEFAULT ticks with min & max also displayed
            nticks - if ticks not specified, can automatically determine a
               number of ticks to be evenly spaced in log space
        """
        for plot in self.plots:
            plot.set_autoscale(False)
            plot.set_zlim(zmin, zmax, *args, **kwargs)

    def set_lim(self, lim):
        r"""Set the x- and y-limits of all plots.

        set_xlim on all plots is called with the parameters passed to this
        function, and then set_ylim is called.

        Parameters
        ----------
        lim : tuple of floats
            (xmin, xmax, ymin, ymax)
        """
        for plot in self.plots:
            plot.set_xlim(*lim[:2])
            plot.set_ylim(*lim[2:])

    def autoscale(self):
        r"""Turn on autoscaling on all plots.

        This has the same effect as:

        >>> for p in pc: p.set_autoscale(True)

        By default, all plots are autoscaled until the colorbar is set
        manually.  This turns autoscaling back on.  The colors may not be
        updated unless _redraw_image is called on the plots, which should occur
        with a change in the width or saving of images.
        """
        for plot in self.plots:
            plot.set_autoscale(True)

    def set_width(self, width, unit):
        r"""Change the width of all image plots.

        This function changes all the widths of the image plots (but notably
        not any phase plots or profile plots) to be a given physical extent.

        Parameters
        ----------
        width : float
            The numeric value of the new width.
        unit : string
            The unit corresponding to the width given.
        """
        for plot in self.plots:
            plot.set_width(width, unit)

    def set_cmap(self, cmap):
        r"""Change the colormap of all plots.

        This function will update the colormap on all plots for which a
        colormap makes sense.  The colors may not be updated unless
        _redraw_image is called on the plots, which should occur with a change
        in the width or saving of images.

        Parameters
        ----------
        cmap : string
            An acceptable colormap.  See either raven.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        """
        for plot in self.plots:
            plot.set_cmap(cmap)

    def switch_field(self, field):
        r"""Change the displayed of all image plots.

        All images that display a field -- slices, cutting planes, projections
        -- will be switched to display the specified field.  For projections,
        this will re-generate the projection, if it is unable to load the
        projected field off-disk.

        Parameters
        ----------
        field : string
            Any field that can be generated or read from disk.
        """
        for plot in self.plots:
            plot.switch_z(field)
    switch_z = switch_field

    def _add_plot(self, plot):
        r"""This function adds a plot to the plot collection.

        This function is typically used internally to add a plot on to the
        current list of plots.  However, if you choose to manually create a
        plot, this can be used to add it to a collection for convenient
        modification.

        Parameters
        ----------
        plot : `yt.raven.RavenPlot`
            A plot, which will be appended to the list of plots handled by this
            plot collection.

        Returns
        -------
        plot : `yt.raven.RavenPlot`
            The plot handed to the function is passed back through.  This is
            unnecessary, but is done for historical reasons.
        """
        self.plots.append(plot)
        return plot

    def add_slice(self, field, axis, coord=None, center=None,
                 use_colorbar=True, figure = None, axes = None, fig_size=None,
                 periodic = True, data_source = None, field_parameters = None):
        r"""Create a slice, from that a slice plot, and add it to the current
        collection.

        This function will generate a `yt.lagos.AMRSliceBase` from the given
        parameters.  This slice then gets passed to a `yt.raven.SlicePlot`, and
        the resultant plot is added to the current collection.  Various
        parameters allow control of the way the slice is displayed, as well as
        how the slice is generated.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 for x, y, z.
        coord : float, optional
            The coordinate to place the slice at, along the slicing axis.
        center : array_like, optional
            The center to be used for things like radius and radial velocity.
            Defaults to the center of the plot collection.
        use_colorbar : bool, optional
            Whether we should leave room for and create a colorbar.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        fig_size : tuple of floats
            This parameter can act as a proxy for the manual creation of a
            figure.  By specifying it, you can create plots with an arbitrarily
            large or small size.  It is in inches, defaulting to 100 dpi.
        periodic : boolean, optional
            By default, the slices are assumed to be periodic, and they will
            wrap around the edges.
        data_source : `yt.lagos.AMRSliceBase`, optional
            If you would like to use an existing slice, you may specify it
            here, in which case a new slice will not be created.
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.

        Returns
        -------
        plot : `yt.raven.SlicePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.lagos.AMRSliceBase : This is the type created by this function and 
                                passed to the plot created here.

        Notes
        -----
        This is the primary mechanism for creating slice plots, and generating
        slice plots along multiple axes was the original purpose of the
        PlotCollection.

        Note that all plots can be modified.  See `callback_list` for more
        information.

        Examples
        --------

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_slice("Density", 0)
        """
        if center == None:
            center = self.c
        if coord == None:
            coord = center[axis]
        if data_source is None:
            if field_parameters == None: field_parameters = {}
            data_source = self.pf.hierarchy.slice(axis, coord, field,
                            center=center, **field_parameters)
        p = self._add_plot(PlotTypes.SlicePlot(
                         data_source, field, use_colorbar=use_colorbar,
                         axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        mylog.info("Added slice of %s at %s = %s with 'center' = %s", field,
                    axis_names[axis], coord, list(center))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_particles(self, axis, width, p_size=1.0, col='k', stride=1.0,
                      data_source=None, figure=None, axes=None):
        r"""Create a plot of a thick slab of particles.

        This function will generate a `yt.lagos.AMRRegionBase` from the given
        parameters, and all particles which are within that region will be
        plotted.

        Parameters
        ----------
        axis : int
            The axis along which to create the thick slab.  Can be 0, 1, or 2
            for x, y, z.
        width : float
            The width of the thick slab, in code units, from which particles
            will be plotted.
        p_size : float, optional
            The size of the points to be used to represent the particles, in
            pixels.
        col : color, optional
            Specified in matplotlib color specifications, the color that
            particles should be.
        stride : float, optional
            The stride through the particles to plot.  Used to plot every
            fifth, every tenth, etc.  Note that the sorted order of particles
            may result in a biased selection of particles.
        data_source : `yt.lagos.AMRData`, optional
            If specified, this will be the data source used for obtaining
            particles.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.

        Returns
        -------
        plot : `yt.raven.ParticlePlot`
            The plot that has been added to the PlotCollection.

        Notes
        -----
        This plot type can be very expensive, and does not necessarily produce
        the best visual results.  Plotting a large number of particles can be
        very tricky, and often it's much better to instead use a slice or a
        (thin) projection of deposited density, like particle_density_pyx.

        Examples
        --------

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_particles(0, 1.0)
        """
        LE = self.pf["DomainLeftEdge"].copy()
        RE = self.pf["DomainRightEdge"].copy()
        LE[axis] = self.c[axis] - width/2.0
        RE[axis] = self.c[axis] + width/2.0
        if data_source is None: data_source = self.pf.h.region(self.c, LE, RE)
        data_source.axis = axis
        p = self._add_plot(PlotTypes.ParticlePlot(data_source, axis,
                                        width, p_size, col, stride, figure,
                                        axes))
        p["Axis"] = lagos.axis_names[axis]
        return p

    def add_cutting_plane(self, field, normal,
                          center=None, use_colorbar=True,
                          figure = None, axes = None, fig_size=None, obj=None,
                           field_parameters = None):
        r"""Create a cutting plane, from that a plot, and add it to the current
        collection.

        A cutting plane is an oblique slice through the simulation volume,
        oriented by a specified normal vector that is perpendicular to the
        image plane.  This function will generate a
        `yt.lagos.AMRCuttingPlaneBase` from the given parameters.  This cutting
        plane then gets passed to a `yt.raven.CuttingPlanePlot`, and the
        resultant plot is added to the current collection.  Various parameters
        allow control of the way the slice is displayed, as well as how the
        plane is generated.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        normal : array_like
            The vector that defines the desired plane.  For instance, the
            angular momentum of a sphere.
        center : array_like, optional
            The center to be used for things like radius and radial velocity.
            Defaults to the center of the plot collection.
        use_colorbar : bool, optional
            Whether we should leave room for and create a colorbar.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        fig_size : tuple of floats
            This parameter can act as a proxy for the manual creation of a
            figure.  By specifying it, you can create plots with an arbitrarily
            large or small size.  It is in inches, defaulting to 100 dpi.
        obj : `AMRCuttingPlaneBase`, optional
            If you would like to use an existing cutting plane, you may specify
            it here, in which case a new cutting plane will not be created.
        field_parameters : dict, optional
            This set of parameters will be passed to the cutting plane upon
            creation, which can be used for passing variables to derived
            fields.

        Returns
        -------
        plot : `yt.raven.CuttingPlanePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.lagos.AMRCuttingPlaneBase : This is the type created by this function.

        Notes
        -----
        This is the primary mechanism for creating cutting plane plots.  Note
        that they are somewhat slow, but useful to orient the image in an
        arbitrary direction.

        Note that all plots can be modified.  See `callback_list` for more
        information.

        Examples
        --------

        Here's a simple mechanism for getting the angular momentum of a
        collapsing cloud and generating a cutting plane aligned with the
        angular momentum vector.

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> v, c = pf.h.find_max("Density")
        >>> sp = pf.h.sphere(c, 1000.0/pf['au'])
        >>> L = sp.quantities["AngularMomentumVector"]()
        >>> pc = PlotCollection(pf)
        >>> p = pc.add_cutting_plane("Density", L)
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

    def add_fixed_res_cutting_plane(self, field, normal, width, res=512,
             center=None, use_colorbar=True, figure = None, axes = None,
             fig_size=None, obj=None, **kwargs):
        """
        Generate a fixed resolution, interpolated cutting plane of
        *field* with *normal*, centered at *center* (defaults to
        PlotCollection center) with *use_colorbar* specifying whether
        the plot is naked or not and optionally providing pre-existing
        Matplotlib *figure* and *axes* objects.  *fig_size* in
        (height_inches, width_inches).  If so desired, *obj* is a
        pre-existing cutting plane object.
        """
        if center == None:
            center = self.c
        if not obj:
            data = self.pf.hierarchy.fixed_res_cutting \
                 (normal, center, width, res, **kwargs)
            #data = frc[field]
        else:
            data = obj
        p = self._add_plot(PlotTypes.FixedResolutionPlot(data, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        mylog.info("Added fixed-res plane of %s with 'center' = %s and "
                   "normal = %s", field, list(center), list(normal))
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
                           lazy_reader=True, id=None,
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
                            fields[0], non_zero = x_log,
                            lazy_reader=lazy_reader)[0]
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
                               lazy_reader=True, id=None,
                               axes = None, figure = None,
                               fractional=False):
        """
        Given a *data_source*, and *fields*, automatically generate a 2D
        profile and plot it.  *id* is used internally to add onto the prefix,
        and will be automatically generated if not given. Remainder of
        arguments are identical to :meth:`add_profile_object`.
        """
        if x_bounds is None:
            x_min, x_max = data_source.quantities["Extrema"](
                                    fields[0], non_zero = x_log,
                                    lazy_reader=lazy_reader)[0]
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = data_source.quantities["Extrema"](
                                    fields[1], non_zero = y_log,
                                    lazy_reader=lazy_reader)[0]
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
            p.switch_z(fields[2], weight=weight, accumulation=accumulation, fractional=fractional)
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

    def add_scatter_source(self, data_source, fields, id=None,
                    axes = None, figure = None, plot_options = None):
        """
        Given a *data_source*, and *fields*, plot a scatter plot.
        *plot_options* are sent to the scatter command.
        """
        if id is None: id = self._get_new_id()
        sp = PlotTypes.ScatterPlot(data_source, fields, id,
                                   plot_options = plot_options,
                                   figure=figure, axes=axes)
        p = self._add_plot(sp)
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

    def add_ray(self, start_point, end_point, field, axes = None,
                figure = None, **kwargs):
        data_source = self.pf.h.ray(start_point, end_point, field)
        p = self._add_plot(PlotTypes.LineQueryPlot(data_source,
                ['t', field], self._get_new_id(),
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
    def save_book(self, filename, info = None):
        """
        This will save out a single PDF, where each page is a plot object.  The
        *info* keyword can be a dictionary composed of the keys and values
        "Author", "Title", "Subject", "Keywords", "Creator", "Producer" ad
        "CreationDate".  Any keywords not filled in will be blank.  The default
        is to use the current settings in Matplotlib for filling them in.
        """
        from matplotlib.backends.backend_pdf import PdfPages
        outfile = PdfPages(filename)
        for plot in self.plots:
            plot.save_to_pdf(outfile)
        if info is not None:
            outfile._file.writeObject(outfile._file.infoObject, info)
        outfile.close()

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
    if colorbar is None:
        fudge_x = fudge_y = 1.0
    elif colorbar.lower() == 'vertical':
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
    if colorbar is None:
        pass
    elif colorbar.lower() == 'horizontal':
        for i in range(nx):
            # left, bottom, width, height
            # Here we want 0.10 on each side of the colorbar
            # We want it to be 0.05 tall
            # And we want a buffer of 0.15
            ax = fig.add_axes([wf*(i+0.10)*fudge_x, hf*fudge_y*0.20,
                               wf*(1-0.20)*fudge_x, hf*fudge_y*0.05])
            cbars.append(ax)
    elif colorbar.lower() == 'vertical':
        for j in range(ny):
            ax = fig.add_axes([wf*(nx+0.05)*fudge_x, hf*fudge_y*(ny-(j+0.95)),
                               wf*fudge_x*0.05, hf*fudge_y*0.90])
            ax.clear()
            cbars.append(ax)
    return fig, tr, cbars
