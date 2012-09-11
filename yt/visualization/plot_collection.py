""" 
All of the base-level stuff for plotting.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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

from matplotlib import figure
import shutil
import tempfile
import numpy as np
import os

from yt.funcs import *

from yt.config import ytcfg
from yt.data_objects.profiles import \
    BinnedProfile1D, \
    BinnedProfile2D
from yt.utilities.definitions import \
    axis_names, inv_axis_names, x_dict, y_dict
from .plot_types import \
    FixedResolutionPlot, \
    PCSlicePlot, \
    PCSlicePlotNaturalNeighbor, \
    PCProjectionPlot, \
    PCProjectionPlotNaturalNeighbor, \
    CuttingPlanePlot, \
    ParticlePlot, \
    ProfilePlot, \
    Profile1DPlot, \
    PhasePlot, \
    LineQueryPlot, \
    ScatterPlot
from yt.utilities.minimal_representation import \
    MinimalImageCollectionData

# No better place to put this
def concatenate_pdfs(output_fn, input_fns):
    from pyPdf import PdfFileWriter, PdfFileReader
    outfile = PdfFileWriter()
    for fn in input_fns:
        infile = PdfFileReader(open(fn, 'rb'))
        outfile.addPage(infile.getPage(0))
    outfile.write(open(output_fn, "wb"))

class ImageCollection(object):
    def __init__(self, pf, name):
        self.pf = pf
        self.name = name
        self.images = []
        self.image_metadata = []

    def add_image(self, fn, descr):
        self.image_metadata.append(descr)
        self.images.append((os.path.basename(fn), np.fromfile(fn, dtype='c')))

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
        >>> pc.add_slice("Density", 'x')
        >>> pc.save()

        """
        self.plots = []
        self.pf = pf
        if center == None:
            v,self.c = pf.h.find_max("Density") # @todo: ensure no caching
        elif center == "center" or center == "c":
            self.c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
        else:
            self.c = np.array(center, dtype='float64')
        mylog.info("Created plot collection with default plot-center = %s",
                    list(self.c))

    def __iter__(self):
        for p in self.plots:
            yield p

    @property
    def _mrep(self):
        ic = ImageCollection(self.pf, "Plot Collection with center %s" % self.c)
        dd = tempfile.mkdtemp()
        fns = self.save(os.path.join(dd, "temp"))
        for fn, p in zip(fns, self.plots):
            ic.add_image(fn, p._pretty_name())
        shutil.rmtree(dd)
        return MinimalImageCollectionData(ic)

    def hub_upload(self):
        self._mrep.upload()

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
        if ytcfg.getboolean("yt", "__withinreason"):
            from yt.gui.reason.bottle_mods import PayloadHandler
            import base64
            ph = PayloadHandler()
            for f in fn:
                if not f.endswith('png'): continue
                img_data = base64.b64encode(open(f,'rb').read())
                payload = {'type':'png_string',
                           'image_data':img_data}
                ph.add_payload(payload)
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

        * ``ticks`` - a list of floating point numbers at which to put ticks
        * ``minmaxtick`` - display DEFAULT ticks with min & max also displayed
        * ``nticks`` - if ticks not specified, can automatically determine a
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
            The unit in which the given width is expressed.
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
            An acceptable colormap.  See either yt.visualization.color_maps or
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
        plot : `yt.visualization.plot_types.RavenPlot`
            A plot, which will be appended to the list of plots handled by this
            plot collection.

        Returns
        -------
        plot : `yt.visualization.plot_types.RavenPlot`
            The plot handed to the function is passed back through.  This is
            unnecessary, but is done for historical reasons.
        """
        self.plots.append(plot)
        return plot

    def add_slice(self, field, axis, coord=None, center=None,
                 use_colorbar=True, figure = None, axes = None, fig_size=None,
                 periodic = True, obj = None, field_parameters = None):
        r"""Create a slice, from that a slice plot, and add it to the current
        collection.

        This function will generate a `yt.data_objects.api.YTSliceBase` from the given
        parameters.  This slice then gets passed to a `yt.visualization.plot_types.PCSlicePlot`, and
        the resultant plot is added to the current collection.  Various
        parameters allow control of the way the slice is displayed, as well as
        how the slice is generated.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 or x, y, z.
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
        obj : `yt.data_objects.api.YTSliceBase`, optional
            If you would like to use an existing slice, you may specify it
            here, in which case a new slice will not be created.
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.

        Returns
        -------
        plot : `yt.visualization.plot_types.PCSlicePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTSliceBase : This is the type created by this function and
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
        >>> p = pc.add_slice("Density", 'x')
        """
        axis = fix_axis(axis)
        if center == None:
            center = self.c
        if coord == None:
            coord = center[axis]
        if obj is None:
            if field_parameters is None: field_parameters = {}
            obj = self.pf.hierarchy.slice(axis, coord, field,
                            center=center, **field_parameters)
        p = self._add_plot(PCSlicePlot(
                         obj, field, use_colorbar=use_colorbar,
                         axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        mylog.info("Added slice of %s at %s = %s with 'center' = %s", field,
                    axis_names[axis], coord, list(center))
        p["Axis"] = axis_names[axis]
        return p

    def add_particles(self, axis, width, p_size=1.0, col='k', stride=1.0,
                      data_source=None, figure=None, axes=None):
        r"""Create a plot of a thick slab of particles.

        This function will generate a `yt.data_objects.api.YTRegionBase` from the given
        parameters, and all particles which are within that region will be
        plotted.

        Parameters
        ----------
        axis : int
            The axis along which to create the thick slab.  Can be 0, 1, or 2
            or x, y, z.
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
        data_source : `yt.data_objects.api.YTDataContainer`, optional
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
        plot : `yt.visualization.plot_types.ParticlePlot`
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
        axis = fix_axis(axis)
        LE = self.pf.domain_left_edge.copy()
        RE = self.pf.domain_right_edge.copy()
        LE[axis] = self.c[axis] - width/2.0
        RE[axis] = self.c[axis] + width/2.0
        if data_source is None: data_source = self.pf.h.region(self.c, LE, RE)
        data_source.axis = axis
        p = self._add_plot(ParticlePlot(data_source, axis,
                                        width, p_size, col, stride, figure,
                                        axes))
        p["Axis"] = axis_names[axis]
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
        `yt.data_objects.api.YTCuttingPlaneBase` from the given parameters.  This cutting
        plane then gets passed to a `yt.visualization.plot_types.CuttingPlanePlot`, and the
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
        obj : `YTCuttingPlaneBase`, optional
            If you would like to use an existing cutting plane, you may specify
            it here, in which case a new cutting plane will not be created.
        field_parameters : dict, optional
            This set of parameters will be passed to the cutting plane upon
            creation, which can be used for passing variables to derived
            fields.

        Returns
        -------
        plot : `yt.visualization.plot_types.CuttingPlanePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTCuttingPlaneBase : This is the type created by this function.

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
            if field_parameters is None: field_parameters = {}
            cp = self.pf.hierarchy.cutting(normal, center, field,
                    **field_parameters)
        else:
            cp = obj
        p = self._add_plot(CuttingPlanePlot(cp, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        mylog.info("Added plane of %s with 'center' = %s and normal = %s", field,
                    list(center), list(normal))
        p["Axis"] = "CuttingPlane"
        return p

    def add_fixed_res_cutting_plane(self, field, normal, width, res=512,
             center=None, use_colorbar=True, figure = None, axes = None,
             fig_size=None, obj=None, field_parameters = None):
        r"""Create a fixed resolution cutting plane, from that a plot, and add
        it to the current collection.

        A cutting plane is an oblique slice through the simulation volume,
        oriented by a specified normal vector that is perpendicular to the
        image plane.  This function will slice through, but instead of
        retaining all the data necessary to rescale the cutting plane at any
        width, it only retains the pixels for a single width.  This function
        will generate a `yt.data_objects.api.YTFixedResCuttingPlaneBase` from the given
        parameters.  This image buffer then gets passed to a
        `yt.visualization.plot_types.FixedResolutionPlot`, and the resultant plot is added to the
        current collection.  Various parameters allow control of the way the
        slice is displayed, as well as how the plane is generated.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        normal : array_like
            The vector that defines the desired plane.  For instance, the
            angular momentum of a sphere.
        width : float
            The width, in code units, of the image plane.
        res : int
            The returned image buffer must be square; this number is how many
            pixels on a side it will have.
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
        obj : `YTCuttingPlaneBase`, optional
            If you would like to use an existing cutting plane, you may specify
            it here, in which case a new cutting plane will not be created.
        field_parameters : dict, optional
            This set of parameters will be passed to the cutting plane upon
            creation, which can be used for passing variables to derived
            fields.

        Returns
        -------
        plot : `yt.visualization.plot_types.FixedResolutionPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTFixedResCuttingPlaneBase : This is the type created by this
                                               function.

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
        >>> p = pc.add_fixed_res_cutting_plane("Density", L, 1000.0/pf['au'])
        """
        if center == None:
            center = self.c
        if not obj:
            if field_parameters is None: field_parameters = {}
            data = self.pf.hierarchy.fixed_res_cutting \
                 (normal, center, width, res, **field_parameters)
            #data = frc[field]
        else:
            data = obj
        p = self._add_plot(FixedResolutionPlot(data, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        mylog.info("Added fixed-res plane of %s with 'center' = %s and "
                   "normal = %s", field, list(center), list(normal))
        p["Axis"] = "CuttingPlane"
        return p

    def add_projection(self, field, axis,  weight_field=None,
                       data_source = None,
                       center=None, use_colorbar=True,
                       figure = None, axes = None, fig_size=None,
                       periodic = True, obj = None, field_parameters = None):
        r"""Create a projection, from that a projection plot, and add it to the
        current collection.

        This function will generate a `yt.data_objects.api.YTOverlapProjBase` from the given
        parameters.  This projection then gets passed to a
        `yt.visualization.plot_types.PCProjectionPlot`, and the resultant plot is added to the
        current collection.  Various parameters allow control of the way the
        slice is displayed, as well as how the slice is generated.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 or x, y, z.
        data_source : `yt.data_objects.api.YTDataContainer`
            This is a data source respecting the `YTDataContainer` protocol (i.e., it
            has grids and so forth) that will be used as input to the
            projection.
        weight_field : string
            If specified, this will be the weighting field and the resultant
            projection will be a line-of-sight average, defined as sum( f_i *
            w_i * dl ) / sum( w_i * dl )
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
        obj : `yt.data_objects.api.YTOverlapProjBase`, optional
            If you would like to use an existing projection, you may specify it
            here, in which case a new projection will not be created.  If this
            option is specified the options data_source, weight_field and
            field_parameters will be ignored.
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.

        Returns
        -------
        plot : `yt.visualization.plot_types.PCProjectionPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTOverlapProjBase : This is the type created by this function and
                               passed to the plot created here.

        Notes
        -----
        This is the primary mechanism for creating projection plots, and
        generating projection plots along multiple axes was the original
        purpose of the PlotCollection.

        Note that all plots can be modified.  See `callback_list` for more
        information.

        Examples
        --------

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_projection("Density", 'x', "Density")
        """
        axis = fix_axis(axis)
        if field_parameters is None: field_parameters = {}
        if center == None:
            center = self.c
        if obj is None:
            obj = self.pf.hierarchy.proj(field, axis, weight_field,
                                         source = data_source, center=center,
                                         **field_parameters)
        p = self._add_plot(PCProjectionPlot(obj, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        p["Axis"] = axis_names[axis]
        return p

    def add_thin_projection(self, field, axis, thickness,
                       weight_field=None, center=None, use_colorbar=True,
                       figure = None, axes = None, fig_size=None,
                       periodic = True, field_parameters = None):
        r"""Create a projection through a thin slice of the domain, from that a
        projection plot, and add it to the current collection.

        This function will generate a rectangular prism region and supply it to
        a`yt.data_objects.api.YTOverlapProjBase` from the given parameters.  This projection
        then gets passed to a `yt.visualization.plot_types.PCProjectionPlot`, and the resultant plot
        is added to the current collection.  Various parameters allow control
        of the way the slice is displayed, as well as how the slice is
        generated.  The center is used as the center of the thin projection.

        Parameters
        ----------
        field : string
            The initial field to slice and display.
        axis : int
            The axis along which to slice.  Can be 0, 1, or 2 or x, y, z.
        thickness : float
            In 'code units' this is the thickness of the region to be
            projected through.
        weight_field : string
            If specified, this will be the weighting field and the resultant
            projection will be a line-of-sight average, defined as sum( f_i *
            w_i * dl ) / sum( w_i * dl )
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
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.

        Returns
        -------
        plot : `yt.visualization.plot_types.PCProjectionPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTOverlapProjBase : This is the type created by this function and
                               passed to the plot created here.

        Notes
        -----
        This is the primary mechanism for creating projection plots, and
        generating projection plots along multiple axes was the original
        purpose of the PlotCollection.

        Note that all plots can be modified.  See `callback_list` for more
        information.

        Examples
        --------

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_thin_projection("Density", 0, 0.1, "Density")
        """
        axis = fix_axis(axis)
        if field_parameters is None: field_parameters = {}
        if center == None:
            center = self.c
        LE = self.pf.domain_left_edge.copy()
        RE = self.pf.domain_right_edge.copy()
        LE[axis] = RE[axis] = center[axis]
        LE[axis] -= thickness/2.0
        RE[axis] += thickness/2.0
        region = self.pf.h.region(center, LE, RE)
        obj = self.pf.hierarchy.proj(field, axis, weight_field,
                                     source = region, center=center,
                                     **field_parameters)
        p = self._add_plot(PCProjectionPlot(obj, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size, periodic=periodic))
        p["Axis"] = axis_names[axis]
        return p

    def add_profile_object(self, data_source, fields,
                           weight="CellMassMsun", accumulation=False,
                           x_bins=128, x_log=True, x_bounds=None,
                           id=None, figure=None, axes=None):
        r"""From an existing object, create a 1D, binned profile.

        This function will accept an existing `YTDataContainer` source and from that,
        it will generate a `Binned1DProfile`, based on the specified options.
        This is useful if you have extracted a region, or if you wish to bin
        some set of massages data -- or even if you wish to bin anything other
        than a sphere.  The profile will be 1D, which means while it can have
        an arbitrary number of fields, those fields will all be binned based on
        a single field.

        Parameters
        ----------
        data_source : `yt.data_objects.api.YTDataContainer`
            This is a data source respecting the `YTDataContainer` protocol (i.e., it
            has grids and so forth) that will be used as input to the profile
            generation.
        fields : list of strings
            The first element of this list is the field by which we will bin;
            all subsequent fields will be binned and their profiles added to
            the underlying `BinnedProfile1D`.
        weight : string, default "CellMassMsun"
            The weighting field for an average.  This defaults to mass-weighted
            averaging.
        accumulation : boolean, optional
            If true, from the low-value to the high-value the values in all
            binned fields will be accumulated.  This is useful for instance
            when adding an unweighted CellMassMsun to a radial plot, as it will
            show mass interior to that radius.
        x_bins : int, optional
            How many bins should there be in the independent variable?
        x_log : boolean, optional
            Should the bin edges be log-spaced?
        x_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.

        Returns
        -------
        plot : `yt.visualization.plot_types.ProfilePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.profiles.BinnedProfile1D : This is the object that does the
                                   transformation of raw data into a 1D
                                   profile.

        Examples
        --------

        >>> reg = pf.h.region([0.1, 0.2, 0.3], [0.0, 0.1, 0.2],
                              [0.2, 0.3, 0.4])
        >>> pc.add_profile_object(reg, ["Density", "Temperature"])
        """
        if x_bounds is None:
            x_min, x_max = data_source.quantities["Extrema"](
                            fields[0], non_zero = x_log)[0]
        else:
            x_min, x_max = x_bounds
        profile = BinnedProfile1D(data_source,
                                  x_bins, fields[0], x_min, x_max, x_log)
        if len(fields) > 1:
            profile.add_fields(fields[1:], weight=weight, accumulation=accumulation)
        if id is None: id = self._get_new_id()
        p = self._add_plot(Profile1DPlot(profile, fields, id,
                                                   axes=axes, figure=figure))
        return p

    def add_profile_sphere(self, radius, unit, fields, center = None,
                           weight="CellMassMsun", accumulation=False,
                           x_bins=128, x_log=True, x_bounds=None,
                           id=None, figure=None, axes=None):
        r"""From a description of a sphere, create a 1D, binned profile.

        This function will accept the radius of a sphere, and from that it will
        generate a `Binned1DProfile`, based on the specified options.  The
        profile will be 1D, which means while it can have an arbitrary number
        of fields, those fields will all be binned based on a single field.

        All subsequent parameters beyond "unit" will be passed verbatim to
        add_profile_object.

        Parameters
        ----------
        radius : float
            The radius of the sphere to generate.
        unit : string
            The unit in which the given radius is expressed.
        fields : list of strings
            The first element of this list is the field by which we will bin;
            all subsequent fields will be binned and their profiles added to
            the underlying `BinnedProfile1D`.
        center : array_like, optional
            The center to be used for things like radius and radial velocity.
            Defaults to the center of the plot collection.
        weight : string, default "CellMassMsun"
            The weighting field for an average.  This defaults to mass-weighted
            averaging.
        accumulation : boolean, optional
            If true, from the low-value to the high-value the values in all
            binned fields will be accumulated.  This is useful for instance
            when adding an unweighted CellMassMsun to a radial plot, as it will
            show mass interior to that radius.
        x_bins : int, optional
            How many bins should there be in the independent variable?
        x_log : boolean, optional
            Should the bin edges be log-spaced?
        x_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.

        Returns
        -------
        plot : `yt.visualization.plot_types.ProfilePlot`
            The plot that has been added to the PlotCollection.  Note that the
            underlying sphere may be accessed as .data.data_source

        See Also
        --------
        yt.data_objects.profiles.BinnedProfile1D : This is the object that does the
                                   transformation of raw data into a 1D
                                   profile.
        yt.data_objects.api.YTSphereBase : This is the object auto-generated by this
                                 function.

        Examples
        --------

        >>> pc.add_profile_sphere(1.0, 'kpc', ["Density", "Electron_Fraction"])
        """
        if center is None:
            center = self.c
        r = radius/self.pf[unit]
        sphere = self.pf.hierarchy.sphere(center, r)
        p = self.add_profile_object(sphere, fields, weight, accumulation,
                           x_bins, x_log, x_bounds, id,
                           figure=figure, axes=axes)
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p

    def add_phase_object(self, data_source, fields, cmap=None,
                               weight="CellMassMsun", accumulation=False,
                               x_bins=128, x_log=True, x_bounds=None,
                               y_bins=128, y_log=True, y_bounds=None,
                               id=None, axes = None, figure = None,
                               fractional=False):
        r"""From an existing object, create a 2D, binned profile.

        This function will accept an existing `YTDataContainer` source and from that,
        it will generate a `Binned2DProfile`, based on the specified options.
        This is useful if you have extracted a region, or if you wish to bin
        some set of massages data -- or even if you wish to bin anything other
        than a sphere.  The profile will be 2D, which means while it can have
        an arbitrary number of fields, those fields will all be binned based on
        two fields.

        Parameters
        ----------
        data_source : `yt.data_objects.api.YTDataContainer`
            This is a data source respecting the `YTDataContainer` protocol (i.e., it
            has grids and so forth) that will be used as input to the profile
            generation.
        fields : list of strings
            The first element of this list is the field by which we will bin
            into the x-axis, the second is the field by which we will bin onto
            the y-axis.  All subsequent fields will be binned and their
            profiles added to the underlying `BinnedProfile2D`.
        cmap : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        weight : string, default "CellMassMsun"
            The weighting field for an average.  This defaults to mass-weighted
            averaging.
        accumulation : list of booleans, optional
            If true, from the low-value to the high-value the values in all
            binned fields will be accumulated.  This is useful for instance
            when adding an unweighted CellMassMsun to a radial plot, as it will
            show mass interior to that radius.  The first value is for the
            x-axis, the second value for the y-axis.  Note that accumulation
            will only be along each row or column.
        x_bins : int, optional
            How many bins should there be in the x-axis variable?
        x_log : boolean, optional
            Should the bin edges be log-spaced?
        x_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        y_bins : int, optional
            How many bins should there be in the y-axis variable?
        y_log : boolean, optional
            Should the bin edges be log-spaced?
        y_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        fractional : boolean
            If true, the plot will be normalized to the sum of all the binned
            values.

        Returns
        -------
        plot : `yt.visualization.plot_types.PlotTypes.PhasePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.profiles.BinnedProfile2D : This is the object that does the
                                   transformation of raw data into a 1D
                                   profile.
        
        Examples
        --------
        This will show the mass-distribution in the Density-Temperature plane.
        
        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> reg = pf.h.region([0.1, 0.2, 0.3], [0.0, 0.1, 0.2],
        ...                   [0.2, 0.3, 0.4])
        >>> pc.add_phase_object(reg, ["Density", "Temperature", "CellMassMsun"],
        ...                     weight = None)
        """
        if x_bounds is None:
            x_min, x_max = data_source.quantities["Extrema"](
                                    fields[0], non_zero = x_log)[0]
        else:
            x_min, x_max = x_bounds
        if y_bounds is None:
            y_min, y_max = data_source.quantities["Extrema"](
                                    fields[1], non_zero = y_log)[0]
        else:
            y_min, y_max = y_bounds
        profile = BinnedProfile2D(data_source,
                                  x_bins, fields[0], x_min, x_max, x_log,
                                  y_bins, fields[1], y_min, y_max, y_log)
        # This will add all the fields to the profile object
        if len(fields)>2:
            profile.add_fields(fields[2:], weight=weight,
                    accumulation=accumulation, fractional=fractional)

        if id is None: id = self._get_new_id()
        p = self._add_plot(PhasePlot(profile, fields, 
                                               id, cmap=cmap,
                                               figure=figure, axes=axes))
        return p

    def add_phase_sphere(self, radius, unit, fields, center = None, cmap=None,
                         weight="CellMassMsun", accumulation=False,
                         x_bins=128, x_log=True, x_bounds=None,
                         y_bins=128, y_log=True, y_bounds=None,
                         id=None, axes = None, figure = None,
                         fractional=False):
        r"""From a description of a sphere, create a 2D, binned profile.

        This function will accept the radius of a sphere, and from that it will
        generate a `Binned1DProfile`, based on the specified options.  The
        profile will be 2D, which means while it can have an arbitrary number
        of fields, those fields will all be binned based on two fields.

        All subsequent parameters beyond "unit" will be passed verbatim to
        add_profile_object.

        Parameters
        ----------
        radius : float
            The radius of the sphere to generate.
        unit : string
            The unit in which the given radius is expressed.
        fields : list of strings
            The first element of this list is the field by which we will bin
            into the y-axis, the second is the field by which we will bin onto
            the y-axis.  All subsequent fields will be binned and their
            profiles added to the underlying `BinnedProfile2D`.
        center : array_like, optional
            The center to be used for things like radius and radial velocity.
            Defaults to the center of the plot collection.
        cmap : string, optional
            An acceptable colormap.  See either yt.visualization.color_maps or
            http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps .
        weight : string, default "CellMassMsun"
            The weighting field for an average.  This defaults to mass-weighted
            averaging.
        accumulation : list of booleans, optional
            If true, from the low-value to the high-value the values in all
            binned fields will be accumulated.  This is useful for instance
            when adding an unweighted CellMassMsun to a radial plot, as it will
            show mass interior to that radius.  The first value is for the
            x-axis, the second value for the y-axis.  Note that accumulation
            will only be along each row or column.
        x_bins : int, optional
            How many bins should there be in the x-axis variable?
        x_log : boolean, optional
            Should the bin edges be log-spaced?
        x_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        y_bins : int, optional
            How many bins should there be in the y-axis variable?
        y_log : boolean, optional
            Should the bin edges be log-spaced?
        y_bounds : tuple of floats, optional
            If specified, the boundary values for the binning.  If unspecified,
            the min/max from the data_source will be used.  (Non-zero min/max
            in case of log-spacing.)
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        fractional : boolean
            If true, the plot will be normalized to the sum of all the binned
            values.

        Returns
        -------
        plot : `yt.visualization.plot_types.PhasePlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.profiles.BinnedProfile2D : This is the object that does the
                                   transformation of raw data into a 1D
                                   profile.

        Examples
        --------

        This will show the mass-distribution in the Density-Temperature plane.

        >>> pc.add_phase_sphere(1.0, 'kpc',
                ["Density", "Temperature", "CellMassMsun"], weight = None)
        """

        if center is None: center = self.c
        r = radius/self.pf[unit]
        data_source = self.pf.hierarchy.sphere(center, r)
        p = self.add_phase_object(data_source, fields, cmap,
                             weight, accumulation,
                             x_bins, x_log, x_bounds,
                             y_bins, y_log, y_bounds,
                             id, axes=axes, figure=figure, fractional=fractional)
        p["Width"] = radius
        p["Unit"] = unit
        p["Axis"] = None
        return p

    def add_scatter_source(self, data_source, fields, id=None,
                     figure = None, axes = None, plot_options = None):
        r"""Given a data source, make a scatter plot from that data source.

        This is a very simple plot: you give it an instance of `YTDataContainer`, two
        field names, and it will plot them on an axis

        Parameters
        ----------
        data_source : `yt.data_objects.api.YTDataContainer`
            This will be the data source from which field values will be
            obtained.
        fields : tuple of strings
            The first of these will be the x-field, and the second the y-field.
        id : int, optional
            If specified, this will be the "semi-unique id" of the resultant
            plot.  This should not be set.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        plot_options : dict
            These options will be given to `matplotlib.axes.Axes.scatter`
        
        Returns
        -------
        plot : `yt.visualization.plot_types.ScatterPlot`
            The plot that has been added to the PlotCollection.

        Notes
        -----
        This is a simpler way of making a phase plot, but note that because
        pixels are deposited in order, the color may be a biased sample.

        Examples
        --------

        >>> reg = pf.h.region([0.1, 0.2, 0.3], [0.0, 0.1, 0.2],
                              [0.2, 0.3, 0.4])
        >>> pc.add_scatter_plot(reg, ["Density", "Temperature"],
        >>>                     plot_options = {'color':'b'})
        """
        if id is None: id = self._get_new_id()
        if plot_options is None: plot_options = {}
        sp = ScatterPlot(data_source, fields, id,
                                   plot_options = plot_options,
                                   figure=figure, axes=axes)
        p = self._add_plot(sp)
        return p

    def add_fixed_resolution_plot(self, frb, field, use_colorbar=True,
                      figure = None, axes = None, fig_size=None):
        r"""Create a fixed resolution image from an existing buffer.

        This accepts a `FixedResolutionBuffer` and will make a plot from that
        buffer.

        Parameters
        ----------
        frb : `yt.visualization.plot_types.FixedResolutionBuffer`
            The buffer from which fields will be pulled.
        field : string
            The initial field to display.
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

        Returns
        -------
        plot : `yt.visualization.plot_types.FixedResolutionPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.extensions.image_writer.write_image : A faster, colorbarless way to
                                                 write out FRBs.

        Examples
        --------

        Here's a simple mechanism for getting the angular momentum of a
        collapsing cloud and generating a cutting plane aligned with the
        angular momentum vector.

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> proj = pf.h.proj("Density", 0)
        >>> frb = FixedResolutionBuffer(proj, (0.2, 0.3, 0.4, 0.5), (512, 512))
        >>> p = pc.add_fixed_resolution_plot(frb, "Density")
        """
        p = self._add_plot(FixedResolutionPlot(frb, field,
                         use_colorbar=use_colorbar, axes=axes, figure=figure,
                         size=fig_size))
        p["Axis"] = "na"
        return p

    def add_ortho_ray(self, axis, coords, field, figure = None,
                      axes = None, field_parameters = None,
                      plot_options = None):
        r"""Create a ray parallel to some axis, from that a line plot, and add
        it to the current collection.

        This function will generate a `yt.data_objects.api.YTOrthoRayBase` from the given
        parameters.  This ray then gets passed to a `yt.visualization.plot_types.LineQueryPLot`, and
        the resultant plot is added to the current collection.  Various
        parameters allow control of the way the line plot is displayed, as well as
        how the ray is generated.

        Parameters
        ----------
        axis : int
            The axis along which to cast the ray.  Can be 0, 1, or 2 or x, y,
            z.
        coords : tuple of floats
            The coordinates to place the ray at.  Note that the axes are in the
            form of x_dict[axis] and y_dict[axis] for some axis.
        field : string
            The initial field to slice and display.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.
        plot_options : dict
            These options will be given to `matplotlib.axes.Axes.plot`

        Returns
        -------
        plot : `yt.visualization.plot_types.LineQueryPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTOrthoRayBase : This is the type created by this function and
                                   passed to the plot created here.

        Examples
        --------

        This will cast a ray from (0.0, 0.5, 0.5) to (1.0, 0.5, 0.5) and plot
        it.

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_ortho_ray(0, (0.5, 0.5), "Density")
        """
        axis = fix_axis(axis)
        if field_parameters is None: field_parameters = {}
        if plot_options is None: plot_options = {}
        data_source = self.pf.h.ortho_ray(axis, coords, field,
                        **field_parameters)
        p = self._add_plot(LineQueryPlot(data_source,
                [axis_names[axis], field], self._get_new_id(),
                figure=figure, axes=axes, plot_options=plot_options))
        return p

    def add_ray(self, start_point, end_point, field, figure = None,
                axes = None, field_parameters = None, plot_options = None):
        r"""Create a ray between two points, from that a line plot, and add
        it to the current collection.

        This function will generate a `yt.data_objects.api.YTRayBase` from the given
        parameters.  This ray then gets passed to a `yt.visualization.plot_types.LineQueryPLot`, and
        the resultant plot is added to the current collection.  Various
        parameters allow control of the way the line plot is displayed, as well as
        how the ray is generated.

        Parameters
        ----------
        start_point : array_like
            The starting point of the ray.
        end_point : array_like
            The ending point of the ray.
        field : string
            The initial field to slice and display.
        figure : `matplotlib.figure.Figure`, optional
            The figure onto which the axes will be placed.  Typically not used
            unless *axes* is also specified.
        axes : `matplotlib.axes.Axes`, optional
            The axes object which will be used to create the image plot.
            Typically used for things like multiplots and the like.
        field_parameters : dict, optional
            This set of parameters will be passed to the slice upon creation,
            which can be used for passing variables to derived fields.
        plot_options : dict
            These options will be given to `matplotlib.axes.Axes.plot`

        Returns
        -------
        plot : `yt.visualization.plot_types.LineQueryPlot`
            The plot that has been added to the PlotCollection.

        See Also
        --------
        yt.data_objects.api.YTRayBase : This is the type created by this function and
                              passed to the plot created here.

        Examples
        --------

        This will cast a ray from (0.1, 0.2, 0.3) to (0.9, 0.7, 0.4) and plot
        it.

        >>> pf = load("RD0005-mine/RedshiftOutput0005")
        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> p = pc.add_ray((0.1, 0.2, 0.3), (0.9, 0.7, 0.4), "Density")
        """
        if field_parameters is None: field_parameters = {}
        if plot_options is None: plot_options = {}
        data_source = self.pf.h.ray(start_point, end_point, field,
                                    **field_parameters)
        p = self._add_plot(LineQueryPlot(data_source,
                ['t', field], self._get_new_id(),
                figure=figure, axes=axes, plot_options=plot_options))
        return p

    def _get_new_id(self):
        self.__id_counter += 1
        return self.__id_counter-1

    @rootonly
    def save_book(self, filename, author = None, title = None, keywords = None,
                  subject = None, creator = None, producer = None,
                  creation_data = None):
        r"""Save a multipage PDF of all the current plots, rather than
        individual image files.

        This function will utilize the matplotlib PDF backend to create a new
        PDF, and for every plot that the PlotCollection currently has, it will
        render a new page into that PDF.  The pages will be in the order of the
        current plots.

        Parameters
        ----------
        filename : string
            The name of the PDF file to generate.  Note that it will be
            overwritten, and '.pdf' will not be appended.
        author : string, optional
            The string to place in the metadata value of the PDF for 'author'.
        title : string, optional
            The string to place in the metadata value of the PDF for 'title'.
        keywords : string, optional
            The string to place in the metadata value of the PDF for 'keywords'.
        subject : string, optional
            The string to place in the metadata value of the PDF for 'subject'.
        creator : string, optional
            The string to place in the metadata value of the PDF for 'creator'.
        producer : string, optional
            The string to place in the metadata value of the PDF for 'producer'.
        creation_date : string, optional
            The string to place in the metadata value of the PDF for
            'creation_date'.

        Returns
        -------
        Nothing

        Examples
        --------
        This will set up a new PlotCollection, add some plots, and then save it
        as a PDF.

        >>> pc = PlotCollection(pf, [0.5, 0.5, 0.5])
        >>> pc.add_projection("Density", 'x')
        >>> pc.add_projection("Density", 'y')
        >>> pc.add_projection("Density", 'z')
        >>> pc.set_width(0.5, 'pc')
        >>> dd = pf.h.all_data()
        >>> pc.add_phase_object(dd, ["Density", "Temperature", "CellMassMsun"],
        ...                     weight = None)
        >>> pc.save_book("my_plots.pdf", author="Matthew Turk", 
        ...              title="Fun plots")
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
        if 'axes' not in kwargs: kwargs['axes'] = self.pylab.gca()
        if 'figure' not in kwargs: kwargs['figure'] = self.pylab.gcf()
        retval = func(self, *args, **kwargs)
        retval._redraw_image()
        retval._fig_num = new_fig.number
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
        r"""Redraw all affiliated plots.

        To ensure that any interactive windows are up to date, this function
        can be called to redraw all images into them.
        """
        for plot in self.plots:
            plot._redraw_image()
        self.pylab.draw()

    def clear_figures(self):
        r"""Clear all interactive figures affiliated with this collection.

        Because reusing figures between plot collections can be tricky,
        occasionally they must be manually cleared to re-obtain empty figures
        for future plotting.  This will clear all figures.
        """
        for plot in self.plots:
            self.pylab.figure(plot._fig_num)
            self.pylab.clf()

    def interactive_zoom(self):
        r"""Enter an interactive zooming session for all plots.

        Use this to enter an interactive session where zoom factors
        can be entered that are used to set the widths of the plot
        collection.  This function has no arguments, but will return
        the final width upon exit.

        Caution: Tested and works with TkAgg and MacOSX backends.  Threaded
        backends like Qt4Agg are likely to fail.

        Controls:
        Any numeric value: Zooms by this factor
        0: Exit interactive zoom session
        -1: Reset to width of 1.0 in code units
        empty: zooms by previously entered factor.

        Returns:
        width: (float) The final width of the plot collection
        """
        print 'Enter Zoom Factor, 0 to exit, -1 to reset to width=1.0'
        zfactor = 1.0
        while(True):
            new_zoom = raw_input('zoom:')
            if new_zoom is not '':
                try:
                    zfactor = float(new_zoom)
                except:
                    print 'Please enter a valid number, or 0 to exit'
                    continue
            else:
                print 'Using previous zoom value of %e' % zfactor
            if zfactor == 0.0:
                break
            elif zfactor == -1.0:
                self.set_width(1.0,'1')
            else:
                self.set_width(self.plots[0].__dict__['width']/zfactor,'1')
        print 'Returning final width of %e' % self.plots[0].width
        return self.plots[0].width

class PlotCollectionIPython(PlotCollection):
    def save(self, basename = None):
        r"""Shows all the plots hanging off this plot collection in the IPython
        web notebook.

        This function will instruct the IPython web notebook to show its
        images.

        Examples
        --------

        >>> pc.save()
        """
        from ._mpl_imports import FigureCanvasAgg
        from IPython.zmq.pylab.backend_inline import \
            send_figure
        if basename is None: basename = str(self.pf)
        for plot in self.plots:
            canvas = FigureCanvasAgg(plot._figure)
            send_figure(plot._figure)

def get_multi_plot(nx, ny, colorbar = 'vertical', bw = 4, dpi=300,
                   cbar_padding = 0.4):
    r"""Construct a multiple axes plot object, with or without a colorbar, into
    which multiple plots may be inserted.

    This will create a set of :class:`matplotlib.axes.Axes`, all lined up into
    a grid, which are then returned to the user and which can be used to plot
    multiple plots on a single figure.

    Parameters
    ----------
    nx : int
        Number of axes to create along the x-direction
    ny : int
        Number of axes to create along the y-direction
    colorbar : {'vertical', 'horizontal', None}, optional
        Should Axes objects for colorbars be allocated, and if so, should they
        correspond to the horizontal or vertical set of axes?
    bw : number
        The base height/width of an axes object inside the figure, in inches
    dpi : number
        The dots per inch fed into the Figure instantiation

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
        The figure created inside which the axes reside
    tr : list of list of :class:`matplotlib.axes.Axes` objects
        This is a list, where the inner list is along the x-axis and the outer
        is along the y-axis
    cbars : list of :class:`matplotlib.axes.Axes` objects
        Each of these is an axes onto which a colorbar can be placed.

    Notes
    -----
    This is a simple implementation for a common use case.  Viewing the source
    can be instructure, and is encouraged to see how to generate more
    complicated or more specific sets of multiplots for your own purposes.
    """
    hf, wf = 1.0/ny, 1.0/nx
    fudge_x = fudge_y = 1.0
    if colorbar is None:
        fudge_x = fudge_y = 1.0
    elif colorbar.lower() == 'vertical':
        fudge_x = nx/(cbar_padding+nx)
        fudge_y = 1.0
    elif colorbar.lower() == 'horizontal':
        fudge_x = 1.0
        fudge_y = ny/(cbar_padding+ny)
    fig = figure.Figure((bw*nx/fudge_x, bw*ny/fudge_y), dpi=dpi)
    from _mpl_imports import FigureCanvasAgg
    fig.set_canvas(FigureCanvasAgg(fig))
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

def _MPLFixImage(data_source, image_obj, field, cbar, cls):
    nx, ny = image_obj.get_size()
    def f(axes):
        x0, x1 = axes.get_xlim()
        y0, y1 = axes.get_ylim()
        frb = cls(data_source, (x0, x1, y0, y1), (nx, ny))
        image_obj.set_data(frb[field])
        mi, ma = frb[field].min(), frb[field].max()
        cbar.norm.autoscale((mi, ma))
        image_obj.set_extent([x0, x1, y0, y1])
        cbar.update_bruteforce(image_obj)
    return f

def matplotlib_widget(data_source, field, npix):
    r"""Create a widget from a data_source that uses the Matplotlib interaction
    method to pan, zoom, and so on.

    This is a simple way to take a yt data source, for instance a projection or
    a slice, and to create a matplotlib view into it that you can pan and zoom.
    It uses the matplotlib interaction engine to manage input and display.

    Parameters
    ----------
    data_source : :class:`yt.data_objects.data_containers.YTOverlapProjBase` or :class:`yt.data_objects.data_containers.YTSliceBase`
        This is the source to be pixelized, which can be a projection or a
        slice.  
    field : string
        The field that you want to display in the window.
    npix : int
        The number of pixels on a side you want the image to be.

    Examples
    --------

    >>> pf = load("DD0030/DD0030")
    >>> p = pf.h.proj("Density", "z")
    >>> matplotlib_widget(p, "Density", 1024)

    """
    import pylab
    import matplotlib.colors
    from .fixed_resolution import FixedResolutionBuffer, \
            ObliqueFixedResolutionBuffer
    pf = data_source.pf
    if getattr(data_source, "axis", 4) < 3:
        cls = FixedResolutionBuffer
        ax = data_source.axis
        extent = [pf.domain_left_edge[x_dict[ax]],
                  pf.domain_right_edge[x_dict[ax]],
                  pf.domain_left_edge[y_dict[ax]],
                  pf.domain_right_edge[y_dict[ax]]]
    else:
        cls = ObliqueFixedResolutionBuffer
        extent = [0.0, 1.0, 0.0, 1.0]
    take_log = pf.field_info[field].take_log
    if take_log:
        norm = matplotlib.colors.LogNorm()
    else:
        norm = matplotlib.colors.Normalize()
    ax = pylab.figure().gca()
    ax.autoscale(False)
    axi = ax.imshow(np.random.random((npix, npix)),
                    extent = extent, norm = norm,
                    origin = 'lower')
    cb = pylab.colorbar(axi, norm = norm)
    showme = _MPLFixImage(data_source, axi, field, cb, cls)
    ax.callbacks.connect("xlim_changed", showme)
    ax.callbacks.connect("ylim_changed", showme)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    return ax
