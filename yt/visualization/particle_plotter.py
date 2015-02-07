"""
This is a simple mechanism for interfacing with Particle Scatter plots



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import os
import matplotlib
import numpy as np

from yt.data_objects.image_array import ImageArray
from .base_plot_types import ImagePlotMPL
from yt.units.yt_array import array_like_field
from .plot_container import \
    ImagePlotContainer, \
    log_transform, linear_transform, get_log_minorticks
from yt.utilities.logger import ytLogger as mylog
from yt.funcs import \
    ensure_list, \
    get_image_suffix
from .profile_plotter import \
    invalidate_plot
from yt.utilities.lib.api import add_points_to_greyscale_image


class ParticleSplatter(object):
    x_log = False
    y_log = False

    def __init__(self, data_source, x_field, y_field,
                 bounds, z_fields=None,
                 x_bins=256, y_bins=256):
        self.data_source = data_source
        self.ds = data_source.ds
        self.x_field = data_source._determine_fields(x_field)[0]
        self.y_field = data_source._determine_fields(y_field)[0]
        self.x_bins = x_bins
        self.y_bins = y_bins
        self.bounds = bounds
        self.data = {}
        self.x = array_like_field(data_source,
                                  np.linspace(bounds[0], bounds[1], x_bins),
                                  x_field)
        self.y = array_like_field(data_source,
                                  np.linspace(bounds[2], bounds[3], y_bins),
                                  y_field)

        if z_fields is None:
            self._splat_particle_field('particle_ones')
        else:
            for f in z_fields:
                self._splat_particle_field(f)

    def keys(self):
        return self.data.keys()

    def __delitem__(self, item):
        del self.data[item]

    def __getitem__(self, item):
        return self.data[item]

    def _splat_particle_field(self, item):

        mylog.info("Splatting (%s) onto a %d by %d mesh in %s by %s space" %
                (item, self.x_bins, self.y_bins, self.x_field, self.y_field))

        bounds = []
        for b in self.bounds:
            if hasattr(b, "in_units"):
                b = float(b.in_units("code_length"))
            bounds.append(b)

        x_data = self.data_source[self.x_field]
        y_data = self.data_source[self.y_field]
        data = self.data_source[item]

        px = (x_data.d - self.bounds[0]) / (self.bounds[1] - self.bounds[0])
        py = (y_data.d - self.bounds[2]) / (self.bounds[3] - self.bounds[2])

        locs1 = np.logical_and(px > 0.0, px < 1.0)
        locs2 = np.logical_and(py > 0.0, py < 1.0)
        locs = np.logical_and(locs1, locs2)
        locs = np.where(locs)

        buff = np.zeros((self.x_bins, self.y_bins))
        add_points_to_greyscale_image(buff, px[locs], py[locs], data[locs])
        ia = ImageArray(buff, input_units=data.units)

        self.data[item] = ia
        return self.data[item]

    def __setitem__(self, item, val):
        self.data[item] = val

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def __iter__(self):
        return sorted(self.items())


class ParticlePlot(ImagePlotContainer):
    r"""
    Create a 2d profile (phase) plot from a data source or from 
    profile object created with 
    `yt.data_objects.profiles.create_profile`.

    Given a data object (all_data, region, sphere, etc.), an x field,
    y field, and z field (or fields), this will create a two-dimensional 
    profile of the average (or total) value of the z field in bins of the 
    x and y fields.

    Parameters
    ----------
    data_source : YTSelectionContainer Object
        The data object to be profiled, such as all_data, region, or 
        sphere.
    x_field : str
        The x binning field for the profile.
    y_field : str
        The y binning field for the profile.
    z_fields : str or list
        The field or fields to be profiled.
    weight_field : str
        The weight field for calculating weighted averages.  If None, 
        the profile values are the sum of the field values within the bin.
        Otherwise, the values are a weighted average.
        Default : "cell_mass".
    x_bins : int
        The number of bins in x field for the profile.
        Default: 128.
    y_bins : int
        The number of bins in y field for the profile.
        Default: 128.
    fontsize: int
        Font size for all text in the plot.
        Default: 18.
    figure_size : int
        Size in inches of the image.
        Default: 8 (8x8)

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = ds.all_data()
    >>> plot = PhasePlot(ad, "density", "temperature", ["cell_mass"],
    ...                  weight_field=None)
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_cmap("cell_mass", "jet")
    >>> plot.set_zlim("cell_mass", 1e8, 1e13)
    >>> plot.set_title("cell_mass", "This is a phase plot")

    """
    x_log = None
    y_log = None
    plot_title = None
    _plot_valid = False
    _plot_type = 'Particle'

    def __init__(self, data_source, x_field, y_field, bounds,
                 z_fields=None, x_bins=256, y_bins=256,
                 fontsize=18, figure_size=8.0):

        self.x_field = data_source._determine_fields(x_field)[0]
        self.y_field = data_source._determine_fields(y_field)[0]
        if z_fields is None:
            self._draw_colorbars = False
        else:
            self._draw_colorbars = True

        splatter = ParticleSplatter(data_source,
                                    x_field,
                                    y_field,
                                    bounds,
                                    z_fields,
                                    x_bins,
                                    y_bins)

        type(self)._initialize_instance(self,
                                        data_source,
                                        splatter,
                                        fontsize,
                                        figure_size)

    @classmethod
    def _initialize_instance(cls, obj, data_source, splatter, fontsize,
                             figure_size):
        obj.plot_title = {}
        obj.z_log = {}
        obj.z_title = {}
        obj._initfinished = False
        obj.x_log = None
        obj.y_log = None
        obj._plot_text = {}
        obj._text_xpos = {}
        obj._text_ypos = {}
        obj._text_kwargs = {}
        obj.splatter = splatter
        super(ParticlePlot, obj).__init__(data_source, figure_size, fontsize)
        obj._setup_plots()
        obj._initfinished = True
        return obj

    def _get_field_title(self, field_z, profile):
        ds = profile.data_source.ds
        field_x = profile.x_field
        field_y = profile.y_field
        xf, yf, zf = profile.data_source._determine_fields(
            [field_x, field_y, field_z])
        xfi = ds._get_field_info(*xf)
        yfi = ds._get_field_info(*yf)
        zfi = ds._get_field_info(*zf)
        x_unit = profile.x.units
        y_unit = profile.y.units
        z_unit = profile.data[field_z].units
        x_label, y_label, z_label = self._get_axes_labels(field_z)
        x_title = x_label or self._get_field_label(field_x, xfi, x_unit)
        y_title = y_label or self._get_field_label(field_y, yfi, y_unit)
        z_title = z_label or self._get_field_label(field_z, zfi, z_unit)
        return (x_title, y_title, z_title)

    def _get_field_label(self, field, field_info, field_unit):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple): field = field[1]
        if field_name is None:
            field_name = r'$\rm{'+field+r'}$'
            field_name = r'$\rm{'+field.replace('_','\ ').title()+r'}$'
            label = field_name
        elif field_name.find('$') == -1:
            field_name = field_name.replace(' ','\ ')
            field_name = r'$\rm{'+field_name+r'}$'
            label = field_name
        elif field_unit is None or field_unit is '':
            label = field_name
        else:
            label = field_name+r'$\ \ ('+field_unit+r')$'
        return label

    def _get_field_log(self, field_z, profile):
        ds = profile.data_source.ds
        zf, = profile.data_source._determine_fields([field_z])
        zfi = ds._get_field_info(*zf)
        if self.x_log is None:
            x_log = profile.x_log
        else:
            x_log = self.x_log
        if self.y_log is None:
            y_log = profile.y_log
        else:
            y_log = self.y_log
        if field_z in self.z_log:
            z_log = self.z_log[field_z]
        else:
            z_log = zfi.take_log
        scales = {True: 'log', False: 'linear'}
        return scales[x_log], scales[y_log], scales[z_log]

    def _setup_plots(self):
        for f, data in self.splatter.items():
            fig = None
            axes = None
            cax = None
            if self._draw_colorbars is False:
                draw_colorbar = False
            else:
                draw_colorbar = True
            draw_axes = True
            zlim = (None, None)
            if f in self.plots:
                draw_colorbar = self.plots[f]._draw_colorbar
                draw_axes = self.plots[f]._draw_axes
                zlim = (self.plots[f].zmin, self.plots[f].zmax)
                if self.plots[f].figure is not None:
                    fig = self.plots[f].figure
                    axes = self.plots[f].axes
                    cax = self.plots[f].cax

            x_scale, y_scale, z_scale = self._get_field_log(f, self.splatter)
            x_title, y_title, z_title = self._get_field_title(f, self.splatter)

            f = self.splatter.data_source._determine_fields(f)[0]

            if zlim == (None, None):
                if z_scale == 'log':
                    positive_values = data[data > 0.0]
                    if len(positive_values) == 0:
                        mylog.warning("Profiled field %s has no positive "
                                      "values.  Max = %f." %
                                      (f, np.nanmax(data)))
                        mylog.warning("Switching to linear colorbar scaling.")
                        zmin = np.nanmin(data)
                        z_scale = 'linear'
                        self._field_transform[f] = linear_transform
                    else:
                        zmin = positive_values.min()
                        self._field_transform[f] = log_transform
                else:
                    zmin = np.nanmin(data)
                    self._field_transform[f] = linear_transform
                zlim = [zmin, np.nanmax(data)]

            font_size = self._font_properties.get_size()

            self.plots[f] = ParticlePlotMPL(data, self.splatter.bounds,
                                            x_scale, y_scale, z_scale,
                                            self._colormaps[f], zlim,
                                            self.figure_size, font_size,
                                            fig, axes, cax)

            self.plots[f]._toggle_axes(draw_axes)
            self.plots[f]._toggle_colorbar(draw_colorbar)

            self.plots[f].axes.xaxis.set_label_text(x_title)
            self.plots[f].axes.yaxis.set_label_text(y_title)
            self.plots[f].cax.yaxis.set_label_text(z_title)

            if f in self._plot_text:
                self.plots[f].axes.text(self._text_xpos[f], self._text_ypos[f],
                                        self._plot_text[f],
                                        fontproperties=self._font_properties,
                                        **self._text_kwargs[f])

            if f in self.plot_title:
                self.plots[f].axes.set_title(self.plot_title[f])

            # x-y axes minorticks
            if f not in self._minorticks:
                self._minorticks[f] = True
            if self._minorticks[f] is True:
                self.plots[f].axes.minorticks_on()
            else:
                self.plots[f].axes.minorticks_off()

            # colorbar minorticks
            if f not in self._cbar_minorticks:
                self._cbar_minorticks[f] = True
            if self._cbar_minorticks[f] is True:
                if self._field_transform[f] == linear_transform:
                    self.plots[f].cax.minorticks_on()
                else:
                    vmin = np.float64( self.plots[f].cb.norm.vmin )
                    vmax = np.float64( self.plots[f].cb.norm.vmax )
                    mticks = self.plots[f].image.norm( get_log_minorticks(vmin, vmax) )
                    self.plots[f].cax.yaxis.set_ticks(mticks, minor=True)
            else:
                self.plots[f].cax.minorticks_off()

        self._set_font_properties()

        self._plot_valid = True

    def annotate_text(self, xpos=0.0, ypos=0.0, text=None, **text_kwargs):
        r"""
        Allow the user to insert text onto the plot
        The x-position and y-position must be given as well as the text string. 
        Add *text* tp plot at location *xpos*, *ypos* in plot coordinates
        (see example below).

        Parameters
        ----------
        field: str or tuple
          The name of the field to add text to. 
        xpos: float
          Position on plot in x-coordinates.
        ypos: float
          Position on plot in y-coordinates.
        text: str
          The text to insert onto the plot.
        text_kwargs: dict
          Dictionary of text keyword arguments to be passed to matplotlib

        >>>  plot.annotate_text(1e-15, 5e4, "Hello YT")

        """
        for f in self.data_source._determine_fields(self.plots.keys()):
            if self.plots[f].figure is not None and text is not None:
                self.plots[f].axes.text(xpos, ypos, text,
                                        fontproperties=self._font_properties,
                                        **text_kwargs)
            self._plot_text[f] = text
            self._text_xpos[f] = xpos
            self._text_ypos[f] = ypos
            self._text_kwargs[f] = text_kwargs
        return self

    def save(self, name=None, mpl_kwargs=None):
        r"""
        Saves a 2d profile plot.

        Parameters
        ----------
        name : str
            The output file keyword.
        mpl_kwargs : dict
           A dict of keyword arguments to be passed to matplotlib.

        >>> plot.save(mpl_kwargs={'bbox_inches':'tight'})

        """
        names = []
        if not self._plot_valid:
            self._setup_plots()
        if mpl_kwargs is None:
            mpl_kwargs = {}
        if name is None:
            name = str(self.splatter.ds)
        name = os.path.expanduser(name)
        xfn = self.x_field
        yfn = self.y_field
        if isinstance(xfn, tuple):
            xfn = xfn[1]
        if isinstance(yfn, tuple):
            yfn = yfn[1]
        for f in self.splatter.field_data:
            _f = f
            if isinstance(f, tuple):
                _f = _f[1]
            middle = "2d-Profile_%s_%s_%s" % (xfn, yfn, _f)
            splitname = os.path.split(name)
            if splitname[0] != '' and not os.path.isdir(splitname[0]):
                os.makedirs(splitname[0])
            if os.path.isdir(name) and name != str(self.splatter.ds):
                prefix = name + (os.sep if name[-1] != os.sep else '')
                prefix += str(self.splatter.ds)
            else:
                prefix = name
            suffix = get_image_suffix(name)
            if suffix != '':
                for k, v in self.plots.iteritems():
                    names.append(v.save(name, mpl_kwargs))
                return names
            fn = "%s_%s%s" % (prefix, middle, '.png')
            names.append(fn)
            self.plots[f].save(fn, mpl_kwargs)
        return names

    @invalidate_plot
    def set_title(self, field, title):
        """Set a title for the plot.

        Parameters
        ----------
        field : str
            The z field of the plot to add the title.
        title : str
            The title to add.

        Examples
        --------

        >>> plot.set_title("cell_mass", "This is a phase plot")

        """
        self.plot_title[self.data_source._determine_fields(field)[0]] = title
        return self

    @invalidate_plot
    def reset_plot(self):
        self.plots = {}
        return self

    @invalidate_plot
    def set_log(self, field, log):
        """set a field to log or linear.

        Parameters
        ----------
        field : string
            the field to set a transform
        log : boolean
            Log on/off.
        """
        if field == "all":
            self.x_log = log
            self.y_log = log
            for field in self.splatter.field_data:
                self.z_log[field] = log
        else:
            if field == self.x_field[1]:
                self.x_log = log
            elif field == self.y_field[1]:
                self.y_log = log
            elif field in self.splatter.keys():
                self.z_log[field] = log
            else:
                raise KeyError("Field %s not in phase plot!" % (field))
        return self

    @invalidate_plot
    def set_unit(self, field, unit):
        """Sets a new unit for the requested field

        Parameters
        ----------
        field : string
           The name of the field that is to be changed.

        new_unit : string or Unit object
           The name of the new unit.
        """
        fields = [fd[1] for fd in self.splatter.data]
        if field == self.splatter.x_field[1]:
            self.splatter.set_x_unit(unit)
        elif field == self.splatter.y_field[1]:
            self.profile.set_y_unit(unit)
        elif field in fields:
            self.profile.set_field_unit(field, unit)
            self.plots[field].zmin, self.plots[field].zmax = (None, None)
        else:
            raise KeyError("Field %s not in phase plot!" % (field))
        return self

    @invalidate_plot
    def set_xlim(self, xmin=None, xmax=None):
        """Sets the limits of the x bin field

        Parameters
        ----------

        xmin : float or None
          The new x minimum.  Defaults to None, which leaves the xmin
          unchanged.

        xmax : float or None
          The new x maximum.  Defaults to None, which leaves the xmax
          unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> pp = yt.PhasePlot(ds.all_data(), 'density', 'temperature', 'cell_mass')
        >>> pp.set_xlim(1e-29, 1e-24)
        >>> pp.save()

        """
        p = self.profile
        if xmin is None:
            xmin = p.x_bins.min()
        if xmax is None:
            xmax = p.x_bins.max()
        units = {p.x_field: str(p.x.units),
                 p.y_field: str(p.y.units)}
        zunits = dict((field, str(p.field_units[field])) for field in p.field_units)
        extrema = {p.x_field: ((xmin, str(p.x.units)), (xmax, str(p.x.units))),
                   p.y_field: ((p.y_bins.min(), str(p.y.units)),
                               (p.y_bins.max(), str(p.y.units)))}
        self.profile = create_profile(
            p.data_source,
            [p.x_field, p.y_field],
            p.field_map.values(),
            n_bins=[len(p.x_bins)-2, len(p.y_bins)-2],
            weight_field=p.weight_field,
            units=units,
            extrema=extrema)
        for field in zunits:
            self.profile.set_field_unit(field, zunits[field])
        return self

    @invalidate_plot
    def set_ylim(self, ymin=None, ymax=None):
        """Sets the plot limits for the y bin field.

        Parameters
        ----------

        ymin : float or None
          The new y minimum.  Defaults to None, which leaves the ymin
          unchanged.

        ymax : float or None
          The new y maximum.  Defaults to None, which leaves the ymax
          unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load('IsolatedGalaxy/galaxy0030/galaxy0030')
        >>> pp = yt.PhasePlot(ds.all_data(), 'density', 'temperature', 'cell_mass')
        >>> pp.set_ylim(1e4, 1e6)
        >>> pp.save()

        """
        p = self.profile
        if ymin is None:
            ymin = p.y_bins.min()
        if ymax is None:
            ymax = p.y_bins.max()
        units = {p.x_field: str(p.x.units),
                 p.y_field: str(p.y.units)}
        zunits = dict((field, str(p.field_units[field])) for field in p.field_units)
        extrema = {p.x_field: ((p.x_bins.min(), str(p.x.units)),
                               (p.x_bins.max(), str(p.x.units))),
                   p.y_field: ((ymin, str(p.y.units)), (ymax, str(p.y.units)))}
        self.profile = create_profile(
            p.data_source,
            [p.x_field, p.y_field],
            p.field_map.values(),
            n_bins=[len(p.x_bins), len(p.y_bins)],
            weight_field=p.weight_field,
            units=units,
            extrema=extrema)
        for field in zunits:
            self.profile.set_field_unit(field, zunits[field])
        return self

    def run_callbacks(self, *args):
        raise NotImplementedError
    def setup_callbacks(self, *args):
        raise NotImplementedError


class ParticlePlotMPL(ImagePlotMPL):
    """A container for a single matplotlib figure and axes for a ParticlePlot"""
    def __init__(self, data, bounds,
                 x_scale, y_scale, z_scale, cmap,
                 zlim, figure_size, fontsize, figure, axes, cax):
        self._initfinished = False
        self._draw_colorbar = True
        self._draw_axes = True
        self._figure_size = figure_size

        # Compute layout
        fontscale = float(fontsize) / 18.0
        if fontscale < 1.0:
            fontscale = np.sqrt(fontscale)

        self._cb_size = 0.0375*figure_size
        self._ax_text_size = [1.1*fontscale, 0.9*fontscale]
        self._top_buff_size = 0.30*fontscale
        self._aspect = 1.0

        size, axrect, caxrect = self._get_best_layout()

        super(ParticlePlotMPL, self).__init__(size, axrect, caxrect, zlim,
                                              figure, axes, cax)

        self._init_image(data, bounds, x_scale, y_scale, z_scale,
                         zlim, cmap)

        self._initfinished = True

    def _init_image(self, image_data, bounds,
                    x_scale, y_scale, z_scale, zlim, cmap):
        """Store output of imshow in image variable"""
        if (z_scale == 'log'):
            norm = matplotlib.colors.LogNorm(zlim[0], zlim[1])
        elif (z_scale == 'linear'):
            norm = matplotlib.colors.Normalize(zlim[0], zlim[1])
        self.image = None
        self.cb = None
        self.image = self.axes.imshow(np.array(image_data),
                                      extent=bounds,
                                      aspect='auto',
                                      norm=norm,
                                      cmap=cmap,
                                      origin='lower')
        self.axes.set_xscale(x_scale)
        self.axes.set_yscale(y_scale)
        self.cb = self.figure.colorbar(self.image, self.cax)
        if z_scale == 'linear':
            self.cb.formatter.set_scientific(True)
            self.cb.formatter.set_powerlimits((-2,3))
            self.cb.update_ticks()
