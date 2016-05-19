"""
This is a simple mechanism for interfacing with Particle plots



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.visualization.fixed_resolution import \
    ParticleImageBuffer
from yt.funcs import \
    ensure_list, \
    fix_axis
from yt.units.yt_array import YTArray
from .plot_window import \
    get_axes_unit, \
    get_window_parameters, \
    PWViewerMPL
from yt.data_objects.profiles import \
    create_profile
from yt.visualization.profile_plotter import \
    PhasePlot


class ParticleAxisAlignedDummyDataSource(object):
    _type_name = 'Particle'
    _dimensionality = 2
    _con_args = ('center', 'axis', 'width', 'fields', 'weight_field')
    _tds_attrs = ()
    _key_fields = []

    def __init__(self, center, ds, axis, width, fields,
                 weight_field=None,
                 field_parameters=None, data_source=None):
        self.center = center
        self.ds = ds
        self.axis = axis
        self.width = width
        self.weight_field = weight_field

        if field_parameters is None:
            self.field_parameters = {}
        else:
            self.field_parameters = field_parameters

        LE = center - 0.5*YTArray(width)
        RE = center + 0.5*YTArray(width)
        self.dd = ds.region(center, LE, RE, fields,
                            field_parameters=field_parameters,
                            data_source=data_source)

        fields = self.dd._determine_fields(fields)
        self.fields = fields

    def _determine_fields(self, *args):
        return self.dd._determine_fields(*args)

    def get_field_parameter(self, name, default=None):
        """
        This is typically only used by derived field functions, but
        it returns parameters used to generate fields.
        """
        if name in self.field_parameters:
            return self.field_parameters[name]
        else:
            return default


class ParticleProjectionPlot(PWViewerMPL):
    r"""Creates a particle plot from a dataset

    Given a ds object, an axis to slice along, and a field name
    string, this will return a PWViewerMPL object containing
    the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : `Dataset`
         This is the dataset object corresponding to the
         simulation output to be plotted.
    axis : int or one of 'x', 'y', 'z'
         An int corresponding to the axis to slice along (0=x, 1=y, 2=z)
         or the axis name itself
    fields : string, list or None
         If a string or list, the name of the particle field(s) to be used 
         one the colorbar. If None, the particle positions will be indicated
         using a fixed color, instead. Default is None.
    color : 'b', 'g', 'r', 'c', 'm', 'y', 'k', or 'w'
         One the matplotlib-recognized color strings.
         The color that will indicate the particle locations
         on the mesh. This argument is ignored if z_fields is
         not None. Default is 'b'.
    center : A sequence of floats, a string, or a tuple.
         The coordinate of the center of the image. If set to 'c', 'center' or
         left blank, the plot is centered on the middle of the domain. If set to
         'max' or 'm', the center will be located at the maximum of the
         ('gas', 'density') field. Centering on the max or min of a specific
         field is supported by providing a tuple such as ("min","temperature") or
         ("max","dark_matter_density"). Units can be specified by passing in *center*
         as a tuple containing a coordinate and string unit name or by passing
         in a YTArray. If a list or unitless array is supplied, code units are
         assumed.
    width : tuple or a float.
         Width can have four different formats to support windows with variable
         x and y widths.  They are:

         ==================================     =======================
         format                                 example
         ==================================     =======================
         (float, string)                        (10,'kpc')
         ((float, string), (float, string))     ((10,'kpc'),(15,'kpc'))
         float                                  0.2
         (float, float)                         (0.2, 0.3)
         ==================================     =======================

         For example, (10, 'kpc') requests a plot window that is 10 kiloparsecs
         wide in the x and y directions, ((10,'kpc'),(15,'kpc')) requests a
         window that is 10 kiloparsecs wide along the x axis and 15
         kiloparsecs wide along the y axis.  In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
         are assumed. Defaults to the entire domain.
    weight_field : string
         The name of the weighting field.  Set to None for no weight.
    axes_unit : A string
         The name of the unit for the tick labels on the x and y axes.
         Defaults to None, which automatically picks an appropriate unit.
         If axes_unit is '1', 'u', or 'unitary', it will not display the
         units, and only show the axes name.
    origin : string or length 1, 2, or 3 sequence of strings
         The location of the origin of the plot coordinate system.  This is
         represented by '-' separated string or a tuple of strings.  In the
         first index the y-location is given by 'lower', 'upper', or 'center'.
         The second index is the x-location, given as 'left', 'right', or
         'center'.  Finally, the whether the origin is applied in 'domain'
         space, plot 'window' space or 'native' simulation coordinate system
         is given. For example, both 'upper-right-domain' and ['upper',
         'right', 'domain'] both place the origin in the upper right hand
         corner of domain space. If x or y are not given, a value is inffered.
         For instance, 'left-domain' corresponds to the lower-left hand corner
         of the simulation domain, 'center-domain' corresponds to the center
         of the simulation domain, or 'center-window' for the center of the
         plot window. Further examples:

         ==================================     ============================
         format                                 example
         ==================================     ============================
         '{space}'                              'domain'
         '{xloc}-{space}'                       'left-window'
         '{yloc}-{space}'                       'upper-domain'
         '{yloc}-{xloc}-{space}'                'lower-right-window'
         ('{space}',)                           ('window',)
         ('{xloc}', '{space}')                  ('right', 'domain')
         ('{yloc}', '{space}')                  ('lower', 'window')
         ('{yloc}', '{xloc}', '{space}')        ('lower', 'right', 'window')
         ==================================     ============================
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    field_parameters : dictionary
         A dictionary of field parameters than can be accessed by derived
         fields.
    data_source : YTSelectionContainer Object
         Object to be used for data selection.  Defaults to a region covering
         the entire simulation.

    Examples
    --------

    This will save an image the the file
    'galaxy0030_Particle_z_particle_mass.png'

    >>> from yt import load
    >>> ds = load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> p = yt.ParticleProjectionPlot(ds, 2, 'particle_mass')
    >>> p.save()

    """
    _plot_type = 'Particle'
    _frb_generator = ParticleImageBuffer

    def __init__(self, ds, axis, fields=None, color='b', center='c', width=None,
                 depth=(1, '1'), weight_field=None, axes_unit=None,
                 origin='center-window', fontsize=18, field_parameters=None,
                 window_size=8.0, aspect=None, data_source=None):
        # this will handle time series data and controllers
        ts = self._initialize_dataset(ds)
        self.ts = ts
        ds = self.ds = ts[0]
        axis = fix_axis(axis, ds)
        (bounds, center, display_center) = \
            get_window_parameters(axis, center, width, ds)
        if field_parameters is None:
            field_parameters = {}

        if axes_unit is None:
            axes_unit = get_axes_unit(width, ds)

        # if no fields are passed in, we simply mark the x and
        # y fields using a given color. Use the 'particle_ones'
        # field to do this. We also turn off the colorbar in
        # this case.
        self._use_cbar = True
        splat_color = None
        if fields is None:
            fields = ['particle_ones']
            weight_field = 'particle_ones'
            self._use_cbar = False
            splat_color = color

        depth = ds.coordinates.sanitize_depth(depth)

        x_coord = ds.coordinates.x_axis[axis]
        y_coord = ds.coordinates.y_axis[axis]

        width = np.zeros_like(center)
        width[x_coord] = bounds[1] - bounds[0]
        width[y_coord] = bounds[3] - bounds[2]
        width[axis] = depth[0].in_units(width[x_coord].units)

        ParticleSource = ParticleAxisAlignedDummyDataSource(center, ds, axis,
                                        width, fields, weight_field,
                                        field_parameters=field_parameters,
                                        data_source=data_source)

        PWViewerMPL.__init__(self, ParticleSource, bounds, origin=origin,
                             fontsize=fontsize, fields=fields,
                             window_size=window_size, aspect=aspect,
                             splat_color=splat_color)

        self.set_axes_unit(axes_unit)

        if self._use_cbar is False:
            self.hide_colorbar()


class ParticlePhasePlot(PhasePlot):
    r"""
    Create a 2d particle phase plot from a data source or from
    a `yt.data_objects.profiles.ParticleProfile` object.

    Given a data object (all_data, region, sphere, etc.), an x field,
    y field, and z field (or fields), this will create a particle plot 
    by depositing the particles onto a two-dimensional mesh, using either
    nearest grid point or cloud-in-cell deposition.

    Parameters
    ----------
    data_source : YTSelectionContainer Object
        The data object to be profiled, such as all_data, region, or 
        sphere.
    x_field : str
        The x field for the mesh.
    y_field : str
        The y field for the mesh.
    z_fields : None, str, or list
        If None, particles will be splatted onto the mesh,
        but no colormap will be used.
        If str or list, the name of the field or fields to
        be displayed on the colorbar.
        Default: None.
    color : 'b', 'g', 'r', 'c', 'm', 'y', 'k', or 'w'
        One the matplotlib-recognized color strings.
        The color that will indicate the particle locations
        on the mesh. This argument is ignored if z_fields is
        not None.
        Default : 'b'
    x_bins : int
        The number of bins in x field for the mesh.
        Default: 800.
    y_bins : int
        The number of bins in y field for the mesh.
        Default: 800.
    weight_field : str
        The field to weight by. Default: None.
    deposition : str
        Either 'ngp' or 'cic'. Controls what type of
        interpolation will be used to deposit the
        particle z_fields onto the mesh.
        Default: 'ngp'
    fontsize: int
        Font size for all text in the plot.
        Default: 18.
    figure_size : int
        Size in inches of the image.
        Default: 8 (8x8)

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> plot = ParticlePhasePlot(ad, "particle_position_x,
                                 "particle_position_y", ["particle_mass"],
    ...                          x_bins=800, y_bins=800)
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_log('particle_mass', True)
    >>> plot.set_unit('particle_position_x', 'Mpc')
    >>> plot.set_unit('particle_velocity_z', 'km/s')
    >>> plot.set_unit('particle_mass', 'Msun')

    """
    _plot_type = 'ParticlePhase'

    def __init__(self, data_source, x_field, y_field, z_fields=None,
                 color='b', x_bins=800, y_bins=800, weight_field=None,
                 deposition='ngp', fontsize=18, figure_size=8.0):

        # if no z_fields are passed in, use a constant color
        if z_fields is None:
            self.use_cbar = False
            self.splat_color = color
            z_fields = ['particle_ones']

        profile = create_profile(
            data_source,
            [x_field, y_field],
            ensure_list(z_fields),
            n_bins=[x_bins, y_bins],
            weight_field=weight_field,
            deposition=deposition)

        type(self)._initialize_instance(self, data_source, profile, fontsize,
                                        figure_size)


def ParticlePlot(ds, x_field, y_field, z_fields=None, color='b', *args, **
                 kwargs):
    r"""
    A factory function for
    :class:`yt.visualization.particle_plots.ParticleProjectionPlot`
    and :class:`yt.visualization.profile_plotter.ParticlePhasePlot` objects.
    This essentially allows for a single entry point to both types of particle
    plots, the distinction being determined by the fields passed in.

    If the x_field and y_field combination corresponds to a valid, right-handed
    spatial plot, an 'ParticleProjectionPlot` will be returned. This plot
    object can be updated using one of the many helper functions defined in
    PlotWindow.

    If the x_field and y_field combo do not correspond to a valid
    'ParticleProjectionPlot`, then a `ParticlePhasePlot`. This object can be
    modified by its own set of  helper functions defined in PhasePlot.

    Parameters
    ----------

    ds : :class:`yt.data_objects.api.Dataset`
        This is the dataset object corresponding to the
        simulation output to be plotted.
    x_field : string
        This is the particle field that will be plotted on the x-axis.
    y_field : string
        This is the particle field that will be plotted on the y-axis.
    z_fields : string, list, or None.
        If None, particles will be splatted onto the plot, but no colormap
        will be used. The particle color will instead be determined by
        the 'color' argument. If str or list, the name of the field or fields
        to be displayed on the colorbar.
        Default: None.
    color : 'b', 'g', 'r', 'c', 'm', 'y', 'k', or 'w'
         One the matplotlib-recognized color strings.
         The color that will indicate the particle locations
         on the plot. This argument is ignored if z_fields is
         not None. Default is 'b'.

    Examples
    --------

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y',
    ...                     'particle_mass', width=(0.5, 0.5))
    >>> p.set_unit('particle_mass', 'Msun')
    >>> p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_velocity_z',
    ...                     color='g')

    """

    ad = ds.all_data()
    x_field = ad._determine_fields(x_field)[0]
    y_field = ad._determine_fields(y_field)[0]

    direction = 3
    # try potential axes for a ParticleProjectionPlot:
    for axis in [0, 1, 2]:
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        ax_field_template = 'particle_position_%s'
        xf = ax_field_template % ds.coordinates.axis_name[xax]
        yf = ax_field_template % ds.coordinates.axis_name[yax]
        if (x_field[1], y_field[1]) == (xf, yf):
            direction = axis
            break

    if direction < 3:
        # Make a ParticleProjectionPlot
        return ParticleProjectionPlot(ds, direction, z_fields, color,
                                      *args, **kwargs)

    # Does not correspond to any valid PlotWindow-style plot,
    # use ParticlePhasePlot instead
    else:
        return ParticlePhasePlot(ad, x_field, y_field,
                                 z_fields, color, *args, **kwargs)
