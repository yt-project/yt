import numpy as np

from yt._maintenance.deprecation import issue_deprecation_warning
from yt.data_objects.profiles import create_profile
from yt.data_objects.static_output import Dataset
from yt.funcs import fix_axis, iter_fields, mylog
from yt.units.yt_array import YTArray
from yt.utilities.orientation import Orientation
from yt.visualization.fixed_resolution import ParticleImageBuffer
from yt.visualization.profile_plotter import PhasePlot

from .plot_window import (
    NormalPlot,
    PWViewerMPL,
    get_axes_unit,
    get_oblique_window_parameters,
    get_window_parameters,
)


class ParticleDummyDataSource:
    _type_name = "Particle"
    _dimensionality = 2
    _con_args = ("center", "axis", "width", "fields", "weight_field")
    _tds_attrs = ()
    _key_fields: list[str] = []

    def __init__(
        self,
        center,
        ds,
        width,
        fields,
        dd,
        *,
        weight_field=None,
        field_parameters=None,
        deposition="ngp",
        density=False,
    ):
        self.center = center
        self.ds = ds
        self.width = width
        self.dd = dd

        if weight_field is not None:
            weight_field = self._determine_fields(weight_field)[0]
        self.weight_field = weight_field

        self.deposition = deposition
        self.density = density

        if field_parameters is None:
            self.field_parameters = {}
        else:
            self.field_parameters = field_parameters

        fields = self._determine_fields(fields)
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


class ParticleAxisAlignedDummyDataSource(ParticleDummyDataSource):
    def __init__(
        self,
        center,
        ds,
        axis,
        width,
        fields,
        *,
        weight_field=None,
        field_parameters=None,
        data_source=None,
        deposition="ngp",
        density=False,
    ):
        self.axis = axis

        LE = center - 0.5 * YTArray(width)
        RE = center + 0.5 * YTArray(width)
        for ax in range(3):
            if not ds.periodicity[ax]:
                LE[ax] = max(LE[ax], ds.domain_left_edge[ax])
                RE[ax] = min(RE[ax], ds.domain_right_edge[ax])

        dd = ds.region(
            center,
            LE,
            RE,
            fields,
            field_parameters=field_parameters,
            data_source=data_source,
        )

        super().__init__(
            center,
            ds,
            width,
            fields,
            dd,
            weight_field=weight_field,
            field_parameters=field_parameters,
            deposition=deposition,
            density=density,
        )


class ParticleOffAxisDummyDataSource(ParticleDummyDataSource):
    def __init__(
        self,
        center,
        ds,
        normal_vector,
        width,
        fields,
        *,
        weight_field=None,
        field_parameters=None,
        data_source=None,
        deposition="ngp",
        density=False,
        north_vector=None,
    ):
        self.axis = None  # always true for oblique data objects
        normal = np.array(normal_vector)
        normal = normal / np.linalg.norm(normal)

        # If north_vector is None, we set the default here.
        # This is chosen so that if normal_vector is one of the
        # cartesian coordinate axes, the projection will match
        # the corresponding on-axis projection.
        if north_vector is None:
            vecs = np.identity(3)
            t = np.cross(vecs, normal).sum(axis=1)
            ax = t.argmax()
            east_vector = np.cross(vecs[ax, :], normal).ravel()
            north = np.cross(normal, east_vector).ravel()
        else:
            north = np.array(north_vector)
            north = north / np.linalg.norm(north)
        self.normal_vector = normal
        self.north_vector = north

        if data_source is None:
            dd = ds.all_data()
        else:
            dd = data_source

        self.orienter = Orientation(normal_vector, north_vector=north_vector)

        super().__init__(
            center,
            ds,
            width,
            fields,
            dd,
            weight_field=weight_field,
            field_parameters=field_parameters,
            deposition=deposition,
            density=density,
        )


class ParticleProjectionPlot(PWViewerMPL, NormalPlot):
    r"""Creates a particle plot from a dataset

    Given a ds object, a normal to project along, and a field name
    string, this will return a PWViewerMPL object containing
    the plot.

    The plot can be updated using one of the many helper functions
    defined in PlotWindow.

    Parameters
    ----------
    ds : `Dataset`
         This is the dataset object corresponding to the
         simulation output to be plotted.
    normal : int, str, or 3-element sequence of floats
        This specifies the normal vector to the projection.
        Valid int values are 0, 1 and 2. Corresponding str values depend on the
        geometry of the dataset and are generally given by `ds.coordinates.axis_order`.
        E.g. in cartesian they are 'x', 'y' and 'z'.
        An arbitrary normal vector may be specified as a 3-element sequence of floats.
    fields : string, list or None
         If a string or list, the name of the particle field(s) to be used
         on the colorbar. The color shown will correspond to the sum of the
         given field along the line of sight. If None, the particle positions
         will be indicated using a fixed color, instead. Default is None.
    color : 'b', 'g', 'r', 'c', 'm', 'y', 'k', or 'w'
         One the matplotlib-recognized color strings.
         The color that will indicate the particle locations
         on the mesh. This argument is ignored if z_fields is
         not None. Default is 'b'.
    center : 'center', 'c', 'left', 'l', 'right', 'r', id of a global extremum, or array-like
        The coordinate of the selection's center.
        Defaults to the 'center', i.e. center of the domain.

        Centering on the min or max of a field is supported by passing a tuple
        such as ('min', ('gas', 'density')) or ('max', ('gas', 'temperature'). A
        single string may also be used (e.g. "min_density" or
        "max_temperature"), though it's not as flexible and does not allow to
        select an exact field/particle type. With this syntax, the first field
        matching the provided name is selected.
        'max' or 'm' can be used as a shortcut for ('max', ('gas', 'density'))
        'min' can be used as a shortcut for ('min', ('gas', 'density'))

        One can also select an exact point as a 3 element coordinate sequence,
        e.g. [0.5, 0.5, 0]
        Units can be specified by passing in *center* as a tuple containing a
        3-element coordinate sequence and string unit name, e.g. ([0, 0.5, 0.5], "cm"),
        or by passing in a YTArray. Code units are assumed if unspecified.

        The domain edges along the selected *axis* can be selected with
        'left'/'l' and 'right'/'r' respectively.
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
         window that is 10 kiloparsecs wide along the x-axis and 15
         kiloparsecs wide along the y-axis. In the other two examples, code
         units are assumed, for example (0.2, 0.3) requests a plot that has an
         x width of 0.2 and a y width of 0.3 in code units.  If units are
         provided the resulting plot axis labels will use the supplied units.
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
         are assumed. Defaults to the entire domain.
    weight_field : string
         The name of the weighting field.  Set to None for no weight.
         If given, the plot will show a weighted average along the line of
         sight of the fields given in the ``fields`` argument.
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
         'center'.  Finally, whether the origin is applied in 'domain'
         space, plot 'window' space or 'native' simulation coordinate system
         is given. For example, both 'upper-right-domain' and ['upper',
         'right', 'domain'] both place the origin in the upper right hand
         corner of domain space. If x or y are not given, a value is inferred.
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
    window_size : float
        The size of the window on the longest axis (in units of inches),
        including the margins but not the colorbar.
    aspect : float
         The aspect ratio of the plot.  Set to None for 1.
    data_source : YTSelectionContainer object
         The object to be used for data selection.  Defaults to a region covering
         the entire simulation.
    deposition : string
        Controls the order of the interpolation of the particles onto the
        mesh. "ngp" is 0th-order "nearest-grid-point" method (the default),
        "cic" is 1st-order "cloud-in-cell".
    density : boolean
        If True, the quantity to be projected will be divided by the area of
        the cells, to make a projected density of the quantity. The plot
        name and units will also reflect this. Default: False
    north_vector : a sequence of floats
        A vector defining the 'up' direction in off-axis particle projection plots;
        not used if the plot is on-axis. This option sets the orientation of the
        projected plane.  If not set, an arbitrary grid-aligned north-vector is
        chosen.

    Examples
    --------

    This will save an image to the file
    'galaxy0030_Particle_z_particle_mass.png'

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> p = yt.ParticleProjectionPlot(ds, 2, "particle_mass")
    >>> p.save()

    """
    _plot_type = "Particle"
    _frb_generator = ParticleImageBuffer

    def __init__(
        self,
        ds,
        normal=None,
        fields=None,
        color="b",
        center="center",
        width=None,
        depth=(1, "1"),
        weight_field=None,
        axes_unit=None,
        origin="center-window",
        fontsize=18,
        field_parameters=None,
        window_size=8.0,
        aspect=None,
        data_source=None,
        deposition="ngp",
        density=False,
        *,
        north_vector=None,
        axis=None,
    ):
        if axis is not None:
            issue_deprecation_warning(
                "The 'axis' argument is a deprecated alias for the 'normal' argument. ",
                stacklevel=3,
                since="4.2",
            )
            normal = axis
        if normal is None:
            raise TypeError(
                "ParticleProjectionPlot() missing 1 required positional argument: 'normal'"
            )
        # this will handle time series data and controllers
        ts = self._initialize_dataset(ds)
        self.ts = ts
        ds = self.ds = ts[0]
        normal = self.sanitize_normal_vector(ds, normal)
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
            fields = [("all", "particle_ones")]
            weight_field = ("all", "particle_ones")
            self._use_cbar = False
            splat_color = color

        if isinstance(normal, str):
            axis = fix_axis(normal, ds)
            (bounds, center, display_center) = get_window_parameters(
                axis, center, width, ds
            )
            x_coord = ds.coordinates.x_axis[axis]
            y_coord = ds.coordinates.y_axis[axis]

            depth = ds.coordinates.sanitize_depth(depth)

            width = np.zeros_like(center)
            width[x_coord] = bounds[1] - bounds[0]
            width[y_coord] = bounds[3] - bounds[2]
            width[axis] = depth[0].in_units(width[x_coord].units)

            ParticleSource = ParticleAxisAlignedDummyDataSource(
                center,
                ds,
                axis,
                width,
                fields,
                weight_field=weight_field,
                field_parameters=field_parameters,
                data_source=data_source,
                deposition=deposition,
                density=density,
            )

            oblique = False
            plt_origin = origin
            periodic = True

        else:
            (bounds, center_rot) = get_oblique_window_parameters(
                normal, center, width, ds, depth=depth
            )

            width = ds.coordinates.sanitize_width(normal, width, depth)

            ParticleSource = ParticleOffAxisDummyDataSource(
                center_rot,
                ds,
                normal,
                width,
                fields,
                weight_field=weight_field,
                field_parameters=field_parameters,
                data_source=None,
                deposition=deposition,
                density=density,
                north_vector=north_vector,
            )

            oblique = True
            periodic = False
            if origin != "center-window":
                mylog.warning(
                    "The 'origin' keyword is ignored for off-axis "
                    "particle projections, it is always 'center-window'"
                )
            plt_origin = "center-window"

        self.projected = weight_field is None

        PWViewerMPL.__init__(
            self,
            ParticleSource,
            bounds,
            origin=plt_origin,
            fontsize=fontsize,
            fields=fields,
            window_size=window_size,
            aspect=aspect,
            splat_color=splat_color,
            geometry=ds.geometry,
            periodic=periodic,
            oblique=oblique,
        )

        self.set_axes_unit(axes_unit)

        if not self._use_cbar:
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
    data_source : YTSelectionContainer or Dataset
        The data object to be profiled, such as all_data, region, or
        sphere. If data_source is a Dataset, data_source.all_data() will be used.
    x_field : str
        The x field for the mesh.
    y_field : str
        The y field for the mesh.
    z_fields : None, str, or list
        If None, particles will be splatted onto the mesh,
        but no colormap will be used.
        If str or list, the name of the field or fields to
        be displayed on the colorbar. The displayed values will
        correspond to the sum of the field or fields along the
        line of sight.
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
        The field to weight by. If given, the plot will show a weighted
        average along the line of sight of the fields given in the
        ``z_fields`` argument. Default: None.
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
    shading : str
        This argument is directly passed down to matplotlib.axes.Axes.pcolormesh
        see
        https://matplotlib.org/3.3.1/gallery/images_contours_and_fields/pcolormesh_grids.html#sphx-glr-gallery-images-contours-and-fields-pcolormesh-grids-py  # noqa
        Default: 'nearest'

    Examples
    --------

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> ad = ds.all_data()
    >>> plot = ParticlePhasePlot(
    ...     ad,
    ...     "particle_position_x",
    ...     "particle_position_y",
    ...     ["particle_mass"],
    ...     x_bins=800,
    ...     y_bins=800,
    ... )
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_log("particle_mass", True)
    >>> plot.set_unit("particle_position_x", "Mpc")
    >>> plot.set_unit("particle_velocity_z", "km/s")
    >>> plot.set_unit("particle_mass", "Msun")

    """
    _plot_type = "ParticlePhase"

    def __init__(
        self,
        data_source,
        x_field,
        y_field,
        z_fields=None,
        color="b",
        x_bins=800,
        y_bins=800,
        weight_field=None,
        deposition="ngp",
        fontsize=18,
        figure_size=8.0,
        shading="nearest",
    ):
        if isinstance(data_source, Dataset):
            data_source = data_source.all_data()
        # if no z_fields are passed in, use a constant color
        if z_fields is None:
            self.use_cbar = False
            self.splat_color = color
            z_fields = [("all", "particle_ones")]

        profile = create_profile(
            data_source,
            [x_field, y_field],
            list(iter_fields(z_fields)),
            n_bins=[x_bins, y_bins],
            weight_field=weight_field,
            deposition=deposition,
        )

        type(self)._initialize_instance(
            self, data_source, profile, fontsize, figure_size, shading
        )


def ParticlePlot(ds, x_field, y_field, z_fields=None, color="b", *args, **kwargs):
    r"""
    A factory function for
    :class:`yt.visualization.particle_plots.ParticleProjectionPlot`
    and :class:`yt.visualization.profile_plotter.ParticlePhasePlot` objects.
    This essentially allows for a single entry point to both types of particle
    plots, the distinction being determined by the fields passed in.

    If the x_field and y_field combination corresponds to a valid, right-handed
    spatial plot, an ``ParticleProjectionPlot`` will be returned. This plot
    object can be updated using one of the many helper functions defined in
    ``PlotWindow``.

    If the x_field and y_field combo do not correspond to a valid
    ``ParticleProjectionPlot``, then a ``ParticlePhasePlot``. This object can be
    modified by its own set of  helper functions defined in PhasePlot. We note
    below which arguments are only accepted by ``ParticleProjectionPlot`` and
    which arguments are only accepted by ``ParticlePhasePlot``.

    Parameters
    ----------

    ds : :class:`yt.data_objects.static_output.Dataset`
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
    weight_field : string
         The name of the weighting field.  Set to None for no weight.
    fontsize : integer
         The size of the fonts for the axis, colorbar, and tick labels.
    data_source : YTSelectionContainer Object
         Object to be used for data selection.  Defaults to a region covering
         the entire simulation.
    center : 'center', 'c', 'left', 'l', 'right', 'r', id of a global extremum, or array-like
        The coordinate of the selection's center.
        Defaults to the 'center', i.e. center of the domain.

        Centering on the min or max of a field is supported by passing a tuple
        such as ('min', ('gas', 'density')) or ('max', ('gas', 'temperature'). A
        single string may also be used (e.g. "min_density" or
        "max_temperature"), though it's not as flexible and does not allow to
        select an exact field/particle type. With this syntax, the first field
        matching the provided name is selected.
        'max' or 'm' can be used as a shortcut for ('max', ('gas', 'density'))
        'min' can be used as a shortcut for ('min', ('gas', 'density'))

        One can also select an exact point as a 3 element coordinate sequence,
        e.g. [0.5, 0.5, 0]
        Units can be specified by passing in *center* as a tuple containing a
        3-element coordinate sequence and string unit name, e.g. ([0, 0.5, 0.5], "cm"),
        or by passing in a YTArray. Code units are assumed if unspecified.

        The domain edges along the selected *axis* can be selected with
        'left'/'l' and 'right'/'r' respectively.

        This argument is only accepted by ``ParticleProjectionPlot``.
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
         This argument is only accepted by ``ParticleProjectionPlot``.
    depth : A tuple or a float
         A tuple containing the depth to project through and the string
         key of the unit: (width, 'unit').  If set to a float, code units
         are assumed. Defaults to the entire domain. This argument is only
         accepted by ``ParticleProjectionPlot``.
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
         corner of domain space. If x or y are not given, a value is inferred.
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

         This argument is only accepted by ``ParticleProjectionPlot``.
    window_size : float
         The size of the window on the longest axis (in units of inches),
         including the margins but not the colorbar. This argument is only
         accepted by ``ParticleProjectionPlot``.
    aspect : float
         The aspect ratio of the plot.  Set to None for 1. This argument is
         only accepted by ``ParticleProjectionPlot``.
    x_bins : int
        The number of bins in x field for the mesh. Defaults to 800. This
        argument is only accepted by ``ParticlePhasePlot``.
    y_bins : int
        The number of bins in y field for the mesh. Defaults to 800. This
        argument is only accepted by ``ParticlePhasePlot``.
    deposition : str
        Either 'ngp' or 'cic'. Controls what type of interpolation will be
        used to deposit the particle z_fields onto the mesh. Defaults to 'ngp'.
    figure_size : int
        Size in inches of the image. Defaults to 8 (product an 8x8 inch figure).
        This argument is only accepted by ``ParticlePhasePlot``.

    Examples
    --------

    >>> from yt import load
    >>> ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> p = yt.ParticlePlot(
    ...     ds,
    ...     "particle_position_x",
    ...     "particle_position_y",
    ...     "particle_mass",
    ...     width=(0.5, 0.5),
    ... )
    >>> p.set_unit("particle_mass", "Msun")
    >>> p = yt.ParticlePlot(ds, "particle_position_x", "particle_velocity_z", color="g")

    """
    dd = kwargs.get("data_source", None)
    if dd is None:
        dd = ds.all_data()
    x_field = dd._determine_fields(x_field)[0]
    y_field = dd._determine_fields(y_field)[0]

    direction = 3
    # try potential axes for a ParticleProjectionPlot:
    for axis in [0, 1, 2]:
        xax = ds.coordinates.x_axis[axis]
        yax = ds.coordinates.y_axis[axis]
        ax_field_template = "particle_position_%s"
        xf = ax_field_template % ds.coordinates.axis_name[xax]
        yf = ax_field_template % ds.coordinates.axis_name[yax]
        if (x_field[1], y_field[1]) in [(xf, yf), (yf, xf)]:
            direction = axis
            break

    if direction < 3:
        # Make a ParticleProjectionPlot
        return ParticleProjectionPlot(ds, direction, z_fields, color, *args, **kwargs)

    # Does not correspond to any valid PlotWindow-style plot,
    # use ParticlePhasePlot instead
    else:
        return ParticlePhasePlot(dd, x_field, y_field, z_fields, color, *args, **kwargs)
