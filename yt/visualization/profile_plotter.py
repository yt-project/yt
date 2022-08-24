import base64
import builtins
import os
from functools import wraps
from typing import Any, Dict, Iterable, Optional, Tuple, Union

import matplotlib
import numpy as np
from more_itertools.more import always_iterable, unzip

from yt.data_objects.profiles import create_profile, sanitize_field_tuple_keys
from yt.data_objects.static_output import Dataset
from yt.frontends.ytdata.data_structures import YTProfileDataset
from yt.funcs import iter_fields, matplotlib_style_context
from yt.utilities.exceptions import YTNotInsideNotebook
from yt.visualization._handlers import ColorbarHandler, NormHandler
from yt.visualization.base_plot_types import ImagePlotMPL, PlotMPL

from ..data_objects.selection_objects.data_selection_objects import YTSelectionContainer
from ._commons import validate_image_name
from .plot_container import (
    BaseLinePlot,
    ImagePlotContainer,
    invalidate_plot,
    validate_plot,
)


def invalidate_profile(f):
    @wraps(f)
    def newfunc(*args, **kwargs):
        rv = f(*args, **kwargs)
        args[0]._profile_valid = False
        return rv

    return newfunc


def sanitize_label(labels, nprofiles):
    labels = list(always_iterable(labels)) or [None]

    if len(labels) == 1:
        labels = labels * nprofiles

    if len(labels) != nprofiles:
        raise ValueError(
            f"Number of labels {len(labels)} must match number of profiles {nprofiles}"
        )

    invalid_data = [
        (label, type(label))
        for label in labels
        if label is not None and not isinstance(label, str)
    ]
    if invalid_data:
        invalid_labels, types = unzip(invalid_data)
        raise TypeError(
            "All labels must be None or a string, "
            f"received {invalid_labels} with type {types}"
        )

    return labels


def data_object_or_all_data(data_source):
    if isinstance(data_source, Dataset):
        data_source = data_source.all_data()

    if not isinstance(data_source, YTSelectionContainer):
        raise RuntimeError("data_source must be a yt selection data object")

    return data_source


class ProfilePlot(BaseLinePlot):
    r"""
    Create a 1d profile plot from a data source or from a list
    of profile objects.

    Given a data object (all_data, region, sphere, etc.), an x field,
    and a y field (or fields), this will create a one-dimensional profile
    of the average (or total) value of the y field in bins of the x field.

    This can be used to create profiles from given fields or to plot
    multiple profiles created from
    `yt.data_objects.profiles.create_profile`.

    Parameters
    ----------
    data_source : YTSelectionContainer Object
        The data object to be profiled, such as all_data, region, or
        sphere. If a dataset is passed in instead, an all_data data object
        is generated internally from the dataset.
    x_field : str
        The binning field for the profile.
    y_fields : str or list
        The field or fields to be profiled.
    weight_field : str
        The weight field for calculating weighted averages. If None,
        the profile values are the sum of the field values within the bin.
        Otherwise, the values are a weighted average.
        Default : ("gas", "mass")
    n_bins : int
        The number of bins in the profile.
        Default: 64.
    accumulation : bool
        If True, the profile values for a bin N are the cumulative sum of
        all the values from bin 0 to N.
        Default: False.
    fractional : If True the profile values are divided by the sum of all
        the profile data such that the profile represents a probability
        distribution function.
    label : str or list of strings
        If a string, the label to be put on the line plotted.  If a list,
        this should be a list of labels for each profile to be overplotted.
        Default: None.
    plot_spec : dict or list of dicts
        A dictionary or list of dictionaries containing plot keyword
        arguments.  For example, dict(color="red", linestyle=":").
        Default: None.
    x_log : bool
        Whether the x_axis should be plotted with a logarithmic
        scaling (True), or linear scaling (False).
        Default: True.
    y_log : dict or bool
        A dictionary containing field:boolean pairs, setting the logarithmic
        property for that field. May be overridden after instantiation using
        set_log
        A single boolean can be passed to signify all fields should use
        logarithmic (True) or linear scaling (False).
        Default: True.

    Examples
    --------

    This creates profiles of a single dataset.

    >>> import yt
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = ds.all_data()
    >>> plot = yt.ProfilePlot(
    ...     ad,
    ...     ("gas", "density"),
    ...     [("gas", "temperature"), ("gas", "velocity_x")],
    ...     weight_field=("gas", "mass"),
    ...     plot_spec=dict(color="red", linestyle="--"),
    ... )
    >>> plot.save()

    This creates profiles from a time series object.

    >>> es = yt.load_simulation("AMRCosmology.enzo", "Enzo")
    >>> es.get_time_series()

    >>> profiles = []
    >>> labels = []
    >>> plot_specs = []
    >>> for ds in es[-4:]:
    ...     ad = ds.all_data()
    ...     profiles.append(
    ...         create_profile(
    ...             ad,
    ...             [("gas", "density")],
    ...             fields=[("gas", "temperature"), ("gas", "velocity_x")],
    ...         )
    ...     )
    ...     labels.append(ds.current_redshift)
    ...     plot_specs.append(dict(linestyle="--", alpha=0.7))

    >>> plot = yt.ProfilePlot.from_profiles(
    ...     profiles, labels=labels, plot_specs=plot_specs
    ... )
    >>> plot.save()

    Use set_line_property to change line properties of one or all profiles.

    """
    _default_figure_size = (10.0, 8.0)
    _default_font_size = 18.0

    x_log = None
    y_log = None
    x_title = None
    y_title = None
    _plot_valid = False

    def __init__(
        self,
        data_source,
        x_field,
        y_fields,
        weight_field=("gas", "mass"),
        n_bins=64,
        accumulation=False,
        fractional=False,
        label=None,
        plot_spec=None,
        x_log=True,
        y_log=True,
    ):
        data_source = data_object_or_all_data(data_source)
        y_fields = list(iter_fields(y_fields))
        logs = {x_field: bool(x_log)}
        if isinstance(y_log, bool):
            y_log = {y_field: y_log for y_field in y_fields}

        if isinstance(data_source.ds, YTProfileDataset):
            profiles = [data_source.ds.profile]
        else:
            profiles = [
                create_profile(
                    data_source,
                    [x_field],
                    n_bins=[n_bins],
                    fields=y_fields,
                    weight_field=weight_field,
                    accumulation=accumulation,
                    fractional=fractional,
                    logs=logs,
                )
            ]

        if plot_spec is None:
            plot_spec = [dict() for p in profiles]
        if not isinstance(plot_spec, list):
            plot_spec = [plot_spec.copy() for p in profiles]

        ProfilePlot._initialize_instance(
            self, data_source, profiles, label, plot_spec, y_log
        )

    @classmethod
    def _initialize_instance(
        cls,
        obj,
        data_source,
        profiles,
        labels,
        plot_specs,
        y_log,
    ):
        obj._plot_title = {}
        obj._plot_text = {}
        obj._text_xpos = {}
        obj._text_ypos = {}
        obj._text_kwargs = {}

        super(ProfilePlot, obj).__init__(data_source)
        obj.profiles = list(always_iterable(profiles))
        obj.x_log = None
        obj.y_log = sanitize_field_tuple_keys(y_log, data_source) or {}
        obj.y_title = {}
        obj.x_title = None
        obj.label = sanitize_label(labels, len(obj.profiles))
        if plot_specs is None:
            plot_specs = [dict() for p in obj.profiles]
        obj.plot_spec = plot_specs
        obj._xlim = (None, None)
        obj._setup_plots()
        return obj

    def _get_axrect(self):
        return (0.1, 0.1, 0.8, 0.8)

    @validate_plot
    def save(
        self,
        name: Optional[str] = None,
        suffix: Optional[str] = None,
        mpl_kwargs: Optional[Dict[str, Any]] = None,
    ):
        r"""
        Saves a 1d profile plot.

        Parameters
        ----------
        name : str, optional
            The output file keyword.
        suffix : string, optional
            Specify the image type by its suffix. If not specified, the output
            type will be inferred from the filename. Defaults to '.png'.
        mpl_kwargs : dict, optional
            A dict of keyword arguments to be passed to matplotlib.
        """
        if not self._plot_valid:
            self._setup_plots()

        # Mypy is hardly convinced that we have a `profiles` attribute
        # at this stage, so we're lasily going to deactivate it locally
        unique = set(self.plots.values())
        iters: Iterable[Tuple[Union[int, Tuple[str, str]], PlotMPL]]
        if len(unique) < len(self.plots):
            iters = enumerate(sorted(unique))
        else:
            iters = self.plots.items()

        if name is None:
            if len(self.profiles) == 1:  # type: ignore
                name = str(self.profiles[0].ds)  # type: ignore
            else:
                name = "Multi-data"

        name = validate_image_name(name, suffix)
        prefix, suffix = os.path.splitext(name)

        xfn = self.profiles[0].x_field  # type: ignore
        if isinstance(xfn, tuple):
            xfn = xfn[1]

        names = []
        for uid, plot in iters:
            if isinstance(uid, tuple):
                uid = uid[1]  # type: ignore
            uid_name = f"{prefix}_1d-Profile_{xfn}_{uid}{suffix}"
            names.append(uid_name)
            with matplotlib_style_context():
                plot.save(uid_name, mpl_kwargs=mpl_kwargs)
        return names

    @validate_plot
    def show(self):
        r"""This will send any existing plots to the IPython notebook.

        If yt is being run from within an IPython session, and it is able to
        determine this, this function will send any existing plots to the
        notebook for display.

        If yt can't determine if it's inside an IPython session, it will raise
        YTNotInsideNotebook.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> pp = ProfilePlot(ds.all_data(), ("gas", "density"), ("gas", "temperature"))
        >>> pp.show()

        """
        if "__IPYTHON__" in dir(builtins):
            from IPython.display import display

            display(self)
        else:
            raise YTNotInsideNotebook

    @validate_plot
    def _repr_html_(self):
        """Return an html representation of the plot object. Will display as a
        png for each WindowPlotMPL instance in self.plots"""
        ret = ""
        unique = set(self.plots.values())
        if len(unique) < len(self.plots):
            iters = sorted(unique)
        else:
            iters = self.plots.values()
        for plot in iters:
            with matplotlib_style_context():
                img = plot._repr_png_()
            img = base64.b64encode(img).decode()
            ret += (
                r'<img style="max-width:100%;max-height:100%;" '
                r'src="data:image/png;base64,{}"><br>'.format(img)
            )
        return ret

    def _setup_plots(self):
        if self._plot_valid:
            return
        for f, p in self.plots.items():
            p.axes.cla()
            if f in self._plot_text:
                p.axes.text(
                    self._text_xpos[f],
                    self._text_ypos[f],
                    self._plot_text[f],
                    fontproperties=self._font_properties,
                    **self._text_kwargs[f],
                )
        self._set_font_properties()

        for i, profile in enumerate(self.profiles):
            for field, field_data in profile.items():
                plot = self._get_plot_instance(field)
                plot.axes.plot(
                    np.array(profile.x),
                    np.array(field_data),
                    label=self.label[i],
                    **self.plot_spec[i],
                )

        for profile in self.profiles:
            for fname in profile.keys():
                axes = self.plots[fname].axes
                xscale, yscale = self._get_field_log(fname, profile)
                xtitle, ytitle = self._get_field_title(fname, profile)

                axes.set_xscale(xscale)
                axes.set_yscale(yscale)

                axes.set_ylabel(ytitle)
                axes.set_xlabel(xtitle)

                pnh = self.plots[fname].norm_handler

                axes.set_ylim(pnh.vmin, pnh.vmax)
                axes.set_xlim(*self._xlim)

                if fname in self._plot_title:
                    axes.set_title(self._plot_title[fname])

                if any(self.label):
                    axes.legend(loc="best")
        self._set_font_properties()
        self._plot_valid = True

    @classmethod
    def from_profiles(cls, profiles, labels=None, plot_specs=None, y_log=None):
        r"""
        Instantiate a ProfilePlot object from a list of profiles
        created with :func:`~yt.data_objects.profiles.create_profile`.

        Parameters
        ----------
        profiles : a profile or list of profiles
            A single profile or list of profile objects created with
            :func:`~yt.data_objects.profiles.create_profile`.
        labels : list of strings
            A list of labels for each profile to be overplotted.
            Default: None.
        plot_specs : list of dicts
            A list of dictionaries containing plot keyword
            arguments.  For example, [dict(color="red", linestyle=":")].
            Default: None.

        Examples
        --------

        >>> from yt import load_simulation
        >>> es = load_simulation("AMRCosmology.enzo", "Enzo")
        >>> es.get_time_series()

        >>> profiles = []
        >>> labels = []
        >>> plot_specs = []
        >>> for ds in es[-4:]:
        ...     ad = ds.all_data()
        ...     profiles.append(
        ...         create_profile(
        ...             ad,
        ...             [("gas", "density")],
        ...             fields=[("gas", "temperature"), ("gas", "velocity_x")],
        ...         )
        ...     )
        ...     labels.append(ds.current_redshift)
        ...     plot_specs.append(dict(linestyle="--", alpha=0.7))
        >>> plot = ProfilePlot.from_profiles(
        ...     profiles, labels=labels, plot_specs=plot_specs
        ... )
        >>> plot.save()

        """
        if labels is not None and len(profiles) != len(labels):
            raise RuntimeError("Profiles list and labels list must be the same size.")
        if plot_specs is not None and len(plot_specs) != len(profiles):
            raise RuntimeError(
                "Profiles list and plot_specs list must be the same size."
            )
        obj = cls.__new__(cls)
        profiles = list(always_iterable(profiles))
        return cls._initialize_instance(
            obj, profiles[0].data_source, profiles, labels, plot_specs, y_log
        )

    @invalidate_plot
    def set_line_property(self, property, value, index=None):
        r"""
        Set properties for one or all lines to be plotted.

        Parameters
        ----------
        property : str
            The line property to be set.
        value : str, int, float
            The value to set for the line property.
        index : int
            The index of the profile in the list of profiles to be
            changed.  If None, change all plotted lines.
            Default : None.

        Examples
        --------

        Change all the lines in a plot
        plot.set_line_property("linestyle", "-")

        Change a single line.
        plot.set_line_property("linewidth", 4, index=0)

        """
        if index is None:
            specs = self.plot_spec
        else:
            specs = [self.plot_spec[index]]
        for spec in specs:
            spec[property] = value
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
            for field in list(self.profiles[0].field_data.keys()):
                self.y_log[field] = log
        else:
            (field,) = self.profiles[0].data_source._determine_fields([field])
            if field == self.profiles[0].x_field:
                self.x_log = log
            elif field in self.profiles[0].field_data:
                self.y_log[field] = log
            else:
                raise KeyError(f"Field {field} not in profile plot!")
        return self

    @invalidate_plot
    def set_ylabel(self, field, label):
        """Sets a new ylabel for the specified fields

        Parameters
        ----------
        field : string
           The name of the field that is to be changed.

        label : string
           The label to be placed on the y-axis
        """
        if field == "all":
            for field in self.profiles[0].field_data:
                self.y_title[field] = label
        else:
            (field,) = self.profiles[0].data_source._determine_fields([field])
            if field in self.profiles[0].field_data:
                self.y_title[field] = label
            else:
                raise KeyError(f"Field {field} not in profile plot!")

        return self

    @invalidate_plot
    def set_xlabel(self, label):
        """Sets a new xlabel for all profiles

        Parameters
        ----------
        label : string
           The label to be placed on the x-axis
        """
        self.x_title = label

        return self

    @invalidate_plot
    def set_unit(self, field, unit):
        """Sets a new unit for the requested field

        Parameters
        ----------
        field : string
           The name of the field that is to be changed.

        unit : string or Unit object
           The name of the new unit.
        """
        fd = self.profiles[0].data_source._determine_fields(field)[0]
        for profile in self.profiles:
            if fd == profile.x_field:
                profile.set_x_unit(unit)
            elif fd[1] in self.profiles[0].field_map:
                profile.set_field_unit(field, unit)
            else:
                raise KeyError(f"Field {field} not in profile plot!")
        return self

    @invalidate_plot
    def set_xlim(self, xmin=None, xmax=None):
        """Sets the limits of the bin field

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
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> pp = yt.ProfilePlot(
        ...     ds.all_data(), ("gas", "density"), ("gas", "temperature")
        ... )
        >>> pp.set_xlim(1e-29, 1e-24)
        >>> pp.save()

        """
        self._xlim = (xmin, xmax)
        for i, p in enumerate(self.profiles):
            if xmin is None:
                xmi = p.x_bins.min()
            else:
                xmi = xmin
            if xmax is None:
                xma = p.x_bins.max()
            else:
                xma = xmax
            extrema = {p.x_field: ((xmi, str(p.x.units)), (xma, str(p.x.units)))}
            units = {p.x_field: str(p.x.units)}
            if self.x_log is None:
                logs = None
            else:
                logs = {p.x_field: self.x_log}
            for field in p.field_map.values():
                units[field] = str(p.field_data[field].units)
            self.profiles[i] = create_profile(
                p.data_source,
                p.x_field,
                n_bins=len(p.x_bins) - 1,
                fields=list(p.field_map.values()),
                weight_field=p.weight_field,
                accumulation=p.accumulation,
                fractional=p.fractional,
                logs=logs,
                extrema=extrema,
                units=units,
            )
        return self

    @invalidate_plot
    def set_ylim(self, field, ymin=None, ymax=None):
        """Sets the plot limits for the specified field we are binning.

        Parameters
        ----------

        field : string or field tuple

        The field that we want to adjust the plot limits for.

        ymin : float or None
          The new y minimum.  Defaults to None, which leaves the ymin
          unchanged.

        ymax : float or None
          The new y maximum.  Defaults to None, which leaves the ymax
          unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> pp = yt.ProfilePlot(
        ...     ds.all_data(),
        ...     ("gas", "density"),
        ...     [("gas", "temperature"), ("gas", "velocity_x")],
        ... )
        >>> pp.set_ylim(("gas", "temperature"), 1e4, 1e6)
        >>> pp.save()

        """
        fields = list(self.plots.keys()) if field == "all" else field
        for profile in self.profiles:
            for field in profile.data_source._determine_fields(fields):
                if field in profile.field_map:
                    field = profile.field_map[field]
                pnh = self.plots[field].norm_handler
                pnh.vmin = ymin
                pnh.vmax = ymax
                # Continue on to the next profile.
                break
        return self

    def _set_font_properties(self):
        for f in self.plots:
            self.plots[f]._set_font_properties(self._font_properties, self._font_color)

    def _get_field_log(self, field_y, profile):
        yfi = profile.field_info[field_y]
        if self.x_log is None:
            x_log = profile.x_log
        else:
            x_log = self.x_log
        y_log = self.y_log.get(field_y, yfi.take_log)
        scales = {True: "log", False: "linear"}
        return scales[x_log], scales[y_log]

    def _get_field_label(self, field, field_info, field_unit, fractional=False):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple):
            field = field[1]
        if field_name is None:
            field_name = field_info.get_latex_display_name()
        elif field_name.find("$") == -1:
            field_name = field_name.replace(" ", r"\ ")
            field_name = r"$\rm{" + field_name + r"}$"
        if fractional:
            label = field_name + r"$\rm{\ Probability\ Density}$"
        elif field_unit is None or field_unit == "":
            label = field_name
        elif r"\frac" in field_unit:
            label = field_name + r"$\ \ \left(" + field_unit + r"\right)$"
        else:
            label = field_name + r"$\ \ (" + field_unit + r")$"
        return label

    def _get_field_title(self, field_y, profile):
        field_x = profile.x_field
        xfi = profile.field_info[field_x]
        yfi = profile.field_info[field_y]
        x_unit = profile.x.units
        y_unit = profile.field_units[field_y]
        fractional = profile.fractional
        x_title = self.x_title or self._get_field_label(field_x, xfi, x_unit)
        y_title = self.y_title.get(field_y, None) or self._get_field_label(
            field_y, yfi, y_unit, fractional
        )

        return (x_title, y_title)

    @invalidate_plot
    def annotate_title(self, title, field="all"):
        r"""Set a title for the plot.

        Parameters
        ----------
        title : str
          The title to add.
        field : str or list of str
          The field name for which title needs to be set.

        Examples
        --------
        >>> # To set title for all the fields:
        >>> plot.annotate_title("This is a Profile Plot")

        >>> # To set title for specific fields:
        >>> plot.annotate_title("Profile Plot for Temperature", ("gas", "temperature"))

        >>> # Setting same plot title for both the given fields
        >>> plot.annotate_title(
        ...     "Profile Plot: Temperature-Dark Matter Density",
        ...     [("gas", "temperature"), ("deposit", "dark_matter_density")],
        ... )

        """
        fields = list(self.plots.keys()) if field == "all" else field
        for profile in self.profiles:
            for field in profile.data_source._determine_fields(fields):
                if field in profile.field_map:
                    field = profile.field_map[field]
                self._plot_title[field] = title
        return self

    @invalidate_plot
    def annotate_text(self, xpos=0.0, ypos=0.0, text=None, field="all", **text_kwargs):
        r"""Allow the user to insert text onto the plot

        The x-position and y-position must be given as well as the text string.
        Add *text* to plot at location *xpos*, *ypos* in plot coordinates for
        the given fields or by default for all fields.
        (see example below).

        Parameters
        ----------
        xpos : float
          Position on plot in x-coordinates.
        ypos : float
          Position on plot in y-coordinates.
        text : str
          The text to insert onto the plot.
        field : str or tuple
          The name of the field to add text to.
        **text_kwargs : dict
          Extra keyword arguments will be passed to matplotlib text instance

        >>> import yt
        >>> from yt.units import kpc
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> my_galaxy = ds.disk(ds.domain_center, [0.0, 0.0, 1.0], 10 * kpc, 3 * kpc)
        >>> plot = yt.ProfilePlot(
        ...     my_galaxy, ("gas", "density"), [("gas", "temperature")]
        ... )

        >>> # Annotate text for all the fields
        >>> plot.annotate_text(1e-26, 1e5, "This is annotated text in the plot area.")
        >>> plot.save()

        >>> # Annotate text for a given field
        >>> plot.annotate_text(1e-26, 1e5, "Annotated text", ("gas", "temperature"))
        >>> plot.save()

        >>> # Annotate text for multiple fields
        >>> fields = [("gas", "temperature"), ("gas", "density")]
        >>> plot.annotate_text(1e-26, 1e5, "Annotated text", fields)
        >>> plot.save()

        """
        fields = list(self.plots.keys()) if field == "all" else field
        for profile in self.profiles:
            for field in profile.data_source._determine_fields(fields):
                if field in profile.field_map:
                    field = profile.field_map[field]
                self._plot_text[field] = text
                self._text_xpos[field] = xpos
                self._text_ypos[field] = ypos
                self._text_kwargs[field] = text_kwargs
        return self


class PhasePlot(ImagePlotContainer):
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
        sphere. If a dataset is passed in instead, an all_data data object
        is generated internally from the dataset.
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
        Default : ("gas", "mass")
    x_bins : int
        The number of bins in x field for the profile.
        Default: 128.
    y_bins : int
        The number of bins in y field for the profile.
        Default: 128.
    accumulation : bool or list of bools
        If True, the profile values for a bin n are the cumulative sum of
        all the values from bin 0 to n.  If -True, the sum is reversed so
        that the value for bin n is the cumulative sum from bin N (total bins)
        to n.  A list of values can be given to control the summation in each
        dimension independently.
        Default: False.
    fractional : If True the profile values are divided by the sum of all
        the profile data such that the profile represents a probability
        distribution function.
    fontsize : int
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
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> ad = ds.all_data()
    >>> plot = yt.PhasePlot(
    ...     ad,
    ...     ("gas", "density"),
    ...     ("gas", "temperature"),
    ...     [("gas", "mass")],
    ...     weight_field=None,
    ... )
    >>> plot.save()

    >>> # Change plot properties.
    >>> plot.set_cmap(("gas", "mass"), "jet")
    >>> plot.set_zlim(("gas", "mass"), 1e8, 1e13)
    >>> plot.annotate_title("This is a phase plot")

    """
    x_log = None
    y_log = None
    plot_title = None
    _plot_valid = False
    _profile_valid = False
    _plot_type = "Phase"
    _xlim = (None, None)
    _ylim = (None, None)

    def __init__(
        self,
        data_source,
        x_field,
        y_field,
        z_fields,
        weight_field=("gas", "mass"),
        x_bins=128,
        y_bins=128,
        accumulation=False,
        fractional=False,
        fontsize=18,
        figure_size=8.0,
        shading="nearest",
    ):

        data_source = data_object_or_all_data(data_source)

        if isinstance(z_fields, tuple):
            z_fields = [z_fields]
        z_fields = list(always_iterable(z_fields))

        if isinstance(data_source.ds, YTProfileDataset):
            profile = data_source.ds.profile
        else:
            profile = create_profile(
                data_source,
                [x_field, y_field],
                z_fields,
                n_bins=[x_bins, y_bins],
                weight_field=weight_field,
                accumulation=accumulation,
                fractional=fractional,
            )

        type(self)._initialize_instance(
            self, data_source, profile, fontsize, figure_size, shading
        )

    @classmethod
    def _initialize_instance(
        cls, obj, data_source, profile, fontsize, figure_size, shading
    ):
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
        obj._profile = profile
        obj._shading = shading
        obj._profile_valid = True
        obj._xlim = (None, None)
        obj._ylim = (None, None)
        super(PhasePlot, obj).__init__(data_source, figure_size, fontsize)
        obj._setup_plots()
        obj._initfinished = True
        return obj

    def _get_field_title(self, field_z, profile):
        field_x = profile.x_field
        field_y = profile.y_field
        xfi = profile.field_info[field_x]
        yfi = profile.field_info[field_y]
        zfi = profile.field_info[field_z]
        x_unit = profile.x.units
        y_unit = profile.y.units
        z_unit = profile.field_units[field_z]
        fractional = profile.fractional
        x_label, y_label, z_label = self._get_axes_labels(field_z)
        x_title = x_label or self._get_field_label(field_x, xfi, x_unit)
        y_title = y_label or self._get_field_label(field_y, yfi, y_unit)
        z_title = z_label or self._get_field_label(field_z, zfi, z_unit, fractional)
        return (x_title, y_title, z_title)

    def _get_field_label(self, field, field_info, field_unit, fractional=False):
        field_unit = field_unit.latex_representation()
        field_name = field_info.display_name
        if isinstance(field, tuple):
            field = field[1]
        if field_name is None:
            field_name = field_info.get_latex_display_name()
        elif field_name.find("$") == -1:
            field_name = field_name.replace(" ", r"\ ")
            field_name = r"$\rm{" + field_name + r"}$"
        if fractional:
            label = field_name + r"$\rm{\ Probability\ Density}$"
        elif field_unit is None or field_unit == "":
            label = field_name
        elif r"\frac" in field_unit:
            label = field_name + r"$\ \ \left(" + field_unit + r"\right)$"
        else:
            label = field_name + r"$\ \ (" + field_unit + r")$"
        return label

    def _get_field_log(self, field_z, profile):
        zfi = profile.field_info[field_z]
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
        scales = {True: "log", False: "linear"}
        return scales[x_log], scales[y_log], scales[z_log]

    @property
    def profile(self):
        if not self._profile_valid:
            self._recreate_profile()
        return self._profile

    @property
    def fields(self):
        return list(self.plots.keys())

    def _setup_plots(self):
        if self._plot_valid:
            return
        for f, data in self.profile.items():
            fig = None
            axes = None
            cax = None
            draw_axes = True
            xlim = self._xlim
            ylim = self._ylim

            if f in self.plots:
                pnh = self.plots[f].norm_handler
                cbh = self.plots[f].colorbar_handler
                draw_axes = self.plots[f]._draw_axes
                if self.plots[f].figure is not None:
                    fig = self.plots[f].figure
                    axes = self.plots[f].axes
                    cax = self.plots[f].cax
            else:
                pnh, cbh = self._get_default_handlers(
                    field=f, default_display_units=self.profile[f].units
                )

            x_scale, y_scale, z_scale = self._get_field_log(f, self.profile)
            x_title, y_title, z_title = self._get_field_title(f, self.profile)

            font_size = self._font_properties.get_size()
            f = self.profile.data_source._determine_fields(f)[0]

            # if this is a Particle Phase Plot AND if we using a single color,
            # override the colorbar here.
            splat_color = getattr(self, "splat_color", None)
            if splat_color is not None:
                cbh.cmap = matplotlib.colors.ListedColormap(splat_color, "dummy")

            masked_data = data.copy()
            masked_data[~self.profile.used] = np.nan
            self.plots[f] = PhasePlotMPL(
                self.profile.x,
                self.profile.y,
                masked_data,
                x_scale,
                y_scale,
                self.figure_size,
                font_size,
                fig,
                axes,
                cax,
                shading=self._shading,
                norm_handler=pnh,
                colorbar_handler=cbh,
            )

            self.plots[f]._toggle_axes(draw_axes)
            self.plots[f]._toggle_colorbar(cbh.draw_cbar)

            self.plots[f].axes.xaxis.set_label_text(x_title)
            self.plots[f].axes.yaxis.set_label_text(y_title)
            self.plots[f].cax.yaxis.set_label_text(z_title)

            self.plots[f].axes.set_xlim(xlim)
            self.plots[f].axes.set_ylim(ylim)

            if f in self._plot_text:
                self.plots[f].axes.text(
                    self._text_xpos[f],
                    self._text_ypos[f],
                    self._plot_text[f],
                    fontproperties=self._font_properties,
                    **self._text_kwargs[f],
                )

            if f in self.plot_title:
                self.plots[f].axes.set_title(self.plot_title[f])

            # x-y axes minorticks
            if f not in self._minorticks:
                self._minorticks[f] = True
            if self._minorticks[f]:
                self.plots[f].axes.minorticks_on()
            else:
                self.plots[f].axes.minorticks_off()

        self._set_font_properties()

        # if this is a particle plot with one color only, hide the cbar here
        if hasattr(self, "use_cbar") and not self.use_cbar:
            self.plots[f].hide_colorbar()

        self._plot_valid = True

    @classmethod
    def from_profile(cls, profile, fontsize=18, figure_size=8.0, shading="nearest"):
        r"""
        Instantiate a PhasePlot object from a profile object created
        with :func:`~yt.data_objects.profiles.create_profile`.

        Parameters
        ----------
        profile : An instance of :class:`~yt.data_objects.profiles.ProfileND`
             A single profile object.
        fontsize : float
             The fontsize to use, in points.
        figure_size : float
             The figure size to use, in inches.
        shading : str
            This argument is directly passed down to matplotlib.axes.Axes.pcolormesh
            see
            https://matplotlib.org/3.3.1/gallery/images_contours_and_fields/pcolormesh_grids.html#sphx-glr-gallery-images-contours-and-fields-pcolormesh-grids-py  # noqa
            Default: 'nearest'

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> extrema = {
        ...     ("gas", "density"): (1e-31, 1e-24),
        ...     ("gas", "temperature"): (1e1, 1e8),
        ...     ("gas", "mass"): (1e-6, 1e-1),
        ... }
        >>> profile = yt.create_profile(
        ...     ds.all_data(),
        ...     [("gas", "density"), ("gas", "temperature")],
        ...     fields=[("gas", "mass")],
        ...     extrema=extrema,
        ...     fractional=True,
        ... )
        >>> ph = yt.PhasePlot.from_profile(profile)
        >>> ph.save()
        """
        obj = cls.__new__(cls)
        data_source = profile.data_source
        return cls._initialize_instance(
            obj, data_source, profile, fontsize, figure_size, shading
        )

    def annotate_text(self, xpos=0.0, ypos=0.0, text=None, **text_kwargs):
        r"""
        Allow the user to insert text onto the plot
        The x-position and y-position must be given as well as the text string.
        Add *text* tp plot at location *xpos*, *ypos* in plot coordinates
        (see example below).

        Parameters
        ----------
        xpos : float
          Position on plot in x-coordinates.
        ypos : float
          Position on plot in y-coordinates.
        text : str
          The text to insert onto the plot.
        **text_kwargs : dict
          Extra keyword arguments will be passed to matplotlib text instance

        >>> plot.annotate_text(1e-15, 5e4, "Hello YT")

        """
        for f in self.data_source._determine_fields(list(self.plots.keys())):
            if self.plots[f].figure is not None and text is not None:
                self.plots[f].axes.text(
                    xpos,
                    ypos,
                    text,
                    fontproperties=self._font_properties,
                    **text_kwargs,
                )
            self._plot_text[f] = text
            self._text_xpos[f] = xpos
            self._text_ypos[f] = ypos
            self._text_kwargs[f] = text_kwargs
        return self

    @validate_plot
    def save(
        self, name: Optional[str] = None, suffix: Optional[str] = None, mpl_kwargs=None
    ):
        r"""
        Saves a 2d profile plot.

        Parameters
        ----------
        name : str, optional
            The output file keyword.
        suffix : string, optional
           Specify the image type by its suffix. If not specified, the output
           type will be inferred from the filename. Defaults to '.png'.
        mpl_kwargs : dict, optional
           A dict of keyword arguments to be passed to matplotlib.

        >>> plot.save(mpl_kwargs={"bbox_inches": "tight"})

        """
        names = []
        if not self._plot_valid:
            self._setup_plots()
        if mpl_kwargs is None:
            mpl_kwargs = {}
        if name is None:
            name = str(self.profile.ds)
        name = os.path.expanduser(name)
        xfn = self.profile.x_field
        yfn = self.profile.y_field
        if isinstance(xfn, tuple):
            xfn = xfn[1]
        if isinstance(yfn, tuple):
            yfn = yfn[1]
        for f in self.profile.field_data:
            _f = f
            if isinstance(f, tuple):
                _f = _f[1]
            middle = f"2d-Profile_{xfn}_{yfn}_{_f}"
            splitname = os.path.split(name)
            if splitname[0] != "" and not os.path.isdir(splitname[0]):
                os.makedirs(splitname[0])
            if os.path.isdir(name) and name != str(self.profile.ds):
                name = name + (os.sep if name[-1] != os.sep else "")
                name += str(self.profile.ds)

            new_name = validate_image_name(name, suffix)
            if new_name == name:
                for v in self.plots.values():
                    out_name = v.save(name, mpl_kwargs)
                    names.append(out_name)
                return names

            name = new_name
            prefix, suffix = os.path.splitext(name)
            name = f"{prefix}_{middle}{suffix}"

            names.append(name)
            self.plots[f].save(name, mpl_kwargs)
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

        >>> plot.set_title(("gas", "mass"), "This is a phase plot")
        """
        self.plot_title[self.data_source._determine_fields(field)[0]] = title
        return self

    @invalidate_plot
    def annotate_title(self, title):
        """Set a title for the plot.

        Parameters
        ----------
        title : str
            The title to add.

        Examples
        --------

        >>> plot.annotate_title("This is a phase plot")

        """
        for f in self._profile.field_data:
            if isinstance(f, tuple):
                f = f[1]
            self.plot_title[self.data_source._determine_fields(f)[0]] = title
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
        p = self._profile
        if field == "all":
            self.x_log = log
            self.y_log = log
            for field in p.field_data:
                self.z_log[field] = log
            self._profile_valid = False
        else:
            (field,) = self.profile.data_source._determine_fields([field])
            if field == p.x_field:
                self.x_log = log
                self._profile_valid = False
            elif field == p.y_field:
                self.y_log = log
                self._profile_valid = False
            elif field in p.field_data:
                super().set_log(field, log)
            else:
                raise KeyError(f"Field {field} not in phase plot!")
        return self

    @invalidate_plot
    def set_unit(self, field, unit):
        """Sets a new unit for the requested field

        Parameters
        ----------
        field : string
           The name of the field that is to be changed.

        unit : string or Unit object
           The name of the new unit.
        """
        fd = self.data_source._determine_fields(field)[0]
        if fd == self.profile.x_field:
            self.profile.set_x_unit(unit)
        elif fd == self.profile.y_field:
            self.profile.set_y_unit(unit)
        elif fd in self.profile.field_data.keys():
            self.profile.set_field_unit(field, unit)
            self.plots[field].norm_handler.display_units = unit
        else:
            raise KeyError(f"Field {field} not in phase plot!")
        return self

    @invalidate_plot
    @invalidate_profile
    def set_xlim(self, xmin=None, xmax=None):
        """Sets the limits of the x bin field

        Parameters
        ----------

        xmin : float or None
          The new x minimum in the current x-axis units.  Defaults to None,
          which leaves the xmin unchanged.

        xmax : float or None
          The new x maximum in the current x-axis units.  Defaults to None,
          which leaves the xmax unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> pp = yt.PhasePlot(ds.all_data(), "density", "temperature", ("gas", "mass"))
        >>> pp.set_xlim(1e-29, 1e-24)
        >>> pp.save()

        """
        p = self._profile
        if xmin is None:
            xmin = p.x_bins.min()
        elif not hasattr(xmin, "units"):
            xmin = self.ds.quan(xmin, p.x_bins.units)
        if xmax is None:
            xmax = p.x_bins.max()
        elif not hasattr(xmax, "units"):
            xmax = self.ds.quan(xmax, p.x_bins.units)
        self._xlim = (xmin, xmax)
        return self

    @invalidate_plot
    @invalidate_profile
    def set_ylim(self, ymin=None, ymax=None):
        """Sets the plot limits for the y bin field.

        Parameters
        ----------

        ymin : float or None
          The new y minimum in the current y-axis units.  Defaults to None,
          which leaves the ymin unchanged.

        ymax : float or None
          The new y maximum in the current y-axis units.  Defaults to None,
          which leaves the ymax unchanged.

        Examples
        --------

        >>> import yt
        >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
        >>> pp = yt.PhasePlot(
        ...     ds.all_data(),
        ...     ("gas", "density"),
        ...     ("gas", "temperature"),
        ...     ("gas", "mass"),
        ... )
        >>> pp.set_ylim(1e4, 1e6)
        >>> pp.save()

        """
        p = self._profile
        if ymin is None:
            ymin = p.y_bins.min()
        elif not hasattr(ymin, "units"):
            ymin = self.ds.quan(ymin, p.y_bins.units)
        if ymax is None:
            ymax = p.y_bins.max()
        elif not hasattr(ymax, "units"):
            ymax = self.ds.quan(ymax, p.y_bins.units)
        self._ylim = (ymin, ymax)
        return self

    def _recreate_profile(self):
        p = self._profile
        units = {p.x_field: str(p.x.units), p.y_field: str(p.y.units)}
        zunits = {field: str(p.field_units[field]) for field in p.field_units}
        extrema = {p.x_field: self._xlim, p.y_field: self._ylim}
        if self.x_log is not None or self.y_log is not None:
            logs = {}
        else:
            logs = None
        if self.x_log is not None:
            logs[p.x_field] = self.x_log
        if self.y_log is not None:
            logs[p.y_field] = self.y_log
        deposition = getattr(p, "deposition", None)
        additional_kwargs = {
            "accumulation": p.accumulation,
            "fractional": p.fractional,
            "deposition": deposition,
        }
        self._profile = create_profile(
            p.data_source,
            [p.x_field, p.y_field],
            list(p.field_map.values()),
            n_bins=[len(p.x_bins) - 1, len(p.y_bins) - 1],
            weight_field=p.weight_field,
            units=units,
            extrema=extrema,
            logs=logs,
            **additional_kwargs,
        )
        for field in zunits:
            self._profile.set_field_unit(field, zunits[field])
        self._profile_valid = True


class PhasePlotMPL(ImagePlotMPL):
    """A container for a single matplotlib figure and axes for a PhasePlot"""

    def __init__(
        self,
        x_data,
        y_data,
        data,
        x_scale,
        y_scale,
        figure_size,
        fontsize,
        figure,
        axes,
        cax,
        shading="nearest",
        *,
        norm_handler: NormHandler,
        colorbar_handler: ColorbarHandler,
    ):
        self._initfinished = False
        self._shading = shading
        self._setup_layout_constraints(figure_size, fontsize)

        # this line is added purely to prevent exact image comparison tests
        # to fail, but eventually we should embrace the change and
        # use similar values for PhasePlotMPL and WindowPlotMPL
        self._ax_text_size[0] *= 1.1 / 1.2  # TODO: remove this

        super().__init__(
            figure=figure,
            axes=axes,
            cax=cax,
            norm_handler=norm_handler,
            colorbar_handler=colorbar_handler,
        )

        self._init_image(x_data, y_data, data, x_scale, y_scale)

        self._initfinished = True

    def _init_image(
        self,
        x_data,
        y_data,
        image_data,
        x_scale,
        y_scale,
    ):
        """Store output of imshow in image variable"""
        norm = self.norm_handler.get_norm(image_data)
        self.image = None
        self.cb = None

        self.image = self.axes.pcolormesh(
            np.array(x_data),
            np.array(y_data),
            np.array(image_data.T),
            norm=norm,
            cmap=self.colorbar_handler.cmap,
            shading=self._shading,
        )

        self._set_axes(norm)
        self.axes.set_xscale(x_scale)
        self.axes.set_yscale(y_scale)
