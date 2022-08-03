import abc
import itertools
import uuid
from collections import defaultdict
from contextlib import contextmanager

import numpy as np
from more_itertools import always_iterable
from unyt.exceptions import UnitConversionError, UnitParseError

import yt.geometry
from yt.data_objects.data_containers import YTDataContainer
from yt.data_objects.derived_quantities import DerivedQuantityCollection
from yt.data_objects.field_data import YTFieldData
from yt.fields.field_exceptions import NeedsGridType
from yt.funcs import fix_axis, is_sequence, iter_fields, validate_width_tuple
from yt.geometry.selection_routines import compose_selector
from yt.units import YTArray
from yt.utilities.exceptions import (
    GenerationInProgress,
    YTBooleanObjectError,
    YTBooleanObjectsWrongDataset,
    YTDataSelectorNotImplemented,
    YTDimensionalityError,
    YTFieldUnitError,
    YTFieldUnitParseError,
)
from yt.utilities.lib.marching_cubes import march_cubes_grid, march_cubes_grid_flux
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    ParallelAnalysisInterface,
)


class YTSelectionContainer(YTDataContainer, ParallelAnalysisInterface, abc.ABC):
    _locked = False
    _sort_by = None
    _selector = None
    _current_chunk = None
    _data_source = None
    _dimensionality: int
    _max_level = None
    _min_level = None
    _derived_quantity_chunking = "io"

    def __init__(self, ds, field_parameters, data_source=None):
        ParallelAnalysisInterface.__init__(self)
        super().__init__(ds, field_parameters)
        self._data_source = data_source
        if data_source is not None:
            if data_source.ds != self.ds:
                raise RuntimeError(
                    "Attempted to construct a DataContainer with a data_source "
                    "from a different Dataset",
                    ds,
                    data_source.ds,
                )
            if data_source._dimensionality < self._dimensionality:
                raise RuntimeError(
                    "Attempted to construct a DataContainer with a data_source "
                    "of lower dimensionality (%u vs %u)"
                    % (data_source._dimensionality, self._dimensionality)
                )
            self.field_parameters.update(data_source.field_parameters)
        self.quantities = DerivedQuantityCollection(self)

    @property
    def selector(self):
        if self._selector is not None:
            return self._selector
        s_module = getattr(self, "_selector_module", yt.geometry.selection_routines)
        sclass = getattr(s_module, f"{self._type_name}_selector", None)
        if sclass is None:
            raise YTDataSelectorNotImplemented(self._type_name)

        if self._data_source is not None:
            self._selector = compose_selector(
                self, self._data_source.selector, sclass(self)
            )
        else:
            self._selector = sclass(self)
        return self._selector

    def chunks(self, fields, chunking_style, **kwargs):
        # This is an iterator that will yield the necessary chunks.
        self.get_data()  # Ensure we have built ourselves
        if fields is None:
            fields = []
        # chunk_ind can be supplied in the keyword arguments.  If it's a
        # scalar, that'll be the only chunk that gets returned; if it's a list,
        # those are the ones that will be.
        chunk_ind = kwargs.pop("chunk_ind", None)
        if chunk_ind is not None:
            chunk_ind = list(always_iterable(chunk_ind))
        for ci, chunk in enumerate(self.index._chunk(self, chunking_style, **kwargs)):
            if chunk_ind is not None and ci not in chunk_ind:
                continue
            with self._chunked_read(chunk):
                self.get_data(fields)
                # NOTE: we yield before releasing the context
                yield self

    def _identify_dependencies(self, fields_to_get, spatial=False):
        inspected = 0
        fields_to_get = fields_to_get[:]
        for field in itertools.cycle(fields_to_get):
            if inspected >= len(fields_to_get):
                break
            inspected += 1
            fi = self.ds._get_field_info(*field)
            fd = self.ds.field_dependencies.get(
                field, None
            ) or self.ds.field_dependencies.get(field[1], None)
            # This is long overdue.  Any time we *can't* find a field
            # dependency -- for instance, if the derived field has been added
            # after dataset instantiation -- let's just try to
            # recalculate it.
            if fd is None:
                try:
                    fd = fi.get_dependencies(ds=self.ds)
                    self.ds.field_dependencies[field] = fd
                except Exception:
                    continue
            requested = self._determine_fields(list(set(fd.requested)))
            deps = [d for d in requested if d not in fields_to_get]
            fields_to_get += deps
        return sorted(fields_to_get)

    def get_data(self, fields=None):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        if fields is None:
            return
        nfields = []
        apply_fields = defaultdict(list)
        for field in self._determine_fields(fields):
            # We need to create the field on the raw particle types
            # for particles types (when the field is not directly
            # defined for the derived particle type only)
            finfo = self.ds.field_info[field]

            if (
                field[0] in self.ds.filtered_particle_types
                and finfo._inherited_particle_filter
            ):
                f = self.ds.known_filters[field[0]]
                apply_fields[field[0]].append((f.filtered_type, field[1]))
            else:
                nfields.append(field)
        for filter_type in apply_fields:
            f = self.ds.known_filters[filter_type]
            with f.apply(self):
                self.get_data(apply_fields[filter_type])
        fields = nfields
        if len(fields) == 0:
            return
        # Now we collect all our fields
        # Here is where we need to perform a validation step, so that if we
        # have a field requested that we actually *can't* yet get, we put it
        # off until the end.  This prevents double-reading fields that will
        # need to be used in spatial fields later on.
        fields_to_get = []
        # This will be pre-populated with spatial fields
        fields_to_generate = []
        for field in self._determine_fields(fields):
            if field in self.field_data:
                continue
            finfo = self.ds._get_field_info(*field)
            try:
                finfo.check_available(self)
            except NeedsGridType:
                fields_to_generate.append(field)
                continue
            fields_to_get.append(field)
        if len(fields_to_get) == 0 and len(fields_to_generate) == 0:
            return
        elif self._locked:
            raise GenerationInProgress(fields)
        # Track which ones we want in the end
        ofields = set(list(self.field_data.keys()) + fields_to_get + fields_to_generate)
        # At this point, we want to figure out *all* our dependencies.
        fields_to_get = self._identify_dependencies(fields_to_get, self._spatial)
        # We now split up into readers for the types of fields
        fluids, particles = [], []
        finfos = {}
        for ftype, fname in fields_to_get:
            finfo = self.ds._get_field_info(ftype, fname)
            finfos[ftype, fname] = finfo
            if finfo.sampling_type == "particle":
                particles.append((ftype, fname))
            elif (ftype, fname) not in fluids:
                fluids.append((ftype, fname))
        # The _read method will figure out which fields it needs to get from
        # disk, and return a dict of those fields along with the fields that
        # need to be generated.
        read_fluids, gen_fluids = self.index._read_fluid_fields(
            fluids, self, self._current_chunk
        )
        for f, v in read_fluids.items():
            self.field_data[f] = self.ds.arr(v, units=finfos[f].units)
            self.field_data[f].convert_to_units(finfos[f].output_units)

        read_particles, gen_particles = self.index._read_particle_fields(
            particles, self, self._current_chunk
        )

        for f, v in read_particles.items():
            self.field_data[f] = self.ds.arr(v, units=finfos[f].units)
            self.field_data[f].convert_to_units(finfos[f].output_units)

        fields_to_generate += gen_fluids + gen_particles
        self._generate_fields(fields_to_generate)
        for field in list(self.field_data.keys()):
            if field not in ofields:
                self.field_data.pop(field)

    def _generate_fields(self, fields_to_generate):
        index = 0
        with self._field_lock():
            # At this point, we assume that any fields that are necessary to
            # *generate* a field are in fact already available to us.  Note
            # that we do not make any assumption about whether or not the
            # fields have a spatial requirement.  This will be checked inside
            # _generate_field, at which point additional dependencies may
            # actually be noted.
            while any(f not in self.field_data for f in fields_to_generate):
                field = fields_to_generate[index % len(fields_to_generate)]
                index += 1
                if field in self.field_data:
                    continue
                fi = self.ds._get_field_info(*field)
                try:
                    fd = self._generate_field(field)
                    if hasattr(fd, "units"):
                        fd.units.registry = self.ds.unit_registry
                    if fd is None:
                        raise RuntimeError
                    if fi.units is None:
                        # first time calling a field with units='auto', so we
                        # infer the units from the units of the data we get back
                        # from the field function and use these units for future
                        # field accesses
                        units = getattr(fd, "units", "")
                        if units == "":
                            sunits = ""
                            dimensions = 1
                        else:
                            sunits = str(
                                units.get_base_equivalent(self.ds.unit_system.name)
                            )
                            dimensions = units.dimensions

                        if fi.dimensions is None:
                            mylog.warning(
                                "Field %s was added without specifying units or dimensions, "
                                "auto setting units to %s",
                                fi.name,
                                sunits,
                            )
                        elif fi.dimensions != dimensions:
                            raise YTDimensionalityError(fi.dimensions, dimensions)
                        fi.units = sunits
                        fi.dimensions = dimensions
                        self.field_data[field] = self.ds.arr(fd, units)
                    if fi.output_units is None:
                        fi.output_units = fi.units

                    try:
                        fd.convert_to_units(fi.units)
                    except AttributeError:
                        # If the field returns an ndarray, coerce to a
                        # dimensionless YTArray and verify that field is
                        # supposed to be unitless
                        fd = self.ds.arr(fd, "")
                        if fi.units != "":
                            raise YTFieldUnitError(fi, fd.units) from None
                    except UnitConversionError as e:
                        raise YTFieldUnitError(fi, fd.units) from e
                    except UnitParseError as e:
                        raise YTFieldUnitParseError(fi) from e
                    self.field_data[field] = fd
                except GenerationInProgress as gip:
                    for f in gip.fields:
                        if f not in fields_to_generate:
                            fields_to_generate.append(f)

    def __or__(self, other):
        if not isinstance(other, YTSelectionContainer):
            raise YTBooleanObjectError(other)
        if self.ds is not other.ds:
            raise YTBooleanObjectsWrongDataset()
        # Should maybe do something with field parameters here
        from yt.data_objects.selection_objects.boolean_operations import (
            YTBooleanContainer,
        )

        return YTBooleanContainer("OR", self, other, ds=self.ds)

    def __invert__(self):
        # ~obj
        asel = yt.geometry.selection_routines.AlwaysSelector(self.ds)
        from yt.data_objects.selection_objects.boolean_operations import (
            YTBooleanContainer,
        )

        return YTBooleanContainer("NOT", self, asel, ds=self.ds)

    def __xor__(self, other):
        if not isinstance(other, YTSelectionContainer):
            raise YTBooleanObjectError(other)
        if self.ds is not other.ds:
            raise YTBooleanObjectsWrongDataset()
        from yt.data_objects.selection_objects.boolean_operations import (
            YTBooleanContainer,
        )

        return YTBooleanContainer("XOR", self, other, ds=self.ds)

    def __and__(self, other):
        if not isinstance(other, YTSelectionContainer):
            raise YTBooleanObjectError(other)
        if self.ds is not other.ds:
            raise YTBooleanObjectsWrongDataset()
        from yt.data_objects.selection_objects.boolean_operations import (
            YTBooleanContainer,
        )

        return YTBooleanContainer("AND", self, other, ds=self.ds)

    def __add__(self, other):
        return self.__or__(other)

    def __sub__(self, other):
        if not isinstance(other, YTSelectionContainer):
            raise YTBooleanObjectError(other)
        if self.ds is not other.ds:
            raise YTBooleanObjectsWrongDataset()
        from yt.data_objects.selection_objects.boolean_operations import (
            YTBooleanContainer,
        )

        return YTBooleanContainer("NEG", self, other, ds=self.ds)

    @contextmanager
    def _field_lock(self):
        self._locked = True
        yield
        self._locked = False

    @contextmanager
    def _ds_hold(self, new_ds):
        """
        This contextmanager is used to take a data object and preserve its
        attributes but allow the dataset that underlies it to be swapped out.
        This is typically only used internally, and differences in unit systems
        may present interesting possibilities.
        """
        old_ds = self.ds
        old_index = self._index
        self.ds = new_ds
        self._index = new_ds.index
        old_chunk_info = self._chunk_info
        old_chunk = self._current_chunk
        old_size = self.size
        self._chunk_info = None
        self._current_chunk = None
        self.size = None
        self._index._identify_base_chunk(self)
        with self._chunked_read(None):
            yield
        self._index = old_index
        self.ds = old_ds
        self._chunk_info = old_chunk_info
        self._current_chunk = old_chunk
        self.size = old_size

    @contextmanager
    def _chunked_read(self, chunk):
        # There are several items that need to be swapped out
        # field_data, size, shape
        obj_field_data = []
        if hasattr(chunk, "objs"):
            for obj in chunk.objs:
                obj_field_data.append(obj.field_data)
                obj.field_data = YTFieldData()
        old_field_data, self.field_data = self.field_data, YTFieldData()
        old_chunk, self._current_chunk = self._current_chunk, chunk
        old_locked, self._locked = self._locked, False
        yield
        self.field_data = old_field_data
        self._current_chunk = old_chunk
        self._locked = old_locked
        if hasattr(chunk, "objs"):
            for obj in chunk.objs:
                obj.field_data = obj_field_data.pop(0)

    @contextmanager
    def _activate_cache(self):
        cache = self._field_cache or {}
        old_fields = {}
        for field in (f for f in cache if f in self.field_data):
            old_fields[field] = self.field_data[field]
        self.field_data.update(cache)
        yield
        for field in cache:
            self.field_data.pop(field)
            if field in old_fields:
                self.field_data[field] = old_fields.pop(field)
        self._field_cache = None

    def _initialize_cache(self, cache):
        # Wipe out what came before
        self._field_cache = {}
        self._field_cache.update(cache)

    @property
    def icoords(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.icoords

    @property
    def fcoords(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fcoords

    @property
    def ires(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.ires

    @property
    def fwidth(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fwidth

    @property
    def fcoords_vertex(self):
        if self._current_chunk is None:
            self.index._identify_base_chunk(self)
        return self._current_chunk.fcoords_vertex

    @property
    def max_level(self):
        if self._max_level is None:
            try:
                return self.ds.max_level
            except AttributeError:
                return None
        return self._max_level

    @max_level.setter
    def max_level(self, value):
        if self._selector is not None:
            del self._selector
            self._selector = None
        self._current_chunk = None
        self.size = None
        self.shape = None
        self.field_data.clear()
        self._max_level = value

    @property
    def min_level(self):
        if self._min_level is None:
            try:
                return 0
            except AttributeError:
                return None
        return self._min_level

    @min_level.setter
    def min_level(self, value):
        if self._selector is not None:
            del self._selector
            self._selector = None
        self.field_data.clear()
        self.size = None
        self.shape = None
        self._current_chunk = None
        self._min_level = value


class YTSelectionContainer0D(YTSelectionContainer):
    _spatial = False
    _dimensionality = 0

    def __init__(self, ds, field_parameters=None, data_source=None):
        super().__init__(ds, field_parameters, data_source)


class YTSelectionContainer1D(YTSelectionContainer):
    _spatial = False
    _dimensionality = 1

    def __init__(self, ds, field_parameters=None, data_source=None):
        super().__init__(ds, field_parameters, data_source)
        self._grids = None
        self._sortkey = None
        self._sorted = {}


class YTSelectionContainer2D(YTSelectionContainer):
    _key_fields = ["px", "py", "pdx", "pdy"]
    _dimensionality = 2
    """
    Prepares the YTSelectionContainer2D, normal to *axis*.  If *axis* is 4, we are not
    aligned with any axis.
    """
    _spatial = False

    def __init__(self, axis, ds, field_parameters=None, data_source=None):
        super().__init__(ds, field_parameters, data_source)
        # We need the ds, which will exist by now, for fix_axis.
        self.axis = fix_axis(axis, self.ds)
        self.set_field_parameter("axis", axis)

    def _convert_field_name(self, field):
        return field

    def _get_pw(self, fields, center, width, origin, plot_type):
        from yt.visualization.fixed_resolution import FixedResolutionBuffer as frb
        from yt.visualization.plot_window import PWViewerMPL, get_window_parameters

        axis = self.axis
        skip = self._key_fields
        skip += list(set(frb._exclude_fields).difference(set(self._key_fields)))
        self.fields = [k for k in self.field_data if k not in skip]
        if fields is not None:
            self.fields = list(iter_fields(fields)) + self.fields
        if len(self.fields) == 0:
            raise ValueError("No fields found to plot in get_pw")
        (bounds, center, display_center) = get_window_parameters(
            axis, center, width, self.ds
        )
        pw = PWViewerMPL(
            self,
            bounds,
            fields=self.fields,
            origin=origin,
            frb_generator=frb,
            plot_type=plot_type,
            geometry=self.ds.geometry,
        )
        pw._setup_plots()
        return pw

    def to_frb(self, width, resolution, center=None, height=None, periodic=False):
        r"""This function returns a FixedResolutionBuffer generated from this
        object.

        A FixedResolutionBuffer is an object that accepts a variable-resolution
        2D object and transforms it into an NxM bitmap that can be plotted,
        examined or processed.  This is a convenience function to return an FRB
        directly from an existing 2D data object.

        Parameters
        ----------
        width : width specifier
            This can either be a floating point value, in the native domain
            units of the simulation, or a tuple of the (value, unit) style.
            This will be the width of the FRB.
        height : height specifier
            This will be the physical height of the FRB, by default it is equal
            to width.  Note that this will not make any corrections to
            resolution for the aspect ratio.
        resolution : int or tuple of ints
            The number of pixels on a side of the final FRB.  If iterable, this
            will be the width then the height.
        center : array-like of floats, optional
            The center of the FRB.  If not specified, defaults to the center of
            the current object.
        periodic : bool
            Should the returned Fixed Resolution Buffer be periodic?  (default:
            False).

        Returns
        -------
        frb : :class:`~yt.visualization.fixed_resolution.FixedResolutionBuffer`
            A fixed resolution buffer, which can be queried for fields.

        Examples
        --------

        >>> proj = ds.proj(("gas", "density"), 0)
        >>> frb = proj.to_frb((100.0, "kpc"), 1024)
        >>> write_image(np.log10(frb[("gas", "density")]), "density_100kpc.png")
        """

        if (self.ds.geometry == "cylindrical" and self.axis == 1) or (
            self.ds.geometry == "polar" and self.axis == 2
        ):
            if center is not None and center != (0.0, 0.0):
                raise NotImplementedError(
                    "Currently we only support images centered at R=0. "
                    + "We plan to generalize this in the near future"
                )
            from yt.visualization.fixed_resolution import (
                CylindricalFixedResolutionBuffer,
            )

            validate_width_tuple(width)
            if is_sequence(resolution):
                resolution = max(resolution)
            frb = CylindricalFixedResolutionBuffer(self, width, resolution)
            return frb

        if center is None:
            center = self.center
            if center is None:
                center = (self.ds.domain_right_edge + self.ds.domain_left_edge) / 2.0
        elif is_sequence(center) and not isinstance(center, YTArray):
            center = self.ds.arr(center, "code_length")
        if is_sequence(width):
            w, u = width
            if isinstance(w, tuple) and isinstance(u, tuple):
                height = u
                w, u = w
            width = self.ds.quan(w, units=u)
        elif not isinstance(width, YTArray):
            width = self.ds.quan(width, "code_length")
        if height is None:
            height = width
        elif is_sequence(height):
            h, u = height
            height = self.ds.quan(h, units=u)
        elif not isinstance(height, YTArray):
            height = self.ds.quan(height, "code_length")
        if not is_sequence(resolution):
            resolution = (resolution, resolution)
        from yt.visualization.fixed_resolution import FixedResolutionBuffer

        xax = self.ds.coordinates.x_axis[self.axis]
        yax = self.ds.coordinates.y_axis[self.axis]
        bounds = (
            center[xax] - width * 0.5,
            center[xax] + width * 0.5,
            center[yax] - height * 0.5,
            center[yax] + height * 0.5,
        )
        frb = FixedResolutionBuffer(self, bounds, resolution, periodic=periodic)
        return frb


class YTSelectionContainer3D(YTSelectionContainer):
    """
    Returns an instance of YTSelectionContainer3D, or prepares one.  Usually only
    used as a base class.  Note that *center* is supplied, but only used
    for fields and quantities that require it.
    """

    _key_fields = ["x", "y", "z", "dx", "dy", "dz"]
    _spatial = False
    _num_ghost_zones = 0
    _dimensionality = 3

    def __init__(self, center, ds, field_parameters=None, data_source=None):
        super().__init__(ds, field_parameters, data_source)
        self._set_center(center)
        self.coords = None
        self._grids = None

    def cut_region(self, field_cuts, field_parameters=None, locals=None):
        """
        Return a YTCutRegion, where the a cell is identified as being inside
        the cut region based on the value of one or more fields.  Note that in
        previous versions of yt the name 'grid' was used to represent the data
        object used to construct the field cut, as of yt 3.0, this has been
        changed to 'obj'.

        Parameters
        ----------
        field_cuts : list of strings
           A list of conditionals that will be evaluated. In the namespace
           available, these conditionals will have access to 'obj' which is a
           data object of unknown shape, and they must generate a boolean array.
           For instance, conditionals = ["obj[('gas', 'temperature')] < 1e3"]
        field_parameters : dictionary
           A dictionary of field parameters to be used when applying the field
           cuts.
        locals : dictionary
            A dictionary of local variables to use when defining the cut region.

        Examples
        --------
        To find the total mass of hot gas with temperature greater than 10^6 K
        in your volume:

        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.cut_region(["obj[('gas', 'temperature')] > 1e6"])
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        if locals is None:
            locals = {}
        cr = self.ds.cut_region(
            self, field_cuts, field_parameters=field_parameters, locals=locals
        )
        return cr

    def _build_operator_cut(self, operation, field, value, units=None):
        """
        Given an operation (>, >=, etc.), a field and a value,
        return the cut_region implementing it.

        This is only meant to be used internally.

        Examples
        --------
        >>> ds._build_operator_cut(">", ("gas", "density"), 1e-24)
        ... # is equivalent to
        ... ds.cut_region(['obj[("gas", "density")] > 1e-24'])
        """
        ftype, fname = self._determine_fields(field)[0]
        if units is None:
            field_cuts = f'obj["{ftype}", "{fname}"] {operation} {value}'
        else:
            field_cuts = (
                f'obj["{ftype}", "{fname}"].in_units("{units}") {operation} {value}'
            )
        return self.cut_region(field_cuts)

    def _build_function_cut(self, function, field, units=None, **kwargs):
        """
        Given a function (np.abs, np.all) and a field,
        return the cut_region implementing it.

        This is only meant to be used internally.

        Examples
        --------
        >>> ds._build_function_cut("np.isnan", ("gas", "density"), locals={"np": np})
        ... # is equivalent to
        ... ds.cut_region(['np.isnan(obj[("gas", "density")])'], locals={"np": np})
        """
        ftype, fname = self._determine_fields(field)[0]
        if units is None:
            field_cuts = f'{function}(obj["{ftype}", "{fname}"])'
        else:
            field_cuts = f'{function}(obj["{ftype}", "{fname}"].in_units("{units}"))'
        return self.cut_region(field_cuts, **kwargs)

    def exclude_above(self, field, value, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field is above a given value are masked.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field above the given value masked.

        Examples
        --------
        To find the total mass of hot gas with temperature colder than 10^6 K
        in your volume:

        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_above(("gas", "temperature"), 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))

        """
        return self._build_operator_cut("<=", field, value, units)

    def include_above(self, field, value, units=None):
        """
        This function will return a YTCutRegion where only the regions
        whose field is above a given value are included.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field above the given value masked.

        Examples
        --------
        To find the total mass of hot gas with temperature warmer than 10^6 K
        in your volume:

        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.include_above(("gas", "temperature"), 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """

        return self._build_operator_cut(">", field, value, units)

    def exclude_equal(self, field, value, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field are equal to given value are masked.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field equal to the given value masked.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_equal(("gas", "temperature"), 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        return self._build_operator_cut("!=", field, value, units)

    def include_equal(self, field, value, units=None):
        """
        This function will return a YTCutRegion where only the regions
        whose field are equal to given value are included.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field equal to the given value included.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.include_equal(("gas", "temperature"), 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        return self._build_operator_cut("==", field, value, units)

    def exclude_inside(self, field, min_value, max_value, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field are inside the interval from min_value to max_value.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        min_value : float
            The minimum value inside the interval to be excluded.
        max_value : float
            The maximum value inside the interval to be excluded.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field inside the given interval excluded.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_inside(("gas", "temperature"), 1e5, 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        ftype, fname = self._determine_fields(field)[0]
        if units is None:
            field_cuts = (
                f'(obj["{ftype}", "{fname}"] <= {min_value}) | '
                f'(obj["{ftype}", "{fname}"] >= {max_value})'
            )
        else:
            field_cuts = (
                f'(obj["{ftype}", "{fname}"].in_units("{units}") <= {min_value}) | '
                f'(obj["{ftype}", "{fname}"].in_units("{units}") >= {max_value})'
            )
        cr = self.cut_region(field_cuts)
        return cr

    def include_inside(self, field, min_value, max_value, units=None):
        """
        This function will return a YTCutRegion where only the regions
        whose field are inside the interval from min_value to max_value are
        included.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        min_value : float
            The minimum value inside the interval to be excluded.
        max_value : float
            The maximum value inside the interval to be excluded.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field inside the given interval excluded.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.include_inside(("gas", "temperature"), 1e5, 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        ftype, fname = self._determine_fields(field)[0]
        if units is None:
            field_cuts = (
                f'(obj["{ftype}", "{fname}"] > {min_value}) & '
                f'(obj["{ftype}", "{fname}"] < {max_value})'
            )
        else:
            field_cuts = (
                f'(obj["{ftype}", "{fname}"].in_units("{units}") > {min_value}) & '
                f'(obj["{ftype}", "{fname}"].in_units("{units}") < {max_value})'
            )
        cr = self.cut_region(field_cuts)
        return cr

    def exclude_outside(self, field, min_value, max_value, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field are outside the interval from min_value to max_value.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        min_value : float
            The minimum value inside the interval to be excluded.
        max_value : float
            The maximum value inside the interval to be excluded.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field outside the given interval excluded.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_outside(("gas", "temperature"), 1e5, 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        cr = self.exclude_below(field, min_value, units)
        cr = cr.exclude_above(field, max_value, units)
        return cr

    def include_outside(self, field, min_value, max_value, units=None):
        """
        This function will return a YTCutRegion where only the regions
        whose field are outside the interval from min_value to max_value are
        included.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        min_value : float
            The minimum value inside the interval to be excluded.
        max_value : float
            The maximum value inside the interval to be excluded.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field outside the given interval excluded.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_outside(("gas", "temperature"), 1e5, 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        cr = self.exclude_inside(field, min_value, max_value, units)
        return cr

    def exclude_below(self, field, value, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field is below a given value are masked.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the field below the given value masked.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_below(("gas", "temperature"), 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        return self._build_operator_cut(">=", field, value, units)

    def exclude_nan(self, field, units=None):
        """
        This function will return a YTCutRegion where all of the regions
        whose field is NaN are masked.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with the NaN entries of the field masked.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.exclude_nan(("gas", "temperature"))
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        return self._build_function_cut("~np.isnan", field, units, locals={"np": np})

    def include_below(self, field, value, units=None):
        """
        This function will return a YTCutRegion where only the regions
        whose field is below a given value are included.

        Parameters
        ----------
        field : string
            The field in which the conditional will be applied.
        value : float
            The minimum value that will not be masked in the output
            YTCutRegion.
        units : string or None
            The units of the value threshold. None will use the default units
            given in the field.

        Returns
        -------
        cut_region : YTCutRegion
            The YTCutRegion with only regions with the field below the given
            value included.

        Examples
        --------
        >>> ds = yt.load("RedshiftOutput0005")
        >>> ad = ds.all_data()
        >>> cr = ad.include_below(("gas", "temperature"), 1e5, 1e6)
        >>> print(cr.quantities.total_quantity(("gas", "cell_mass")).in_units("Msun"))
        """
        return self._build_operator_cut("<", field, value, units)

    def extract_isocontours(
        self, field, value, filename=None, rescale=False, sample_values=None
    ):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and returns the vertices of the
        Triangles in that isocontour.

        This function simply returns the vertices of all the triangles
        calculated by the `marching cubes
        <https://en.wikipedia.org/wiki/Marching_cubes>`_ algorithm; for more
        complex operations, such as identifying connected sets of cells above a
        given threshold, see the extract_connected_sets function.  This is more
        useful for calculating, for instance, total isocontour area, or
        visualizing in an external program (such as `MeshLab
        <http://www.meshlab.net>`_.)

        Parameters
        ----------
        field : string
            Any field that can be obtained in a data object.  This is the field
            which will be isocontoured.
        value : float
            The value at which the isocontour should be calculated.
        filename : string, optional
            If supplied, this file will be filled with the vertices in .obj
            format.  Suitable for loading into meshlab.
        rescale : bool, optional
            If true, the vertices will be rescaled within their min/max.
        sample_values : string, optional
            Any field whose value should be extracted at the center of each
            triangle.

        Returns
        -------
        verts : array of floats
            The array of vertices, x,y,z.  Taken in threes, these are the
            triangle vertices.
        samples : array of floats
            If `sample_values` is specified, this will be returned and will
            contain the values of the field specified at the center of each
            triangle.

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        output the vertices to "triangles.obj" after rescaling them.

        >>> dd = ds.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     ("gas", "density"), weight=("gas", "cell_mass")
        ... )
        >>> verts = dd.extract_isocontours(
        ...     ("gas", "density"), rho, "triangles.obj", True
        ... )
        """
        from yt.data_objects.static_output import ParticleDataset
        from yt.frontends.stream.data_structures import StreamParticlesDataset

        verts = []
        samples = []
        if isinstance(self.ds, (ParticleDataset, StreamParticlesDataset)):
            raise NotImplementedError
        for block, mask in self.blocks:
            my_verts = self._extract_isocontours_from_grid(
                block, mask, field, value, sample_values
            )
            if sample_values is not None:
                my_verts, svals = my_verts
                samples.append(svals)
            verts.append(my_verts)
        verts = np.concatenate(verts).transpose()
        verts = self.comm.par_combine_object(verts, op="cat", datatype="array")
        verts = verts.transpose()
        if sample_values is not None:
            samples = np.concatenate(samples)
            samples = self.comm.par_combine_object(samples, op="cat", datatype="array")
        if rescale:
            mi = np.min(verts, axis=0)
            ma = np.max(verts, axis=0)
            verts = (verts - mi) / (ma - mi).max()
        if filename is not None and self.comm.rank == 0:
            if hasattr(filename, "write"):
                f = filename
            else:
                f = open(filename, "w")
            for v1 in verts:
                f.write(f"v {v1[0]:0.16e} {v1[1]:0.16e} {v1[2]:0.16e}\n")
            for i in range(len(verts) // 3):
                f.write(f"f {i * 3 + 1} {i * 3 + 2} {i * 3 + 3}\n")
            if not hasattr(filename, "write"):
                f.close()
        if sample_values is not None:
            return verts, samples
        return verts

    def _extract_isocontours_from_grid(
        self, grid, mask, field, value, sample_values=None
    ):
        vc_fields = [field]
        if sample_values is not None:
            vc_fields.append(sample_values)

        vc_data = grid.get_vertex_centered_data(vc_fields, no_ghost=False)
        try:
            svals = vc_data[sample_values]
        except KeyError:
            svals = None

        my_verts = march_cubes_grid(
            value, vc_data[field], mask, grid.LeftEdge, grid.dds, svals
        )
        return my_verts

    def calculate_isocontour_flux(
        self, field, value, field_x, field_y, field_z, fluxing_field=None
    ):
        r"""This identifies isocontours on a cell-by-cell basis, with no
        consideration of global connectedness, and calculates the flux over
        those contours.

        This function will conduct `marching cubes
        <https://en.wikipedia.org/wiki/Marching_cubes>`_ on all the cells in a
        given data container (grid-by-grid), and then for each identified
        triangular segment of an isocontour in a given cell, calculate the
        gradient (i.e., normal) in the isocontoured field, interpolate the local
        value of the "fluxing" field, the area of the triangle, and then return:

        area * local_flux_value * (n dot v)

        Where area, local_value, and the vector v are interpolated at the barycenter
        (weighted by the vertex values) of the triangle.  Note that this
        specifically allows for the field fluxing across the surface to be
        *different* from the field being contoured.  If the fluxing_field is
        not specified, it is assumed to be 1.0 everywhere, and the raw flux
        with no local-weighting is returned.

        Additionally, the returned flux is defined as flux *into* the surface,
        not flux *out of* the surface.

        Parameters
        ----------
        field : string
            Any field that can be obtained in a data object.  This is the field
            which will be isocontoured and used as the "local_value" in the
            flux equation.
        value : float
            The value at which the isocontour should be calculated.
        field_x : string
            The x-component field
        field_y : string
            The y-component field
        field_z : string
            The z-component field
        fluxing_field : string, optional
            The field whose passage over the surface is of interest.  If not
            specified, assumed to be 1.0 everywhere.

        Returns
        -------
        flux : float
            The summed flux.  Note that it is not currently scaled; this is
            simply the code-unit area times the fields.

        Examples
        --------
        This will create a data object, find a nice value in the center, and
        calculate the metal flux over it.

        >>> dd = ds.all_data()
        >>> rho = dd.quantities["WeightedAverageQuantity"](
        ...     ("gas", "density"), weight=("gas", "cell_mass")
        ... )
        >>> flux = dd.calculate_isocontour_flux(
        ...     ("gas", "density"),
        ...     rho,
        ...     ("gas", "velocity_x"),
        ...     ("gas", "velocity_y"),
        ...     ("gas", "velocity_z"),
        ...     ("gas", "metallicity"),
        ... )
        """
        flux = 0.0
        for block, mask in self.blocks:
            flux += self._calculate_flux_in_grid(
                block, mask, field, value, field_x, field_y, field_z, fluxing_field
            )
        flux = self.comm.mpi_allreduce(flux, op="sum")
        return flux

    def _calculate_flux_in_grid(
        self, grid, mask, field, value, field_x, field_y, field_z, fluxing_field=None
    ):

        vc_fields = [field, field_x, field_y, field_z]
        if fluxing_field is not None:
            vc_fields.append(fluxing_field)

        vc_data = grid.get_vertex_centered_data(vc_fields)

        if fluxing_field is None:
            ff = np.ones_like(vc_data[field], dtype="float64")
        else:
            ff = vc_data[fluxing_field]

        return march_cubes_grid_flux(
            value,
            vc_data[field],
            vc_data[field_x],
            vc_data[field_y],
            vc_data[field_z],
            ff,
            mask,
            grid.LeftEdge,
            grid.dds,
        )

    def extract_connected_sets(
        self, field, num_levels, min_val, max_val, log_space=True, cumulative=True
    ):
        """
        This function will create a set of contour objects, defined
        by having connected cell structures, which can then be
        studied and used to 'paint' their source grids, thus enabling
        them to be plotted.

        Note that this function *can* return a connected set object that has no
        member values.
        """
        if log_space:
            cons = np.logspace(np.log10(min_val), np.log10(max_val), num_levels + 1)
        else:
            cons = np.linspace(min_val, max_val, num_levels + 1)
        contours = {}
        for level in range(num_levels):
            contours[level] = {}
            if cumulative:
                mv = max_val
            else:
                mv = cons[level + 1]
            from yt.data_objects.level_sets.api import identify_contours
            from yt.data_objects.level_sets.clump_handling import add_contour_field

            nj, cids = identify_contours(self, field, cons[level], mv)
            unique_contours = set()
            for sl_list in cids.values():
                for _sl, ff in sl_list:
                    unique_contours.update(np.unique(ff))
            contour_key = uuid.uuid4().hex
            # In case we're a cut region already...
            base_object = getattr(self, "base_object", self)
            add_contour_field(base_object.ds, contour_key)
            for cid in sorted(unique_contours):
                if cid == -1:
                    continue
                contours[level][cid] = base_object.cut_region(
                    [f"obj['contours_{contour_key}'] == {cid}"],
                    {f"contour_slices_{contour_key}": cids},
                )
        return cons, contours

    def _get_bbox(self):
        """
        Return the bounding box for this data container.
        This generic version will return the bounds of the entire domain.
        """
        return self.ds.domain_left_edge, self.ds.domain_right_edge

    def get_bbox(self):
        """
        Return the bounding box for this data container.
        """
        if self.ds.geometry != "cartesian":
            raise NotImplementedError(
                "get_bbox is currently only implemented for cartesian geometries!"
            )
        le, re = self._get_bbox()
        le.convert_to_units("code_length")
        re.convert_to_units("code_length")
        return le, re

    def volume(self):
        """
        Return the volume of the data container.
        This is found by adding up the volume of the cells with centers
        in the container, rather than using the geometric shape of
        the container, so this may vary very slightly
        from what might be expected from the geometric volume.
        """
        return self.quantities.total_quantity(("index", "cell_volume"))
