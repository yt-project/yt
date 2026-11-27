"""
AMRVAC data structures



"""

import os
import struct
import warnings
import weakref
from pathlib import Path

import numpy as np
from more_itertools import always_iterable

from yt.config import ytcfg
from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.index_subobjects.stretched_grid import StretchedGrid
from yt.data_objects.static_output import Dataset, _unsupported_object
from yt.funcs import mylog, setdefaultattr
from yt.geometry.api import Geometry
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.physical_constants import boltzmann_constant_cgs as kb_cgs

from .datfile_utils import get_header, get_tree_info
from .fields import AMRVACFieldInfo
from .io import read_amrvac_namelist


def _parse_geometry(geometry_tag: str) -> Geometry:
    """Translate AMRVAC's geometry tag to yt's format.

    Parameters
    ----------
    geometry_tag : str
        A geometry tag as read from AMRVAC's datfile from v5.

    Returns
    -------
    geometry_yt : Geometry
        An enum member of the yt.geometry.geometry_enum.Geometry class

    Examples
    --------
    >>> _parse_geometry("Polar_2.5D")
     <Geometry.POLAR: 'polar'>
    >>> _parse_geometry("Cartesian_2.5D")
    <Geometry.CARTESIAN: 'cartesian'>

    """
    geometry_str, _, _dimension_str = geometry_tag.partition("_")
    return Geometry(geometry_str.lower())


# Stretching types as written in AMRVAC parfiles (&meshlist:stretch_dim)
_NO_STRETCH = ("", "none")
_UNIFORM_STRETCH = ("uni", "uniform")
_SYMMETRIC_STRETCH = ("symm", "symmetric")


def _parse_stretch_dim(meshlist, ndim: int) -> list[str]:
    """Return the stretching type ("none" or "uni") for each dimension.

    Accepts a scalar (``stretch_dim = 'uni'``, i.e. first dimension only, following
    Fortran semantics), an array definition (``stretch_dim = 'uni', '', ''``) or
    indexed assignments (``stretch_dim(2) = 'uni'``, parsed as ``[None, 'uni']``).
    """
    stretch_dim = ["none"] * ndim
    raw = meshlist.get("stretch_dim")
    if raw is None:
        return stretch_dim

    for idim, value in enumerate(always_iterable(raw)):
        if value is None:
            continue
        if not isinstance(value, str):
            raise ValueError(
                f"Expected a string for stretch_dim({idim + 1}), got {value!r}."
            )
        value = value.strip().lower()
        if value in _NO_STRETCH:
            continue
        if idim >= ndim:
            raise ValueError(
                f"stretch_dim({idim + 1}) = {value!r} is set, "
                f"but the dataset only has {ndim} dimension(s)."
            )
        if value in _UNIFORM_STRETCH:
            stretch_dim[idim] = "uni"
        elif value in _SYMMETRIC_STRETCH:
            raise NotImplementedError(
                "Symmetric stretching is not currently supported for AMRVAC data."
            )
        else:
            raise ValueError(f"Unknown stretch_dim({idim + 1}) = {value!r}.")
    return stretch_dim


def _parse_qstretch_baselevel(meshlist, ndim: int) -> list[float | None]:
    """Return the user-supplied qstretch_baselevel for each dimension (or None).

    A scalar applies to the first dimension only, following Fortran semantics.
    """
    qstretch_baselevel: list[float | None] = [None] * ndim
    raw = meshlist.get("qstretch_baselevel")
    if raw is None:
        return qstretch_baselevel

    for idim, value in enumerate(always_iterable(raw)):
        if idim >= ndim or value is None:
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(
                f"Expected a number for qstretch_baselevel({idim + 1}), got {value!r}."
            )
        qstretch_baselevel[idim] = float(value)
    return qstretch_baselevel


def _stretched_cell_edges(xmin, dxfirst, qstretch, k):
    """Left edges of cells with integer indices ``k`` on a uniformly stretched axis.

    Cell ``k`` has width ``dxfirst * qstretch**k``, so its left edge is the geometric
    sum ``xmin + dxfirst * (1 - qstretch**k) / (1 - qstretch)``, as in AMRVAC.
    All arguments may be arrays that broadcast together.
    """
    return xmin + dxfirst * (1.0 - qstretch**k) / (1.0 - qstretch)


class AMRVACGrid(AMRGridPatch):
    """A class to populate AMRVACHierarchy.grids, setting parent/children relations."""

    _id_offset = 0

    def __init__(self, id, index, level, morton_index=None):
        # <level> should use yt's convention (start from 0)
        super().__init__(id, filename=index.index_filename, index=index)
        self._setup_amrvac_attributes(level, morton_index)

    def _setup_amrvac_attributes(self, level, morton_index):
        self.Parent = None
        self.Children = []
        self.Level = level
        # AMRVAC's (1-based) block index along each dimension, at this grid's level
        self._morton_index = morton_index

    def get_global_startindex(self):
        """Refresh and retrieve the starting index for each dimension at current level.

        Returns
        -------
        self.start_index : int
        """
        if self._morton_index is None:
            start_index = (self.LeftEdge - self.ds.domain_left_edge) / self.dds
            self.start_index = np.rint(start_index).astype("int64").ravel()
            return self.start_index

        # exact, and independent of cell widths (valid for stretched grids)
        ndim = self.ds.dimensionality
        start_index = np.zeros(3, dtype="int64")
        start_index[:ndim] = (np.asarray(self._morton_index) - 1) * np.asarray(
            self.ActiveDimensions[:ndim]
        )
        self.start_index = start_index
        return self.start_index

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False, smoothed=False):
        if smoothed:
            warnings.warn(
                "ghost-zones interpolation/smoothing is not "
                "currently supported for AMRVAC data.",
                category=RuntimeWarning,
                stacklevel=2,
            )
            smoothed = False
        return super().retrieve_ghost_zones(
            n_zones, fields, all_levels=all_levels, smoothed=smoothed
        )


class AMRVACStretchedGrid(AMRVACGrid, StretchedGrid):
    """An AMRVAC grid with non-uniform cell widths along at least one dimension."""

    def __init__(self, id, index, level, morton_index, cell_widths):
        # skip AMRVACGrid.__init__: StretchedGrid needs the cell widths
        StretchedGrid.__init__(
            self, id, cell_widths, filename=index.index_filename, index=index
        )
        self._setup_amrvac_attributes(level, morton_index)

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False, smoothed=False):
        raise NotImplementedError(
            "Ghost zones are not currently supported for stretched AMRVAC grids."
        )


class AMRVACHierarchy(GridIndex):
    grid = AMRVACGrid

    def __init__(self, ds, dataset_type="amrvac"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # the index file *is* the datfile
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64

        ndim = self.dataset.dimensionality
        meshlist = {}
        if self.dataset.namelist is not None:
            meshlist = self.dataset.namelist.get("meshlist", {})
        self.stretch_dim = _parse_stretch_dim(meshlist, ndim)
        self._qstretch_baselevel = _parse_qstretch_baselevel(meshlist, ndim)
        self.is_stretched = any(s != "none" for s in self.stretch_dim)
        if self.is_stretched:
            self.grid = AMRVACStretchedGrid

        super().__init__(ds, dataset_type)

        if self.is_stretched:
            # these rely on uniform cell widths
            for name in ("covering_grid", "smoothed_covering_grid"):
                setattr(self.dataset, name, _unsupported_object(self.dataset, name))

    def _detect_output_fields(self):
        """Parse field names from the header, as stored in self.dataset.parameters"""
        self.field_list = [
            (self.dataset_type, f) for f in self.dataset.parameters["w_names"]
        ]

    def _count_grids(self):
        """Set self.num_grids from datfile header."""
        self.num_grids = self.dataset.parameters["nleafs"]

    def _parse_index(self):
        """Populate self.grid_* attributes from tree info from datfile header."""
        with open(self.index_filename, "rb") as istream:
            vaclevels, morton_indices, block_offsets = get_tree_info(istream)
            assert (
                len(vaclevels)
                == len(morton_indices)
                == len(block_offsets)
                == self.num_grids
            )

        self.block_offsets = block_offsets
        # YT uses 0-based grid indexing:
        # lowest level = 0, while AMRVAC uses 1 for lowest level
        ytlevels = np.array(vaclevels, dtype="int32") - 1
        self.grid_levels.flat[:] = ytlevels
        self.min_level = np.min(ytlevels)
        self.max_level = np.max(ytlevels)
        assert self.max_level == self.dataset.parameters["levmax"] - 1

        if self.is_stretched:
            self._setup_stretching()

        block_nx = self.dataset.parameters["block_nx"]
        ndim = self.dataset.dimensionality

        self.grids = np.empty(self.num_grids, dtype="object")
        for igrid, (ytlevel, morton_index) in enumerate(
            zip(ytlevels, morton_indices, strict=True)
        ):
            # edges and dimensions are filled in a dimensionality-agnostic way
            cell_widths = []
            for idim in range(ndim):
                k = (morton_index[idim] - 1) * block_nx[idim] + np.arange(
                    block_nx[idim] + 1
                )
                edges = self._cell_edges(idim, ytlevel, k)
                # neighbouring grids compute their shared edge identically
                self.grid_left_edge[igrid, idim] = edges[0]
                self.grid_right_edge[igrid, idim] = edges[-1]
                cell_widths.append(np.diff(edges))
            self.grid_dimensions[igrid, :ndim] = block_nx

            if not self.is_stretched:
                self.grids[igrid] = self.grid(igrid, self, ytlevel, morton_index)
                continue

            for idim in range(ndim, 3):
                # unused dimensions: a single cell spanning the domain
                width = (
                    self.grid_right_edge[igrid, idim].d
                    - self.grid_left_edge[igrid, idim].d
                )
                cell_widths.append(np.array([width]))
            self.grids[igrid] = self.grid(
                igrid, self, ytlevel, morton_index, cell_widths=cell_widths
            )

    def _setup_stretching(self):
        """Compute AMRVAC's qstretch and dxfirst for every level and dimension.

        This mirrors AMRVAC's own setup: arrays are indexed by AMRVAC level
        (1 = coarsest), and only stretched dimensions get non-zero values.
        """
        ndim = self.dataset.dimensionality
        xmin = self.dataset.parameters["xmin"]
        xmax = self.dataset.parameters["xmax"]
        domain_nx = self.dataset.parameters["domain_nx"]
        nlevels = self.max_level + 2  # AMRVAC levels 0 to levmax

        self._qstretch = np.zeros((nlevels, ndim), dtype="float64")
        self._dxfirst = np.zeros((nlevels, ndim), dtype="float64")
        for idim in range(ndim):
            if self.stretch_dim[idim] == "none":
                continue
            q = self._qstretch_baselevel[idim]
            if q is None or q < 1.0:
                # same fallback as AMRVAC
                if xmin[idim] <= 0.0:
                    raise ValueError(
                        f"qstretch_baselevel({idim + 1}) must be set (> 1) in the "
                        "parfile when the stretched dimension has xprobmin <= 0."
                    )
                q = (xmax[idim] / xmin[idim]) ** (1.0 / domain_nx[idim])
            if q == 1.0:
                raise ValueError(
                    f"qstretch_baselevel({idim + 1}) = 1 means no stretching; "
                    f"use stretch_dim({idim + 1}) = 'none' instead."
                )
            q_lev = self._qstretch[:, idim]
            dx_lev = self._dxfirst[:, idim]
            q_lev[1] = q
            dx_lev[1] = (
                (xmax[idim] - xmin[idim]) * (1.0 - q) / (1.0 - q ** domain_nx[idim])
            )
            q_lev[0] = q**2
            dx_lev[0] = dx_lev[1] * (1.0 + q)
            for ilev in range(2, nlevels):
                q_lev[ilev] = np.sqrt(q_lev[ilev - 1])
                dx_lev[ilev] = dx_lev[ilev - 1] / (1.0 + q_lev[ilev])

    def _cell_edges(self, idim, ytlevel, k):
        """Left edge of cells with global integer index ``k`` at yt level ``ytlevel``.

        ``ytlevel`` and ``k`` may be arrays that broadcast together. Cell widths are
        differences of consecutive edges.
        """
        xmin = self.dataset.parameters["xmin"][idim]
        if self.stretch_dim[idim] == "none":
            dx0 = (
                self.dataset.parameters["xmax"][idim] - xmin
            ) / self.dataset.parameters["domain_nx"][idim]
            return xmin + k * (dx0 / self.dataset.refine_by ** np.asarray(ytlevel))
        amrvac_level = np.asarray(ytlevel) + 1
        return _stretched_cell_edges(
            xmin,
            self._dxfirst[amrvac_level, idim],
            self._qstretch[amrvac_level, idim],
            k,
        )

    def _icoords_to_fcoords(self, icoords, ires, axes=None):
        """Cell centers and widths from integer coordinates, for stretched grids."""
        if not self.is_stretched:
            return super()._icoords_to_fcoords(icoords, ires, axes)
        if axes is None:
            axes = (0, 1, 2)
        coords = np.empty(icoords.shape, dtype="float64")
        widths = np.empty(icoords.shape, dtype="float64")
        ires = np.asarray(ires).ravel()
        for i, ax in enumerate(axes):
            if ax >= self.dataset.dimensionality:
                # unused dimension: a single cell spanning the domain
                le = self.dataset.domain_left_edge.d[ax]
                re = self.dataset.domain_right_edge.d[ax]
                coords[:, i] = 0.5 * (le + re)
                widths[:, i] = re - le
                continue
            left = self._cell_edges(ax, ires, icoords[:, i])
            right = self._cell_edges(ax, ires, icoords[:, i] + 1)
            coords[:, i] = 0.5 * (left + right)
            widths[:, i] = right - left
        return coords, widths

    def _populate_grid_objects(self):
        # required method
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()


class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="amrvac",
        units_override=None,
        unit_system="cgs",
        geometry_override=None,
        parfiles=None,
        default_species_fields=None,
    ):
        """Instantiate AMRVACDataset.

        Parameters
        ----------
        filename : str
            Path to a datfile.

        dataset_type : str, optional
            This should always be 'amrvac'.

        units_override : dict, optional
            A dictionary of physical normalisation factors to interpret on disk data.

        unit_system : str, optional
            Either "cgs" (default), "mks" or "code"

        geometry_override : str, optional
            A geometry flag formatted either according to either
            AMRVAC or yt standards.
            When this parameter is passed along with v5 or more newer datfiles,
            will precede over their internal "geometry" tag.

        parfiles : str or list, optional
            One or more parfiles to be passed to
            yt.frontends.amrvac.read_amrvac_parfiles()

        """
        # note: geometry_override and parfiles are specific to this frontend

        self._geometry_override = geometry_override
        self._parfiles = []

        if parfiles is not None:
            parfiles = list(always_iterable(parfiles))
            ppf = Path(parfiles[0])
            if not ppf.is_absolute() and Path(filename).resolve().is_relative_to(
                ytcfg["yt", "test_data_dir"]
            ):
                mylog.debug(
                    "Looks like %s is relative to your test_data_dir. Assuming this is also true for parfiles.",
                    filename,
                )
                parfiles = [Path(ytcfg["yt", "test_data_dir"]) / pf for pf in parfiles]
            self._parfiles = parfiles

        super().__init__(
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

        namelist = None
        namelist_gamma = None
        c_adiab = None
        e_is_internal = None

        if parfiles is not None:
            namelist = read_amrvac_namelist(parfiles)
            if "hd_list" in namelist:
                c_adiab = namelist["hd_list"].get("hd_adiab", 1.0)
                namelist_gamma = namelist["hd_list"].get("hd_gamma")
            elif "mhd_list" in namelist:
                c_adiab = namelist["mhd_list"].get("mhd_adiab", 1.0)
                namelist_gamma = namelist["mhd_list"].get("mhd_gamma")

            if "eos_list" in namelist:
                namelist_gamma = namelist["eos_list"].get("gamma")

            if namelist_gamma is not None and self.gamma != namelist_gamma:
                mylog.error(
                    "Inconsistent values in gamma: datfile %s, parfiles %s",
                    self.gamma,
                    namelist_gamma,
                )

            if "method_list" in namelist:
                e_is_internal = namelist["method_list"].get("solve_internal_e", False)

        if c_adiab is not None:
            # this complicated unit is required for the adiabatic equation
            # of state to make physical sense
            c_adiab *= (
                self.mass_unit ** (1 - self.gamma)
                * self.length_unit ** (2 + 3 * (self.gamma - 1))
                / self.time_unit**2
            )

        self.namelist = namelist
        self._c_adiab = c_adiab
        self._e_is_internal = e_is_internal

        self.fluid_types += ("amrvac",)
        # refinement factor between a grid and its subgrid
        self.refine_by = 2

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        """At load time, check whether data is recognized as AMRVAC formatted."""
        validation = False
        if filename.endswith(".dat"):
            try:
                with open(filename, mode="rb") as istream:
                    fmt = "=i"
                    [datfile_version] = struct.unpack(
                        fmt, istream.read(struct.calcsize(fmt))
                    )
                    if 3 <= datfile_version < 6:
                        fmt = "=ii"
                        offset_tree, offset_blocks = struct.unpack(
                            fmt, istream.read(struct.calcsize(fmt))
                        )
                        istream.seek(0, 2)
                        file_size = istream.tell()
                        validation = (
                            offset_tree < file_size and offset_blocks < file_size
                        )
            except Exception:
                pass
        return validation

    def _parse_parameter_file(self):
        """Parse input datfile's header. Apply geometry_override if specified."""
        # required method

        # populate self.parameters with header data
        with open(self.parameter_filename, "rb") as istream:
            self.parameters.update(get_header(istream))

        self.current_time = self.parameters["time"]
        self.dimensionality = self.parameters["ndim"]

        # force 3D for this definition
        dd = np.ones(3, dtype="int64")
        dd[: self.dimensionality] = self.parameters["domain_nx"]
        self.domain_dimensions = dd

        if self.parameters.get("staggered", False):
            mylog.warning(
                "'staggered' flag was found, but is currently ignored (unsupported)"
            )

        # parse geometry
        # by order of decreasing priority, we use
        # - geometry_override
        # - "geometry" parameter from datfile
        # - if all fails, default to "cartesian"
        self.geometry = Geometry.CARTESIAN

        amrvac_geom = self.parameters.get("geometry", None)
        if amrvac_geom is not None:
            self.geometry = _parse_geometry(amrvac_geom)
        elif self.parameters["datfile_version"] > 4:
            mylog.error(
                "No 'geometry' flag found in datfile with version %d >4.",
                self.parameters["datfile_version"],
            )

        if self._geometry_override is not None:
            try:
                new_geometry = _parse_geometry(self._geometry_override)
                if new_geometry == self.geometry:
                    mylog.info("geometry_override is identical to datfile parameter.")
                else:
                    self.geometry = new_geometry
                    mylog.warning(
                        "Overriding geometry, this may lead to surprising results."
                    )
            except ValueError:
                mylog.error(
                    "Unable to parse geometry_override '%s' (will be ignored).",
                    self._geometry_override,
                )

        # parse peridiocity
        periodicity = self.parameters.get("periodic", ())
        missing_dim = 3 - len(periodicity)
        self._periodicity = (*periodicity, *(missing_dim * (False,)))

        self.gamma = self.parameters.get("gamma", 5.0 / 3.0)

        # parse domain edges
        dle = np.zeros(3)
        dre = np.ones(3)
        dle[: self.dimensionality] = self.parameters["xmin"]
        dre[: self.dimensionality] = self.parameters["xmax"]
        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # defaulting to non-cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_matter = 0.0
        self.omega_lambda = 0.0
        self.hubble_constant = 0.0

    # units stuff ======================================================================
    def _set_code_unit_attributes(self):
        """Reproduce how AMRVAC internally set up physical normalisation factors."""
        # This gets called later than Dataset._override_code_units()
        # This is the reason why it uses setdefaultattr: it will only fill in the gaps
        # left by the "override", instead of overriding them again.

        # note: yt sets hydrogen mass equal to proton mass, amrvac doesn't.
        mp_cgs = self.quan(1.672621898e-24, "g")  # This value is taken from AstroPy

        # get self.length_unit if overrides are supplied, otherwise use default
        length_unit = getattr(self, "length_unit", self.quan(1, "cm"))
        namelist = read_amrvac_namelist(self._parfiles)
        if namelist.get("eos_list", {}):
            He_abundance = namelist.get("eos_list", {}).get("he_abundance", 0.1)
        else:
            He_abundance = namelist.get("mhd_list", {}).get("he_abundance", 0.1)

        # 1. calculations for mass, density, numberdensity
        if "mass_unit" in self.units_override:
            # in this case unit_mass is supplied (and has been set as attribute)
            mass_unit = self.mass_unit
            density_unit = mass_unit / length_unit**3
            nd_unit = density_unit / ((1.0 + 4.0 * He_abundance) * mp_cgs)
        else:
            # other case: numberdensity is supplied.
            # Fall back to one (default) if no overrides supplied
            try:
                nd_unit = self.quan(self.units_override["numberdensity_unit"])
            except KeyError:
                nd_unit = self.quan(
                    1.0, self.__class__.default_units["numberdensity_unit"]
                )
            density_unit = (1.0 + 4.0 * He_abundance) * mp_cgs * nd_unit
            mass_unit = density_unit * length_unit**3

        # 2. calculations for velocity
        if "time_unit" in self.units_override:
            # in this case time was supplied
            velocity_unit = length_unit / self.time_unit
        else:
            # other case: velocity was supplied.
            # Fall back to None if no overrides supplied
            velocity_unit = getattr(self, "velocity_unit", None)

        # 3. calculations for pressure and temperature
        if velocity_unit is None:
            # velocity and time not given, see if temperature is given.
            # Fall back to one (default) if not
            temperature_unit = getattr(self, "temperature_unit", self.quan(1, "K"))
            pressure_unit = (
                (2.0 + 3.0 * He_abundance) * nd_unit * kb_cgs * temperature_unit
            ).in_cgs()
            velocity_unit = (np.sqrt(pressure_unit / density_unit)).in_cgs()
        else:
            # velocity is not zero if either time was given OR velocity was given
            pressure_unit = (density_unit * velocity_unit**2).in_cgs()
            temperature_unit = (
                pressure_unit / ((2.0 + 3.0 * He_abundance) * nd_unit * kb_cgs)
            ).in_cgs()

        # 4. calculations for magnetic unit and time
        time_unit = getattr(
            self, "time_unit", length_unit / velocity_unit
        )  # if time given use it, else calculate
        magnetic_unit = (np.sqrt(4 * np.pi * pressure_unit)).to("gauss")

        setdefaultattr(self, "mass_unit", mass_unit)
        setdefaultattr(self, "density_unit", density_unit)

        setdefaultattr(self, "length_unit", length_unit)
        setdefaultattr(self, "velocity_unit", velocity_unit)
        setdefaultattr(self, "time_unit", time_unit)

        setdefaultattr(self, "temperature_unit", temperature_unit)
        setdefaultattr(self, "pressure_unit", pressure_unit)
        setdefaultattr(self, "magnetic_unit", magnetic_unit)

    allowed_unit_combinations = [
        {"numberdensity_unit", "temperature_unit", "length_unit"},
        {"mass_unit", "temperature_unit", "length_unit"},
        {"mass_unit", "time_unit", "length_unit"},
        {"numberdensity_unit", "velocity_unit", "length_unit"},
        {"mass_unit", "velocity_unit", "length_unit"},
    ]

    default_units = {
        "length_unit": "cm",
        "time_unit": "s",
        "mass_unit": "g",
        "velocity_unit": "cm/s",
        "magnetic_unit": "gauss",
        "temperature_unit": "K",
        # this is the one difference with Dataset.default_units:
        # we accept numberdensity_unit as a valid override
        "numberdensity_unit": "cm**-3",
    }

    @classmethod
    def _validate_units_override_keys(cls, units_override):
        """Check that keys in units_override are consistent with AMRVAC's internal
        normalisations factors.
        """
        # YT supports overriding other normalisations, this method ensures consistency
        # between supplied 'units_override' items and those used by AMRVAC.

        # AMRVAC's normalisations/units have 3 degrees of freedom.
        # Moreover, if temperature unit is specified then velocity unit will be
        # calculated accordingly, and vice-versa.
        # We replicate this by allowing a finite set of combinations.

        # there are only three degrees of freedom, so explicitly check for this
        if len(units_override) > 3:
            raise ValueError(
                "More than 3 degrees of freedom were specified "
                f"in units_override ({len(units_override)} given)"
            )
        # temperature and velocity cannot both be specified
        if "temperature_unit" in units_override and "velocity_unit" in units_override:
            raise ValueError(
                "Either temperature or velocity is allowed in units_override, not both."
            )
        # check if provided overrides are allowed
        suo = set(units_override)
        for allowed_combo in cls.allowed_unit_combinations:
            if suo.issubset(allowed_combo):
                break
        else:
            raise ValueError(
                f"Combination {suo} passed to units_override "
                "is not consistent with AMRVAC.\n"
                f"Allowed combinations are {cls.allowed_unit_combinations}"
            )

        # syntax for mixing super with classmethod is weird...
        super(cls, cls)._validate_units_override_keys(units_override)
