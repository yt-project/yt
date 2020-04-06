"""
AMRVAC data structures



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import stat
import weakref
import struct

import numpy as np

from yt.data_objects.grid_patch import \
   AMRGridPatch
from yt.geometry.grid_geometry_handler import \
   GridIndex
from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.static_output import \
   Dataset
from yt.utilities.physical_constants import \
    boltzmann_constant_cgs as kb_cgs

from .fields import AMRVACFieldInfo
from .datfile_utils import get_header, get_tree_info
from . import read_amrvac_namelist

ALLOWED_UNIT_COMBINATIONS = [{'numberdensity_unit', 'temperature_unit', 'length_unit'},
                             {'mass_unit', 'temperature_unit', 'length_unit'},
                             {'mass_unit', 'time_unit', 'length_unit'},
                             {'numberdensity_unit', 'velocity_unit', 'length_unit'},
                             {'mass_unit', 'velocity_unit', 'length_unit'}]


class AMRVACGrid(AMRGridPatch):
    """A class to populate AMRVACHierarchy.grids, setting parent/children relations."""
    _id_offset = 0

    def __init__(self, id, index, level):
        #<level> should use yt's convention (start from 0)
        super(AMRVACGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "AMRVACGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

    def get_global_startindex(self):
        """Refresh and retrieve the starting index for each dimension at current level.

        Returns
        -------
        self.start_index : int
        """
        start_index = (self.LeftEdge - self.ds.domain_left_edge)/self.dds
        self.start_index = np.rint(start_index).astype('int64').ravel()
        return self.start_index


class AMRVACHierarchy(GridIndex):
    grid = AMRVACGrid
    def __init__(self, ds, dataset_type="amrvac"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # the index file *is* the datfile
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64

        super(AMRVACHierarchy, self).__init__(ds, dataset_type)

    def _detect_output_fields(self):
        """Parse field names from datfile header, which is stored in self.dataset.parameters"""
        # required method
        self.field_list = [(self.dataset_type, f) for f in self.dataset.parameters["w_names"]]

    def _count_grids(self):
        """Set self.num_grids from datfile header."""
        # required method
        self.num_grids = self.dataset.parameters['nleafs']

    def _parse_index(self):
        """Populate self.grid_* attributes from tree info from datfile header."""
        # required method
        with open(self.index_filename, "rb") as istream:
            vaclevels, morton_indices, block_offsets = get_tree_info(istream)
            assert len(vaclevels) == len(morton_indices) == len(block_offsets) == self.num_grids

        self.block_offsets = block_offsets
        # YT uses 0-based grid indexing, lowest level = 0 (AMRVAC uses 1 for lowest level)
        ytlevels = np.array(vaclevels, dtype="int32") - 1
        self.grid_levels.flat[:] = ytlevels
        self.min_level = np.min(ytlevels)
        self.max_level = np.max(ytlevels)
        assert self.max_level == self.dataset.parameters["levmax"] - 1

        # some aliases for left/right edges computation in the coming loop
        domain_width = self.dataset.parameters["xmax"] - self.dataset.parameters["xmin"]
        block_nx = self.dataset.parameters["block_nx"]
        xmin = self.dataset.parameters["xmin"]
        dx0 = domain_width / self.dataset.parameters["domain_nx"] # dx at coarsest grid level (YT level 0)
        dim = self.dataset.dimensionality

        self.grids = np.empty(self.num_grids, dtype='object')
        for igrid, (ytlevel, morton_index) in enumerate(zip(ytlevels, morton_indices)):
            dx = dx0 / self.dataset.refine_by**ytlevel
            left_edge = xmin + (morton_index-1) * block_nx * dx

            # edges and dimensions are filled in a dimensionality-agnostic way
            self.grid_left_edge[igrid, :dim] = left_edge
            self.grid_right_edge[igrid, :dim] = left_edge + block_nx * dx
            self.grid_dimensions[igrid, :dim] = block_nx
            self.grids[igrid] = self.grid(igrid, self, ytlevels[igrid])

    def _populate_grid_objects(self):
        # required method
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()



class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(self, filename, dataset_type='amrvac',
                units_override=None, unit_system="cgs",
                geometry_override=None,
                parfiles=None):
        """Instanciate AMRVACDataset.

        Parameters
        ----------
        filename : str
            Path to a datfile.

        dataset_type : str, optional
            This should always be 'amrvac'.

        units_override : dict, optional
            A dictionnary of physical normalisation factors to interpret on disk data.

        unit_system : str, optional
            Either "cgs" (default), "mks" or "code"

        geometry_override : str, optional
            A geometry flag formatted either according to either AMRVAC's or yt's standards.
            When this parameter is passed along with v5 or more newer datfiles, will precede over
            their internal "geometry" tag.

        parfiles : str or list, optional
            One or more parfiles to be passed to yt.frontends.amrvac.read_amrvac_parfiles()

        """
        # note: geometry_override and parfiles are specific to this frontend

        self._geometry_override = geometry_override
        super(AMRVACDataset, self).__init__(filename, dataset_type,
                                            units_override=units_override, unit_system=unit_system)

        self._parfiles = parfiles

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

            if namelist_gamma is not None and self.gamma != namelist_gamma:
                mylog.error("Inconsistent values in gamma: datfile {}, parfiles {}".format(self.gamma, namelist_gamma))

            if "method_list" in namelist:
                e_is_internal = namelist["method_list"].get("solve_internal_e", False)

        if c_adiab is not None:
            # this complicated unit is required for the adiabatic equation of state to make physical sense
            c_adiab *= self.mass_unit**(1-self.gamma) * self.length_unit**(2+3*(self.gamma-1)) / self.time_unit**2

        self.namelist = namelist
        self._c_adiab = c_adiab
        self._e_is_internal = e_is_internal

        self.fluid_types += ('amrvac',)
        # refinement factor between a grid and its subgrid
        self.refine_by = 2

    @classmethod
    def _is_valid(self, *args, **kwargs):
        """At load time, check whether data is recognized as AMRVAC formatted."""
        # required class method
        validation = False
        if args[0].endswith(".dat"):
            try:
                with open(args[0], mode="rb") as istream:
                    fmt = "=i"
                    [datfile_version] = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
                    if 3 <= datfile_version < 6:
                        fmt = "=ii"
                        offset_tree, offset_blocks = struct.unpack(fmt, istream.read(struct.calcsize(fmt)))
                        istream.seek(0,2)
                        file_size = istream.tell()
                        validation = offset_tree < file_size and offset_blocks < file_size
            except Exception: pass
        return validation

    def _parse_geometry(self, geometry_tag):
        """Translate AMRVAC's geometry tag to yt's format.

        Parameters
        ----------
        geometry_tag : str
            A geometry tag as read from AMRVAC's datfile from v5.

        Returns
        -------
        geometry_yt : str
            Lower case geometry tag among "cartesian", "polar", "cylindrical", "spherical"

        Raises
        ------
        ValueError
            In case the tag is not properly formatted or recognized.

        Examples
        --------
        >>> print(self._parse_geometry("Polar_2.5D"))
        "polar"
        >>> print(self._parse_geometry("Cartesian_2.5D"))

        """
        # frontend specific method
        geometry_yt = geometry_tag.split("_")[0].lower()
        if geometry_yt not in ("cartesian", "polar", "cylindrical", "spherical"):
            raise ValueError
        return geometry_yt

    def _parse_parameter_file(self):
        """Parse input datfile's header. Apply geometry_override if specified."""
        # required method
        self.unique_identifier = int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # populate self.parameters with header data
        with open(self.parameter_filename, 'rb') as istream:
            self.parameters.update(get_header(istream))

        self.current_time = self.parameters['time']
        self.dimensionality = self.parameters['ndim']

        # force 3D for this definition
        dd = np.ones(3, dtype="int64")
        dd[:self.dimensionality] = self.parameters['domain_nx']
        self.domain_dimensions = dd

        if self.parameters.get("staggered", False):
            mylog.warning("'staggered' flag was found, but is currently ignored (unsupported)")

        # parse geometry
        # by order of decreasing priority, we use
        # - geometry_override
        # - "geometry" parameter from datfile
        # - if all fails, default to "cartesian"
        self.geometry = None
        amrvac_geom = self.parameters.get("geometry", None)
        if amrvac_geom is not None:
            self.geometry = self._parse_geometry(amrvac_geom)
        elif self.parameters["datfile_version"] > 4:
            # py38: walrus here
            mylog.error("No 'geometry' flag found in datfile with version %d >4." % self.parameters["datfile_version"])

        if self._geometry_override is not None:
            # py38: walrus here
            try:
                new_geometry = self._parse_geometry(self._geometry_override)
                if new_geometry == self.geometry:
                    mylog.info("geometry_override is identical to datfile parameter.")
                else:
                    self.geometry = new_geometry
                    mylog.warning("Overriding geometry, this may lead to surprising results.")
            except ValueError:
                mylog.error("Unable to parse geometry_override '%s' (will be ignored)." % self._geometry_override)

        if self.geometry is None:
            mylog.warning("No geometry parameter supplied or found, defaulting to cartesian.")
            self.geometry = "cartesian"

        # parse peridiocity
        per = self.parameters.get("periodic", np.array([False, False, False]))
        missing_dim = 3 - len(per)
        self.periodicity = np.append(per, [False]*missing_dim)

        self.gamma = self.parameters.get("gamma", 5.0/3.0)

        # parse domain edges
        dle = np.zeros(3)
        dre = np.ones(3)
        dle[:self.dimensionality] = self.parameters['xmin']
        dre[:self.dimensionality] = self.parameters['xmax']
        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # defaulting to non-cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_matter = 0.0
        self.omega_lambda = 0.0
        self.hubble_constant = 0.0

    # units stuff ===============================================================================
    def _set_code_unit_attributes(self):
        """Reproduce how AMRVAC internally set up physical normalisation factors."""
        # required method
        # devnote: this method is never defined in the parent abstract class Dataset
        # but it is called in Dataset.set_code_units(), which is part of Dataset.__init__()
        # so it must be defined here.

        # devnote: this gets called later than Dataset._override_code_units()
        # This is the reason why it uses setdefaultattr: it will only fill in the gaps left
        # by the "override", instead of overriding them again.
        # For the same reason, self.units_override is set, as well as corresponding *_unit instance attributes
        # which may include up to 3 of the following items: length, time, mass, velocity, number_density, temperature

        # note: yt sets hydrogen mass equal to proton mass, amrvac doesn't.
        mp_cgs = self.quan(1.672621898e-24, 'g')  # This value is taken from AstroPy
        He_abundance = 0.1  # hardcoded parameter in AMRVAC

        # get self.length_unit if overrides are supplied, otherwise use default
        length_unit = getattr(self, 'length_unit', self.quan(1, 'cm'))

        # 1. calculations for mass, density, numberdensity
        if 'mass_unit' in self.units_override:
            # in this case unit_mass is supplied (and has been set as attribute)
            mass_unit = self.mass_unit
            density_unit = mass_unit / length_unit**3
            numberdensity_unit = density_unit / ((1.0 + 4.0 * He_abundance) * mp_cgs)
        else:
            # other case: numberdensity is supplied. Fall back to one (default) if no overrides supplied
            numberdensity_override = self.units_override.get('numberdensity_unit', (1, 'cm**-3'))
            if 'numberdensity_unit' in self.units_override: # print similar warning as yt when overriding numberdensity
                mylog.info("Overriding numberdensity_unit: %g %s.", *numberdensity_override)
            numberdensity_unit = self.quan(*numberdensity_override)  # numberdensity is never set as attribute
            density_unit = (1.0 + 4.0 * He_abundance) * mp_cgs * numberdensity_unit
            mass_unit = density_unit * length_unit**3

        # 2. calculations for velocity
        if 'time_unit' in self.units_override:
            # in this case time was supplied
            velocity_unit = length_unit / self.time_unit
        else:
            # other case: velocity was supplied. Fall back to None if no overrides supplied
            velocity_unit = getattr(self, 'velocity_unit', None)

        # 3. calculations for pressure and temperature
        if velocity_unit is None:
            # velocity and time not given, see if temperature is given. Fall back to one (default) if not
            temperature_unit = getattr(self, 'temperature_unit', self.quan(1, 'K'))
            pressure_unit = ((2.0 + 3.0 * He_abundance) * numberdensity_unit * kb_cgs * temperature_unit).in_cgs()
            velocity_unit = (np.sqrt(pressure_unit / density_unit)).in_cgs()
        else:
            # velocity is not zero if either time was given OR velocity was given
            pressure_unit = (density_unit * velocity_unit ** 2).in_cgs()
            temperature_unit = (pressure_unit / ((2.0 + 3.0 * He_abundance) * numberdensity_unit * kb_cgs)).in_cgs()

        # 4. calculations for magnetic unit and time
        time_unit = getattr(self, 'time_unit', length_unit / velocity_unit)  # if time given use it, else calculate
        magnetic_unit = (np.sqrt(4 * np.pi * pressure_unit)).to('gauss')

        setdefaultattr(self, 'mass_unit', mass_unit)
        setdefaultattr(self, 'density_unit', density_unit)
        setdefaultattr(self, 'numberdensity_unit', numberdensity_unit)

        setdefaultattr(self, 'length_unit', length_unit)
        setdefaultattr(self, 'velocity_unit', velocity_unit)
        setdefaultattr(self, 'time_unit', time_unit)

        setdefaultattr(self, 'temperature_unit', temperature_unit)
        setdefaultattr(self, 'pressure_unit', pressure_unit)
        setdefaultattr(self, 'magnetic_unit', magnetic_unit)


    def _override_code_units(self):
        """Add a check step to the base class' method (Dataset)."""
        self._check_override_consistency()
        super(AMRVACDataset, self)._override_code_units()


    def _check_override_consistency(self):
        """Check that keys in units_override are consistent with respect to AMRVAC's internal way to
        set up normalisations factors.

        """
        # frontend specific method
        # YT supports overriding other normalisations, this method ensures consistency between
        # supplied 'units_override' items and those used by AMRVAC.

        # AMRVAC's normalisations/units have 3 degrees of freedom.
        # Moreover, if temperature unit is specified then velocity unit will be calculated
        # accordingly, and vice-versa.
        # Currently we replicate this by allowing a finite set of combinations in units_override
        if not self.units_override:
            return
        overrides = set(self.units_override)

        # there are only three degrees of freedom, so explicitly check for this
        if len(overrides) > 3:
            raise ValueError('More than 3 degrees of freedom were specified '
                             'in units_override ({} given)'.format(len(overrides)))
        # temperature and velocity cannot both be specified
        if 'temperature_unit' in overrides and 'velocity_unit' in overrides:
            raise ValueError('Either temperature or velocity is allowed in units_override, not both.')
        # check if provided overrides are allowed
        for allowed_combo in ALLOWED_UNIT_COMBINATIONS:
            if overrides.issubset(allowed_combo):
                break
        else:
            raise ValueError('Combination {} passed to units_override is not consistent with AMRVAC. \n'
                             'Allowed combinations are {}'.format(overrides, ALLOWED_UNIT_COMBINATIONS))
