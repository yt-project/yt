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
import numpy as np
import weakref

from yt.data_objects.grid_patch import \
   AMRGridPatch
from yt.geometry.grid_geometry_handler import \
   GridIndex
from yt.funcs import \
    mylog, \
    setdefaultattr
from yt.data_objects.static_output import \
   Dataset

from .fields import AMRVACFieldInfo
from .datfile_utils import get_header, get_tree_info



class AMRVACGrid(AMRGridPatch):
    """devnote : a patch represent part of a block. The hierarchy/index is a collection of patches"""
    _id_offset = 0

    def __init__(self, id, index, level):
        #note: the <level> here should use yt's convention (start from 0)
        super(AMRVACGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level

    def __repr__(self):
        return "AMRVACGrid_%04i (%s)" % (self.id, self.ActiveDimensions)

    def get_global_startindex(self):
        """
        Return the integer starting index for each dimension at the current
        level.

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
        self.field_list = [(self.dataset_type, f) for f in self.dataset.parameters["w_names"]]

    def _count_grids(self):
        """Set self.num_grids using datfile header"""
        self.num_grids = self.dataset.parameters['nleafs']

    def _get_patch(self, ytlevel, morton_index):
        # Width of the domain, used to correctly calculate fractions
        domain_width = self.dataset.parameters["xmax"] - self.dataset.parameters["xmin"]
        block_nx = self.dataset.parameters["block_nx"]

        # dx at coarsest grid level (YT level 0)
        dx0 = domain_width / self.dataset.parameters["domain_nx"]
        # dx at ytlevel
        dx = dx0 * (1./self.dataset.refine_by)**ytlevel
        l_edge = self.dataset.parameters["xmin"] + (morton_index-1) * block_nx * dx
        r_edge = l_edge + block_nx * dx

        # force return arrays to 3D
        missing_dim = 3 - self.dataset.dimensionality
        l_edge = np.append(l_edge, [0]*missing_dim)
        r_edge = np.append(r_edge, [1]*missing_dim)
        block_nx = np.append(block_nx, [1]*missing_dim)

        patch = {
            "left_edge": l_edge,
            "right_edge": r_edge,
            "width": block_nx, # number of cells along an axis
        }
        return patch

    def _parse_index(self):
        with open(self.index_filename, "rb") as istream:
            vaclevels, idxs, block_offsets = get_tree_info(istream)
            assert len(vaclevels) == len(idxs) == len(block_offsets) == self.num_grids

        self.block_offsets = block_offsets
        # YT uses 0-based grid indexing, lowest level = 0 (AMRVAC uses 1 for lowest level)
        ytlevels = np.array(vaclevels, dtype="int32") - 1
        self.grid_levels.flat[:] = ytlevels
        self.max_level = np.max(ytlevels)
        assert self.max_level == self.dataset.parameters["levmax"] - 1

        self.grids = np.empty(self.num_grids, dtype='object')

        for igrid, (lvl, idx) in enumerate(zip(ytlevels, idxs)):
            # devnote : idx is the index on the Morton Curve
            # maybe it ought to be properly translated to yt indexing first...
            patch = self._get_patch(lvl, idx)
            self.grid_left_edge[igrid, :] = patch["left_edge"]
            self.grid_right_edge[igrid, :] = patch["right_edge"]
            self.grid_dimensions[igrid, :] = patch["width"]

            self.grids[igrid] = self.grid(igrid, self, ytlevels[igrid])

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()



class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(self, filename, dataset_type='amrvac', units_override=None, unit_system="code", geometry_override=None):
        self._geometry_override = geometry_override
        super(AMRVACDataset, self).__init__(filename, dataset_type,
                                            units_override=units_override, unit_system=unit_system)
        self.fluid_types += ('amrvac',)
        # refinement factor between a grid and its subgrid
        self.refine_by = 2

    @classmethod
    def _is_valid(self, *args, **kwargs):
        validation = False
        try:
            with open(args[0], "rb") as istream:
                get_header(istream)
            validation = True
        finally:
            return validation

    def parse_geometry(self, geometry_string):
        geom = geometry_string.split("_")[0].lower()
        if geom not in ("cartesian", "polar", "cylindrical", "spherical"):
            raise ValueError
        return geom

    def _parse_parameter_file(self):
        self.unique_identifier = int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        # populate self.parameters with header data
        with open(self.parameter_filename, 'rb') as istream:
            self.parameters.update(get_header(istream))

        self.current_time = self.parameters['time']
        self.dimensionality = self.parameters['ndim']

        # force 3D for this definition
        self.domain_dimensions = np.ones(3, dtype="int64")
        self.domain_dimensions[:self.dimensionality] = self.parameters['domain_nx']

        # the following parameters may not be present in the datfile,
        # dependending on format version
        if self.parameters["datfile_version"] < 5:
            mylog.warning("This data format does not contain geometry or periodicity info")
        if self.parameters.get("staggered", False):
            mylog.warning("'staggered' flag was found, but is currently ignored (unsupported)")

        # parse geometry
        # by order of descending priority, we use
        # - geometry_override 
        # - "geometry" parameter from datfile
        # - if all fails, default ("cartesian")
        geom_candidates = {"param": None, "override": None}
        amrvac_geom = self.parameters.get("geometry", None)
        if amrvac_geom is None:
            mylog.warning("Could not find a 'geometry' parameter in source file.")
        else:
            geom_candidates.update({"param": self.parse_geometry(amrvac_geom)})

        if self._geometry_override is not None:
            try:
                geom_candidates.update({"override": self.parse_geometry(self._geometry_override)})
            except ValueError:
                mylog.error("Unknown value for geometry_override (will be ignored).")
    
        if geom_candidates["override"] is not None:
            mylog.warning("Using override geometry, this may lead to surprising results for inappropriate values.")
            self.geometry = geom_candidates["override"]
        elif geom_candidates["param"] is not None:
            mylog.info("Using parameter geometry")
            self.geometry = geom_candidates["param"]
        else:
            mylog.warning("No geometry parameter supplied or found, defaulting to cartesian.")
            self.geometry = "cartesian"

        # parse peridicity
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
        self.current_redshift        = 0.0
        self.omega_matter            = 0.0
        self.omega_lambda            = 0.0
        self.hubble_constant         = 0.0

    # units stuff ===============================================================================
    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.

        mp  = 1.672621777e-24    # cgs units (g)
        kB  = 1.3806488e-16      # cgs units (erg / K)
        mu0 = 4 * np.pi
        He_abundance = 0.1       # hardcoded in AMRVAC

        unit_system = self.parameters.get("unit_system", "cgs")
        # devnote : leaving this idle for now, waiting for datfile format 6 to be released
        #if self.parameters["datfile_version"] < 6:
        #    mylog.warning("This datfile format does not contain unit system. Defaulting to cgs.")

        if unit_system not in ("cgs", "si"):
            mylog.warning("Unknown unit_system parameter '%s'" % self.parameters["unit_system"])

        if unit_system == "si":
            mp *= 1e-3
            kB *= 1e-7
            mu0 = 1.2566370614e-6

        # Obtain unit normalisations from .dat file
        # default values mock AMRVAC itself
        unit_length = self.parameters.get("unit_length", 1)
        unit_numberdensity = self.parameters.get("unit_numberdensity", 1)
        unit_velocity = self.parameters.get("unit_velocity", 0)
        unit_temperature = self.parameters.get("unit_temperature", 1)

        unit_density = (1.0 + 4.0*He_abundance) * mp * unit_numberdensity

        # numberdensity, length and temperature defined in mod_usr.t file.
        if unit_velocity == 0:
            unit_pressure = (2.0 + 3.0*He_abundance) * unit_numberdensity * kB * unit_temperature
            unit_velocity = np.sqrt(unit_pressure / unit_density)
        # numberdensity, length and velocity defined in mod_usr.t file.
        else:
            unit_pressure    = unit_density * unit_velocity**2
            unit_temperature = unit_pressure / ((2.0 + 3.0*He_abundance) * unit_numberdensity * kB)

        unit_magneticfield = np.sqrt(mu0 * unit_pressure)
        unit_time          = unit_length / unit_velocity

        # TOREVIEW @Niels: this seems really off, I don't think density and numberdensity share a unit
        unit_mass          = unit_numberdensity * unit_length**3


        # Set unit attributes
        if unit_system == "cgs":
            setdefaultattr(self, "length_unit", self.quan(unit_length, "cm"))
            setdefaultattr(self, "numberdensity_unit", self.quan(unit_numberdensity, "cm**-3"))
            setdefaultattr(self, "velocity_unit", self.quan(unit_velocity, "cm/s"))
            setdefaultattr(self, "temperature_unit", self.quan(unit_temperature, "K"))
            setdefaultattr(self, "density_unit", self.quan(unit_density, "g*cm**-3"))
            setdefaultattr(self, "pressure_unit", self.quan(unit_pressure, "dyn*cm**-2"))
            setdefaultattr(self, "magnetic_unit", self.quan(unit_magneticfield, "gauss"))
            setdefaultattr(self, "mass_unit", self.quan(unit_mass, "g"))
            setdefaultattr(self, "time_unit", self.quan(unit_time, "s"))
        elif unit_system == "si":
            setdefaultattr(self, "length_unit", self.quan(unit_length, "m"))
            setdefaultattr(self, "numberdensity_unit", self.quan(unit_numberdensity, "m**-3"))
            setdefaultattr(self, "velocity_unit", self.quan(unit_velocity, "m/s"))
            setdefaultattr(self, "temperature_unit", self.quan(unit_temperature, "K"))
            setdefaultattr(self, "density_unit", self.quan(unit_density, "kg*m**-3"))
            setdefaultattr(self, "pressure_unit", self.quan(unit_pressure, "pa"))
            setdefaultattr(self, "mass_unit", self.quan(unit_mass, "kg"))
            setdefaultattr(self, "time_unit", self.quan(unit_time, "s"))

    def set_code_unit(self):
        super(AMRVACDataset, self).set_code_units()
