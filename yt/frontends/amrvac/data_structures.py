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
from yt.geometry.oct_geometry_handler import \
    OctreeIndex
from yt.geometry.geometry_handler import \
    YTDataChunk
from yt.data_objects.static_output import \
    Dataset
from yt.data_objects.octree_subset import \
    OctreeSubset

from .fields import AMRVACFieldInfo
from .datreader import get_header, get_block_info, get_block_byte_limits


class AMRVACGrid(AMRGridPatch):
    """devnote : a patch represent part of a block. The hierarchy/index is a collection of patches"""
    _id_offset = 0

    def __init__(self, id, index, level, block_idx):
        super(AMRVACGrid, self).__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        # used to keep track of block index in the AMRVAC Morton curve
        # (useful in reading the data itself in io.py)
        self.block_idx = block_idx

    def __repr__(self):
        return "AMRVACGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class AMRVACHierarchy(GridIndex):
    grid = AMRVACGrid

    def __init__(self, ds, dataset_type="amrvac"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # the index file *is* the datfile
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64

        # init everything to make it clear what's in there
        self.field_list = []
        self.num_grids = None

        super(AMRVACHierarchy, self).__init__(ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [(self.dataset_type, f) for f in self.dataset.parameters["w_names"]]

    def _count_grids(self):
        """Set self.num_grids using datfile header"""
        self.num_grids = self.dataset.parameters['nleafs']

    def _create_patch(self, lvl, idx):
        grid_difference = 2**(self.dataset.parameters["levmax"] - lvl)
        max_idx = idx * grid_difference
        min_idx = max_idx - grid_difference

        # inner indices of block
        idx0 = min_idx * self.dataset.parameters["block_nx"]
        # outer indices of block    TODO: these depend on AMR level, right?
        if lvl == self.dataset.parameters["levmax"]:
            idx1 = idx0 + self.dataset.parameters["block_nx"]
        else:
            idx1 = idx0 + (self.dataset.parameters["block_nx"] * grid_difference)

        # Outer index of domain, taking AMR into account
        domain_end_idx = self.dataset.parameters["block_nx"] * 2**self.dataset.parameters["levmax"]
        # Width of the domain, used to correctly calculate fractions
        domain_width   = self.dataset.parameters["xmax"] - self.dataset.parameters["xmin"]

        # TOREVIEW @Niels
        # So idx0 / domain_end_idx gives the "fraction" (between 0 and 1) of the current block
        # position. Multiply this by domain_width to take the width of the domain into account,
        # as this can vary from one.
        l_edge = idx0 / domain_end_idx * domain_width
        r_edge = idx1 / domain_end_idx * domain_width
        block_nx = self.dataset.parameters["block_nx"]
        dim = self.dataset.dimensionality
        if dim < 3:
            d = 3-dim
            l_edge = np.append(l_edge, [0]*d)
            r_edge = np.append(r_edge, [1]*d)
            block_nx = np.append(block_nx, [1]*d)
        patch = {
            "left_edge":  l_edge,
            "right_edge": r_edge,
            # TOREVIEW: is this what "width" is supposed to represent ?
            "width": block_nx,
            "block_idx": idx
        }
        return patch

    def _add_patch(self, igrid, patch):
        # for idim, left_edge in enumerate(patch['left_edge']):
        #     self.grid_left_edge[igrid, idim] = left_edge
        # for idim, right_edge in enumerate(patch['right_edge']):
        #     self.grid_right_edge[igrid, idim] = right_edge
        # for idim, width  in enumerate(patch['width']):
        #     self.grid_dimensions[igrid, idim] = width
        self.grid_left_edge[igrid, :] = patch["left_edge"]
        self.grid_right_edge[igrid, :] = patch["right_edge"]
        self.grid_dimensions[igrid, :] = patch["width"]

        # TOREVIEW: is it intended that this attribute changes at each call ?
        self.grid_block_idx = patch["block_idx"]

    def _parse_index(self):
        with open(self.index_filename, "rb") as istream:
            lvls, idxs = get_block_info(istream)
            assert len(lvls) == len(idxs) == self.num_grids

        for igrid, (lvl, idx) in enumerate(zip(lvls, idxs)):
            # devnote : idx is the index on the Morton Curve
            # maybe it ought to be properly translated to yt indexing first...
            patch = self._create_patch(lvl, idx)
            self._add_patch(igrid, patch)

        # YT uses 0-based grid indexing, lowest level = 0 (AMRVAC uses 1 for lowest level)
        ytlevels = np.array(lvls, dtype="int32") - 1

        self.grid_levels = ytlevels.reshape(self.num_grids, 1)
        self.max_level = self.dataset.parameters["levmax"] - 1
        assert self.max_level == max(ytlevels)

        self.grids = np.empty(self.num_grids, dtype='object')
        for i in range(self.num_grids):
            self.grids[i] = self.grid(i, self, self.grid_levels[i, 0], self.grid_block_idx)

        with open(self.index_filename, "rb") as istream:
            self.block_offsets, self.block_shapes = get_block_byte_limits(istream)

    def _populate_grid_objects(self):
        for g in self.grids:
            g._prepare_grid()
            g._setup_dx()



class AMRVACDataset(Dataset):
    _index_class = AMRVACHierarchy
    _field_info_class = AMRVACFieldInfo

    def __init__(self, filename, dataset_type='amrvac',
                 storage_filename=None,
                 units_override=None,):
        self.fluid_types += ('amrvac',) #devnote: input 'gas', 'dust' here ?
        super(AMRVACDataset, self).__init__(filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename
        # refinement factor between a grid and its subgrid
        self.refine_by = 2
        self.cosmological_simulation = False

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.

        # TODO @Niels: Default to cgs. This will be implemented in upcoming .dat file extensions.

        unit_system = self.parameters.get("unit_system", "cgs")

        # TODO @Niels: Units in AMRVAC are defined by specifying unit_numberdensity and unit_length, together with
        #              EITHER unit_temperature or unit_velocity. This will again be implemented in .dat file changes.

        # @Niels: can we use astropy here?
        # @clm : no, but we *could* use unyt
        mp  = 1.672621777e-24    # cgs units (g)
        kB  = 1.3806488e-16      # cgs units (erg / K)
        mu0 = 4 * np.pi
        He_abundance = 0.1       # hardcoded in AMRVAC

        if unit_system == "cgs":
            pass
        elif unit_system == "si":
            mp *= 1e-3
            kB *= 1e-7
            mu0 = 1.2566370614e-6
        else:
            raise RuntimeError("AMRVAC data file contains an "
                               "unknown unit_system parameter {}".format(self.parameters["unit_system"]))

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
        else:
            setdefaultattr(self, "length_unit", self.quan(unit_length, "m"))
            setdefaultattr(self, "numberdensity_unit", self.quan(unit_numberdensity, "m**-3"))
            setdefaultattr(self, "velocity_unit", self.quan(unit_velocity, "m/s"))
            setdefaultattr(self, "temperature_unit", self.quan(unit_temperature, "K"))
            setdefaultattr(self, "density_unit", self.quan(unit_density, "kg*m**-3"))
            setdefaultattr(self, "pressure_unit", self.quan(unit_pressure, "pa"))
            setdefaultattr(self, "mass_unit", self.quan(unit_mass, "kg"))
            setdefaultattr(self, "time_unit", self.quan(unit_time, "s"))


    # TODO: uncomment this when using "setdefaultattr" above, it overrides default yt code units (see flash)
    def set_code_units(self):
        super(AMRVACDataset, self).set_code_units()


    def _parse_parameter_file(self):
        # This needs to set up the following items.  Note that these are all
        # assumed to be in code units; domain_left_edge and domain_right_edge
        # will be converted to YTArray automatically at a later time.
        # This includes the cosmological parameters.
        #
        #   self.unique_identifier      <= unique identifier for the dataset
        #                                  being read (e.g., UUID or ST_CTIME)
        #   self.parameters             <= full of code-specific items of use
        #   self.domain_left_edge       <= array of float64                         OK
        #   self.domain_right_edge      <= array of float64                         OK
        #   self.dimensionality         <= int                                      OK
        #   self.domain_dimensions      <= array of int64                           OK
        #   self.periodicity            <= three-element tuple of booleans          OK
        #   self.current_time           <= simulation time in code units            OK
        #
        # We also set up cosmological information.  Set these to zero if
        # non-cosmological.
        #
        #   self.cosmological_simulation    <= int, 0 or 1                          OK
        #   self.current_redshift           <= float                                OK
        #   self.omega_lambda               <= float                                OK
        #   self.omega_matter               <= float                                OK
        #   self.hubble_constant            <= float                                OK
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        with open(self.parameter_filename, 'rb') as df:
            self.parameters = get_header(df)

        self.current_time      = self.parameters['time']
        self.dimensionality    = self.parameters['ndim'] #devnote, warining : ndir != ndim
        # for domain_dimensions, found similar thing in Flash frontend
        self.domain_dimensions = self.parameters['block_nx'] * 2**self.parameters['levmax']
        self.gamma             = self.parameters['gamma']

        dle = np.zeros(3)
        dre = np.ones(3)
        for idim in range(self.dimensionality):
            dle[idim] = self.parameters['xmin'][idim]
            dre[idim] = self.parameters['xmax'][idim]

        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # TODO @Niels: this must also be included in .dat file.
        self.periodicity = tuple([False, False, False])

        #devnote: these could be made optional if needed
        self.cosmological_simulation = 0
        self.current_redshift        = 0.0
        self.omega_matter            = 0.0
        self.omega_lambda            = 0.0
        self.hubble_constant         = 0.0

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        validation = False
        try:
            with open(args[0], "rb") as istream:
                h = get_header(istream)
            validation = True
        finally:
            return validation
