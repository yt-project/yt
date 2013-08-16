"""
ARTIO-specific data structures

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt-project.org/
License:
  Copyright (C) 2010-2011 Matthew Turk.  All Rights Reserved.

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
import numpy as np
import stat
import weakref
import cStringIO

from .definitions import yt_to_art, art_to_yt, ARTIOconstants
from _artio_caller import \
    artio_is_valid, artio_fileset, ARTIOOctreeContainer, \
    ARTIORootMeshContainer
import _artio_caller
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from .fields import ARTIOFieldInfo, KnownARTIOFields, b2t

from yt.funcs import *
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput
from yt.data_objects.octree_subset import \
    OctreeSubset
from yt.data_objects.data_containers import \
    YTFieldData

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

class ARTIOOctreeSubset(OctreeSubset):
    _domain_offset = 0
    domain_id = 2
    _con_args = ("base_region", "sfc_start", "sfc_end", "pf")
    _type_name = 'octree_subset'

    def __init__(self, base_region, sfc_start, sfc_end, pf):
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.sfc_start = sfc_start
        self.sfc_end = sfc_end
        self.pf = pf
        self.hierarchy = self.pf.hierarchy
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = 'all'
        self._current_fluid_type = self.pf.default_fluid_type
        self.base_region = base_region
        self.base_selector = base_region.selector

    _oct_handler = None

    @property
    def oct_handler(self):
        if self._oct_handler is None: 
            self._oct_handler = ARTIOOctreeContainer(
                self.pf.domain_dimensions/2, # Octs, not cells
                self.pf.domain_left_edge, self.pf.domain_right_edge,
                self.sfc_start, self.sfc_end, self.pf._handle)
        return self._oct_handler

    @property
    def min_ind(self):
        return self.sfc_start

    @property
    def max_ind(self):
        return self.sfc_end

    def fill(self, fields, selector):
        # Here we get a copy of the file, which we skip through and read the
        # bits we want.
        handle = self.oct_handler.artio_handle
        field_indices = [handle.parameters["grid_variable_labels"].index(
                        yt_to_art[f]) for (ft, f) in fields]
        cell_count = selector.count_oct_cells(
            self.oct_handler, self.domain_id)
        self.data_size = cell_count
        levels, cell_inds, file_inds = self.oct_handler.file_index_octs(
            selector, self.domain_id, cell_count)
        tr = [np.zeros(cell_count, dtype="float64") for field in fields]
        self.oct_handler.fill_sfc(levels, cell_inds, file_inds, field_indices, tr)
        tr = dict((field, v) for field, v in zip(fields, tr))
        return tr

    def fill_particles(self, fields):
        art_fields = []
        for s, f in fields:
            fn = yt_to_art[f]
            for i in self.pf.particle_type_map[s]:
                if fn in self.pf.particle_variables[i]:
                    art_fields.append((i, fn))
        species_data = self.oct_handler.fill_sfc_particles(art_fields)
        tr = defaultdict(dict)
        # Now we need to sum things up and then fill
        for s, f in fields:
            count = 0
            fn = yt_to_art[f]
            dt = "float64" # default
            for i in self.pf.particle_type_map[s]:
                if (i, fn) not in species_data: continue
                # No vector fields in ARTIO
                count += species_data[i, fn].size
                dt = species_data[i, fn].dtype
            tr[s][f] = np.zeros(count, dtype=dt)
            cp = 0
            for i in self.pf.particle_type_map[s]:
                if (i, fn) not in species_data: continue
                v = species_data.pop((i, fn))
                tr[s][f][cp:cp+v.size] = v
                cp += v.size
        return tr

# We create something of a fake octree here.  This is primarily to enable us to
# reuse code for things like __getitem__ and the like.  We will also create a
# new oct_handler type that is functionally equivalent, except that it will
# only manage the root mesh.
class ARTIORootMeshSubset(ARTIOOctreeSubset):
    _num_zones = 1
    _type_name = 'sfc_subset'
    _selector_module = _artio_caller
    domain_id = 1

    @property
    def oct_handler(self):
        if self._oct_handler is None: 
            self._oct_handler = ARTIORootMeshContainer(
                self.pf.domain_dimensions, # Cells, not octs
                self.pf.domain_left_edge, self.pf.domain_right_edge,
                self.pf._handle,
                self.sfc_start, self.sfc_end)
        return self._oct_handler

    def fill(self, fields, selector):
        # We know how big these will be.
        handle = self.pf._handle
        field_indices = [handle.parameters["grid_variable_labels"].index(
                        yt_to_art[f]) for (ft, f) in fields]
        tr = self.oct_handler.fill_sfc(selector, field_indices)
        self.data_size = tr[0].size
        tr = dict((field, v) for field, v in zip(fields, tr))
        return tr

class ARTIOChunk(object):

    def __init__(self, pf, sfc_start, sfc_end):
        self.pf = pf
        self.sfc_start = sfc_start
        self.sfc_end = sfc_end

    _data_size = None

    @property
    def data_size(self):
        if self._data_size is None:
            mylog.error("ARTIOChunk.data_size called before fill")
            raise RuntimeError
        return self._data_size

    _fcoords = None
    def select_fcoords(self, dobj):
        if self._fcoords is None:
            mylog.error("ARTIOChunk.fcoords called before fill")
            raise RuntimeError
        return self._fcoords

    _ires = None
    def select_ires(self, dobj):
        if self._ires is None:
            raise RuntimeError("ARTIOChunk.select_ires called before fill")
        return self._ires

    def select_fwidth(self, dobj):
        if self._ires is None:
            raise RuntimeErorr("ARTIOChunk.fwidth called before fill")
        return np.array([2.**-self._ires, 2.**-self._ires,
                         2.**-self._ires]).transpose()

    def select_icoords(self, dobj):
        if self._fcoords is None or self._ires is None:
            raise RuntimeError("ARTIOChunk.icoords called before fill")
        return (int)(self._fcoords/2**-self._ires)

    def fill(self, fields, selector):
        art_fields = [yt_to_art[f[1]] for f in fields]
        (self._fcoords, self._ires, artdata) = \
            self.pf._handle.read_grid_chunk(selector,
                                            self.sfc_start,
                                            self.sfc_end, art_fields)
        data = {}
        for i, f in enumerate(fields):
            data[f] = artdata[i]
        self._data_size = len(self._fcoords)
        return data

    def fill_particles(self, field_data, fields, selector):
        art_fields = {}
        for s, f in fields:
            for i in range(self.pf.num_species):
                if s == "all" or self.pf.particle_species[i] == yt_to_art[s]:
                    if yt_to_art[f] in self.pf.particle_variables[i]:
                        art_fields[(i, yt_to_art[f])] = 1

        species_data = self.pf._handle.read_particle_chunk(
            selector, self.sfc_start, self.sfc_end, art_fields.keys())

        for s, f in fields:
            af = yt_to_art[f]
            np = sum(len(species_data[(i, af)])
                     for i in range(self.pf.num_species)
                     if s == "all"
                     or self.pf.particle_species[i] == yt_to_art[s])

            cp = len(field_data[(s, f)])
            field_data[(s, f)].resize(cp + np)
            for i in range(self.pf.num_species):
                if s == "all" or self.pf.particle_species[i] == yt_to_art[s]:
                    np = len(species_data[(i, yt_to_art[f])])
                    field_data[(s, f)][cp:cp+np] = \
                        species_data[(i, yt_to_art[f])]
                    cp += np


class ARTIOGeometryHandler(GeometryHandler):

    def __init__(self, pf, data_style='artio'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)

        self.max_level = pf.max_level
        self.float_type = np.float64
        super(ARTIOGeometryHandler, self).__init__(pf, data_style)

    def _setup_geometry(self):
        mylog.debug("Initializing Geometry Handler empty for now.")

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return  1.0/(2**self.max_level)

    def convert(self, unit):
        return self.parameter_file.conversion_factors[unit]

    def find_max(self, field, finest_levels=3):
        """
        Returns (value, center) of location of maximum for a given field.
        """
        if (field, finest_levels) in self._max_locations:
            return self._max_locations[(field, finest_levels)]
        mv, pos = self.find_max_cell_location(field, finest_levels)
        self._max_locations[(field, finest_levels)] = (mv, pos)
        return mv, pos

    def find_max_cell_location(self, field, finest_levels=3):
        source = self.all_data()
        if finest_levels is not False:
            source.min_level = self.max_level - finest_levels
        mylog.debug("Searching for maximum value of %s", field)
        max_val, maxi, mx, my, mz = \
            source.quantities["MaxLocation"](field)
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f",
                   max_val, mx, my, mz)
        self.pf.parameters["Max%sValue" % (field)] = max_val
        self.pf.parameters["Max%sPos" % (field)] = "%s" % ((mx, my, mz),)
        return max_val, np.array((mx, my, mz), dtype='float64')

    def _detect_fields(self):
        self.fluid_field_list = self._detect_fluid_fields()
        self.particle_field_list = self._detect_particle_fields()
        self.field_list = self.fluid_field_list + self.particle_field_list
        mylog.debug("Detected fields:", (self.field_list,))

    def _detect_fluid_fields(self):
        return [art_to_yt[f] for f in yt_to_art.values() if f in
                self.pf.artio_parameters["grid_variable_labels"]]

    def _detect_particle_fields(self):
        fields = set()
        for ptype in self.pf.particle_types:
            if ptype == "all": continue
            for f in yt_to_art.values():
                if all(f in self.pf.particle_variables[i]
                       for i in range(self.pf.num_species)
                       if art_to_yt[self.pf.particle_species[i]] == ptype):
                    fields.add((ptype, art_to_yt[f]))
        return list(fields)

    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ARTIOGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            try:
                all_data = all(dobj.left_edge == self.pf.domain_left_edge) and\
                    all(dobj.right_edge == self.pf.domain_right_edge)
            except:
                all_data = False
            base_region = getattr(dobj, "base_region", dobj)
            sfc_start = getattr(dobj, "sfc_start", None)
            sfc_end = getattr(dobj, "sfc_end", None)
            domain = getattr(dobj, "domain_id", 0)
            if all_data:
                mylog.debug("Selecting entire artio domain")
                list_sfc_ranges = self.pf._handle.root_sfc_ranges_all()
            elif sfc_start is not None and sfc_end is not None:
                mylog.debug("Restricting to %s .. %s", sfc_start, sfc_end)
                list_sfc_ranges = [(sfc_start, sfc_end)]
            else:
                mylog.debug("Running selector on artio base grid")
                list_sfc_ranges = self.pf._handle.root_sfc_ranges(
                    dobj.selector)
            ci = []
            if domain != 2:
                ci += [ARTIORootMeshSubset(base_region, start, end, self.pf)
                        for (start, end) in list_sfc_ranges]
            if domain != 1:
                ci += [ARTIOOctreeSubset(base_region, start, end, self.pf)
                       for (start, end) in list_sfc_ranges]
            dobj._chunk_info = ci
            mylog.info("Created %d chunks for ARTIO" % len(list_sfc_ranges))
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _data_size(self, dobj, dobjs):
        size = 0
        for d in dobjs:
            size += d.data_size
        return size

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz):
        if ngz > 0:
            raise NotImplementedError
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        # These are ARTIOChunk objects
        for i,og in enumerate(sobjs):
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g], None)

    def _chunk_io(self, dobj, cache = True):
        # _current_chunk is made from identify_base_chunk
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for chunk in oobjs:
            yield YTDataChunk(dobj, "io", [chunk], None,
                              cache = cache)

    def _read_fluid_fields(self, fields, dobj, chunk=None):
        if len(fields) == 0:
            return {}, []
        if chunk is None:
            self._identify_base_chunk(dobj)
        fields_to_return = {}
        fields_to_read, fields_to_generate = self._split_fields(fields)
        if len(fields_to_read) == 0:
            return {}, fields_to_generate
        fields_to_return = self.io._read_fluid_selection(self._chunk_io(dobj),
                                                         dobj.selector,
                                                         fields_to_read)
        for field in fields_to_read:
            ftype, fname = field
            conv_factor = self.pf.field_info[fname]._convert_function(self)
            np.multiply(fields_to_return[field], conv_factor,
                        fields_to_return[field])
        return fields_to_return, fields_to_generate


class ARTIOStaticOutput(StaticOutput):
    _handle = None
    _hierarchy_class = ARTIOGeometryHandler
    _fieldinfo_fallback = ARTIOFieldInfo
    _fieldinfo_known = KnownARTIOFields
    _particle_mass_name = "particle_mass"
    _particle_coordinates_name = "Coordinates"

    def __init__(self, filename, data_style='artio',
                 storage_filename=None):
        if self._handle is not None:
            return
        self._filename = filename
        self._fileset_prefix = filename[:-4]
        self._handle = artio_fileset(self._fileset_prefix)
        self.artio_parameters = self._handle.parameters
        # Here we want to initiate a traceback, if the reader is not built.
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

    def _set_units(self):
        """
        Generates the conversion to physical units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        for unit in mpc_conversion.keys():
            self.units[unit] = self.parameters['unit_l']\
                * mpc_conversion[unit] / mpc_conversion["cm"]
        if self.cosmological_simulation:
            for unit in mpc_conversion:
                self.units["%sh" % unit] = self.units[unit] * \
                    self.hubble_constant
                self.units["%shcm" % unit] = \
                    (self.units["%sh" % unit] /
                        (1 + self.current_redshift))
                self.units["%scm" % unit] = \
                    self.units[unit] / (1 + self.current_redshift)

        for unit in sec_conversion.keys():
            self.time_units[unit] = self.parameters['unit_t']\
                / sec_conversion[unit]

        constants = ARTIOconstants()
        mb = constants.XH*constants.mH + constants.XHe * constants.mHe

        self.parameters['unit_d'] = self.parameters['unit_m']\
            / self.parameters['unit_l']**3.0
        self.parameters['unit_v'] = self.parameters['unit_l']\
            / self.parameters['unit_t']
        self.parameters['unit_E'] = self.parameters['unit_m']\
            * self.parameters['unit_v']**2.0
        self.parameters['unit_T'] = self.parameters['unit_v']**2.0*mb\
            / constants.k
        self.parameters['unit_rhoE'] = self.parameters['unit_E']\
            / self.parameters['unit_l']**3.0
        self.parameters['unit_nden'] = self.parameters['unit_d'] / mb
        self.parameters['Gamma'] = constants.gamma

        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge -
                                       self.domain_left_edge).max()
        self.conversion_factors["Density"] = self.parameters['unit_d']
        self.conversion_factors["x-velocity"] = self.parameters['unit_v']
        self.conversion_factors["y-velocity"] = self.parameters['unit_v']
        self.conversion_factors["z-velocity"] = self.parameters['unit_v']
        #*cell_gas_internal_energy(cell)/cell_gas_density(cell);
        self.conversion_factors["Temperature"] = \
            self.parameters['unit_T'] * constants.wmu * (constants.gamma-1)

        mylog.info('note artio T uses fixed gamma not variable')

        for ax in 'xyz':
            self.conversion_factors["particle_velocity_%s" % ax] =\
                self.parameters['unit_v']
        for unit in sec_conversion.keys():
            self.time_units[unit] = 1.0 / sec_conversion[unit]
        self.conversion_factors['particle_mass'] = self.parameters['unit_m']
        self.conversion_factors['particle_creation_time'] =\
            self.parameters['unit_t']
        self.conversion_factors['particle_mass_msun'] =\
            self.parameters['unit_m'] / constants.Msun

        #for mult_halo_profiler.py:
        self.parameters['TopGridDimensions'] = 3 * [self._handle.num_grid]
        self.parameters['RefineBy'] = 2
        self.parameters['DomainLeftEdge'] = 3 * [0]
        self.parameters['DomainRightEdge'] = 3 * [self._handle.num_grid]
        self.parameters['TopGridRank'] = 3  # number of dimensions

    def _parse_parameter_file(self):
        # hard-coded -- not provided by headers
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'artio'
        self.parameters["Time"] = 1.  # default unit is 1...

        # read header
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])

        self.num_grid = self._handle.num_grid
        self.domain_dimensions = np.ones(3, dtype='int32') * self.num_grid
        self.domain_left_edge = np.zeros(3, dtype="float64")
        self.domain_right_edge = np.ones(3, dtype='float64')*self.num_grid

        # TODO: detect if grid exists
        self.min_level = 0  # ART has min_level=0
        self.max_level = self.artio_parameters["max_refinement_level"][0]

        # TODO: detect if particles exist
        self.num_species = self.artio_parameters["num_particle_species"][0]
        self.particle_variables = [["PID", "SPECIES"]
                                   for i in range(self.num_species)]
        self.particle_species =\
            self.artio_parameters["particle_species_labels"]
        self.particle_type_map = {}
        for i, s in enumerate(self.particle_species):
            f = art_to_yt[s]
            if f not in self.particle_type_map:
                self.particle_type_map[f] = []
            self.particle_type_map[f].append(i)

        for species in range(self.num_species):
            # Mass would be best as a derived field,
            # but wouldn't detect under 'all'
            if self.artio_parameters["particle_species_labels"][species]\
                    == "N-BODY":
                self.particle_variables[species].append("MASS")

            if self.artio_parameters["num_primary_variables"][species] > 0:
                self.particle_variables[species].extend(
                    self.artio_parameters[
                        "species_%02d_primary_variable_labels"
                        % (species, )])
            if self.artio_parameters["num_secondary_variables"][species] > 0:
                self.particle_variables[species].extend(
                    self.artio_parameters[
                        "species_%02d_secondary_variable_labels"
                        % (species, )])

        self.particle_types = ("all",) + tuple(
            set(art_to_yt[s] for s in
                self.artio_parameters["particle_species_labels"]))
        self.particle_types = tuple(self.particle_types)

        self.current_time = b2t(self.artio_parameters["tl"][0])

        # detect cosmology
        if "abox" in self.artio_parameters:
            abox = self.artio_parameters["abox"][0]
            self.cosmological_simulation = True
            self.omega_lambda = self.artio_parameters["OmegaL"][0]
            self.omega_matter = self.artio_parameters["OmegaM"][0]
            self.hubble_constant = self.artio_parameters["hubble"][0]
            self.current_redshift = 1.0/self.artio_parameters["abox"][0] - 1.0

            self.parameters["initial_redshift"] =\
                1.0 / self.artio_parameters["auni_init"][0] - 1.0
            self.parameters["CosmologyInitialRedshift"] =\
                self.parameters["initial_redshift"]
        else:
            self.cosmological_simulation = False

        #units
        if self.cosmological_simulation:
            self.parameters['unit_m'] = self.artio_parameters["mass_unit"][0]
            self.parameters['unit_t'] =\
                self.artio_parameters["time_unit"][0] * abox**2
            self.parameters['unit_l'] =\
                self.artio_parameters["length_unit"][0] * abox
        else:
            self.parameters['unit_l'] = self.artio_parameters["length_unit"][0]
            self.parameters['unit_t'] = self.artio_parameters["time_unit"][0]
            self.parameters['unit_m'] = self.artio_parameters["mass_unit"][0]

        # hard coded assumption of 3D periodicity (add to parameter file)
        self.periodicity = (True, True, True)

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"):
            return False
        return artio_is_valid(args[0][:-4])
