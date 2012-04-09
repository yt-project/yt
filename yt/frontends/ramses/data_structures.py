"""
RAMSES-specific data structures

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

import numpy as na
import stat
import weakref

from yt.funcs import *
from yt.data_objects.grid_patch import \
      AMRGridPatch
from yt.geometry.oct_geometry_handler import \
      OctreeGeometryHandler
from yt.data_objects.static_output import \
      StaticOutput

from .fields import RAMSESFieldInfo, KnownRAMSESFields
from .definitions import ramses_header
from yt.utilities.definitions import \
    mpc_conversion
from yt.utilities.amr_utils import \
    get_box_grids_level
from yt.utilities.io_handler import \
    io_registry
from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc
import yt.utilities.fortran_utils as fpu
from yt.geometry.oct_container import \
    OctreeContainer

class RAMSESGeometryHandler(OctreeGeometryHandler):

    def __init__(self, pf, data_style='ramses'):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)

        self.float_type = na.float64
        super(RAMSESGeometryHandler, self).__init__(pf, data_style)

    def _initialize_oct_handler(self):
        # We have a standard filename 
        base = os.path.dirname(self.parameter_file.parameter_filename)
        fn = self.parameter_file.parameter_filename
        output_id = os.path.basename(fn).split(".")[0].split("_")[1]
        self.oct_handler = None
        for i in range(self.parameter_file['ncpu']):
            fn = os.path.join(base, "amr_%s.out%05i" % (output_id, i + 1))
            self._read_domain(i + 1, fn)

    def _read_domain(self, domain, domain_fn):
        hvals = {}
        f = open(domain_fn, "rb")
        for header in ramses_header(hvals):
            hvals.update(fpu.read_attrs(f, header))
        # That's the header, now we skip a few.
        hvals['numbl'] = na.array(hvals['numbl']).reshape(
            (hvals['nlevelmax'], hvals['ncpu']))
        fpu.skip(f)
        if hvals['nboundary'] > 0:
            fpu.skip(f, 2)
            ngridbound = fpu.read_vector(f, 'i')
        free_mem = fpu.read_attrs(f, (('free_mem', 5, 'i'), ) )
        ordering = fpu.read_vector(f, 'c')
        fpu.skip(f, 4)
        # Now we're at the tree itself
        if self.oct_handler is None:
            self.oct_handler = OctreeContainer(hvals['nx'],
                self.parameter_file.domain_left_edge,
                self.parameter_file.domain_right_edge)
        # Now we iterate over each level and each CPU.
        mylog.debug("Inspecting domain % 4i", domain)
        for level in range(hvals['nlevelmax']):
            # Easier if do this 1-indexed
            for cpu in range(hvals['nboundary'] + hvals['ncpu']):
                if cpu < hvals['ncpu']:
                    ng = hvals['numbl'][level, cpu]
                else:
                    ng = ngridbound[cpu - hvals['ncpu'] + hvals['nboundary']*level]
                if ng == 0: continue
                ind = fpu.read_vector(f, "I").astype("int64")
                fpu.skip(f, 2)
                pos = na.empty((ng, 3), dtype='float64')
                v1 = fpu.read_vector(f, "d")
                v2 = fpu.read_vector(f, "d")
                v3 = fpu.read_vector(f, "d")
                pos[:,0] = v1
                pos[:,1] = v2
                pos[:,2] = v3
                parents = fpu.read_vector(f, "I")
                fpu.skip(f, 6)
                children = na.empty((ng, 8), dtype='int64')
                for i in range(8):
                    children[:,i] = fpu.read_vector(f, "I")
                cpu_map = na.empty((ng, 8), dtype="int64")
                for i in range(8):
                    cpu_map[:,i] = fpu.read_vector(f, "I")
                rmap = na.empty((ng, 8), dtype="int64")
                for i in range(8):
                    rmap[:,i] = fpu.read_vector(f, "I")
                if cpu + 1 == domain:
                    self.oct_handler.add_ramses(domain, level, ng, pos, ind, cpu_map)

    def _detect_fields(self):
        # TODO: Add additional fields
        self.field_list = [ "Density", "x-velocity", "y-velocity",
	                        "z-velocity", "Pressure", "Metallicity"]
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(RAMSESGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

class RAMSESStaticOutput(StaticOutput):
    _hierarchy_class = RAMSESGeometryHandler
    _fieldinfo_fallback = RAMSESFieldInfo
    _fieldinfo_known = KnownRAMSESFields
    
    def __init__(self, filename, data_style='ramses',
                 storage_filename = None):
        # Here we want to initiate a traceback, if the reader is not built.
        StaticOutput.__init__(self, filename, data_style)
        self.storage_filename = storage_filename

    def __repr__(self):
        return self.basename.rsplit(".", 1)[0]
        
    def _set_units(self):
        """
        Generates the conversion to various physical _units based on the parameter file
        """
        self.units = {}
        self.time_units = {}
        if len(self.parameters) == 0:
            self._parse_parameter_file()
        self._setup_nounits_units()
        self.conversion_factors = defaultdict(lambda: 1.0)
        self.time_units['1'] = 1
        self.units['1'] = 1.0
        self.units['unitary'] = 1.0 / (self.domain_right_edge - self.domain_left_edge).max()
        seconds = self.parameters['unit_t']
        self.time_units['years'] = seconds / (365*3600*24.0)
        self.time_units['days']  = seconds / (3600*24.0)
        self.time_units['Myr'] = self.time_units['years'] / 1.0e6
        self.time_units['Gyr']  = self.time_units['years'] / 1.0e9
        self.conversion_factors["Density"] = self.parameters['unit_d']
        vel_u = self.parameters['unit_l'] / self.parameters['unit_t']
        self.conversion_factors["x-velocity"] = vel_u
        self.conversion_factors["y-velocity"] = vel_u
        self.conversion_factors["z-velocity"] = vel_u

    def _setup_nounits_units(self):
        for unit in mpc_conversion.keys():
            self.units[unit] = self.parameters['unit_l'] * mpc_conversion[unit] / mpc_conversion["cm"]

    def _parse_parameter_file(self):
        # hardcoded for now
        # These should be explicitly obtained from the file, but for now that
        # will wait until a reorganization of the source tree and better
        # generalization.
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = 'ramses'
        self.parameters["Time"] = 1. # default unit is 1...

        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # We now execute the same logic Oliver's code does
        rheader = {}
        f = open(self.parameter_filename)
        def read_rhs(cast):
            line = f.readline()
            p, v = line.split("=")
            rheader[p.strip()] = cast(v)
        for i in range(6): read_rhs(int)
        f.readline()
        for i in range(11): read_rhs(float)
        f.readline()
        read_rhs(str)
        # Now we read the hilber indices
        self.hilbert_indices = {}
        if rheader['ordering type'] == "hilbert":
            f.readline() # header
            for n in range(rheader['ncpu']):
                dom, mi, ma = f.readline().split()
                self.hilbert_indices[int(dom)] = (float(mi), float(ma))
        self.parameters.update(rheader)
        self.current_time = self.parameters['time'] * self.parameters['unit_t']
        self.domain_right_edge = na.ones(3, dtype='float64') \
                                           * rheader['boxlen']
        self.domain_left_edge = na.zeros(3, dtype='float64')
        self.domain_dimensions = na.ones(3, dtype='int32') * 2
        # This is likely not true, but I am not sure how to otherwise
        # distinguish them.
        mylog.warning("No current mechanism of distinguishing cosmological simulations in RAMSES!")
        self.cosmological_simulation = 1
        self.current_redshift = (1.0 / rheader["aexp"]) - 1.0
        self.omega_lambda = rheader["omega_l"]
        self.omega_matter = rheader["omega_m"]
        self.hubble_constant = rheader["H0"]

    @classmethod
    def _is_valid(self, *args, **kwargs):
        if not os.path.basename(args[0]).startswith("info_"): return False
        fn = args[0].replace("info_", "amr_").replace(".txt", ".out00001")
        print fn
        return os.path.exists(fn)

