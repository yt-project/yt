"""
Data structures for a generic SPH/Gadget frontend.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2012 Matthew Turk.  All Rights Reserved.

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

import h5py
import numpy as np
import stat
import weakref
import struct
from itertools import izip

from yt.utilities.fortran_utils import read_record
from yt.funcs import *
from yt.geometry.oct_geometry_handler import \
    OctreeGeometryHandler
from yt.geometry.oct_container import \
    ParticleOctreeContainer
from yt.geometry.geometry_handler import \
    GeometryHandler, YTDataChunk
from yt.data_objects.static_output import \
    StaticOutput
from yt.data_objects.octree_subset import \
    OctreeSubset
from yt.utilities.definitions import \
    mpc_conversion, sec_conversion
from .fields import \
    OWLSFieldInfo, \
    KnownOWLSFields, \
    GadgetFieldInfo, \
    KnownGadgetFields, \
    TipsyFieldInfo, \
    KnownTipsyFields

from yt.data_objects.field_info_container import \
    FieldInfoContainer, NullFunc

class ParticleDomainFile(object):
    def __init__(self, pf, io, domain_filename, domain_id):
        self.pf = pf
        self.io = weakref.proxy(io)
        self.domain_filename = domain_filename
        self.domain_id = domain_id
        self.total_particles = self.io._count_particles(self)

    def select(self, selector):
        pass

    def count(self, selector):
        pass

    def _calculate_offsets(self, fields):
        pass

class ParticleDomainSubset(OctreeSubset):
    pass

class ParticleGeometryHandler(OctreeGeometryHandler):

    def __init__(self, pf, data_style):
        self.data_style = data_style
        self.parameter_file = weakref.proxy(pf)
        # for now, the hierarchy file is the parameter file!
        self.hierarchy_filename = self.parameter_file.parameter_filename
        self.directory = os.path.dirname(self.hierarchy_filename)
        self.float_type = np.float64
        super(ParticleGeometryHandler, self).__init__(pf, data_style)
        
    def _initialize_oct_handler(self):
        self._setup_data_io()
        template = self.parameter_file.domain_template
        ndoms = self.parameter_file.domain_count
        cls = self.parameter_file._domain_class
        self.domains = [cls(self.parameter_file, self.io, template % {'num':i}, i)
                        for i in range(ndoms)]
        total_particles = sum(sum(d.total_particles.values())
                              for d in self.domains)
        self.oct_handler = ParticleOctreeContainer(
            self.parameter_file.domain_dimensions/2,
            self.parameter_file.domain_left_edge,
            self.parameter_file.domain_right_edge)
        self.oct_handler.n_ref = 64
        mylog.info("Allocating for %0.3e particles", total_particles)
        for dom in self.domains:
            self.io._initialize_octree(dom, self.oct_handler)
        self.oct_handler.finalize()
        self.max_level = self.oct_handler.max_level
        tot = self.oct_handler.linearly_count()
        mylog.info("Identified %0.3e octs", tot)

    def _detect_fields(self):
        # TODO: Add additional fields
        pfl = []
        for dom in self.domains:
            fl = self.io._identify_fields(dom)
            dom._calculate_offsets(fl)
            for f in fl:
                if f not in pfl: pfl.append(f)
        self.field_list = pfl
        pf = self.parameter_file
        pf.particle_types = tuple(set(pt for pt, pf in pfl))
    
    def _setup_classes(self):
        dd = self._get_data_reader_dict()
        super(ParticleGeometryHandler, self)._setup_classes(dd)
        self.object_types.sort()

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            mask = dobj.selector.select_octs(self.oct_handler)
            counts = self.oct_handler.count_cells(dobj.selector, mask)
            subsets = [ParticleDomainSubset(d, mask, c)
                       for d, c in zip(self.domains, counts) if c > 0]
            dobj._chunk_info = subsets
            dobj.size = sum(counts)
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, dobj.size)

    def _chunk_spatial(self, dobj, ngz, sort = None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for i,og in enumerate(sobjs):
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            size = og.cell_count
            if size == 0: continue
            yield YTDataChunk(dobj, "spatial", [g], size)

    def _chunk_io(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], subset.cell_count)

class OWLSStaticOutput(StaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _domain_class = ParticleDomainFile
    _fieldinfo_fallback = OWLSFieldInfo
    _fieldinfo_known = KnownOWLSFields

    def __init__(self, filename, data_style="OWLS", root_dimensions = 64):
        self._root_dimensions = root_dimensions
        # Set up the template for domain files
        self.storage_filename = None
        super(OWLSStaticOutput, self).__init__(filename, data_style)

    def __repr__(self):
        return os.path.basename(self.parameter_filename).split(".")[0]

    def _set_units(self):
        self.units = {}
        self.time_units = {}
        self.conversion_factors = {}
        DW = self.domain_right_edge - self.domain_left_edge
        self.units["unitary"] = 1.0 / DW.max()

    def _parse_parameter_file(self):
        handle = h5py.File(self.parameter_filename)
        hvals = {}
        hvals.update(handle["/Header"].attrs)

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # Set standard values
        self.current_time = hvals["Time_GYR"] * sec_conversion["Gyr"]
        self.domain_left_edge = np.zeros(3, "float64")
        self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]
        self.domain_dimensions = np.ones(3, "int32") * self._root_dimensions
        self.cosmological_simulation = 1
        self.periodicity = (True, True, True)
        self.current_redshift = hvals["Redshift"]
        self.omega_lambda = hvals["OmegaLambda"]
        self.omega_matter = hvals["Omega0"]
        self.hubble_constant = hvals["HubbleParam"]
        self.parameters = hvals

        prefix = self.parameter_filename.split(".", 1)[0]
        suffix = self.parameter_filename.rsplit(".", 1)[-1]
        self.domain_template = "%s.%%(num)i.%s" % (prefix, suffix)
        self.domain_count = hvals["NumFilesPerSnapshot"]

        handle.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        try:
            fileh = h5py.File(args[0],'r')
            if "Constants" in fileh["/"].keys() and \
               "Header" in fileh["/"].keys():
                fileh.close()
                return True
            fileh.close()
        except:
            pass
        return False

class GadgetBinaryDomainFile(ParticleDomainFile):
    def __init__(self, pf, io, domain_filename, domain_id):
        with open(domain_filename, "rb") as f:
            self.header = read_record(f, pf._header_spec)
            self._position_offset = f.tell()

        super(GadgetBinaryDomainFile, self).__init__(pf, io,
                domain_filename, domain_id)

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_field_offsets(
                field_list, self.total_particles)

class ParticleStaticOutput(StaticOutput):
    _unit_base = None

    def _set_units(self):
        self.units = {}
        self.time_units = {}
        self.conversion_factors = {}
        self.units['1'] = 1.0
        self.units['unitary'] = (self.domain_right_edge -
                                 self.domain_left_edge).max()
        # Check 
        base = None
        mpch = {}
        mpch.update(mpc_conversion)
        unit_base = self._unit_base or {}
        for unit in mpc_conversion:
            mpch['%sh' % unit] = mpch[unit] / self.hubble_constant
            mpch['%shcm' % unit] = (mpch["%sh" % unit] / 
                    (1 + self.current_redshift))
            mpch['%scm' % unit] = mpch[unit] / (1 + self.current_redshift)
        for unit_registry in [mpch, sec_conversion]:
            for unit in sorted(unit_base):
                if unit in unit_registry:
                    base = unit_registry[unit] / unit_registry['mpc'] 
                    break
            if base is None: continue
            for unit in unit_registry:
                self.units[unit] = unit_registry[unit] / base

class GadgetStaticOutput(ParticleStaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _domain_class = GadgetBinaryDomainFile
    _fieldinfo_fallback = GadgetFieldInfo
    _fieldinfo_known = KnownGadgetFields
    _header_spec = (('Npart', 6, 'i'),
                    ('Massarr', 6, 'd'),
                    ('Time', 1, 'd'),
                    ('Redshift', 1, 'd'),
                    ('FlagSfr', 1, 'i'),
                    ('FlagFeedback', 1, 'i'),
                    ('Nall', 6, 'i'),
                    ('FlagCooling', 1, 'i'),
                    ('NumFiles', 1, 'i'),
                    ('BoxSize', 1, 'd'),
                    ('Omega0', 1, 'd'),
                    ('OmegaLambda', 1, 'd'),
                    ('HubbleParam', 1, 'd'),
                    ('FlagAge', 1, 'i'),
                    ('FlagMEtals', 1, 'i'),
                    ('NallHW', 6, 'i'),
                    ('unused', 16, 'i') )

    def __init__(self, filename, data_style="gadget_binary",
                 additional_fields = (), root_dimensions = 64,
                 unit_base = None):
        self._root_dimensions = root_dimensions
        # Set up the template for domain files
        self.storage_filename = None
        self._unit_base = unit_base
        super(GadgetStaticOutput, self).__init__(filename, data_style)

    def __repr__(self):
        return os.path.basename(self.parameter_filename).split(".")[0]

    def _parse_parameter_file(self):

        # The entries in this header are capitalized and named to match Table 4
        # in the GADGET-2 user guide.

        f = open(self.parameter_filename)
        hvals = read_record(f, self._header_spec)
        for i in hvals:
            if len(hvals[i]) == 1:
                hvals[i] = hvals[i][0]
        
        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # Set standard values

        # This may not be correct.
        self.current_time = hvals["Time"] * sec_conversion["Gyr"]

        self.domain_left_edge = np.zeros(3, "float64")
        self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]
        self.domain_dimensions = np.ones(3, "int32") * self._root_dimensions
        self.periodicity = (True, True, True)

        self.cosmological_simulation = 1

        self.current_redshift = hvals["Redshift"]
        self.omega_lambda = hvals["OmegaLambda"]
        self.omega_matter = hvals["Omega0"]
        self.hubble_constant = hvals["HubbleParam"]
        self.parameters = hvals

        prefix = self.parameter_filename.split(".", 1)[0]
        suffix = self.parameter_filename.rsplit(".", 1)[-1]

        if hvals["NumFiles"] > 1:
            self.domain_template = "%s.%%(num)s" % (prefix)
        else:
            self.domain_template = self.parameter_filename

        self.domain_count = hvals["NumFiles"]

        f.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We do not allow load() of these files.
        return False

class TipsyDomainFile(ParticleDomainFile):

    def _calculate_offsets(self, field_list):
        self.field_offsets = self.io._calculate_particle_offsets(self)

    def __init__(self, pf, io, domain_filename, domain_id):
        # To go above 1 domain, we need to include an indexing step in the
        # IOHandler, rather than simply reading from a single file.
        assert domain_id == 0 
        super(TipsyDomainFile, self).__init__(pf, io,
                domain_filename, domain_id)
        io._create_dtypes(self)


class TipsyStaticOutput(ParticleStaticOutput):
    _hierarchy_class = ParticleGeometryHandler
    _domain_class = TipsyDomainFile
    _fieldinfo_fallback = TipsyFieldInfo
    _fieldinfo_known = KnownTipsyFields
    _header_spec = (('time',    'd'),
                    ('nbodies', 'i'),
                    ('ndim',    'i'),
                    ('nsph',    'i'),
                    ('ndark',   'i'),
                    ('nstar',   'i'),
                    ('dummy',   'i'))

    def __init__(self, filename, data_style="tipsy",
                 root_dimensions = 64, endian = ">",
                 field_dtypes = None,
                 domain_left_edge = None,
                 domain_right_edge = None,
                 unit_base = None):
        self.endian = endian
        self._root_dimensions = root_dimensions
        # Set up the template for domain files
        self.storage_filename = None
        if domain_left_edge is None:
            domain_left_edge = np.zeros(3, "float64") - 0.5
        if domain_right_edge is None:
            domain_right_edge = np.zeros(3, "float64") + 0.5

        self.domain_left_edge = np.array(domain_left_edge, dtype="float64")
        self.domain_right_edge = np.array(domain_right_edge, dtype="float64")

        # My understanding is that dtypes are set on a field by field basis,
        # not on a (particle type, field) basis
        if field_dtypes is None: field_dtypes = {}
        self._field_dtypes = field_dtypes

        unit_base = self._unit_base or {}
        super(TipsyStaticOutput, self).__init__(filename, data_style)

    def __repr__(self):
        return os.path.basename(self.parameter_filename)

    def _set_units(self):
        self.units = {}
        self.time_units = {}
        self.conversion_factors = {}
        DW = self.domain_right_edge - self.domain_left_edge
        self.units["unitary"] = 1.0 / DW.max()

    def _parse_parameter_file(self):

        # The entries in this header are capitalized and named to match Table 4
        # in the GADGET-2 user guide.

        f = open(self.parameter_filename, "rb")
        hh = self.endian + "".join(["%s" % (b) for a,b in self._header_spec])
        hvals = dict([(a, c) for (a, b), c in zip(self._header_spec,
                     struct.unpack(hh, f.read(struct.calcsize(hh))))])
        self._header_offset = f.tell()

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        self.unique_identifier = \
            int(os.stat(self.parameter_filename)[stat.ST_CTIME])
        # Set standard values

        # This may not be correct.
        self.current_time = hvals["time"]

        # NOTE: These are now set in the main initializer.
        #self.domain_left_edge = np.zeros(3, "float64") - 0.5
        #self.domain_right_edge = np.ones(3, "float64") + 0.5
        self.domain_dimensions = np.ones(3, "int32") * self._root_dimensions
        self.periodicity = (True, True, True)

        self.cosmological_simulation = 1

        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0
        self.parameters = hvals

        self.domain_template = self.parameter_filename
        self.domain_count = 1

        f.close()

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # We do not allow load() of these files.
        return False
