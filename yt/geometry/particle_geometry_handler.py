"""
Particle-only geometry handler




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import collections
import numpy as np
import os
import weakref

from yt.funcs import only_on_root
from yt.utilities.logger import ytLogger as mylog
from yt.data_objects.octree_subset import ParticleOctreeSubset
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer, ParticleRegions

class ParticleIndex(Index):
    """The Index subclass for particle datasets"""
    _global_mesh = False

    def __init__(self, ds, dataset_type):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64
        super(ParticleIndex, self).__init__(ds, dataset_type)

    @property
    def index_ptype(self):
        if hasattr(self.dataset, "index_ptype"):
            return self.dataset.index_ptype
        else:
            return "all"

    def _setup_geometry(self):
        mylog.debug("Initializing Particle Geometry Handler.")
        self._initialize_particle_handler()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        ML = self.oct_handler.max_level
        dx = 1.0/(self.dataset.domain_dimensions*2**ML)
        dx = dx * (self.dataset.domain_right_edge -
                   self.dataset.domain_left_edge)
        return dx.min()

    def _get_particle_type_counts(self):
        result = collections.defaultdict(lambda: 0)
        for df in self.data_files:
            for k in df.total_particles.keys():
                result[k] += df.total_particles[k]
        return dict(result)

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _initialize_particle_handler(self):
        self._setup_data_io()
        template = self.dataset.filename_template
        ndoms = self.dataset.file_count
        cls = self.dataset._file_class
        self.data_files = [cls(self.dataset, self.io, template % {'num':i}, i)
                           for i in range(ndoms)]
        index_ptype = self.index_ptype
        if index_ptype == "all":
            self.total_particles = sum(
                    sum(d.total_particles.values()) for d in self.data_files)
        else:
            self.total_particles = sum(
                    d.total_particles[index_ptype] for d in self.data_files)
        ds = self.dataset
        self.oct_handler = ParticleOctreeContainer(
            [1, 1, 1], ds.domain_left_edge, ds.domain_right_edge,
            over_refine = ds.over_refine_factor)
        self.oct_handler.n_ref = ds.n_ref
        only_on_root(mylog.info, "Allocating for %0.3e particles "
                                 "(index particle type '%s')",
                     self.total_particles, index_ptype)
        # No more than 256^3 in the region finder.
        N = min(len(self.data_files), 256) 
        self.regions = ParticleRegions(
                ds.domain_left_edge, ds.domain_right_edge,
                [N, N, N], len(self.data_files))
        self._initialize_indices()
        self.oct_handler.finalize()
        self.max_level = self.oct_handler.max_level
        self.dataset.max_level = self.max_level
        tot = sum(self.oct_handler.recursively_count().values())
        only_on_root(mylog.info, "Identified %0.3e octs", tot)

    def _initialize_indices(self):
        # This will be replaced with a parallel-aware iteration step.
        # Roughly outlined, what we will do is:
        #   * Generate Morton indices on each set of files that belong to
        #     an individual processor
        #   * Create a global, accumulated histogram
        #   * Cut based on estimated load balancing
        #   * Pass particles to specific processors, along with NREF buffer
        #   * Broadcast back a serialized octree to join
        #
        # For now we will do this in serial.
        index_ptype = self.index_ptype
        # Set the index_ptype attribute of self.io dynamically here, so we don't
        # need to assume that the dataset has the attribute.
        self.io.index_ptype = index_ptype
        morton = np.empty(self.total_particles, dtype="uint64")
        ind = 0
        for data_file in self.data_files:
            if index_ptype == "all":
                npart = sum(data_file.total_particles.values())
            else:
                npart = data_file.total_particles[index_ptype]
            morton[ind:ind + npart] = \
                self.io._initialize_index(data_file, self.regions)
            ind += npart
        morton.sort()
        # Now we add them all at once.
        self.oct_handler.add(morton)

    def _detect_output_fields(self):
        # TODO: Add additional fields
        dsl = []
        units = {}
        for dom in self.data_files:
            fl, _units = self.io._identify_fields(dom)
            units.update(_units)
            dom._calculate_offsets(fl)
            for f in fl:
                if f not in dsl: dsl.append(f)
        self.field_list = dsl
        ds = self.dataset
        ds.particle_types = tuple(set(pt for pt, ds in dsl))
        # This is an attribute that means these particle types *actually*
        # exist.  As in, they are real, in the dataset.
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types

    def _identify_base_chunk(self, dobj):
        if getattr(dobj, "_chunk_info", None) is None:
            data_files = getattr(dobj, "data_files", None)
            if data_files is None:
                data_files = [self.data_files[i] for i in
                              self.regions.identify_data_files(dobj.selector)]
            base_region = getattr(dobj, "base_region", dobj)
            oref = self.dataset.over_refine_factor
            subset = [ParticleOctreeSubset(base_region, data_files, 
                        self.dataset, over_refine_factor = oref)]
            dobj._chunk_info = subset
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        # We actually do not really use the data files except as input to the
        # ParticleOctreeSubset.
        # This is where we will perform cutting of the Octree and
        # load-balancing.  That may require a specialized selector object to
        # cut based on some space-filling curve index.
        for i,og in enumerate(sobjs):
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache = True, local_only = False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in oobjs:
            yield YTDataChunk(dobj, "io", [subset], None, cache = cache)

class ParticleDataChunk(YTDataChunk):
    def __init__(self, oct_handler, regions, *args, **kwargs):
        self.oct_handler = oct_handler
        self.regions = regions
        super(ParticleDataChunk, self).__init__(*args, **kwargs)
