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

import numpy as np
import os
import weakref

from yt.funcs import get_pbar, only_on_root
from yt.utilities.logger import ytLogger as mylog
from yt.data_objects.octree_subset import ParticleOctreeSubset
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer, ParticleForest
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.extern.functools32 import lru_cache

from yt.data_objects.data_containers import data_object_registry
from yt.data_objects.octree_subset import ParticleOctreeSubset
from yt.data_objects.particle_container import ParticleContainer

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

    def _setup_geometry(self):
        self.regions = None
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

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _initialize_particle_handler(self):
        self._setup_data_io()
        template = self.dataset.filename_template
        ndoms = self.dataset.file_count
        cls = self.dataset._file_class
        self.data_files = [cls(self.dataset, self.io, template % {'num':i}, i)
                           for i in range(ndoms)]
        N = min(len(4*self.data_files), 256) 
        self.ds.domain_dimensions[:] = N*(1<<self.ds.over_refine_factor)
        self.total_particles = sum(
                sum(d.total_particles.values()) for d in self.data_files)

    def _initialize_coarse_index(self):
        ds = self.dataset
        only_on_root(mylog.info, "Allocating for %0.3e particles",
          self.total_particles)
        # No more than 256^3 in the region finder.
        N = self.ds.domain_dimensions / (1<<self.ds.over_refine_factor)
        self.regions = ParticleForest(
                ds.domain_left_edge, ds.domain_right_edge,
                N, len(self.data_files), ds.over_refine_factor,
                ds.n_ref)
        pb = get_pbar("Initializing coarse index ", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            for pos in self.io._yield_coordinates(data_file):
                self.regions.add_data_file(pos, data_file.file_id, 0, 0.0)
        pb.finish()

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
        if self.regions is None:
            self._initialize_coarse_index()
        if getattr(dobj, "_chunk_info", None) is None:
            data_files = getattr(dobj, "data_files", None)
            buffer_files = getattr(dobj, "buffer_files", None)
            if data_files is None:
                dfi, = np.where(self.regions.identify_data_files(dobj.selector))
                #n_cells = omask.sum()
                data_files = [self.data_files[i] for i in dfi]
                #mylog.debug("Maximum particle count of %s identified", count)
            base_region = getattr(dobj, "base_region", dobj)
            dobj._chunk_info = [ParticleOctreeSubset(dobj, df, self.ds,
                over_refine_factor = self.ds.over_refine_factor)
                for df in data_files]
        dobj._current_chunk = list(self._chunk_all(dobj))[0]

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None,
                       ghost_particles = False):
        ghost_particles = False
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
            with g._expand_data_files(ghost_particles):
                yield YTDataChunk(dobj, "spatial", [g])

    def old_chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None):
        dfi, count, omask = self.regions.identify_data_files(
                                dobj.selector)
        # We actually do not really use the data files except as input to the
        # ParticleOctreeSubset.
        # This is where we will perform cutting of the Octree and
        # load-balancing.  That may require a specialized selector object to
        # cut based on some space-filling curve index.
        for df in (self.data_files[i] for i in dfi):
            if ngz > 0:
                raise NotImplementedError
            else:
                oct_handler = self.regions.construct_forest(
                        df.file_id, dobj.selector, self.io, self.data_files,
                        (dfi, count, omask))
                g = ParticleOctreeSubset(dobj, df, self.ds,
                        over_refine_factor = self.ds.over_refine_factor)
            yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache = True, local_only = False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for container in oobjs:
            yield YTDataChunk(dobj, "io", [container], None, cache = cache)
