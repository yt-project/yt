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

from yt.config import ytcfg
from yt.funcs import get_pbar, only_on_root, ensure_list
from yt.utilities.logger import ytLogger as mylog
from yt.data_objects.octree_subset import ParticleOctreeSubset
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.geometry.particle_oct_container import \
    ParticleOctreeContainer, ParticleBitmap
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.io_handler import io_registry
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, communication_system, MPI

from yt.data_objects.data_containers import data_object_registry
from yt.data_objects.octree_subset import ParticleOctreeSubset, _use_global_octree
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

    @property
    def index_ptype(self):
        if hasattr(self.dataset, "index_ptype"):
            return self.dataset.index_ptype
        else:
            return "all"

    def _setup_geometry(self):
        self.regions = None
        mylog.debug("Initializing Particle Geometry Handler.")
        self._initialize_particle_handler()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        # ML = self.oct_handler.max_level 
        ML = self.regions.index_order1 # was self.oct_handler.max_level
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
        self.data_files = []
        CHUNKSIZE = 64**3
        fi = 0
        for i in range(ndoms):
            start = 0
            end = start + CHUNKSIZE
            while 1:
                df = cls(self.dataset, self.io, template % {'num':i}, fi,
                        (start, end))
                if max(df.total_particles.values()) == 0:
                    break
                fi += 1
                self.data_files.append(df)
                start = end
                end += CHUNKSIZE
        #self.data_files = [cls(self.dataset, self.io, template % {'num':i}, i)
        #                   for i in range(ndoms)]
        self.total_particles = sum(
                sum(d.total_particles.values()) for d in self.data_files)
        # Get index & populate octree
        self._initialize_index()
        # Intialize global octree
        if _use_global_octree:
            index_ptype = self.index_ptype
            if index_ptype == "all":
                self.total_particles = sum(
                        sum(d.total_particles.values()) for d in self.data_files)
            else:
                self.total_particles = sum(
                        d.total_particles[index_ptype] for d in self.data_files)
            ds = self.dataset
            # TODO: Re-insert the usage of index_ptype here
            self.global_oct_handler = ParticleOctreeContainer(
                [1, 1, 1], ds.domain_left_edge, ds.domain_right_edge,
                over_refine = ds.over_refine_factor)
            self.global_oct_handler.n_ref = ds.n_ref
            self.global_oct_handler.add(self.regions.primary_indices(), 
                                        self.regions.index_order1)
            self.global_oct_handler.finalize()
            self.max_level = self.global_oct_handler.max_level
            for i in range(ndoms):
                fmask = self.regions.file_ownership_mask(i)
                self.global_oct_handler.apply_domain(i+1, fmask, 
                                                     self.regions.index_order1)
            tot = sum(self.global_oct_handler.recursively_count().values())
            only_on_root(mylog.info, "Allocating for %0.3e particles "
                                    "(index particle type '%s')",
                        self.total_particles, index_ptype)

    def _index_filename(self,o1,o2):
        import shutil
        # Uses parameter file
        fname_old = os.path.join(self.dataset.fullpath, 
                                 "index{}_{}.ewah".format(o1,o2))
        fname_new = self.index_filename + ".index{}_{}.ewah".format(o1,o2)
        if os.path.isfile(fname_old):
            shutil.move(fname_old,fname_new)
        return fname_new

    def _initialize_indices(self):
        # TODO: THIS IS NOT CURRENTLY USED AND MUST BE MERGED WITH ABOVE
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
        morton = np.empty(self.total_particles, dtype="uint64")
        ind = 0
        for data_file in self.data_files:
            npart = sum(data_file.total_particles.values())
            # TODO: Make this take index_ptype
            morton[ind:ind + npart] = \
                self.io._initialize_index(data_file, self.regions)
            ind += npart
        morton.sort()
        # Now we add them all at once.
        self.oct_handler.add(morton)

    def _initialize_index(self, fname=None, noref=False,
                          order1=None, order2=None, dont_cache=False):
        ds = self.dataset
        only_on_root(mylog.info, "Allocating for %0.3e particles",
                     self.total_particles, global_rootonly = True)
        # No more than 256^3 in the region finder.
        self.regions = ParticleBitmap(
                ds.domain_left_edge, ds.domain_right_edge,
                len(self.data_files), 
                index_order1=order1, index_order2=order2)
        N = 1<<(self.regions.index_order1 + self.ds.over_refine_factor)
        self.ds.domain_dimensions[:] = N
        # Load indices from file if provided
        if fname is None: 
            fname = self._index_filename(self.regions.index_order1,
                                         self.regions.index_order2)
        try:
            rflag = self.regions.load_bitmasks(fname)
            if rflag == 0:
                self._initialize_owners()
                self.regions.save_bitmasks(fname)
            rflag = self.regions.check_bitmasks()
            if rflag == 0:
                raise IOError()
                # dont_cache = True
                # raise IOError()
            #     self._initialize_owners()
            # rflag = self.regions.check_bitmasks()
            # else: # Save pcounts in file?
            #     self._initialize_owners()
        except IOError:
            self.regions.reset_bitmasks()
            self._initialize_coarse_index()
            if not noref:
                self._initialize_refined_index()
                self.regions.set_owners()
            else:
                self._initialize_owners()
            if not dont_cache:
                self.regions.save_bitmasks(fname)
            rflag = self.regions.check_bitmasks()
        # These are now invalid, but I don't know what to replace them with:
        #self.max_level = self.oct_handler.max_level
        #self.dataset.max_level = self.max_level

    def _initialize_coarse_index(self):
        pb = get_pbar("Initializing coarse index ", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            for pos in self.io._yield_coordinates(data_file):
                self.regions._coarse_index_data_file(pos, data_file.file_id)
            self.regions._set_coarse_index_data_file(data_file.file_id)
        pb.finish()
        self.regions.find_collisions_coarse()

    def _initialize_refined_index(self):
        mask = self.regions.masks.sum(axis=1).astype('uint8')
        max_npart = max(sum(d.total_particles.values())
                        for d in self.data_files)
        sub_mi1 = np.zeros(max_npart, "uint64")
        sub_mi2 = np.zeros(max_npart, "uint64")
        pcount = np.zeros(1 << (self.regions.index_order1*3), 'uint32')
        pb = get_pbar("Initializing refined index", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            pcount[:] = 0
            nsub_mi = 0
            for pos in self.io._yield_coordinates(data_file):
                nsub_mi = self.regions._refined_index_data_file(pos, 
                    mask, pcount, sub_mi1, sub_mi2,
                    data_file.file_id, nsub_mi)
            self.regions._set_refined_index_data_file(
                mask, sub_mi1, sub_mi2,
                data_file.file_id, nsub_mi)
        pb.finish()
        self.regions.find_collisions_refined()

    def _initialize_owners(self):
        pcount = np.zeros(1 << (self.regions.index_order1*3), 'uint32')
        pb = get_pbar("Initializing owners", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            pcount[:] = 0
            for pos in self.io._yield_coordinates(data_file):
                self.regions._owners_data_file(pos, pcount, data_file.file_id)
        pb.finish()
        self.regions.set_owners()
            
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
            self._initialize_index()
        # Must check that chunk_info contains the right number of ghost zones
        if getattr(dobj, "_chunk_info", None) is None:
            if isinstance(dobj, (ParticleContainer, ParticleOctreeSubset)):
                dobj._chunk_info = [dobj]
            else:
                dfi, file_masks, addfi = self.regions.identify_file_masks(dobj.selector)
                nfiles = len(file_masks)
                dobj._chunk_info = [None for _ in range(nfiles)]
                for i in range(nfiles):
                    domain_id = i+1
                    dobj._chunk_info[i] = ParticleContainer(
                        dobj, [self.data_files[dfi[i]]],
                        overlap_files = [self.data_files[k] for k in addfi[i]],
                        selector_mask = file_masks[i], domain_id = domain_id)
                    # dobj._chunk_info[i] = ParticleOctreeSubset(
                    #     dobj, [self.data_files[dfi[i]]], 
                    #     overlap_files = [self.data_files[k] for k in addfi[i]],
                    #     selector_mask = file_masks[i], domain_id = domain_id,
                    #     over_refine_factor = self.ds.over_refine_factor)
                # NOTE: One fun thing about the way IO works is that it
                # consolidates things quite nicely.  So we should feel free to
                # create as many objects as part of the chunk as we want, since
                # it'll take the set() of them.  So if we break stuff up like this
                # here, we end up in a situation where we have the ability to break
                # things down further later on for buffer zones and the like.
        dobj._current_chunk, = self._chunk_all(dobj)

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None,
                       ghost_particles = False):
        if ngz == 0 and ghost_particles:
            ngz = 1
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            with og._expand_data_files():
                if ghost_particles: # change to ngz > 0?
                    g = og.retrieve_ghost_zones(ngz)
                else:
                    g = og
                with g._as_spatial():
                    yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache = True, local_only = False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for container in oobjs:
            yield YTDataChunk(dobj, "io", [container], None, cache = cache)
