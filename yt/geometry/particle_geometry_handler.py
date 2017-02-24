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

from yt.funcs import \
    get_pbar, \
    only_on_root
from yt.utilities.logger import ytLogger as mylog
from yt.geometry.geometry_handler import \
    Index, \
    YTDataChunk
from yt.geometry.particle_oct_container import ParticleBitmap
from yt.data_objects.particle_container import ParticleContainer

class ParticleIndex(Index):
    """The Index subclass for particle datasets"""
    _global_mesh = False

    def __init__(self, ds, dataset_type):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self._index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self._index_filename)
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
        return self.ds.arr(0, 'code_length')

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
                try:
                    df = cls(self.dataset, self.io, template % {'num':i}, fi,
                            (start, end))
                except:
                    df = cls(self.dataset, self.io, template % {'num':i}, fi)
                    self.data_files.append(df)
                    break
                if max(df.total_particles.values()) == 0:
                    break
                fi += 1
                self.data_files.append(df)
                start = end
                end += CHUNKSIZE
        self.total_particles = sum(
                sum(d.total_particles.values()) for d in self.data_files)
        # Get index & populate octree
        self._initialize_index()

    def _initialize_index(self):
        ds = self.dataset
        only_on_root(mylog.info, "Allocating for %0.3e particles",
                     self.total_particles, global_rootonly = True)

        # use a trivial morton index for datasets containing a single data file
        if len(self.data_files) == 1:
            order1 = 1
            order2 = 1
        else:
            order1 = ds.index_order[0]
            order2 = ds.index_order[1]

        if order1 == 1 and order2 == 1:
            dont_cache = True
        else:
            dont_cache = False

        self.regions = ParticleBitmap(
            ds.domain_left_edge, ds.domain_right_edge,
            ds.periodicity,
            len(self.data_files), 
            index_order1=order1,
            index_order2=order2)

        # Load indices from file if provided
        if getattr(ds, 'index_filename', None) is None:
            fname = self._index_filename + ".index{}_{}.ewah".format(
                self.regions.index_order1, self.regions.index_order2)
        else:
            fname = self.index_filename

        try:
            rflag = self.regions.load_bitmasks(fname)
            rflag = self.regions.check_bitmasks()
            if rflag == 0:
                raise IOError
        except IOError:
            self.regions.reset_bitmasks()
            self._initialize_coarse_index()
            self._initialize_refined_index()
            if not dont_cache:
                self.regions.save_bitmasks(fname)
            rflag = self.regions.check_bitmasks()

    def _initialize_coarse_index(self):
        pb = get_pbar("Initializing coarse index ", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            for ptype, pos in self.io._yield_coordinates(data_file):
                if hasattr(self.ds, '_sph_ptype') and ptype == self.ds._sph_ptype:
                    hsml = self.io._yield_smoothing_length(data_file, pos.dtype)
                else:
                    hsml = None
                self.regions._coarse_index_data_file(
                    pos, hsml, data_file.file_id)
            self.regions._set_coarse_index_data_file(data_file.file_id)
        pb.finish()
        self.regions.find_collisions_coarse()

    def _initialize_refined_index(self):
        mask = self.regions.masks.sum(axis=1).astype('uint8')
        max_npart = max(sum(d.total_particles.values())
                        for d in self.data_files) * 28
        sub_mi1 = np.zeros(max_npart, "uint64")
        sub_mi2 = np.zeros(max_npart, "uint64")
        pb = get_pbar("Initializing refined index", len(self.data_files))
        for i, data_file in enumerate(self.data_files):
            pb.update(i)
            nsub_mi = 0
            for ptype, pos in self.io._yield_coordinates(data_file):
                if hasattr(self.ds, '_sph_ptype') and ptype == self.ds._sph_ptype:
                    hsml = self.io._yield_smoothing_length(data_file, pos.dtype)
                else:
                    hsml = None
                nsub_mi = self.regions._refined_index_data_file(
                    pos, hsml, mask, sub_mi1, sub_mi2,
                    data_file.file_id, nsub_mi)
            self.regions._set_refined_index_data_file(
                sub_mi1, sub_mi2,
                data_file.file_id, nsub_mi)
        pb.finish()
        self.regions.find_collisions_refined()

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
            if isinstance(dobj, ParticleContainer):
                dobj._chunk_info = [dobj]
            else:
                # TODO: only return files
                dfi, file_masks, addfi = self.regions.identify_file_masks(
                    dobj.selector)
                nfiles = len(file_masks)
                dobj._chunk_info = [None for _ in range(nfiles)]
                for i in range(nfiles):
                    domain_id = i+1
                    dobj._chunk_info[i] = ParticleContainer(
                        dobj, [self.data_files[dfi[i]]],
                        domain_id = domain_id)
                # NOTE: One fun thing about the way IO works is that it
                # consolidates things quite nicely.  So we should feel free to
                # create as many objects as part of the chunk as we want, since
                # it'll take the set() of them.  So if we break stuff up like
                # this here, we end up in a situation where we have the ability
                # to break things down further later on for buffer zones and the
                # like.
        dobj._current_chunk, = self._chunk_all(dobj)

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort = None, preload_fields = None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            with og._expand_data_files():
                if ngz > 0:
                    g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
                else:
                    g = og
                yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache = True, local_only = False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for container in oobjs:
            yield YTDataChunk(dobj, "io", [container], None, cache = cache)
