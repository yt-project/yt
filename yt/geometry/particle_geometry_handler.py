import collections
import errno
import os
import struct
import weakref

import numpy as np

from yt.data_objects.index_subobjects.particle_container import ParticleContainer
from yt.funcs import get_pbar, only_on_root
from yt.geometry.geometry_handler import Index, YTDataChunk
from yt.geometry.particle_oct_container import ParticleBitmap
from yt.utilities.lib.ewah_bool_wrap import BoolArrayCollection
from yt.utilities.lib.fnv_hash import fnv_hash
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_objects


class ParticleIndex(Index):
    """The Index subclass for particle datasets"""

    def __init__(self, ds, dataset_type):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.float_type = np.float64
        super().__init__(ds, dataset_type)
        self._initialize_index()

    def _setup_geometry(self):
        self.regions = None

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return self.ds.arr(0, "code_length")

    def _get_particle_type_counts(self):
        result = collections.defaultdict(lambda: 0)
        for df in self.data_files:
            for k in df.total_particles.keys():
                result[k] += df.total_particles[k]
        return dict(result)

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    @property
    def chunksize(self):
        # This can be overridden in subclasses
        return 64**3

    _data_files = None

    @property
    def data_files(self):
        if self._data_files is not None:
            return self._data_files

        self._setup_filenames()
        return self._data_files

    @data_files.setter
    def data_files(self, value):
        self._data_files = value

    _total_particles = None

    @property
    def total_particles(self):
        if self._total_particles is not None:
            return self._total_particles

        self._total_particles = sum(
            sum(d.total_particles.values()) for d in self.data_files
        )
        return self._total_particles

    def _setup_filenames(self):
        template = self.dataset.filename_template
        ndoms = self.dataset.file_count
        cls = self.dataset._file_class
        self.data_files = []
        fi = 0
        for i in range(int(ndoms)):
            start = 0
            if self.chunksize > 0:
                end = start + self.chunksize
            else:
                end = None
            while True:
                try:
                    _filename = template % {"num": i}
                    df = cls(self.dataset, self.io, _filename, fi, (start, end))
                except FileNotFoundError:
                    mylog.warning(
                        "Failed to load '%s' (missing file or directory)", _filename
                    )
                    break
                if max(df.total_particles.values()) == 0:
                    break
                fi += 1
                self.data_files.append(df)
                if end is None:
                    break
                start = end
                end += self.chunksize

    def _initialize_index(self):
        ds = self.dataset
        only_on_root(
            mylog.info,
            "Allocating for %0.3e particles",
            self.total_particles,
            global_rootonly=True,
        )

        # if we have not yet set domain_left_edge and domain_right_edge then do
        # an I/O pass over the particle coordinates to determine a bounding box
        if self.ds.domain_left_edge is None:
            min_ppos = np.empty(3, dtype="float64")
            min_ppos[:] = np.nan
            max_ppos = np.empty(3, dtype="float64")
            max_ppos[:] = np.nan
            only_on_root(
                mylog.info,
                "Bounding box cannot be inferred from metadata, reading "
                "particle positions to infer bounding box",
            )
            for df in self.data_files:
                for _, ppos in self.io._yield_coordinates(df):
                    min_ppos = np.nanmin(np.vstack([min_ppos, ppos]), axis=0)
                    max_ppos = np.nanmax(np.vstack([max_ppos, ppos]), axis=0)
            only_on_root(
                mylog.info,
                "Load this dataset with bounding_box=[%s, %s] to avoid I/O "
                "overhead from inferring bounding_box." % (min_ppos, max_ppos),
            )
            ds.domain_left_edge = ds.arr(1.05 * min_ppos, "code_length")
            ds.domain_right_edge = ds.arr(1.05 * max_ppos, "code_length")
            ds.domain_width = ds.domain_right_edge - ds.domain_left_edge

        # use a trivial morton index for datasets containing a single chunk
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

        # If we have applied a bounding box then we can't cache the
        # ParticleBitmap because it is domain dependent
        if getattr(ds, "_domain_override", False):
            dont_cache = True

        if not hasattr(self.ds, "_file_hash"):
            self.ds._file_hash = self._generate_hash()

        self.regions = ParticleBitmap(
            ds.domain_left_edge,
            ds.domain_right_edge,
            ds.periodicity,
            self.ds._file_hash,
            len(self.data_files),
            index_order1=order1,
            index_order2=order2,
        )

        # Load Morton index from file if provided
        def _current_fname():
            if getattr(ds, "index_filename", None) is None:
                fname = ds.parameter_filename + ".index{}_{}.ewah".format(
                    self.regions.index_order1, self.regions.index_order2
                )
            else:
                fname = ds.index_filename
            return fname

        fname = _current_fname()

        dont_load = dont_cache and not hasattr(ds, "index_filename")
        try:
            if dont_load:
                raise OSError
            rflag = self.regions.load_bitmasks(fname)
            rflag = self.regions.check_bitmasks()
            self._initialize_frontend_specific()
            if rflag == 0:
                raise OSError
        except (OSError, struct.error):
            self.regions.reset_bitmasks()
            self._initialize_coarse_index()
            self._initialize_refined_index()
            # We now update fname since index_order2 may have changed
            fname = _current_fname()
            wdir = os.path.dirname(fname)
            if not dont_cache and os.access(wdir, os.W_OK):
                # Sometimes os mis-reports whether a directory is writable,
                # So pass if writing the bitmask file fails.
                try:
                    self.regions.save_bitmasks(fname)
                except OSError:
                    pass
            rflag = self.regions.check_bitmasks()

    def _initialize_coarse_index(self):
        max_hsml = 0.0
        pb = get_pbar("Initializing coarse index ", len(self.data_files))
        for i, data_file in parallel_objects(enumerate(self.data_files)):
            pb.update(i + 1)
            for ptype, pos in self.io._yield_coordinates(data_file):
                ds = self.ds
                if hasattr(ds, "_sph_ptypes") and ptype == ds._sph_ptypes[0]:
                    hsml = self.io._get_smoothing_length(
                        data_file, pos.dtype, pos.shape
                    )
                    if hsml is not None and hsml.size > 0.0:
                        max_hsml = max(max_hsml, hsml.max())
                else:
                    hsml = None
                self.regions._coarse_index_data_file(pos, hsml, data_file.file_id)
        pb.finish()
        self.regions.masks = self.comm.mpi_allreduce(self.regions.masks, op="sum")
        self.regions.particle_counts = self.comm.mpi_allreduce(
            self.regions.particle_counts, op="sum"
        )
        for data_file in self.data_files:
            self.regions._set_coarse_index_data_file(data_file.file_id)
        self.regions.find_collisions_coarse()
        if max_hsml > 0.0 and len(self.data_files) > 1:
            # By passing this in, we only allow index_order2 to be increased by
            # two at most, never increased.  One place this becomes particularly
            # useful is in the case of an extremely small section of gas
            # particles embedded in a much much larger domain.  The max
            # smoothing length will be quite small, so based on the larger
            # domain, it will correspond to a very very high index order, which
            # is a large amount of memory!  Having multiple indexes, one for
            # each particle type, would fix this.
            new_order2 = self.regions.update_mi2(max_hsml, ds.index_order[1] + 2)
            mylog.info(
                "Updating index_order2 from %s to %s", ds.index_order[1], new_order2
            )
            self.ds.index_order = (self.ds.index_order[0], new_order2)

    def _initialize_refined_index(self):
        mask = self.regions.masks.sum(axis=1).astype("uint8")
        max_npart = max(sum(d.total_particles.values()) for d in self.data_files) * 28
        sub_mi1 = np.zeros(max_npart, "uint64")
        sub_mi2 = np.zeros(max_npart, "uint64")
        pb = get_pbar("Initializing refined index", len(self.data_files))
        mask_threshold = getattr(self, "_index_mask_threshold", 2)
        count_threshold = getattr(self, "_index_count_threshold", 256)
        mylog.debug(
            "Using estimated thresholds of %s and %s for refinement",
            mask_threshold,
            count_threshold,
        )
        total_refined = 0
        total_coarse_refined = (
            (mask >= 2) & (self.regions.particle_counts > count_threshold)
        ).sum()
        mylog.debug(
            "This should produce roughly %s zones, for %s of the domain",
            total_coarse_refined,
            100 * total_coarse_refined / mask.size,
        )
        storage = {}
        for sto, (i, data_file) in parallel_objects(
            enumerate(self.data_files), storage=storage
        ):
            coll = None
            pb.update(i + 1)
            nsub_mi = 0
            for ptype, pos in self.io._yield_coordinates(data_file):
                if pos.size == 0:
                    continue
                if hasattr(self.ds, "_sph_ptypes") and ptype == self.ds._sph_ptypes[0]:
                    hsml = self.io._get_smoothing_length(
                        data_file, pos.dtype, pos.shape
                    )
                else:
                    hsml = None
                nsub_mi, coll = self.regions._refined_index_data_file(
                    coll,
                    pos,
                    hsml,
                    mask,
                    sub_mi1,
                    sub_mi2,
                    data_file.file_id,
                    nsub_mi,
                    count_threshold=count_threshold,
                    mask_threshold=mask_threshold,
                )
                total_refined += nsub_mi
            sto.result_id = i
            if coll is None:
                coll_str = b""
            else:
                coll_str = coll.dumps()
            sto.result = (data_file.file_id, coll_str)
        pb.finish()
        for i in sorted(storage):
            file_id, coll_str = storage[i]
            coll = BoolArrayCollection()
            coll.loads(coll_str)
            self.regions.bitmasks.append(file_id, coll)
        self.regions.find_collisions_refined()

    def _detect_output_fields(self):
        # TODO: Add additional fields
        dsl = []
        units = {}
        pcounts = self._get_particle_type_counts()
        field_cache = {}
        for dom in self.data_files:
            if dom.filename in field_cache:
                fl, _units = field_cache[dom.filename]
            else:
                fl, _units = self.io._identify_fields(dom)
                field_cache[dom.filename] = fl, _units
            units.update(_units)
            dom._calculate_offsets(fl, pcounts)
            for f in fl:
                if f not in dsl:
                    dsl.append(f)
        self.field_list = dsl
        ds = self.dataset
        ds.particle_types = tuple({pt for pt, ds in dsl})
        # This is an attribute that means these particle types *actually*
        # exist.  As in, they are real, in the dataset.
        ds.field_units.update(units)
        ds.particle_types_raw = ds.particle_types

    def _identify_base_chunk(self, dobj):
        # Must check that chunk_info contains the right number of ghost zones
        if getattr(dobj, "_chunk_info", None) is None:
            if isinstance(dobj, ParticleContainer):
                dobj._chunk_info = [dobj]
            else:
                # TODO: only return files
                if getattr(dobj.selector, "is_all_data", False):
                    nfiles = self.regions.nfiles
                    dfi = np.arange(nfiles)
                else:
                    dfi, file_masks, addfi = self.regions.identify_file_masks(
                        dobj.selector
                    )
                    nfiles = len(file_masks)
                dobj._chunk_info = [None for _ in range(nfiles)]

                # The following was moved here from ParticleContainer in order
                # to make the ParticleContainer object pickleable. By having
                # the base_selector as its own argument, we avoid having to
                # rebuild the index on unpickling a ParticleContainer.
                if hasattr(dobj, "base_selector"):
                    base_selector = dobj.base_selector
                    base_region = dobj.base_region
                else:
                    base_region = dobj
                    base_selector = dobj.selector

                for i in range(nfiles):
                    domain_id = i + 1
                    dobj._chunk_info[i] = ParticleContainer(
                        base_region,
                        base_selector,
                        [self.data_files[dfi[i]]],
                        domain_id=domain_id,
                    )
                # NOTE: One fun thing about the way IO works is that it
                # consolidates things quite nicely.  So we should feel free to
                # create as many objects as part of the chunk as we want, since
                # it'll take the set() of them.  So if we break stuff up like
                # this here, we end up in a situation where we have the ability
                # to break things down further later on for buffer zones and the
                # like.
        (dobj._current_chunk,) = self._chunk_all(dobj)

    def _chunk_all(self, dobj):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        yield YTDataChunk(dobj, "all", oobjs, None)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        sobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for og in sobjs:
            with og._expand_data_files():
                if ngz > 0:
                    g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
                else:
                    g = og
                yield YTDataChunk(dobj, "spatial", [g])

    def _chunk_io(self, dobj, cache=True, local_only=False):
        oobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for container in oobjs:
            yield YTDataChunk(dobj, "io", [container], None, cache=cache)

    def _generate_hash(self):
        # Generate an FNV hash by creating a byte array containing the
        # modification time of as well as the first and last 1 MB of data in
        # every output file
        ret = bytearray()
        for pfile in self.data_files:

            # only look at "real" files, not "fake" files generated by the
            # chunking system
            if pfile.start not in (0, None):
                continue
            try:
                mtime = os.path.getmtime(pfile.filename)
            except OSError as e:
                if e.errno == errno.ENOENT:
                    # this is an in-memory file so we return with a dummy
                    # value
                    return -1
                else:
                    raise
            ret.extend(str(mtime).encode("utf-8"))
            size = os.path.getsize(pfile.filename)
            if size > 1e6:
                size = int(1e6)
            with open(pfile.filename, "rb") as fh:
                # read in first and last 1 MB of data
                data = fh.read(size)
                fh.seek(-size, os.SEEK_END)
                data = fh.read(size)
                ret.extend(data)
            return fnv_hash(ret)

    def _initialize_frontend_specific(self):
        """This is for frontend-specific initialization code

        If there are frontend-specific things that need to be set while
        creating the index, this function forces these operations to happen
        in cases where we are reloading the index from a sidecar file.
        """
        pass
