import os
import re
import string
import sys
import time
import weakref
from collections import defaultdict

import numpy as np
from more_itertools import always_iterable

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.fields.field_info_container import NullFunc
from yt.frontends.enzo.misc import cosmology_get_units
from yt.funcs import get_pbar, iter_fields, setdefaultattr
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py, _libconf as libconf

from .fields import EnzoFieldInfo

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property


class EnzoGrid(AMRGridPatch):
    """
    Class representing a single Enzo Grid instance.
    """

    def __init__(self, id, index):
        """
        Returns an instance of EnzoGrid with *id*, associated with
        *filename* and *index*.
        """
        # All of the field parameters will be passed to us as needed.
        AMRGridPatch.__init__(self, id, filename=None, index=index)
        self._children_ids = []
        self._parent_id = -1
        self.Level = -1

    def set_filename(self, filename):
        """
        Intelligently set the filename.
        """
        if filename is None:
            self.filename = filename
            return
        if self.index._strip_path:
            self.filename = os.path.join(
                self.index.directory, os.path.basename(filename)
            )
        elif filename[0] == os.path.sep:
            self.filename = filename
        else:
            self.filename = os.path.join(self.index.directory, filename)
        return

    @property
    def Parent(self):
        if self._parent_id == -1:
            return None
        return self.index.grids[self._parent_id - self._id_offset]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset] for cid in self._children_ids]

    @property
    def NumberOfActiveParticles(self):
        if not hasattr(self.index, "grid_active_particle_count"):
            return {}
        id = self.id - self._id_offset
        nap = {
            ptype: self.index.grid_active_particle_count[ptype][id]
            for ptype in self.index.grid_active_particle_count
        }
        return nap


class EnzoGridInMemory(EnzoGrid):
    __slots__ = ["proc_num"]

    def set_filename(self, filename):
        pass


class EnzoGridGZ(EnzoGrid):

    __slots__ = ()

    def retrieve_ghost_zones(self, n_zones, fields, all_levels=False, smoothed=False):
        NGZ = self.ds.parameters.get("NumberOfGhostZones", 3)
        if n_zones > NGZ:
            return EnzoGrid.retrieve_ghost_zones(
                self, n_zones, fields, all_levels, smoothed
            )

        # ----- Below is mostly the original code, except we remove the field
        # ----- access section
        # We will attempt this by creating a datacube that is exactly bigger
        # than the grid by nZones*dx in each direction
        nl = self.get_global_startindex() - n_zones
        new_left_edge = nl * self.dds + self.ds.domain_left_edge
        # Something different needs to be done for the root grid, though
        level = self.Level
        kwargs = {
            "dims": self.ActiveDimensions + 2 * n_zones,
            "num_ghost_zones": n_zones,
            "use_pbar": False,
        }
        # This should update the arguments to set the field parameters to be
        # those of this grid.
        kwargs.update(self.field_parameters)
        if smoothed:
            # cube = self.index.smoothed_covering_grid(
            #    level, new_left_edge, new_right_edge, **kwargs)
            cube = self.index.smoothed_covering_grid(level, new_left_edge, **kwargs)
        else:
            cube = self.index.covering_grid(level, new_left_edge, **kwargs)
        # ----- This is EnzoGrid.get_data, duplicated here mostly for
        # ----  efficiency's sake.
        start_zone = NGZ - n_zones
        if start_zone == 0:
            end_zone = None
        else:
            end_zone = -(NGZ - n_zones)
        sl = tuple(slice(start_zone, end_zone) for i in range(3))
        if fields is None:
            return cube
        for field in iter_fields(fields):
            if field in self.field_list:
                conv_factor = 1.0
                if field in self.ds.field_info:
                    conv_factor = self.ds.field_info[field]._convert_function(self)
                if self.ds.field_info[field].sampling_type == "particle":
                    continue
                temp = self.index.io._read_raw_data_set(self, field)
                temp = temp.swapaxes(0, 2)
                cube.field_data[field] = np.multiply(temp, conv_factor, temp)[sl]
        return cube


class EnzoHierarchy(GridIndex):

    _strip_path = False
    grid = EnzoGrid
    _preload_implemented = True

    def __init__(self, ds, dataset_type):

        self.dataset_type = dataset_type
        self.index_filename = os.path.abspath(f"{ds.parameter_filename}.hierarchy")
        if os.path.getsize(self.index_filename) == 0:
            raise OSError(-1, "File empty", self.index_filename)
        self.directory = os.path.dirname(self.index_filename)

        # For some reason, r8 seems to want Float64
        if "CompilerPrecision" in ds and ds["CompilerPrecision"] == "r4":
            self.float_type = "float32"
        else:
            self.float_type = "float64"

        GridIndex.__init__(self, ds, dataset_type)
        # sync it back
        self.dataset.dataset_type = self.dataset_type

    def _count_grids(self):
        self.num_grids = None
        test_grid = test_grid_id = None
        self.num_stars = 0
        for line in rlines(open(self.index_filename, "rb")):
            if (
                line.startswith("BaryonFileName")
                or line.startswith("ParticleFileName")
                or line.startswith("FileName ")
            ):
                test_grid = line.split("=")[-1].strip().rstrip()
            if line.startswith("NumberOfStarParticles"):
                self.num_stars = int(line.split("=")[-1])
            if line.startswith("Grid "):
                if self.num_grids is None:
                    self.num_grids = int(line.split("=")[-1])
                test_grid_id = int(line.split("=")[-1])
                if test_grid is not None:
                    break
        self._guess_dataset_type(self.ds.dimensionality, test_grid, test_grid_id)

    def _guess_dataset_type(self, rank, test_grid, test_grid_id):
        if test_grid[0] != os.path.sep:
            test_grid = os.path.join(self.directory, test_grid)
        if not os.path.exists(test_grid):
            test_grid = os.path.join(self.directory, os.path.basename(test_grid))
            mylog.debug("Your data uses the annoying hardcoded path.")
            self._strip_path = True
        if self.dataset_type is not None:
            return
        if rank == 3:
            mylog.debug("Detected packed HDF5")
            if self.parameters.get("WriteGhostZones", 0) == 1:
                self.dataset_type = "enzo_packed_3d_gz"
                self.grid = EnzoGridGZ
            else:
                self.dataset_type = "enzo_packed_3d"
        elif rank == 2:
            mylog.debug("Detect packed 2D")
            self.dataset_type = "enzo_packed_2d"
        elif rank == 1:
            mylog.debug("Detect packed 1D")
            self.dataset_type = "enzo_packed_1d"
        else:
            raise NotImplementedError

    # Sets are sorted, so that won't work!
    def _parse_index(self):
        def _next_token_line(token, f):
            for line in f:
                if line.startswith(token):
                    return line.split()[2:]

        pattern = r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)\s+$"
        patt = re.compile(pattern)
        f = open(self.index_filename)
        self.grids = [self.grid(1, self)]
        self.grids[0].Level = 0
        si, ei, LE, RE, fn, npart = [], [], [], [], [], []
        pbar = get_pbar("Parsing Hierarchy ", self.num_grids)
        version = self.dataset.parameters.get("VersionNumber", None)
        params = self.dataset.parameters
        if version is None and "Internal" in params:
            version = float(params["Internal"]["Provenance"]["VersionNumber"])
        if version >= 3.0:
            active_particles = True
            nap = {
                ap_type: []
                for ap_type in params["Physics"]["ActiveParticles"][
                    "ActiveParticlesEnabled"
                ]
            }
        else:
            if "AppendActiveParticleType" in self.parameters:
                nap = {}
                active_particles = True
                for type in self.parameters.get("AppendActiveParticleType", []):
                    nap[type] = []
            else:
                nap = None
                active_particles = False
        for grid_id in range(self.num_grids):
            pbar.update(grid_id + 1)
            # We will unroll this list
            si.append(_next_token_line("GridStartIndex", f))
            ei.append(_next_token_line("GridEndIndex", f))
            LE.append(_next_token_line("GridLeftEdge", f))
            RE.append(_next_token_line("GridRightEdge", f))
            nb = int(_next_token_line("NumberOfBaryonFields", f)[0])
            fn.append([None])
            if nb > 0:
                fn[-1] = _next_token_line("BaryonFileName", f)
            npart.append(int(_next_token_line("NumberOfParticles", f)[0]))
            # Below we find out what active particles exist in this grid,
            # and add their counts individually.
            if active_particles:
                ptypes = _next_token_line("PresentParticleTypes", f)
                counts = [int(c) for c in _next_token_line("ParticleTypeCounts", f)]
                for ptype in self.parameters.get("AppendActiveParticleType", []):
                    if ptype in ptypes:
                        nap[ptype].append(counts[ptypes.index(ptype)])
                    else:
                        nap[ptype].append(0)
            if nb == 0 and npart[-1] > 0:
                fn[-1] = _next_token_line("ParticleFileName", f)
            for line in f:
                if len(line) < 2:
                    break
                if line.startswith("Pointer:"):
                    vv = patt.findall(line)[0]
                    self.__pointer_handler(vv)
        pbar.finish()
        self._fill_arrays(ei, si, LE, RE, npart, nap)
        temp_grids = np.empty(self.num_grids, dtype="object")
        temp_grids[:] = self.grids
        self.grids = temp_grids
        self.filenames = fn

    def _initialize_grid_arrays(self):
        super()._initialize_grid_arrays()
        if "AppendActiveParticleType" in self.parameters.keys() and len(
            self.parameters["AppendActiveParticleType"]
        ):
            gac = {
                ptype: np.zeros(self.num_grids, dtype="i4")
                for ptype in self.parameters["AppendActiveParticleType"]
            }
            self.grid_active_particle_count = gac

    def _fill_arrays(self, ei, si, LE, RE, npart, nap):
        self.grid_dimensions.flat[:] = ei
        self.grid_dimensions -= np.array(si, dtype="i4")
        self.grid_dimensions += 1
        self.grid_left_edge.flat[:] = LE
        self.grid_right_edge.flat[:] = RE
        self.grid_particle_count.flat[:] = npart
        if nap is not None:
            for ptype in nap:
                self.grid_active_particle_count[ptype].flat[:] = nap[ptype]

    def __pointer_handler(self, m):
        sgi = int(m[2]) - 1
        if sgi == -1:
            return  # if it's 0, then we're done with that lineage
        # Okay, so, we have a pointer.  We make a new grid, with an id of the length+1
        # (recall, Enzo grids are 1-indexed)
        self.grids.append(self.grid(len(self.grids) + 1, self))
        # We'll just go ahead and make a weakref to cache
        second_grid = self.grids[sgi]  # zero-indexed already
        first_grid = self.grids[int(m[0]) - 1]
        if m[1] == "Next":
            first_grid._children_ids.append(second_grid.id)
            second_grid._parent_id = first_grid.id
            second_grid.Level = first_grid.Level + 1
        elif m[1] == "This":
            if first_grid.Parent is not None:
                first_grid.Parent._children_ids.append(second_grid.id)
                second_grid._parent_id = first_grid._parent_id
            second_grid.Level = first_grid.Level
        self.grid_levels[sgi] = second_grid.Level

    def _rebuild_top_grids(self, level=0):
        mylog.info("Rebuilding grids on level %s", level)
        cmask = self.grid_levels.flat == (level + 1)
        cmsum = cmask.sum()
        mask = np.zeros(self.num_grids, dtype="bool")
        for grid in self.select_grids(level):
            mask[:] = 0
            LE = self.grid_left_edge[grid.id - grid._id_offset]
            RE = self.grid_right_edge[grid.id - grid._id_offset]
            grids, grid_i = self.get_box_grids(LE, RE)
            mask[grid_i] = 1
            grid._children_ids = []
            cgrids = self.grids[(mask * cmask).astype("bool")]
            mylog.info("%s: %s / %s", grid, len(cgrids), cmsum)
            for cgrid in cgrids:
                grid._children_ids.append(cgrid.id)
                cgrid._parent_id = grid.id
        mylog.info("Finished rebuilding")

    def _populate_grid_objects(self):
        for g, f in zip(self.grids, self.filenames):
            g._prepare_grid()
            g._setup_dx()
            g.set_filename(f[0])
        del self.filenames  # No longer needed.
        self.max_level = self.grid_levels.max()

    def _detect_active_particle_fields(self):
        ap_list = self.dataset["AppendActiveParticleType"]
        _fields = {ap: [] for ap in ap_list}
        fields = []
        for ptype in self.dataset["AppendActiveParticleType"]:
            select_grids = self.grid_active_particle_count[ptype].flat
            if not np.any(select_grids):
                current_ptypes = self.dataset.particle_types
                new_ptypes = [p for p in current_ptypes if p != ptype]
                self.dataset.particle_types = new_ptypes
                self.dataset.particle_types_raw = new_ptypes
                continue
            if ptype != "DarkMatter":
                gs = self.grids[select_grids > 0]
                g = gs[0]
                handle = h5py.File(g.filename, "r")
                grid_group = handle[f"/Grid{g.id:08d}"]
                for pname in ["Active Particles", "Particles"]:
                    if pname in grid_group:
                        break
                else:
                    raise RuntimeError("Could not find active particle group in data.")
                node = grid_group[pname]
                for ptype in (str(p) for p in node):
                    if ptype not in _fields:
                        continue
                    for field in (str(f) for f in node[ptype]):
                        _fields[ptype].append(field)
                        if node[ptype][field].ndim > 1:
                            self.io._array_fields[field] = (
                                node[ptype][field].shape[1:],
                            )
                    fields += [(ptype, field) for field in _fields.pop(ptype)]
                handle.close()
        return set(fields)

    def _setup_derived_fields(self):
        super()._setup_derived_fields()
        aps = self.dataset.parameters.get("AppendActiveParticleType", [])
        for fname, field in self.ds.field_info.items():
            if not field.sampling_type == "particle":
                continue
            if isinstance(fname, tuple):
                continue
            if field._function is NullFunc:
                continue
            for apt in aps:
                dd = field._copy_def()
                dd.pop("name")
                self.ds.field_info.add_field((apt, fname), sampling_type="cell", **dd)

    def _detect_output_fields(self):
        self.field_list = []
        # Do this only on the root processor to save disk work.
        if self.comm.rank in (0, None):
            mylog.info("Gathering a field list (this may take a moment.)")
            field_list = set()
            random_sample = self._generate_random_grids()
            for grid in random_sample:
                if not hasattr(grid, "filename"):
                    continue
                try:
                    gf = self.io._read_field_names(grid)
                except self.io._read_exception as e:
                    raise OSError("Grid %s is a bit funky?", grid.id) from e
                mylog.debug("Grid %s has: %s", grid.id, gf)
                field_list = field_list.union(gf)
            if "AppendActiveParticleType" in self.dataset.parameters:
                ap_fields = self._detect_active_particle_fields()
                field_list = list(set(field_list).union(ap_fields))
                if not any(f[0] == "io" for f in field_list):
                    if "io" in self.dataset.particle_types_raw:
                        ptypes_raw = list(self.dataset.particle_types_raw)
                        ptypes_raw.remove("io")
                        self.dataset.particle_types_raw = tuple(ptypes_raw)

                    if "io" in self.dataset.particle_types:
                        ptypes = list(self.dataset.particle_types)
                        ptypes.remove("io")
                        self.dataset.particle_types = tuple(ptypes)
            ptypes = self.dataset.particle_types
            ptypes_raw = self.dataset.particle_types_raw
        else:
            field_list = None
            ptypes = None
            ptypes_raw = None
        self.field_list = list(self.comm.mpi_bcast(field_list))
        self.dataset.particle_types = list(self.comm.mpi_bcast(ptypes))
        self.dataset.particle_types_raw = list(self.comm.mpi_bcast(ptypes_raw))

    def _generate_random_grids(self):
        if self.num_grids > 40:
            starter = np.random.randint(0, 20)
            random_sample = np.mgrid[starter : len(self.grids) - 1 : 20j].astype(
                "int32"
            )
            # We also add in a bit to make sure that some of the grids have
            # particles
            gwp = self.grid_particle_count > 0
            if np.any(gwp) and not np.any(gwp[(random_sample,)]):
                # We just add one grid.  This is not terribly efficient.
                first_grid = np.where(gwp)[0][0]
                random_sample.resize((21,))
                random_sample[-1] = first_grid
                mylog.debug("Added additional grid %s", first_grid)
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = np.mgrid[0 : max(len(self.grids), 1)].astype("int32")
        return self.grids[(random_sample,)]

    def _get_particle_type_counts(self):
        try:
            ret = {}
            for ptype in self.grid_active_particle_count:
                ret[ptype] = self.grid_active_particle_count[ptype].sum()
            return ret
        except AttributeError:
            return super()._get_particle_type_counts()

    def find_particles_by_type(self, ptype, max_num=None, additional_fields=None):
        """
        Returns a structure of arrays with all of the particles'
        positions, velocities, masses, types, IDs, and attributes for
        a particle type **ptype** for a maximum of **max_num**
        particles.  If non-default particle fields are used, provide
        them in **additional_fields**.
        """
        # Not sure whether this routine should be in the general HierarchyType.
        if self.grid_particle_count.sum() == 0:
            mylog.info("Data contains no particles.")
            return None
        if additional_fields is None:
            additional_fields = [
                "metallicity_fraction",
                "creation_time",
                "dynamical_time",
            ]
        pfields = [f for f in self.field_list if f.startswith("particle_")]
        nattr = self.dataset["NumberOfParticleAttributes"]
        if nattr > 0:
            pfields += additional_fields[:nattr]
        # Find where the particles reside and count them
        if max_num is None:
            max_num = 1e100
        total = 0
        pstore = []
        for level in range(self.max_level, -1, -1):
            for grid in self.select_grids(level):
                index = np.where(grid["particle_type"] == ptype)[0]
                total += len(index)
                pstore.append(index)
                if total >= max_num:
                    break
            if total >= max_num:
                break
        result = None
        if total > 0:
            result = {}
            for p in pfields:
                result[p] = np.zeros(total, "float64")
            # Now we retrieve data for each field
            ig = count = 0
            for level in range(self.max_level, -1, -1):
                for grid in self.select_grids(level):
                    nidx = len(pstore[ig])
                    if nidx > 0:
                        for p in pfields:
                            result[p][count : count + nidx] = grid[p][pstore[ig]]
                        count += nidx
                    ig += 1
                    if count >= total:
                        break
                if count >= total:
                    break
            # Crop data if retrieved more than max_num
            if count > max_num:
                for p in pfields:
                    result[p] = result[p][0:max_num]
        return result


class EnzoHierarchyInMemory(EnzoHierarchy):

    grid = EnzoGridInMemory

    @cached_property
    def enzo(self):
        import enzo

        return enzo

    def __init__(self, ds, dataset_type=None):
        self.dataset_type = dataset_type
        self.float_type = "float64"
        self.dataset = weakref.proxy(ds)  # for _obtain_enzo
        self.float_type = self.enzo.hierarchy_information["GridLeftEdge"].dtype
        self.directory = os.getcwd()
        GridIndex.__init__(self, ds, dataset_type)

    def _initialize_data_storage(self):
        pass

    def _count_grids(self):
        self.num_grids = self.enzo.hierarchy_information["GridDimensions"].shape[0]

    def _parse_index(self):
        self._copy_index_structure()
        mylog.debug("Copying reverse tree")
        reverse_tree = self.enzo.hierarchy_information["GridParentIDs"].ravel().tolist()
        # Initial setup:
        mylog.debug("Reconstructing parent-child relationships")
        grids = []
        # We enumerate, so it's 0-indexed id and 1-indexed pid
        self.filenames = ["-1"] * self.num_grids
        for id, pid in enumerate(reverse_tree):
            grids.append(self.grid(id + 1, self))
            grids[-1].Level = self.grid_levels[id, 0]
            if pid > 0:
                grids[-1]._parent_id = pid
                grids[pid - 1]._children_ids.append(grids[-1].id)
        self.max_level = self.grid_levels.max()
        mylog.debug("Preparing grids")
        self.grids = np.empty(len(grids), dtype="object")
        for i, grid in enumerate(grids):
            if (i % 1e4) == 0:
                mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid.filename = "Inline_processor_%07i" % (self.grid_procs[i, 0])
            grid._prepare_grid()
            grid._setup_dx()
            grid.proc_num = self.grid_procs[i, 0]
            self.grids[i] = grid
        mylog.debug("Prepared")

    def _initialize_grid_arrays(self):
        EnzoHierarchy._initialize_grid_arrays(self)
        self.grid_procs = np.zeros((self.num_grids, 1), "int32")

    def _copy_index_structure(self):
        # Dimensions are important!
        self.grid_dimensions[:] = self.enzo.hierarchy_information["GridEndIndices"][:]
        self.grid_dimensions -= self.enzo.hierarchy_information["GridStartIndices"][:]
        self.grid_dimensions += 1
        self.grid_left_edge[:] = self.enzo.hierarchy_information["GridLeftEdge"][:]
        self.grid_right_edge[:] = self.enzo.hierarchy_information["GridRightEdge"][:]
        self.grid_levels[:] = self.enzo.hierarchy_information["GridLevels"][:]
        self.grid_procs = self.enzo.hierarchy_information["GridProcs"].copy()
        self.grid_particle_count[:] = self.enzo.hierarchy_information[
            "GridNumberOfParticles"
        ][:]

    def save_data(self, *args, **kwargs):
        pass

    _cached_field_list = None
    _cached_derived_field_list = None

    def _generate_random_grids(self):
        my_rank = self.comm.rank
        my_grids = self.grids[self.grid_procs.ravel() == my_rank]
        if len(my_grids) > 40:
            starter = np.random.randint(0, 20)
            random_sample = np.mgrid[starter : len(my_grids) - 1 : 20j].astype("int32")
            mylog.debug("Checking grids: %s", random_sample.tolist())
        else:
            random_sample = np.mgrid[0 : max(len(my_grids) - 1, 1)].astype("int32")
        return my_grids[(random_sample,)]

    def _chunk_io(self, dobj, cache=True, local_only=False):
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for g in gobjs:
            gfiles[g.filename].append(g)
        for fn in sorted(gfiles):
            if local_only:
                gobjs = [g for g in gfiles[fn] if g.proc_num == self.comm.rank]
                gfiles[fn] = gobjs
            gs = gfiles[fn]
            count = self._count_selection(dobj, gs)
            yield YTDataChunk(dobj, "io", gs, count, cache=cache)


class EnzoHierarchy1D(EnzoHierarchy):
    def _fill_arrays(self, ei, si, LE, RE, npart, nap):
        self.grid_dimensions[:, :1] = ei
        self.grid_dimensions[:, :1] -= np.array(si, dtype="i4")
        self.grid_dimensions += 1
        self.grid_left_edge[:, :1] = LE
        self.grid_right_edge[:, :1] = RE
        self.grid_particle_count.flat[:] = npart
        self.grid_left_edge[:, 1:] = 0.0
        self.grid_right_edge[:, 1:] = 1.0
        self.grid_dimensions[:, 1:] = 1
        if nap is not None:
            raise NotImplementedError


class EnzoHierarchy2D(EnzoHierarchy):
    def _fill_arrays(self, ei, si, LE, RE, npart, nap):
        self.grid_dimensions[:, :2] = ei
        self.grid_dimensions[:, :2] -= np.array(si, dtype="i4")
        self.grid_dimensions += 1
        self.grid_left_edge[:, :2] = LE
        self.grid_right_edge[:, :2] = RE
        self.grid_particle_count.flat[:] = npart
        self.grid_left_edge[:, 2] = 0.0
        self.grid_right_edge[:, 2] = 1.0
        self.grid_dimensions[:, 2] = 1
        if nap is not None:
            raise NotImplementedError


class EnzoDataset(Dataset):
    """
    Enzo-specific output, set at a fixed time.
    """

    _index_class = EnzoHierarchy
    _field_info_class = EnzoFieldInfo

    def __init__(
        self,
        filename,
        dataset_type=None,
        parameter_override=None,
        conversion_override=None,
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        """
        This class is a stripped down class that simply reads and parses
        *filename* without looking at the index.  *dataset_type* gets passed
        to the index to pre-determine the style of data-output.  However,
        it is not strictly necessary.  Optionally you may specify a
        *parameter_override* dictionary that will override anything in the
        parameter file and a *conversion_override* dictionary that consists
        of {fieldname : conversion_to_cgs} that will override the #DataCGS.
        """
        self.fluid_types += ("enzo",)
        if filename.endswith(".hierarchy"):
            filename = filename[:-10]
        if parameter_override is None:
            parameter_override = {}
        self._parameter_override = parameter_override
        if conversion_override is None:
            conversion_override = {}
        self._conversion_override = conversion_override
        self.storage_filename = storage_filename
        Dataset.__init__(
            self,
            filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

    def _setup_1d(self):
        self._index_class = EnzoHierarchy1D
        self.domain_left_edge = np.concatenate([[self.domain_left_edge], [0.0, 0.0]])
        self.domain_right_edge = np.concatenate([[self.domain_right_edge], [1.0, 1.0]])

    def _setup_2d(self):
        self._index_class = EnzoHierarchy2D
        self.domain_left_edge = np.concatenate([self.domain_left_edge, [0.0]])
        self.domain_right_edge = np.concatenate([self.domain_right_edge, [1.0]])

    def get_parameter(self, parameter, type=None):
        """
        Gets a parameter not in the parameterDict.
        """
        if parameter in self.parameters:
            return self.parameters[parameter]
        for line in open(self.parameter_filename):
            if line.find("#") >= 1:  # Keep the commented lines
                line = line[: line.find("#")]
            line = line.strip().rstrip()
            if len(line) < 2:
                continue
            try:
                param, vals = map(string.strip, map(string.rstrip, line.split("=")))
            except ValueError:
                mylog.error("ValueError: '%s'", line)
            if parameter == param:
                if type is None:
                    t = vals.split()
                else:
                    t = map(type, vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units") and not param.startswith("Temperature"):
                    dataType = param[:-5]
                    self.conversion_factors[dataType] = self.parameters[param]
                return self.parameters[parameter]

        return ""

    @cached_property
    def unique_identifier(self) -> str:
        if "CurrentTimeIdentifier" in self.parameters:
            # enzo2
            return str(self.parameters["CurrentTimeIdentifier"])
        elif "MetaDataDatasetUUID" in self.parameters:
            # enzo2
            return str(self.parameters["MetaDataDatasetUUID"])
        elif "Internal" in self.parameters:
            # enzo3
            return str(
                self.parameters["Internal"]["Provenance"]["CurrentTimeIdentidier"]
            )
        else:
            return super().unique_identifier

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        with open(self.parameter_filename) as f:
            line = f.readline().strip()
            f.seek(0)
            if line == "Internal:":
                self._parse_enzo3_parameter_file(f)
            else:
                self._parse_enzo2_parameter_file(f)

    def _parse_enzo3_parameter_file(self, f):
        self.parameters = p = libconf.load(f)
        sim = p["SimulationControl"]
        internal = p["Internal"]
        phys = p["Physics"]
        self.refine_by = sim["AMR"]["RefineBy"]
        self._periodicity = tuple(
            a == 3 for a in sim["Domain"]["LeftFaceBoundaryCondition"]
        )
        self.dimensionality = sim["Domain"]["TopGridRank"]
        self.domain_dimensions = np.array(
            sim["Domain"]["TopGridDimensions"], dtype="int64"
        )
        self.domain_left_edge = np.array(
            sim["Domain"]["DomainLeftEdge"], dtype="float64"
        )
        self.domain_right_edge = np.array(
            sim["Domain"]["DomainRightEdge"], dtype="float64"
        )
        self.gamma = phys["Hydro"]["Gamma"]
        self.current_time = internal["InitialTime"]
        self.cosmological_simulation = phys["Cosmology"]["ComovingCoordinates"]
        if self.cosmological_simulation == 1:
            cosmo = phys["Cosmology"]
            self.current_redshift = internal["CosmologyCurrentRedshift"]
            self.omega_lambda = cosmo["OmegaLambdaNow"]
            self.omega_matter = cosmo["OmegaMatterNow"]
            self.hubble_constant = cosmo["HubbleConstantNow"]
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0
        self.particle_types = ["DarkMatter"] + phys["ActiveParticles"][
            "ActiveParticlesEnabled"
        ]
        self.particle_types = tuple(self.particle_types)
        self.particle_types_raw = self.particle_types
        if self.dimensionality == 1:
            self._setup_1d()
        elif self.dimensionality == 2:
            self._setup_2d()

    def _parse_enzo2_parameter_file(self, f):
        for line in (l.strip() for l in f):
            if (len(line) < 2) or ("=" not in line):
                continue
            param, vals = (i.strip() for i in line.split("=", 1))
            # First we try to decipher what type of value it is.
            vals = vals.split()
            # Special case approaching.
            if "(do" in vals:
                vals = vals[:1]
            if len(vals) == 0:
                pcast = str  # Assume NULL output
            else:
                v = vals[0]
                # Figure out if it's castable to floating point:
                try:
                    float(v)
                except ValueError:
                    pcast = str
                else:
                    if any("." in v or "e+" in v or "e-" in v for v in vals):
                        pcast = float
                    elif v == "inf":
                        pcast = str
                    else:
                        pcast = int
            # Now we figure out what to do with it.
            if len(vals) == 0:
                vals = ""
            elif len(vals) == 1:
                vals = pcast(vals[0])
            else:
                vals = np.array([pcast(i) for i in vals if i != "-99999"])
            if param.startswith("Append"):
                if param not in self.parameters:
                    self.parameters[param] = []
                self.parameters[param].append(vals)
            else:
                self.parameters[param] = vals
        self.refine_by = self.parameters["RefineBy"]
        _periodicity = tuple(
            always_iterable(self.parameters["LeftFaceBoundaryCondition"] == 3)
        )
        self.dimensionality = self.parameters["TopGridRank"]

        if self.dimensionality > 1:
            self.domain_dimensions = self.parameters["TopGridDimensions"]
            if len(self.domain_dimensions) < 3:
                tmp = self.domain_dimensions.tolist()
                tmp.append(1)
                self.domain_dimensions = np.array(tmp)
                _periodicity += (False,)
            self.domain_left_edge = np.array(
                self.parameters["DomainLeftEdge"], "float64"
            ).copy()
            self.domain_right_edge = np.array(
                self.parameters["DomainRightEdge"], "float64"
            ).copy()
        else:
            self.domain_left_edge = np.array(
                self.parameters["DomainLeftEdge"], "float64"
            )
            self.domain_right_edge = np.array(
                self.parameters["DomainRightEdge"], "float64"
            )
            self.domain_dimensions = np.array(
                [self.parameters["TopGridDimensions"], 1, 1]
            )
            _periodicity += (False, False)
        assert len(_periodicity) == 3
        self._periodicity = _periodicity

        self.gamma = self.parameters["Gamma"]
        if self.parameters["ComovingCoordinates"]:
            self.cosmological_simulation = 1
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.omega_radiation = self.parameters.get(
                "CosmologyOmegaRadiationNow", 0.0
            )
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0
        self.particle_types = []
        self.current_time = self.parameters["InitialTime"]
        if (
            self.parameters["NumberOfParticles"] > 0
            and "AppendActiveParticleType" in self.parameters.keys()
        ):
            # If this is the case, then we know we should have a DarkMatter
            # particle type, and we don't need the "io" type.
            self.parameters["AppendActiveParticleType"].append("DarkMatter")
        else:
            # We do not have an "io" type for Enzo particles if the
            # ActiveParticle machinery is on, as we simply will ignore any of
            # the non-DarkMatter particles in that case.  However, for older
            # datasets, we call this particle type "io".
            self.particle_types = ["io"]
        for ptype in self.parameters.get("AppendActiveParticleType", []):
            self.particle_types.append(ptype)
        self.particle_types = tuple(self.particle_types)
        self.particle_types_raw = self.particle_types

        if self.dimensionality == 1:
            self._setup_1d()
        elif self.dimensionality == 2:
            self._setup_2d()

    def _set_code_unit_attributes(self):
        if self.cosmological_simulation:
            k = cosmology_get_units(
                self.hubble_constant,
                self.omega_matter,
                self.parameters["CosmologyComovingBoxSize"],
                self.parameters["CosmologyInitialRedshift"],
                self.current_redshift,
            )
            # Now some CGS values
            box_size = self.parameters["CosmologyComovingBoxSize"]
            setdefaultattr(self, "length_unit", self.quan(box_size, "Mpccm/h"))
            setdefaultattr(
                self,
                "mass_unit",
                self.quan(k["urho"], "g/cm**3") * (self.length_unit.in_cgs()) ** 3,
            )
            setdefaultattr(self, "time_unit", self.quan(k["utim"], "s"))
            setdefaultattr(self, "velocity_unit", self.quan(k["uvel"], "cm/s"))
        else:
            if "LengthUnits" in self.parameters:
                length_unit = self.parameters["LengthUnits"]
                mass_unit = self.parameters["DensityUnits"] * length_unit**3
                time_unit = self.parameters["TimeUnits"]
            elif "SimulationControl" in self.parameters:
                units = self.parameters["SimulationControl"]["Units"]
                length_unit = units["Length"]
                mass_unit = units["Density"] * length_unit**3
                time_unit = units["Time"]
            else:
                mylog.warning("Setting 1.0 in code units to be 1.0 cm")
                mylog.warning("Setting 1.0 in code units to be 1.0 s")
                length_unit = mass_unit = time_unit = 1.0

            setdefaultattr(self, "length_unit", self.quan(length_unit, "cm"))
            setdefaultattr(self, "mass_unit", self.quan(mass_unit, "g"))
            setdefaultattr(self, "time_unit", self.quan(time_unit, "s"))
            setdefaultattr(self, "velocity_unit", self.length_unit / self.time_unit)

        density_unit = self.mass_unit / self.length_unit**3
        magnetic_unit = np.sqrt(4 * np.pi * density_unit) * self.velocity_unit
        magnetic_unit = np.float64(magnetic_unit.in_cgs())
        setdefaultattr(self, "magnetic_unit", self.quan(magnetic_unit, "gauss"))

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        if str(filename).endswith(".hierarchy"):
            return True
        return os.path.exists(f"{filename}.hierarchy")

    @classmethod
    def _guess_candidates(cls, base, directories, files):
        candidates = [
            _
            for _ in files
            if _.endswith(".hierarchy")
            and os.path.exists(os.path.join(base, _.rsplit(".", 1)[0]))
        ]
        # Typically, Enzo won't have nested outputs.
        return candidates, (len(candidates) == 0)


class EnzoDatasetInMemory(EnzoDataset):
    _index_class = EnzoHierarchyInMemory
    _dataset_type = "enzo_inline"

    def __init__(self, parameter_override=None, conversion_override=None):
        self.fluid_types += ("enzo",)
        if parameter_override is None:
            parameter_override = {}
        self._parameter_override = parameter_override
        if conversion_override is None:
            conversion_override = {}
        self._conversion_override = conversion_override

        Dataset.__init__(self, "InMemoryParameterFile", self._dataset_type)

    def _parse_parameter_file(self):
        enzo = self._obtain_enzo()
        self.basename = "cycle%08i" % (enzo.yt_parameter_file["NumberOfPythonCalls"])
        self.parameters["CurrentTimeIdentifier"] = time.time()
        self.parameters.update(enzo.yt_parameter_file)
        self.conversion_factors.update(enzo.conversion_factors)
        for i in self.parameters:
            if isinstance(self.parameters[i], tuple):
                self.parameters[i] = np.array(self.parameters[i])
            if i.endswith("Units") and not i.startswith("Temperature"):
                dataType = i[:-5]
                self.conversion_factors[dataType] = self.parameters[i]
        self.domain_left_edge = self.parameters["DomainLeftEdge"].copy()
        self.domain_right_edge = self.parameters["DomainRightEdge"].copy()
        for i in self.conversion_factors:
            if isinstance(self.conversion_factors[i], tuple):
                self.conversion_factors[i] = np.array(self.conversion_factors[i])
        for p, v in self._parameter_override.items():
            self.parameters[p] = v
        for p, v in self._conversion_override.items():
            self.conversion_factors[p] = v
        self.refine_by = self.parameters["RefineBy"]
        self._periodicity = tuple(
            always_iterable(self.parameters["LeftFaceBoundaryCondition"] == 3)
        )
        self.dimensionality = self.parameters["TopGridRank"]
        self.domain_dimensions = self.parameters["TopGridDimensions"]
        self.current_time = self.parameters["InitialTime"]
        if self.parameters["ComovingCoordinates"]:
            self.cosmological_simulation = 1
            self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
            self.omega_lambda = self.parameters["CosmologyOmegaLambdaNow"]
            self.omega_matter = self.parameters["CosmologyOmegaMatterNow"]
            self.hubble_constant = self.parameters["CosmologyHubbleConstantNow"]
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0

    def _obtain_enzo(self):
        import enzo

        return enzo

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        return False


# These next two functions are taken from
# http://www.reddit.com/r/Python/comments/6hj75/reverse_file_iterator/c03vms4
# Credit goes to "Brian" on Reddit


def rblocks(f, blocksize=4096):
    """Read file as series of blocks from end of file to start.

    The data itself is in normal order, only the order of the blocks is reversed.
    ie. "hello world" -> ["ld","wor", "lo ", "hel"]
    Note that the file must be opened in binary mode.
    """
    if "b" not in f.mode.lower():
        raise Exception("File must be opened using binary mode.")
    size = os.stat(f.name).st_size
    fullblocks, lastblock = divmod(size, blocksize)

    # The first(end of file) block will be short, since this leaves
    # the rest aligned on a blocksize boundary.  This may be more
    # efficient than having the last (first in file) block be short
    f.seek(-lastblock, 2)
    yield f.read(lastblock).decode("ascii")

    for i in range(fullblocks - 1, -1, -1):
        f.seek(i * blocksize)
        yield f.read(blocksize).decode("ascii")


def rlines(f, keepends=False):
    """Iterate through the lines of a file in reverse order.

    If keepends is true, line endings are kept as part of the line.
    """
    buf = ""
    for block in rblocks(f):
        buf = block + buf
        lines = buf.splitlines(keepends)
        # Return all lines except the first (since may be partial)
        if lines:
            lines.reverse()
            buf = lines.pop()  # Last line becomes end of new first line.
            yield from lines
    yield buf  # First line.
