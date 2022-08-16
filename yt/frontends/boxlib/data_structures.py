import glob
import os
import re
import sys
from collections import namedtuple
from stat import ST_CTIME
from typing import Type

import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.fields.field_info_container import FieldInfoContainer
from yt.funcs import mylog, setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt.utilities.io_handler import io_registry
from yt.utilities.lib.misc_utilities import get_box_grids_level
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_root_only

from .fields import (
    BoxlibFieldInfo,
    CastroFieldInfo,
    MaestroFieldInfo,
    NyxFieldInfo,
    WarpXFieldInfo,
)

if sys.version_info >= (3, 8):
    from functools import cached_property
else:
    from yt._maintenance.backports import cached_property
# This is what we use to find scientific notation that might include d's
# instead of e's.
_scinot_finder = re.compile(r"[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?")
# This is the dimensions in the Cell_H file for each level
# It is different for different dimensionalities, so we make a list
_1dregx = r"-?\d+"
_2dregx = r"-?\d+,-?\d+"
_3dregx = r"-?\d+,-?\d+,-?\d+"
_dim_finder = [
    re.compile(rf"\(\(({ndregx})\) \(({ndregx})\) \({ndregx}\)\)$")
    for ndregx in (_1dregx, _2dregx, _3dregx)
]
# This is the line that prefixes each set of data for a FAB in the FAB file
# It is different for different dimensionalities, so we make a list
_endian_regex = r"^FAB\ \(\(\d+,\ \([\d\ ]+\)\),\((\d+),\ \(([\d\ ]+)\)\)\)"
_header_pattern = [
    re.compile(
        rf"""{_endian_regex}            # match `endianness`
        \(
              \(( {ndregx} )\)          # match `start`
            \ \(( {ndregx} )\)          # match `end`
            \ \(( {ndregx} )\)          # match `centering`
        \)
        \ (-?\d+)                       # match `nc`
        $ # end of line
        """,
        re.VERBOSE,
    )
    for ndregx in (_1dregx, _2dregx, _3dregx)
]


class BoxlibGrid(AMRGridPatch):
    _id_offset = 0
    _offset = -1

    def __init__(self, grid_id, offset, filename=None, index=None):
        super().__init__(grid_id, filename, index)
        self._base_offset = offset
        self._parent_id = []
        self._children_ids = []
        self._pdata = {}

    def _prepare_grid(self):
        super()._prepare_grid()
        my_ind = self.id - self._id_offset
        self.start_index = self.index.grid_start_index[my_ind]

    def get_global_startindex(self):
        return self.start_index

    def _setup_dx(self):
        # has already been read in and stored in index
        self.dds = self.index.ds.arr(self.index.level_dds[self.Level, :], "code_length")
        self.field_data["dx"], self.field_data["dy"], self.field_data["dz"] = self.dds

    @property
    def Parent(self):
        if len(self._parent_id) == 0:
            return None
        return [self.index.grids[pid - self._id_offset] for pid in self._parent_id]

    @property
    def Children(self):
        return [self.index.grids[cid - self._id_offset] for cid in self._children_ids]

    def _get_offset(self, f):
        # This will either seek to the _offset or figure out the correct
        # _offset.
        if self._offset == -1:
            f.seek(self._base_offset, os.SEEK_SET)
            f.readline()
            self._offset = f.tell()
        return self._offset

    # We override here because we can have varying refinement levels
    def select_ires(self, dobj):
        mask = self._get_selector_mask(dobj.selector)
        if mask is None:
            return np.empty(0, dtype="int64")
        coords = np.empty(self._last_count, dtype="int64")
        coords[:] = self.Level + self.ds.level_offsets[self.Level]
        return coords

    # Override this as well, since refine_by can vary
    def _fill_child_mask(self, child, mask, tofill, dlevel=1):
        rf = self.ds.ref_factors[self.Level]
        if dlevel != 1:
            raise NotImplementedError
        gi, cgi = self.get_global_startindex(), child.get_global_startindex()
        startIndex = np.maximum(0, cgi // rf - gi)
        endIndex = np.minimum(
            (cgi + child.ActiveDimensions) // rf - gi, self.ActiveDimensions
        )
        endIndex += startIndex == endIndex
        mask[
            startIndex[0] : endIndex[0],
            startIndex[1] : endIndex[1],
            startIndex[2] : endIndex[2],
        ] = tofill


class BoxLibParticleHeader:
    def __init__(self, ds, directory_name, is_checkpoint, extra_field_names=None):

        self.particle_type = directory_name
        header_filename = os.path.join(ds.output_dir, directory_name, "Header")
        with open(header_filename) as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split("_")[-1]
            known_real_types = {"double": np.float64, "single": np.float32}
            try:
                self.real_type = known_real_types[particle_real_type]
            except KeyError:
                mylog.warning(
                    "yt did not recognize particle real type '%s'. Assuming 'double'.",
                    particle_real_type,
                )
                self.real_type = known_real_types["double"]

            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2 + self.dim
            self.num_real_base = self.dim
            self.num_int_extra = 0  # this should be written by Boxlib, but isn't
            self.num_real_extra = int(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            # Boxlib particles can be written in checkpoint or plotfile mode
            # The base integer fields are only there for checkpoints, but some
            # codes use the checkpoint format for plotting
            if not is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype="int64")
            self.data_map = {}
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())
                self.data_map[level_num] = {}

            pfd = namedtuple(
                "ParticleFileDescriptor", ["file_number", "num_particles", "offset"]
            )

            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.data_map[level_num][grid_num] = pfd(*entry)

        self._generate_particle_fields(extra_field_names)

    def _generate_particle_fields(self, extra_field_names):

        # these are the 'base' integer fields
        self.known_int_fields = [
            (self.particle_type, "particle_id"),
            (self.particle_type, "particle_cpu"),
            (self.particle_type, "particle_cell_x"),
            (self.particle_type, "particle_cell_y"),
            (self.particle_type, "particle_cell_z"),
        ]
        self.known_int_fields = self.known_int_fields[0 : self.num_int_base]

        # these are extra integer fields
        extra_int_fields = [
            "particle_int_comp%d" % i for i in range(self.num_int_extra)
        ]
        self.known_int_fields.extend(
            [(self.particle_type, field) for field in extra_int_fields]
        )

        # these are the base real fields
        self.known_real_fields = [
            (self.particle_type, "particle_position_x"),
            (self.particle_type, "particle_position_y"),
            (self.particle_type, "particle_position_z"),
        ]
        self.known_real_fields = self.known_real_fields[0 : self.num_real_base]

        # these are the extras
        if extra_field_names is not None:
            assert len(extra_field_names) == self.num_real_extra
        else:
            extra_field_names = [
                "particle_real_comp%d" % i for i in range(self.num_real_extra)
            ]

        self.known_real_fields.extend(
            [(self.particle_type, field) for field in extra_field_names]
        )

        self.known_fields = self.known_int_fields + self.known_real_fields

        self.particle_int_dtype = np.dtype(
            [(t[1], self.int_type) for t in self.known_int_fields]
        )

        self.particle_real_dtype = np.dtype(
            [(t[1], self.real_type) for t in self.known_real_fields]
        )


class AMReXParticleHeader:
    def __init__(self, ds, directory_name, is_checkpoint, extra_field_names=None):

        self.particle_type = directory_name
        header_filename = os.path.join(ds.output_dir, directory_name, "Header")
        self.real_component_names = []
        self.int_component_names = []
        with open(header_filename) as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split("_")[-1]

            if particle_real_type == "double":
                self.real_type = np.float64
            elif particle_real_type == "single":
                self.real_type = np.float32
            else:
                raise RuntimeError("yt did not recognize particle real type.")
            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2
            self.num_real_base = self.dim
            self.num_real_extra = int(f.readline().strip())
            for _ in range(self.num_real_extra):
                self.real_component_names.append(f.readline().strip())
            self.num_int_extra = int(f.readline().strip())
            for _ in range(self.num_int_extra):
                self.int_component_names.append(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.is_checkpoint = bool(int(f.readline().strip()))
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            if not self.is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype="int64")
            self.data_map = {}
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())
                self.data_map[level_num] = {}

            pfd = namedtuple(
                "ParticleFileDescriptor", ["file_number", "num_particles", "offset"]
            )

            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.data_map[level_num][grid_num] = pfd(*entry)

        self._generate_particle_fields()

    def _generate_particle_fields(self):

        # these are the 'base' integer fields
        self.known_int_fields = [
            (self.particle_type, "particle_id"),
            (self.particle_type, "particle_cpu"),
        ]
        self.known_int_fields = self.known_int_fields[0 : self.num_int_base]

        self.known_int_fields.extend(
            [
                (self.particle_type, "particle_" + field)
                for field in self.int_component_names
            ]
        )

        # these are the base real fields
        self.known_real_fields = [
            (self.particle_type, "particle_position_x"),
            (self.particle_type, "particle_position_y"),
            (self.particle_type, "particle_position_z"),
        ]
        self.known_real_fields = self.known_real_fields[0 : self.num_real_base]

        self.known_real_fields.extend(
            [
                (self.particle_type, "particle_" + field)
                for field in self.real_component_names
            ]
        )

        self.known_fields = self.known_int_fields + self.known_real_fields

        self.particle_int_dtype = np.dtype(
            [(t[1], self.int_type) for t in self.known_int_fields]
        )

        self.particle_real_dtype = np.dtype(
            [(t[1], self.real_type) for t in self.known_real_fields]
        )


class BoxlibHierarchy(GridIndex):

    grid = BoxlibGrid

    def __init__(self, ds, dataset_type="boxlib_native"):
        self.dataset_type = dataset_type
        self.header_filename = os.path.join(ds.output_dir, "Header")
        self.directory = ds.output_dir
        self.particle_headers = {}

        GridIndex.__init__(self, ds, dataset_type)
        self._cache_endianness(self.grids[-1])

    def _parse_index(self):
        """
        read the global header file for an Boxlib plotfile output.
        """
        self.max_level = self.dataset._max_level
        header_file = open(self.header_filename)

        self.dimensionality = self.dataset.dimensionality
        _our_dim_finder = _dim_finder[self.dimensionality - 1]
        DRE = self.dataset.domain_right_edge  # shortcut
        DLE = self.dataset.domain_left_edge  # shortcut

        # We can now skip to the point in the file we want to start parsing.
        header_file.seek(self.dataset._header_mesh_start)

        dx = []
        for i in range(self.max_level + 1):
            dx.append([float(v) for v in next(header_file).split()])
            # account for non-3d data sets
            if self.dimensionality < 2:
                dx[i].append(DRE[1] - DLE[1])
            if self.dimensionality < 3:
                dx[i].append(DRE[2] - DLE[2])
        self.level_dds = np.array(dx, dtype="float64")
        next(header_file)
        if self.ds.geometry == "cartesian":
            default_ybounds = (0.0, 1.0)
            default_zbounds = (0.0, 1.0)
        elif self.ds.geometry == "cylindrical":
            # Now we check for dimensionality issues
            if self.dimensionality != 2:
                raise RuntimeError("yt needs cylindrical to be 2D")
            self.level_dds[:, 2] = 2 * np.pi
            default_zbounds = (0.0, 2 * np.pi)
        elif self.ds.geometry == "spherical":
            # BoxLib only supports 1D spherical, so ensure
            # the other dimensions have the right extent.
            self.level_dds[:, 1] = np.pi
            self.level_dds[:, 2] = 2 * np.pi
            default_ybounds = (0.0, np.pi)
            default_zbounds = (0.0, 2 * np.pi)
        else:
            raise RuntimeError("Unknown BoxLib coordinate system.")
        if int(next(header_file)) != 0:
            raise RuntimeError("INTERNAL ERROR! This should be a zero.")

        # each level is one group with ngrids on it.
        # each grid has self.dimensionality number of lines of 2 reals
        self.grids = []
        grid_counter = 0
        for level in range(self.max_level + 1):
            vals = next(header_file).split()
            lev, ngrids = int(vals[0]), int(vals[1])
            assert lev == level
            nsteps = int(next(header_file))  # NOQA
            for gi in range(ngrids):
                xlo, xhi = (float(v) for v in next(header_file).split())
                if self.dimensionality > 1:
                    ylo, yhi = (float(v) for v in next(header_file).split())
                else:
                    ylo, yhi = default_ybounds
                if self.dimensionality > 2:
                    zlo, zhi = (float(v) for v in next(header_file).split())
                else:
                    zlo, zhi = default_zbounds
                self.grid_left_edge[grid_counter + gi, :] = [xlo, ylo, zlo]
                self.grid_right_edge[grid_counter + gi, :] = [xhi, yhi, zhi]
            # Now we get to the level header filename, which we open and parse.
            fn = os.path.join(self.dataset.output_dir, next(header_file).strip())
            level_header_file = open(fn + "_H")
            level_dir = os.path.dirname(fn)
            # We skip the first two lines, which contain BoxLib header file
            # version and 'how' the data was written
            next(level_header_file)
            next(level_header_file)
            # Now we get the number of components
            ncomp_this_file = int(next(level_header_file))  # NOQA
            # Skip the next line, which contains the number of ghost zones
            next(level_header_file)
            # To decipher this next line, we expect something like:
            # (8 0
            # where the first is the number of FABs in this level.
            ngrids = int(next(level_header_file).split()[0][1:])
            # Now we can iterate over each and get the indices.
            for gi in range(ngrids):
                # components within it
                start, stop = _our_dim_finder.match(next(level_header_file)).groups()
                # fix for non-3d data
                # note we append '0' to both ends b/c of the '+1' in dims below
                start += ",0" * (3 - self.dimensionality)
                stop += ",0" * (3 - self.dimensionality)
                start = np.array(start.split(","), dtype="int64")
                stop = np.array(stop.split(","), dtype="int64")
                dims = stop - start + 1
                self.grid_dimensions[grid_counter + gi, :] = dims
                self.grid_start_index[grid_counter + gi, :] = start
            # Now we read two more lines.  The first of these is a close
            # parenthesis.
            next(level_header_file)
            # The next is again the number of grids
            next(level_header_file)
            # Now we iterate over grids to find their offsets in each file.
            for gi in range(ngrids):
                # Now we get the data file, at which point we're ready to
                # create the grid.
                dummy, filename, offset = next(level_header_file).split()
                filename = os.path.join(level_dir, filename)
                go = self.grid(grid_counter + gi, int(offset), filename, self)
                go.Level = self.grid_levels[grid_counter + gi, :] = level
                self.grids.append(go)
            grid_counter += ngrids
            # already read the filenames above...
        self.float_type = "float64"

    def _cache_endianness(self, test_grid):
        """
        Cache the endianness and bytes perreal of the grids by using a
        test grid and assuming that all grids have the same
        endianness. This is a pretty safe assumption since Boxlib uses
        one file per processor, and if you're running on a cluster
        with different endian processors, then you're on your own!
        """
        # open the test file & grab the header
        with open(os.path.expanduser(test_grid.filename), "rb") as f:
            header = f.readline().decode("ascii", "ignore")

        bpr, endian, start, stop, centering, nc = (
            _header_pattern[self.dimensionality - 1].search(header).groups()
        )
        # Note that previously we were using a different value for BPR than we
        # use now.  Here is an example set of information directly from BoxLib
        """
        * DOUBLE data
        * FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))((0,0) (63,63) (0,0)) 27  # NOQA: E501
        * FLOAT data
        * FAB ((8, (32 8 23 0 1 9 0 127)),(4, (1 2 3 4)))((0,0) (63,63) (0,0)) 27
        """
        if bpr == endian[0]:
            dtype = f"<f{bpr}"
        elif bpr == endian[-1]:
            dtype = f">f{bpr}"
        else:
            raise ValueError(
                "FAB header is neither big nor little endian. "
                "Perhaps the file is corrupt?"
            )

        mylog.debug("FAB header suggests dtype of %s", dtype)
        self._dtype = np.dtype(dtype)

    def _populate_grid_objects(self):
        mylog.debug("Creating grid objects")
        self.grids = np.array(self.grids, dtype="object")
        self._reconstruct_parent_child()
        for i, grid in enumerate(self.grids):
            if (i % 1e4) == 0:
                mylog.debug("Prepared % 7i / % 7i grids", i, self.num_grids)
            grid._prepare_grid()
            grid._setup_dx()
        mylog.debug("Done creating grid objects")

    def _reconstruct_parent_child(self):
        if self.max_level == 0:
            return
        mask = np.empty(len(self.grids), dtype="int32")
        mylog.debug("First pass; identifying child grids")
        for i, grid in enumerate(self.grids):
            get_box_grids_level(
                self.grid_left_edge[i, :],
                self.grid_right_edge[i, :],
                self.grid_levels[i] + 1,
                self.grid_left_edge,
                self.grid_right_edge,
                self.grid_levels,
                mask,
            )
            ids = np.where(mask.astype("bool"))  # where is a tuple
            grid._children_ids = ids[0] + grid._id_offset
        mylog.debug("Second pass; identifying parents")
        for i, grid in enumerate(self.grids):  # Second pass
            for child in grid.Children:
                child._parent_id.append(i + grid._id_offset)

    def _count_grids(self):
        # We can get everything from the Header file, but note that we're
        # duplicating some work done elsewhere.  In a future where we don't
        # pre-allocate grid arrays, this becomes unnecessary.
        header_file = open(self.header_filename)
        header_file.seek(self.dataset._header_mesh_start)
        # Skip over the level dxs, geometry and the zero:
        [next(header_file) for i in range(self.dataset._max_level + 3)]
        # Now we need to be very careful, as we've seeked, and now we iterate.
        # Does this work?  We are going to count the number of places that we
        # have a three-item line.  The three items would be level, number of
        # grids, and then grid time.
        self.num_grids = 0
        for line in header_file:
            if len(line.split()) != 3:
                continue
            self.num_grids += int(line.split()[1])

    def _initialize_grid_arrays(self):
        super()._initialize_grid_arrays()
        self.grid_start_index = np.zeros((self.num_grids, 3), "int64")

    def _initialize_state_variables(self):
        """override to not re-initialize num_grids in AMRHierarchy.__init__"""
        self._parallel_locking = False
        self._data_file = None
        self._data_mode = None

    def _detect_output_fields(self):
        # This is all done in _parse_header_file
        self.field_list = [("boxlib", f) for f in self.dataset._field_list]
        self.field_indexes = {f[1]: i for i, f in enumerate(self.field_list)}
        # There are times when field_list may change.  We copy it here to
        # avoid that possibility.
        self.field_order = [f for f in self.field_list]

    def _setup_data_io(self):
        self.io = io_registry[self.dataset_type](self.dataset)

    def _determine_particle_output_type(self, directory_name):
        header_filename = os.path.join(self.ds.output_dir, directory_name, "Header")
        with open(header_filename) as f:
            version_string = f.readline().strip()
            if version_string.startswith("Version_Two"):
                return AMReXParticleHeader
            else:
                return BoxLibParticleHeader

    def _read_particles(self, directory_name, is_checkpoint, extra_field_names=None):
        pheader = self._determine_particle_output_type(directory_name)
        self.particle_headers[directory_name] = pheader(
            self.ds, directory_name, is_checkpoint, extra_field_names
        )

        num_parts = self.particle_headers[directory_name].num_particles
        if self.ds._particle_type_counts is None:
            self.ds._particle_type_counts = {}
        self.ds._particle_type_counts[directory_name] = num_parts

        base = os.path.join(self.ds.output_dir, directory_name)
        if len(glob.glob(os.path.join(base, "Level_?", "DATA_????"))) > 0:
            base_particle_fn = os.path.join(base, "Level_%d", "DATA_%.4d")
        elif len(glob.glob(os.path.join(base, "Level_?", "DATA_?????"))) > 0:
            base_particle_fn = os.path.join(base, "Level_%d", "DATA_%.5d")
        else:
            return

        gid = 0
        for lev, data in self.particle_headers[directory_name].data_map.items():
            for pdf in data.values():
                pdict = self.grids[gid]._pdata
                pdict[directory_name] = {}
                pdict[directory_name]["particle_filename"] = base_particle_fn % (
                    lev,
                    pdf.file_number,
                )
                pdict[directory_name]["offset"] = pdf.offset
                pdict[directory_name]["NumberOfParticles"] = pdf.num_particles
                self.grid_particle_count[gid] += pdf.num_particles
                self.grids[gid].NumberOfParticles += pdf.num_particles
                gid += 1

        # add particle fields to field_list
        pfield_list = self.particle_headers[directory_name].known_fields
        self.field_list.extend(pfield_list)


class BoxlibDataset(Dataset):
    """
    This class is a stripped down class that simply reads and parses
    *filename*, without looking at the Boxlib index.
    """

    _index_class = BoxlibHierarchy
    _field_info_class: Type[FieldInfoContainer] = BoxlibFieldInfo
    _output_prefix = None
    _default_cparam_filename = "job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        """
        The paramfile is usually called "inputs"
        and there may be a fortran inputs file usually called "probin"
        plotname here will be a directory name
        as per BoxLib, dataset_type will be Native (implemented here), IEEE (not
        yet implemented) or ASCII (not yet implemented.)
        """
        self.fluid_types += ("boxlib",)
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

        cparam_filename = cparam_filename or self.__class__._default_cparam_filename
        self.cparam_filename = self._lookup_cparam_filepath(
            output_dir, cparam_filename=cparam_filename
        )
        self.fparam_filename = self._localize_check(fparam_filename)
        self.storage_filename = storage_filename

        Dataset.__init__(
            self,
            output_dir,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )

        # These are still used in a few places.
        if "HydroMethod" not in self.parameters.keys():
            self.parameters["HydroMethod"] = "boxlib"
        self.parameters["Time"] = 1.0  # default unit is 1...
        self.parameters["EOSType"] = -1  # default
        self.parameters["gamma"] = self.parameters.get("materials.gamma", 1.6667)

    def _localize_check(self, fn):
        if fn is None:
            return None
        # If the file exists, use it.  If not, set it to None.
        root_dir = os.path.dirname(self.output_dir)
        full_fn = os.path.join(root_dir, fn)
        if os.path.exists(full_fn):
            return full_fn
        return None

    @classmethod
    def _is_valid(cls, filename, *args, cparam_filename=None, **kwargs):
        output_dir = filename
        header_filename = os.path.join(output_dir, "Header")
        # boxlib datasets are always directories, and
        # We *know* it's not boxlib if Header doesn't exist.
        if not os.path.exists(header_filename):
            return False

        if cls is BoxlibDataset:
            # Stop checks here for the boxlib base class.
            # Further checks are performed on subclasses.
            return True

        cparam_filename = cparam_filename or cls._default_cparam_filename
        cparam_filepath = cls._lookup_cparam_filepath(output_dir, cparam_filename)

        if cparam_filepath is None:
            return False

        lines = [line.lower() for line in open(cparam_filepath).readlines()]
        return any(cls._subtype_keyword in line for line in lines)

    @classmethod
    def _lookup_cparam_filepath(cls, output_dir, cparam_filename):
        lookup_table = [
            os.path.abspath(os.path.join(p, cparam_filename))
            for p in (output_dir, os.path.dirname(output_dir))
        ]
        found = [os.path.exists(file) for file in lookup_table]

        if not any(found):
            return None

        return lookup_table[found.index(True)]

    @cached_property
    def unique_identifier(self) -> str:
        hfn = os.path.join(self.output_dir, "Header")
        return str(int(os.stat(hfn)[ST_CTIME]))

    def _parse_parameter_file(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        self._periodicity = (False, False, False)
        self._parse_header_file()
        # Let's read the file
        # the 'inputs' file is now optional
        self._parse_cparams()
        self._parse_fparams()

    def _parse_cparams(self):
        if self.cparam_filename is None:
            return
        for line in (line.split("#")[0].strip() for line in open(self.cparam_filename)):
            try:
                param, vals = (s.strip() for s in line.split("="))
            except ValueError:
                continue
            if param == "amr.ref_ratio":
                vals = self.refine_by = int(vals[0])
            elif param == "Prob.lo_bc":
                vals = tuple(p == "1" for p in vals.split())
                assert len(vals) == self.dimensionality
                periodicity = [False, False, False]  # default to non periodic
                periodicity[: self.dimensionality] = vals  # fill in ndim parsed values
                self._periodicity = tuple(periodicity)
            elif param == "castro.use_comoving":
                vals = self.cosmological_simulation = int(vals)
            else:
                try:
                    vals = _guess_pcast(vals)
                except (IndexError, ValueError):
                    # hitting an empty string or a comment
                    vals = None
            self.parameters[param] = vals

        if getattr(self, "cosmological_simulation", 0) == 1:
            self.omega_lambda = self.parameters["comoving_OmL"]
            self.omega_matter = self.parameters["comoving_OmM"]
            self.hubble_constant = self.parameters["comoving_h"]
            a_file = open(os.path.join(self.output_dir, "comoving_a"))
            line = a_file.readline().strip()
            a_file.close()
            self.current_redshift = 1 / float(line) - 1
        else:
            self.current_redshift = 0.0
            self.omega_lambda = 0.0
            self.omega_matter = 0.0
            self.hubble_constant = 0.0
            self.cosmological_simulation = 0

    def _parse_fparams(self):
        """
        Parses the fortran parameter file for Orion. Most of this will
        be useless, but this is where it keeps mu = mass per
        particle/m_hydrogen.
        """
        if self.fparam_filename is None:
            return
        for line in (l for l in open(self.fparam_filename) if "=" in l):
            param, vals = (v.strip() for v in line.split("="))
            # Now, there are a couple different types of parameters.
            # Some will be where you only have floating point values, others
            # will be where things are specified as string literals.
            # Unfortunately, we're also using Fortran values, which will have
            # things like 1.d-2 which is pathologically difficult to parse if
            # your C library doesn't include 'd' in its locale for strtod.
            # So we'll try to determine this.
            vals = vals.split()
            if any(_scinot_finder.match(v) for v in vals):
                vals = [float(v.replace("D", "e").replace("d", "e")) for v in vals]
            if len(vals) == 1:
                vals = vals[0]
            self.parameters[param] = vals

    def _parse_header_file(self):
        """
        We parse the Boxlib header, which we use as our basis.  Anything in the
        inputs file will override this, but the inputs file is not strictly
        necessary for orientation of the data in space.
        """

        # Note: Python uses a read-ahead buffer, so using next(), which would
        # be my preferred solution, won't work here.  We have to explicitly
        # call readline() if we want to end up with an offset at the very end.
        # Fortunately, elsewhere we don't care about the offset, so we're fine
        # everywhere else using iteration exclusively.
        header_file = open(os.path.join(self.output_dir, "Header"))
        self.orion_version = header_file.readline().rstrip()
        n_fields = int(header_file.readline())

        self._field_list = [header_file.readline().strip() for i in range(n_fields)]

        self.dimensionality = int(header_file.readline())
        self.current_time = float(header_file.readline())
        # This is traditionally a index attribute, so we will set it, but
        # in a slightly hidden variable.
        self._max_level = int(header_file.readline())

        for side, init in zip(["left", "right"], [np.zeros, np.ones]):
            domain_edge = init(3, dtype="float64")
            domain_edge[: self.dimensionality] = header_file.readline().split()
            setattr(self, f"domain_{side}_edge", domain_edge)

        ref_factors = np.array(header_file.readline().split(), dtype="int64")
        if ref_factors.size == 0:
            # We use a default of two, as Nyx doesn't always output this value
            ref_factors = [2] * (self._max_level + 1)
        # We can't vary refinement factors based on dimension, or whatever else
        # they are varied on.  In one curious thing, I found that some Castro 3D
        # data has only two refinement factors, which I don't know how to
        # understand.
        self.ref_factors = ref_factors
        if np.unique(ref_factors).size > 1:
            # We want everything to be a multiple of this.
            self.refine_by = min(ref_factors)
            # Check that they're all multiples of the minimum.
            if not all(
                float(rf) / self.refine_by == int(float(rf) / self.refine_by)
                for rf in ref_factors
            ):
                raise RuntimeError
            base_log = np.log2(self.refine_by)
            self.level_offsets = [0]  # level 0 has to have 0 offset
            lo = 0
            for rf in self.ref_factors:
                lo += int(np.log2(rf) / base_log) - 1
                self.level_offsets.append(lo)
        # assert(np.unique(ref_factors).size == 1)
        else:
            self.refine_by = ref_factors[0]
            self.level_offsets = [0 for l in range(self._max_level + 1)]
        # Now we read the global index space, to get
        index_space = header_file.readline()
        # This will be of the form:
        #  ((0,0,0) (255,255,255) (0,0,0)) ((0,0,0) (511,511,511) (0,0,0))
        # So note that if we split it all up based on spaces, we should be
        # fine, as long as we take the first two entries, which correspond to
        # the root level.  I'm not 100% pleased with this solution.
        root_space = index_space.replace("(", "").replace(")", "").split()[:2]
        start = np.array(root_space[0].split(","), dtype="int64")
        stop = np.array(root_space[1].split(","), dtype="int64")
        dd = np.ones(3, dtype="int64")
        dd[: self.dimensionality] = stop - start + 1
        self.domain_offset[: self.dimensionality] = start
        self.domain_dimensions = dd

        # Skip timesteps per level
        header_file.readline()
        self._header_mesh_start = header_file.tell()
        # Skip the cell size information per level - we'll get this later
        for _ in range(self._max_level + 1):
            header_file.readline()
        # Get the geometry
        next_line = header_file.readline()
        if len(next_line.split()) == 1:
            coordinate_type = int(next_line)
        else:
            coordinate_type = 0

        known_types = {0: "cartesian", 1: "cylindrical", 2: "spherical"}
        try:
            self.geometry = known_types[coordinate_type]
        except KeyError as err:
            raise ValueError(f"Unknown BoxLib coord_type `{coordinate_type}`.") from err

        if self.geometry == "cylindrical":
            dre = self.domain_right_edge
            dre[2] = 2.0 * np.pi
            self.domain_right_edge = dre

    def _set_code_unit_attributes(self):
        setdefaultattr(self, "length_unit", self.quan(1.0, "cm"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "g"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "cm/s"))

    @parallel_root_only
    def print_key_parameters(self):
        for a in [
            "current_time",
            "domain_dimensions",
            "domain_left_edge",
            "domain_right_edge",
        ]:
            if not hasattr(self, a):
                mylog.error("Missing %s in parameter file definition!", a)
                continue
            v = getattr(self, a)
            mylog.info("Parameters: %-25s = %s", a, v)

    def relative_refinement(self, l0, l1):
        offset = self.level_offsets[l1] - self.level_offsets[l0]
        return self.refine_by ** (l1 - l0 + offset)


class OrionHierarchy(BoxlibHierarchy):
    def __init__(self, ds, dataset_type="orion_native"):
        BoxlibHierarchy.__init__(self, ds, dataset_type)
        self._read_particles()
        # self.io = IOHandlerOrion

    def _detect_output_fields(self):
        # This is all done in _parse_header_file
        self.field_list = [("boxlib", f) for f in self.dataset._field_list]
        self.field_indexes = {f[1]: i for i, f in enumerate(self.field_list)}
        # There are times when field_list may change.  We copy it here to
        # avoid that possibility.
        self.field_order = [f for f in self.field_list]

        # look for particle fields
        self.particle_filename = None
        for particle_filename in ["StarParticles", "SinkParticles"]:
            fn = os.path.join(self.ds.output_dir, particle_filename)
            if os.path.exists(fn):
                self.particle_filename = fn

        if self.particle_filename is None:
            return

        pfield_list = [("io", c) for c in self.io.particle_field_index.keys()]
        self.field_list.extend(pfield_list)

    def _read_particles(self):
        """
        Reads in particles and assigns them to grids. Will search for
        Star particles, then sink particles if no star particle file
        is found, and finally will simply note that no particles are
        found if neither works. To add a new Orion particle type,
        simply add it to the if/elif/else block.

        """
        self.grid_particle_count = np.zeros(len(self.grids))

        if self.particle_filename is not None:
            self._read_particle_file(self.particle_filename)

    def _read_particle_file(self, fn):
        """actually reads the orion particle data file itself."""
        if not os.path.exists(fn):
            return
        with open(fn) as f:
            lines = f.readlines()
            self.num_stars = int(lines[0].strip()[0])
            for num, line in enumerate(lines[1:]):
                particle_position_x = float(line.split(" ")[1])
                particle_position_y = float(line.split(" ")[2])
                particle_position_z = float(line.split(" ")[3])
                coord = [particle_position_x, particle_position_y, particle_position_z]
                # for each particle, determine which grids contain it
                # copied from object_finding_mixin.py
                mask = np.ones(self.num_grids)
                for i in range(len(coord)):
                    np.choose(
                        np.greater(self.grid_left_edge.d[:, i], coord[i]),
                        (mask, 0),
                        mask,
                    )
                    np.choose(
                        np.greater(self.grid_right_edge.d[:, i], coord[i]),
                        (0, mask),
                        mask,
                    )
                ind = np.where(mask == 1)
                selected_grids = self.grids[ind]
                # in orion, particles always live on the finest level.
                # so, we want to assign the particle to the finest of
                # the grids we just found
                if len(selected_grids) != 0:
                    grid = sorted(selected_grids, key=lambda grid: grid.Level)[-1]
                    ind = np.where(self.grids == grid)[0][0]
                    self.grid_particle_count[ind] += 1
                    self.grids[ind].NumberOfParticles += 1

                    # store the position in the particle file for fast access.
                    try:
                        self.grids[ind]._particle_line_numbers.append(num + 1)
                    except AttributeError:
                        self.grids[ind]._particle_line_numbers = [num + 1]
        return True


class OrionDataset(BoxlibDataset):

    _index_class = OrionHierarchy
    _subtype_keyword = "hyp."
    _default_cparam_filename = "inputs"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename="probin",
        dataset_type="orion_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        BoxlibDataset.__init__(
            self,
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )


class CastroHierarchy(BoxlibHierarchy):
    def __init__(self, ds, dataset_type="castro_native"):
        super().__init__(ds, dataset_type)

        if "particles" in self.ds.parameters:

            # extra beyond the base real fields that all Boxlib
            # particles have, i.e. the xyz positions
            castro_extra_real_fields = [
                "particle_velocity_x",
                "particle_velocity_y",
                "particle_velocity_z",
            ]

            is_checkpoint = True

            self._read_particles(
                "Tracer",
                is_checkpoint,
                castro_extra_real_fields[0 : self.ds.dimensionality],
            )


class CastroDataset(BoxlibDataset):

    _index_class = CastroHierarchy
    _field_info_class = CastroFieldInfo
    _subtype_keyword = "castro"
    _default_cparam_filename = "job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        super().__init__(
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            storage_filename,
            units_override,
            unit_system,
            default_species_fields=default_species_fields,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        jobinfo_filename = os.path.join(self.output_dir, self.cparam_filename)
        line = ""
        with open(jobinfo_filename) as f:
            while not line.startswith(" Inputs File Parameters"):
                # boundary condition info starts with -x:, etc.
                bcs = ["-x:", "+x:", "-y:", "+y:", "-z:", "+z:"]
                if any(b in line for b in bcs):
                    p, v = line.strip().split(":")
                    self.parameters[p] = v.strip()
                if "git describe" in line or "git hash" in line:
                    # Castro release 17.02 and later
                    #    line format: codename git describe:  the-hash
                    # Castro before release 17.02
                    #    line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()
                line = next(f)

            # runtime parameters that we overrode follow "Inputs File
            # Parameters"
            # skip the "====..." line
            line = next(f)
            for line in f:
                if line.strip() == "" or "fortin parameters" in line:
                    continue
                p, v = line.strip().split("=")
                self.parameters[p] = v.strip()

        # hydro method is set by the base class -- override it here
        self.parameters["HydroMethod"] = "Castro"

        # set the periodicity based on the runtime parameters
        # https://amrex-astro.github.io/Castro/docs/inputs.html?highlight=periodicity
        periodicity = [False, False, False]
        for i, axis in enumerate("xyz"):
            try:
                periodicity[i] = self.parameters[f"-{axis}"] == "interior"
            except KeyError:
                break

        self._periodicity = tuple(periodicity)

        if os.path.isdir(os.path.join(self.output_dir, "Tracer")):
            # we have particles
            self.parameters["particles"] = 1
            self.particle_types = ("Tracer",)
            self.particle_types_raw = self.particle_types


class MaestroDataset(BoxlibDataset):

    _field_info_class = MaestroFieldInfo
    _subtype_keyword = "maestro"
    _default_cparam_filename = "job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        super().__init__(
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            storage_filename,
            units_override,
            unit_system,
            default_species_fields=default_species_fields,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        jobinfo_filename = os.path.join(self.output_dir, self.cparam_filename)

        with open(jobinfo_filename) as f:
            for line in f:
                # get the code git hashes
                if "git hash" in line:
                    # line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()

        with open(jobinfo_filename) as f:
            # get the runtime parameters
            for line in f:
                try:
                    p, v = (_.strip() for _ in line[4:].split("=", 1))
                    if len(v) == 0:
                        self.parameters[p] = ""
                    else:
                        self.parameters[p] = _guess_pcast(v)
                except ValueError:
                    # not a parameter line
                    pass

            # hydro method is set by the base class -- override it here
            self.parameters["HydroMethod"] = "Maestro"

        # set the periodicity based on the integer BC runtime parameters
        periodicity = [False, False, False]
        for i, ax in enumerate("xyz"):
            try:
                periodicity[i] = self.parameters[f"bc{ax}_lo"] != -1
            except KeyError:
                pass

        self._periodicity = tuple(periodicity)


class NyxHierarchy(BoxlibHierarchy):
    def __init__(self, ds, dataset_type="nyx_native"):
        super().__init__(ds, dataset_type)

        if "particles" in self.ds.parameters:
            # extra beyond the base real fields that all Boxlib
            # particles have, i.e. the xyz positions
            nyx_extra_real_fields = [
                "particle_mass",
                "particle_velocity_x",
                "particle_velocity_y",
                "particle_velocity_z",
            ]

            is_checkpoint = False

            self._read_particles(
                "DM",
                is_checkpoint,
                nyx_extra_real_fields[0 : self.ds.dimensionality + 1],
            )


class NyxDataset(BoxlibDataset):

    _index_class = NyxHierarchy
    _field_info_class = NyxFieldInfo
    _subtype_keyword = "nyx"
    _default_cparam_filename = "job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        super().__init__(
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            storage_filename,
            units_override,
            unit_system,
            default_species_fields=default_species_fields,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()

        jobinfo_filename = os.path.join(self.output_dir, self.cparam_filename)

        with open(jobinfo_filename) as f:
            for line in f:
                # get the code git hashes
                if "git hash" in line:
                    # line format: codename git hash:  the-hash
                    fields = line.split(":")
                    self.parameters[fields[0]] = fields[1].strip()

                if line.startswith(" Cosmology Information"):
                    self.cosmological_simulation = 1
                    break
            else:
                self.cosmological_simulation = 0

            if self.cosmological_simulation:
                # note that modern Nyx is always cosmological, but there are some old
                # files without these parameters so we want to special-case them
                for line in f:
                    if "Omega_m (comoving)" in line:
                        self.omega_matter = float(line.split(":")[1])
                    elif "Omega_lambda (comoving)" in line:
                        self.omega_lambda = float(line.split(":")[1])
                    elif "h (comoving)" in line:
                        self.hubble_constant = float(line.split(":")[1])

        # Read in the `comoving_a` file and parse the value. We should fix this
        # in the new Nyx output format...
        a_file = open(os.path.join(self.output_dir, "comoving_a"))
        a_string = a_file.readline().strip()
        a_file.close()

        # Set the scale factor and redshift
        self.cosmological_scale_factor = float(a_string)
        self.parameters["CosmologyCurrentRedshift"] = 1 / float(a_string) - 1

        # alias
        self.current_redshift = self.parameters["CosmologyCurrentRedshift"]
        if os.path.isfile(os.path.join(self.output_dir, "DM/Header")):
            # we have particles
            self.parameters["particles"] = 1
            self.particle_types = ("DM",)
            self.particle_types_raw = self.particle_types

    def _set_code_unit_attributes(self):
        setdefaultattr(self, "mass_unit", self.quan(1.0, "Msun"))
        setdefaultattr(self, "time_unit", self.quan(1.0 / 3.08568025e19, "s"))
        setdefaultattr(
            self, "length_unit", self.quan(1.0 / (1 + self.current_redshift), "Mpc")
        )
        setdefaultattr(self, "velocity_unit", self.length_unit / self.time_unit)


def _guess_pcast(vals):
    # Now we guess some things about the parameter and its type
    # Just in case there are multiple; we'll go
    # back afterward to using vals.
    v = vals.split()[0]
    try:
        float(v.upper().replace("D", "E"))
    except Exception:
        pcast = str
        if v in ("F", "T"):
            pcast = bool
    else:
        syms = (".", "D+", "D-", "E+", "E-", "E", "D")
        if any(sym in v.upper() for sym in syms for v in vals.split()):
            pcast = float
        else:
            pcast = int
    if pcast == bool:
        vals = [value == "T" for value in vals.split()]
    else:
        vals = [pcast(value) for value in vals.split()]
    if len(vals) == 1:
        vals = vals[0]
    return vals


def _read_raw_field_names(raw_file):
    header_files = glob.glob(os.path.join(raw_file, "*_H"))
    return [hf.split(os.sep)[-1][:-2] for hf in header_files]


def _string_to_numpy_array(s):
    return np.array([int(v) for v in s[1:-1].split(",")], dtype=np.int64)


def _line_to_numpy_arrays(line):
    lo_corner = _string_to_numpy_array(line[0][1:])
    hi_corner = _string_to_numpy_array(line[1][:])
    node_type = _string_to_numpy_array(line[2][:-1])
    return lo_corner, hi_corner, node_type


def _get_active_dimensions(box):
    return box[1] - box[2] - box[0] + 1


def _read_header(raw_file, field):
    level_files = glob.glob(os.path.join(raw_file, "Level_*"))
    level_files.sort()

    all_boxes = []
    all_file_names = []
    all_offsets = []

    for level_file in level_files:
        header_file = os.path.join(level_file, field + "_H")
        with open(header_file) as f:

            f.readline()  # version
            f.readline()  # how
            f.readline()  # ncomp

            # nghost_line will be parsed below after the number of dimensions
            # is determined when the boxes are read in
            nghost_line = f.readline().strip().split()

            f.readline()  # num boxes

            # read boxes
            boxes = []
            for line in f:
                clean_line = line.strip().split()
                if clean_line == [")"]:
                    break
                lo_corner, hi_corner, node_type = _line_to_numpy_arrays(clean_line)
                boxes.append((lo_corner, hi_corner, node_type))

            try:
                # nghost_line[0] is a single number
                ng = int(nghost_line[0])
                ndims = len(lo_corner)
                nghost = np.array(ndims * [ng])
            except ValueError:
                # nghost_line[0] is (#,#,#)
                nghost_list = nghost_line[0].strip("()").split(",")
                nghost = np.array(nghost_list, dtype="int64")

            # read the file and offset position for the corresponding box
            file_names = []
            offsets = []
            for line in f:
                if line.startswith("FabOnDisk:"):
                    clean_line = line.strip().split()
                    file_names.append(clean_line[1])
                    offsets.append(int(clean_line[2]))

            all_boxes += boxes
            all_file_names += file_names
            all_offsets += offsets

    return nghost, all_boxes, all_file_names, all_offsets


class WarpXHeader:
    def __init__(self, header_fn):
        self.data = {}
        with open(header_fn) as f:
            self.data["Checkpoint_version"] = int(f.readline().strip().split()[-1])

            self.data["num_levels"] = int(f.readline().strip().split()[-1])
            self.data["istep"] = [int(num) for num in f.readline().strip().split()]
            self.data["nsubsteps"] = [int(num) for num in f.readline().strip().split()]

            self.data["t_new"] = [float(num) for num in f.readline().strip().split()]
            self.data["t_old"] = [float(num) for num in f.readline().strip().split()]
            self.data["dt"] = [float(num) for num in f.readline().strip().split()]

            self.data["moving_window_x"] = float(f.readline().strip().split()[-1])

            #  not all datasets will have is_synchronized
            line = f.readline().strip().split()
            if len(line) == 1:
                self.data["is_synchronized"] = bool(line[-1])
                self.data["prob_lo"] = [
                    float(num) for num in f.readline().strip().split()
                ]
            else:
                self.data["is_synchronized"] = True
                self.data["prob_lo"] = [float(num) for num in line]

            self.data["prob_hi"] = [float(num) for num in f.readline().strip().split()]

            for _ in range(self.data["num_levels"]):
                num_boxes = int(f.readline().strip().split()[0][1:])
                for __ in range(num_boxes):
                    f.readline()
                f.readline()

            i = 0
            line = f.readline()
            while line:
                line = line.strip().split()
                if len(line) == 1:
                    line = f.readline()
                    continue
                self.data["species_%d" % i] = [float(val) for val in line]
                i = i + 1
                line = f.readline()


class WarpXHierarchy(BoxlibHierarchy):
    def __init__(self, ds, dataset_type="boxlib_native"):
        super().__init__(ds, dataset_type)

        is_checkpoint = True
        for ptype in self.ds.particle_types:
            self._read_particles(ptype, is_checkpoint)

        # Additional WarpX particle information (used to set up species)
        self.warpx_header = WarpXHeader(os.path.join(self.ds.output_dir, "WarpXHeader"))

        for key, val in self.warpx_header.data.items():
            if key.startswith("species_"):
                i = int(key.split("_")[-1])
                charge_name = "particle%.1d_charge" % i
                mass_name = "particle%.1d_mass" % i
                self.parameters[charge_name] = val[0]
                self.parameters[mass_name] = val[1]

    def _detect_output_fields(self):
        super()._detect_output_fields()

        # now detect the optional, non-cell-centered fields
        self.raw_file = os.path.join(self.ds.output_dir, "raw_fields")
        self.raw_fields = _read_raw_field_names(os.path.join(self.raw_file, "Level_0"))
        self.field_list += [("raw", f) for f in self.raw_fields]
        self.raw_field_map = {}
        self.ds.nodal_flags = {}
        self.raw_field_nghost = {}
        for field_name in self.raw_fields:
            nghost, boxes, file_names, offsets = _read_header(self.raw_file, field_name)
            self.raw_field_map[field_name] = (boxes, file_names, offsets)
            self.raw_field_nghost[field_name] = nghost
            self.ds.nodal_flags[field_name] = np.array(boxes[0][2])


def _skip_line(line):
    if len(line) == 0:
        return True
    if line[0] == "\n":
        return True
    if line[0] == "=":
        return True
    if line[0] == " ":
        return True


class WarpXDataset(BoxlibDataset):

    _index_class = WarpXHierarchy
    _field_info_class = WarpXFieldInfo
    _subtype_keyword = "warpx"
    _default_cparam_filename = "warpx_job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="mks",
    ):

        self.default_fluid_type = "mesh"
        self.default_field = ("mesh", "density")
        self.fluid_types = ("mesh", "index", "raw")

        super().__init__(
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            storage_filename,
            units_override,
            unit_system,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        jobinfo_filename = os.path.join(self.output_dir, self.cparam_filename)
        with open(jobinfo_filename) as f:
            for line in f.readlines():
                if _skip_line(line):
                    continue
                l = line.strip().split(":")
                if len(l) == 2:
                    self.parameters[l[0].strip()] = l[1].strip()
                l = line.strip().split("=")
                if len(l) == 2:
                    self.parameters[l[0].strip()] = l[1].strip()

        # set the periodicity based on the integer BC runtime parameters
        # https://amrex-codes.github.io/amrex/docs_html/InputsProblemDefinition.html
        periodicity = [False, False, False]
        try:
            is_periodic = self.parameters["geometry.is_periodic"].split()
            periodicity[: len(is_periodic)] = [p == "1" for p in is_periodic]
        except KeyError:
            pass
        self._periodicity = tuple(periodicity)

        particle_types = glob.glob(os.path.join(self.output_dir, "*", "Header"))
        particle_types = [cpt.split(os.sep)[-2] for cpt in particle_types]
        if len(particle_types) > 0:
            self.parameters["particles"] = 1
            self.particle_types = tuple(particle_types)
            self.particle_types_raw = self.particle_types
        else:
            self.particle_types = ()
            self.particle_types_raw = ()

    def _set_code_unit_attributes(self):
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        setdefaultattr(self, "velocity_unit", self.quan(1.0, "m/s"))
        setdefaultattr(self, "magnetic_unit", self.quan(1.0, "T"))


class AMReXHierarchy(BoxlibHierarchy):
    def __init__(self, ds, dataset_type="boxlib_native"):
        super().__init__(ds, dataset_type)

        if "particles" in self.ds.parameters:
            is_checkpoint = True
            for ptype in self.ds.particle_types:
                self._read_particles(ptype, is_checkpoint)


class AMReXDataset(BoxlibDataset):

    _index_class = AMReXHierarchy
    _subtype_keyword = "amrex"
    _default_cparam_filename = "job_info"

    def __init__(
        self,
        output_dir,
        cparam_filename=None,
        fparam_filename=None,
        dataset_type="boxlib_native",
        storage_filename=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):

        super().__init__(
            output_dir,
            cparam_filename,
            fparam_filename,
            dataset_type,
            storage_filename,
            units_override,
            unit_system,
            default_species_fields=default_species_fields,
        )

    def _parse_parameter_file(self):
        super()._parse_parameter_file()
        particle_types = glob.glob(os.path.join(self.output_dir, "*", "Header"))
        particle_types = [cpt.split(os.sep)[-2] for cpt in particle_types]
        if len(particle_types) > 0:
            self.parameters["particles"] = 1
            self.particle_types = tuple(particle_types)
            self.particle_types_raw = self.particle_types
