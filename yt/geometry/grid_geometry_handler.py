import abc
import weakref
from collections import defaultdict
from typing import Optional, Tuple

import numpy as np

from yt.arraytypes import blankRecordArray
from yt.config import ytcfg
from yt.fields.derived_field import ValidateSpatial
from yt.fields.field_detector import FieldDetector
from yt.funcs import ensure_numpy_array, iter_fields
from yt.geometry.geometry_handler import ChunkDataCache, Index, YTDataChunk
from yt.utilities.definitions import MAXLEVEL
from yt.utilities.logger import ytLogger as mylog

from .grid_container import GridTree, MatchPointsToGrids


class GridIndex(Index, abc.ABC):
    """The index class for patch and block AMR datasets."""

    float_type = "float64"
    _preload_implemented = False
    _index_properties = (
        "grid_left_edge",
        "grid_right_edge",
        "grid_levels",
        "grid_particle_count",
        "grid_dimensions",
    )

    def _setup_geometry(self):
        mylog.debug("Counting grids.")
        self._count_grids()

        mylog.debug("Initializing grid arrays.")
        self._initialize_grid_arrays()

        mylog.debug("Parsing index.")
        self._parse_index()

        mylog.debug("Constructing grid objects.")
        self._populate_grid_objects()

        mylog.debug("Re-examining index")
        self._initialize_level_stats()

    @abc.abstractmethod
    def _count_grids(self):
        pass

    @abc.abstractmethod
    def _parse_index(self):
        pass

    @abc.abstractmethod
    def _populate_grid_objects(self):
        pass

    def __del__(self):
        del self.grid_dimensions
        del self.grid_left_edge
        del self.grid_right_edge
        del self.grid_levels
        del self.grid_particle_count
        del self.grids

    @property
    def parameters(self):
        return self.dataset.parameters

    def _detect_output_fields_backup(self):
        # grab fields from backup file as well, if present
        return

    def select_grids(self, level):
        """
        Returns an array of grids at *level*.
        """
        return self.grids[self.grid_levels.flat == level]

    def get_levels(self):
        for level in range(self.max_level + 1):
            yield self.select_grids(level)

    def _initialize_grid_arrays(self):
        mylog.debug("Allocating arrays for %s grids", self.num_grids)
        self.grid_dimensions = np.ones((self.num_grids, 3), "int32")
        self.grid_left_edge = self.ds.arr(
            np.zeros((self.num_grids, 3), self.float_type), "code_length"
        )
        self.grid_right_edge = self.ds.arr(
            np.ones((self.num_grids, 3), self.float_type), "code_length"
        )
        self.grid_levels = np.zeros((self.num_grids, 1), "int32")
        self.grid_particle_count = np.zeros((self.num_grids, 1), "int32")

    def clear_all_data(self):
        """
        This routine clears all the data currently being held onto by the grids
        and the data io handler.
        """
        for g in self.grids:
            g.clear_data()
        self.io.queue.clear()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return self.select_grids(self.grid_levels.max())[0].dds[:].min()

    def _get_particle_type_counts(self):
        return {self.ds.particle_types_raw[0]: self.grid_particle_count.sum()}

    def _initialize_level_stats(self):
        # Now some statistics:
        #   0 = number of grids
        #   1 = number of cells
        #   2 = blank
        desc = {"names": ["numgrids", "numcells", "level"], "formats": ["int64"] * 3}
        self.level_stats = blankRecordArray(desc, MAXLEVEL)
        self.level_stats["level"] = [i for i in range(MAXLEVEL)]
        self.level_stats["numgrids"] = [0 for i in range(MAXLEVEL)]
        self.level_stats["numcells"] = [0 for i in range(MAXLEVEL)]
        for level in range(self.max_level + 1):
            self.level_stats[level]["numgrids"] = np.sum(self.grid_levels == level)
            li = self.grid_levels[:, 0] == level
            self.level_stats[level]["numcells"] = (
                self.grid_dimensions[li, :].prod(axis=1).sum()
            )

    @property
    def grid_corners(self):
        return np.array(
            [
                [
                    self.grid_left_edge[:, 0],
                    self.grid_left_edge[:, 1],
                    self.grid_left_edge[:, 2],
                ],
                [
                    self.grid_right_edge[:, 0],
                    self.grid_left_edge[:, 1],
                    self.grid_left_edge[:, 2],
                ],
                [
                    self.grid_right_edge[:, 0],
                    self.grid_right_edge[:, 1],
                    self.grid_left_edge[:, 2],
                ],
                [
                    self.grid_left_edge[:, 0],
                    self.grid_right_edge[:, 1],
                    self.grid_left_edge[:, 2],
                ],
                [
                    self.grid_left_edge[:, 0],
                    self.grid_left_edge[:, 1],
                    self.grid_right_edge[:, 2],
                ],
                [
                    self.grid_right_edge[:, 0],
                    self.grid_left_edge[:, 1],
                    self.grid_right_edge[:, 2],
                ],
                [
                    self.grid_right_edge[:, 0],
                    self.grid_right_edge[:, 1],
                    self.grid_right_edge[:, 2],
                ],
                [
                    self.grid_left_edge[:, 0],
                    self.grid_right_edge[:, 1],
                    self.grid_right_edge[:, 2],
                ],
            ],
            dtype="float64",
        )

    def lock_grids_to_parents(self):
        r"""This function locks grid edges to their parents.

        This is useful in cases where the grid structure may be somewhat
        irregular, or where setting the left and right edges is a lossy
        process.  It is designed to correct situations where left/right edges
        may be set slightly incorrectly, resulting in discontinuities in images
        and the like.
        """
        mylog.info("Locking grids to parents.")
        for i, g in enumerate(self.grids):
            si = g.get_global_startindex()
            g.LeftEdge = self.ds.domain_left_edge + g.dds * si
            g.RightEdge = g.LeftEdge + g.ActiveDimensions * g.dds
            self.grid_left_edge[i, :] = g.LeftEdge
            self.grid_right_edge[i, :] = g.RightEdge

    def print_stats(self):
        """
        Prints out (stdout) relevant information about the simulation
        """
        header = "{:>3}\t{:>6}\t{:>14}\t{:>14}".format(
            "level", "# grids", "# cells", "# cells^3"
        )
        print(header)
        print(f"{len(header.expandtabs()) * '-'}")
        for level in range(MAXLEVEL):
            if (self.level_stats["numgrids"][level]) == 0:
                continue
            print(
                "% 3i\t% 6i\t% 14i\t% 14i"
                % (
                    level,
                    self.level_stats["numgrids"][level],
                    self.level_stats["numcells"][level],
                    np.ceil(self.level_stats["numcells"][level] ** (1.0 / 3)),
                )
            )
            dx = self.select_grids(level)[0].dds[0]
        print("-" * 46)
        print(
            "   \t% 6i\t% 14i"
            % (self.level_stats["numgrids"].sum(), self.level_stats["numcells"].sum())
        )
        print("\n")
        try:
            print(f"z = {self['CosmologyCurrentRedshift']:0.8f}")
        except Exception:
            pass
        print(
            "t = %0.8e = %0.8e s = %0.8e years"
            % (
                self.ds.current_time.in_units("code_time"),
                self.ds.current_time.in_units("s"),
                self.ds.current_time.in_units("yr"),
            )
        )
        print("\nSmallest Cell:")
        for item in ("Mpc", "pc", "AU", "cm"):
            print(f"\tWidth: {dx.in_units(item):0.3e} {item}")

    def _find_field_values_at_points(self, fields, coords):
        r"""Find the value of fields at a set of coordinates.

        Returns the values [field1, field2,...] of the fields at the given
        (x, y, z) points. Returns a numpy array of field values cross coords
        """
        coords = self.ds.arr(ensure_numpy_array(coords), "code_length")
        grids = self._find_points(coords[:, 0], coords[:, 1], coords[:, 2])[0]
        fields = list(iter_fields(fields))
        mark = np.zeros(3, dtype="int64")
        out = []

        # create point -> grid mapping
        grid_index = {}
        for coord_index, grid in enumerate(grids):
            if grid not in grid_index:
                grid_index[grid] = []
            grid_index[grid].append(coord_index)

        out = []
        for field in fields:
            funit = self.ds._get_field_info(field).units
            out.append(self.ds.arr(np.empty(len(coords)), funit))

        for grid in grid_index:
            cellwidth = (grid.RightEdge - grid.LeftEdge) / grid.ActiveDimensions
            for field_index, field in enumerate(fields):
                for coord_index in grid_index[grid]:
                    mark = (coords[coord_index, :] - grid.LeftEdge) / cellwidth
                    mark = np.array(mark, dtype="int64")
                    out[field_index][coord_index] = grid[field][
                        mark[0], mark[1], mark[2]
                    ]
        if len(fields) == 1:
            return out[0]
        return out

    def _find_points(self, x, y, z):
        """
        Returns the (objects, indices) of leaf grids
        containing a number of (x,y,z) points
        """
        x = ensure_numpy_array(x)
        y = ensure_numpy_array(y)
        z = ensure_numpy_array(z)
        if not len(x) == len(y) == len(z):
            raise ValueError("Arrays of indices must be of the same size")

        grid_tree = self._get_grid_tree()
        pts = MatchPointsToGrids(grid_tree, len(x), x, y, z)
        ind = pts.find_points_in_tree()
        return self.grids[ind], ind

    def _get_grid_tree(self):

        left_edge = self.ds.arr(np.zeros((self.num_grids, 3)), "code_length")
        right_edge = self.ds.arr(np.zeros((self.num_grids, 3)), "code_length")
        level = np.zeros((self.num_grids), dtype="int64")
        parent_ind = np.zeros((self.num_grids), dtype="int64")
        num_children = np.zeros((self.num_grids), dtype="int64")
        dimensions = np.zeros((self.num_grids, 3), dtype="int32")

        for i, grid in enumerate(self.grids):

            left_edge[i, :] = grid.LeftEdge
            right_edge[i, :] = grid.RightEdge
            level[i] = grid.Level
            if grid.Parent is None:
                parent_ind[i] = -1
            else:
                parent_ind[i] = grid.Parent.id - grid.Parent._id_offset
            num_children[i] = np.int64(len(grid.Children))
            dimensions[i, :] = grid.ActiveDimensions

        return GridTree(
            self.num_grids,
            left_edge,
            right_edge,
            dimensions,
            parent_ind,
            level,
            num_children,
        )

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _identify_base_chunk(self, dobj):
        fast_index = None
        if dobj._type_name == "grid":
            dobj._chunk_info = np.empty(1, dtype="object")
            dobj._chunk_info[0] = weakref.proxy(dobj)
        elif getattr(dobj, "_grids", None) is None:
            gi = dobj.selector.select_grids(
                self.grid_left_edge, self.grid_right_edge, self.grid_levels
            )
            if any([g.filename is not None for g in self.grids[gi]]):
                _gsort = _grid_sort_mixed
            else:
                _gsort = _grid_sort_id
            grids = sorted(self.grids[gi], key=_gsort)
            dobj._chunk_info = np.empty(len(grids), dtype="object")
            for i, g in enumerate(grids):
                dobj._chunk_info[i] = g
        # These next two lines, when uncommented, turn "on" the fast index.
        # if dobj._type_name != "grid":
        #    fast_index = self._get_grid_tree()
        if getattr(dobj, "size", None) is None:
            dobj.size = self._count_selection(dobj, fast_index=fast_index)
        if getattr(dobj, "shape", None) is None:
            dobj.shape = (dobj.size,)
        dobj._current_chunk = list(
            self._chunk_all(dobj, cache=False, fast_index=fast_index)
        )[0]

    def _count_selection(self, dobj, grids=None, fast_index=None):
        if fast_index is not None:
            return fast_index.count(dobj.selector)
        if grids is None:
            grids = dobj._chunk_info
        count = sum(g.count(dobj.selector) for g in grids)
        return count

    def _chunk_all(self, dobj, cache=True, fast_index=None):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        fast_index = fast_index or getattr(dobj._current_chunk, "_fast_index", None)
        yield YTDataChunk(dobj, "all", gobjs, dobj.size, cache, fast_index=fast_index)

    def _chunk_spatial(self, dobj, ngz, sort=None, preload_fields=None):
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        if sort in ("+level", "level"):
            giter = sorted(gobjs, key=lambda g: g.Level)
        elif sort == "-level":
            giter = sorted(gobjs, key=lambda g: -g.Level)
        elif sort is None:
            giter = gobjs
        if preload_fields is None:
            preload_fields = []
        preload_fields, _ = self._split_fields(preload_fields)
        if self._preload_implemented and len(preload_fields) > 0 and ngz == 0:
            giter = ChunkDataCache(list(giter), preload_fields, self)
        for og in giter:
            if ngz > 0:
                g = og.retrieve_ghost_zones(ngz, [], smoothed=True)
            else:
                g = og
            size = self._count_selection(dobj, [og])
            if size == 0:
                continue
            # We don't want to cache any of the masks or icoords or fcoords for
            # individual grids.
            yield YTDataChunk(dobj, "spatial", [g], size, cache=False)

    _grid_chunksize = 1000

    def _chunk_io(
        self,
        dobj,
        cache=True,
        local_only=False,
        preload_fields=None,
        chunk_sizing="auto",
    ):
        # local_only is only useful for inline datasets and requires
        # implementation by subclasses.
        if preload_fields is None:
            preload_fields = []
        preload_fields, _ = self._split_fields(preload_fields)
        gfiles = defaultdict(list)
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        fast_index = dobj._current_chunk._fast_index
        for g in gobjs:
            # Force to be a string because sometimes g.filename is None.
            gfiles[str(g.filename)].append(g)
        # We can apply a heuristic here to make sure we aren't loading too
        # many grids all at once.
        if chunk_sizing == "auto":
            chunk_ngrids = len(gobjs)
            if chunk_ngrids > 0:
                nproc = int(ytcfg.get("yt", "internals", "global_parallel_size"))
                chunking_factor = np.ceil(
                    self._grid_chunksize * nproc / chunk_ngrids
                ).astype("int")
                size = max(self._grid_chunksize // chunking_factor, 1)
            else:
                size = self._grid_chunksize
        elif chunk_sizing == "config_file":
            size = ytcfg.get("yt", "chunk_size")
        elif chunk_sizing == "just_one":
            size = 1
        elif chunk_sizing == "old":
            size = self._grid_chunksize
        else:
            raise RuntimeError(
                f"{chunk_sizing} is an invalid value for the 'chunk_sizing' argument."
            )
        for fn in sorted(gfiles):
            gs = gfiles[fn]
            for grids in (gs[pos : pos + size] for pos in range(0, len(gs), size)):
                dc = YTDataChunk(
                    dobj,
                    "io",
                    grids,
                    self._count_selection(dobj, grids),
                    cache=cache,
                    fast_index=fast_index,
                )
                # We allow four full chunks to be included.
                with self.io.preload(dc, preload_fields, 4.0 * size):
                    yield dc

    def _icoords_to_fcoords(
        self,
        icoords: np.ndarray,
        ires: np.ndarray,
        axes: Optional[Tuple[int, ...]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Accepts icoords and ires and returns appropriate fcoords and fwidth.
        Mostly useful for cases where we have irregularly spaced or structured
        grids.
        """
        dds = self.ds.domain_width[(axes,)] / (
            self.ds.domain_dimensions[
                axes,
            ]
            * self.ds.refine_by ** ires[:, None]
        )
        pos = (0.5 + icoords) * dds + self.ds.domain_left_edge[
            axes,
        ]
        return pos, dds

    def _add_mesh_sampling_particle_field(self, deposit_field, ftype, ptype):
        units = self.ds.field_info[ftype, deposit_field].units
        take_log = self.ds.field_info[ftype, deposit_field].take_log
        field_name = f"cell_{ftype}_{deposit_field}"

        def _mesh_sampling_particle_field(field, data):
            pos = data[ptype, "particle_position"]
            field_values = data[ftype, deposit_field]

            if isinstance(data, FieldDetector):
                return np.zeros(pos.shape[0])

            i, j, k = np.floor((pos - data.LeftEdge) / data.dds).astype("int64").T

            # Make sure all particles are within the current grid, otherwise return nan
            maxi, maxj, maxk = field_values.shape

            mask = (i < maxi) & (j < maxj) & (k < maxk)
            mask &= (i >= 0) & (j >= 0) & (k >= 0)

            result = np.full(len(pos), np.nan, dtype="float64")
            if result.shape[0] > 0:
                result[mask] = field_values[i[mask], j[mask], k[mask]]

            return data.ds.arr(result, field_values.units)

        self.ds.add_field(
            (ptype, field_name),
            function=_mesh_sampling_particle_field,
            sampling_type="particle",
            units=units,
            take_log=take_log,
            validators=[ValidateSpatial()],
        )


def _grid_sort_id(g):
    return g.id


def _grid_sort_mixed(g):
    if g.filename is None:
        return str(g.id)
    return g.filename
