import re

import numpy as np

from yt.geometry.selection_routines import GridSelector
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog


class IOHandlerChomboHDF5(BaseIOHandler):
    _dataset_type = "chombo_hdf5"
    _offset_string = "data:offsets=0"
    _data_string = "data:datatype=0"
    _offsets = None

    def __init__(self, ds, *args, **kwargs):
        BaseIOHandler.__init__(self, ds, *args, **kwargs)
        self.ds = ds
        self._handle = ds._handle
        self.dim = self._handle["Chombo_global/"].attrs["SpaceDim"]
        self._read_ghost_info()
        if self._offset_string not in self._handle["level_0"]:
            self._calculate_offsets()

    def _calculate_offsets(self):
        def box_size(corners):
            size = 1
            for idim in range(self.dim):
                size *= corners[idim + self.dim] - corners[idim] + 1
            return size

        self._offsets = {}
        num_comp = self._handle.attrs["num_components"]
        level = 0
        while True:
            lname = "level_%i" % level
            if lname not in self._handle:
                break
            boxes = self._handle["level_0"]["boxes"][()]
            box_sizes = np.array([box_size(box) for box in boxes])

            offsets = np.cumsum(box_sizes * num_comp, dtype="int64")
            offsets -= offsets[0]
            self._offsets[level] = offsets
            level += 1

    def _read_ghost_info(self):
        try:
            self.ghost = tuple(
                self._handle["level_0/data_attributes"].attrs["outputGhost"]
            )
            # pad with zeros if the dataset is low-dimensional
            self.ghost += (3 - self.dim) * (0,)
            self.ghost = np.array(self.ghost)
        except KeyError:
            # assume zero ghosts if outputGhosts not present
            self.ghost = np.zeros(self.dim, "int64")

    _field_dict = None

    @property
    def field_dict(self):
        if self._field_dict is not None:
            return self._field_dict
        field_dict = {}
        for key, val in self._handle.attrs.items():
            if key.startswith("component_"):
                comp_number = int(re.match(r"component_(\d+)", key).groups()[0])
                field_dict[val.decode("utf-8")] = comp_number
        self._field_dict = field_dict
        return self._field_dict

    _particle_field_index = None

    @property
    def particle_field_index(self):
        if self._particle_field_index is not None:
            return self._particle_field_index
        field_dict = {}
        for key, val in self._handle.attrs.items():
            if key.startswith("particle_"):
                comp_number = int(
                    re.match(r"particle_component_(\d+)", key).groups()[0]
                )
                field_dict[val.decode("ascii")] = comp_number
        self._particle_field_index = field_dict
        return self._particle_field_index

    def _read_data(self, grid, field):
        lstring = "level_%i" % grid.Level
        lev = self._handle[lstring]
        dims = grid.ActiveDimensions
        shape = dims + 2 * self.ghost
        boxsize = shape.prod()

        if self._offsets is not None:
            grid_offset = self._offsets[grid.Level][grid._level_id]
        else:
            grid_offset = lev[self._offset_string][grid._level_id]
        start = grid_offset + self.field_dict[field] * boxsize
        stop = start + boxsize
        data = lev[self._data_string][start:stop]
        data_no_ghost = data.reshape(shape, order="F")
        ghost_slice = tuple(slice(g, d + g, None) for g, d in zip(self.ghost, dims))
        ghost_slice = ghost_slice[0 : self.dim]
        return data_no_ghost[ghost_slice]

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        chunks = list(chunks)
        fields.sort(key=lambda a: self.field_dict[a[1]])
        if isinstance(selector, GridSelector):
            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError
            grid = chunks[0].objs[0]
            for ftype, fname in fields:
                rv[ftype, fname] = self._read_data(grid, fname)
            return rv
        if size is None:
            size = sum(g.count(selector) for chunk in chunks for g in chunk.objs)
        for field in fields:
            ftype, fname = field
            fsize = size
            rv[field] = np.empty(fsize, dtype="float64")
        ng = sum(len(c.objs) for c in chunks)
        mylog.debug(
            "Reading %s cells of %s fields in %s grids",
            size,
            [f2 for f1, f2 in fields],
            ng,
        )

        ind = 0
        for chunk in chunks:
            for g in chunk.objs:
                nd = 0
                for field in fields:
                    ftype, fname = field
                    data = self._read_data(g, fname)
                    nd = g.select(selector, data, rv[field], ind)  # caches
                ind += nd
        return rv

    def _read_particle_selection(self, chunks, selector, fields):
        rv = {}
        chunks = list(chunks)

        if isinstance(selector, GridSelector):

            if not (len(chunks) == len(chunks[0].objs) == 1):
                raise RuntimeError

            grid = chunks[0].objs[0]

            for ftype, fname in fields:
                rv[ftype, fname] = self._read_particles(grid, fname)

            return rv

        rv = {f: np.array([]) for f in fields}
        for chunk in chunks:
            for grid in chunk.objs:
                for ftype, fname in fields:
                    data = self._read_particles(grid, fname)
                    rv[ftype, fname] = np.concatenate((data, rv[ftype, fname]))
        return rv

    def _read_particles(self, grid, name):

        field_index = self.particle_field_index[name]
        lev = f"level_{grid.Level}"

        particles_per_grid = self._handle[lev]["particles:offsets"][()]
        items_per_particle = len(self._particle_field_index)

        # compute global offset position
        offsets = items_per_particle * np.cumsum(particles_per_grid)
        offsets = np.append(np.array([0]), offsets)
        offsets = np.array(offsets, dtype=np.int64)

        # convert between the global grid id and the id on this level
        grid_levels = np.array([g.Level for g in self.ds.index.grids])
        grid_ids = np.array([g.id for g in self.ds.index.grids])
        grid_level_offset = grid_ids[np.where(grid_levels == grid.Level)[0][0]]
        lo = grid.id - grid_level_offset
        hi = lo + 1

        # handle the case where this grid has no particles
        if offsets[lo] == offsets[hi]:
            return np.array([], dtype=np.float64)

        data = self._handle[lev]["particles:data"][offsets[lo] : offsets[hi]]
        return np.asarray(
            data[field_index::items_per_particle], dtype=np.float64, order="F"
        )


def parse_orion_sinks(fn):
    r"""
    Orion sink particles are stored in text files. This function
    is for figuring what particle fields are present based on the
    number of entries per line in the \*.sink file.
    """

    # Figure out the format of the particle file
    with open(fn) as f:
        lines = f.readlines()

    try:
        line = lines[1]
    except IndexError:
        # a particle file exists, but there is only one line,
        # so no sinks have been created yet.
        index = {}
        return index

    # The basic fields that all sink particles have
    index = {
        "particle_mass": 0,
        "particle_position_x": 1,
        "particle_position_y": 2,
        "particle_position_z": 3,
        "particle_momentum_x": 4,
        "particle_momentum_y": 5,
        "particle_momentum_z": 6,
        "particle_angmomen_x": 7,
        "particle_angmomen_y": 8,
        "particle_angmomen_z": 9,
        "particle_id": -1,
    }

    if len(line.strip().split()) == 11:
        # these are vanilla sinks, do nothing
        pass

    elif len(line.strip().split()) == 17:
        # these are old-style stars, add stellar model parameters
        index["particle_mlast"] = 10
        index["particle_r"] = 11
        index["particle_mdeut"] = 12
        index["particle_n"] = 13
        index["particle_mdot"] = 14
        index["particle_burnstate"] = 15

    elif len(line.strip().split()) == 18 or len(line.strip().split()) == 19:
        # these are the newer style, add luminosity as well
        index["particle_mlast"] = 10
        index["particle_r"] = 11
        index["particle_mdeut"] = 12
        index["particle_n"] = 13
        index["particle_mdot"] = 14
        index["particle_burnstate"] = 15
        index["particle_luminosity"] = 16
    else:
        # give a warning if none of the above apply:
        mylog.warning("Warning - could not figure out particle output file")
        mylog.warning("These results could be nonsense!")

    return index


class IOHandlerOrion2HDF5(IOHandlerChomboHDF5):
    _dataset_type = "orion_chombo_native"

    _particle_field_index = None

    @property
    def particle_field_index(self):

        fn = self.ds.fullplotdir[:-4] + "sink"

        index = parse_orion_sinks(fn)

        self._particle_field_index = index
        return self._particle_field_index

    def _read_particles(self, grid, field):
        """
        parses the Orion Star Particle text files

        """

        particles = []

        if grid.NumberOfParticles == 0:
            return np.array(particles)

        def read(line, field):
            entry = line.strip().split(" ")[self.particle_field_index[field]]
            return float(entry)

        try:
            lines = self._cached_lines
            for num in grid._particle_line_numbers:
                line = lines[num]
                particles.append(read(line, field))
            return np.array(particles)
        except AttributeError:
            fn = grid.ds.fullplotdir[:-4] + "sink"
            with open(fn) as f:
                lines = f.readlines()
                self._cached_lines = lines
                for num in grid._particle_line_numbers:
                    line = lines[num]
                    particles.append(read(line, field))
            return np.array(particles)
