import os

import numpy as np

from yt.funcs import mylog
from yt.geometry.selection_routines import GridSelector


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


class BoxlibParticleSelectionMixin:
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


class BoxlibReadParticleFileMixin:
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


class BoxlibSetupParticleFieldsMixin:
    def setup_particle_fields(self, ptype):
        def _get_vel(axis):
            def velocity(field, data):
                return (
                    data[(ptype, f"particle_momentum_{axis}")]
                    / data[(ptype, "particle_mass")]
                )

            return velocity

        for ax in "xyz":
            self.add_field(
                (ptype, f"particle_velocity_{ax}"),
                sampling_type="particle",
                function=_get_vel(ax),
                units="code_length/code_time",
            )

        super().setup_particle_fields(ptype)


class IOHandlerParticlesBoxlibMixin:
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
            fn = self.particle_filename
            with open(fn) as f:
                lines = f.readlines()
                self._cached_lines = lines
                for num in grid._particle_line_numbers:
                    line = lines[num]
                    particles.append(read(line, field))
            return np.array(particles)

    _particle_field_index = None

    @property
    def particle_field_index(self):
        if self._particle_field_index:
            return self._particle_field_index

        index = parse_orion_sinks(self.particle_filename)

        self._particle_field_index = index
        return self._particle_field_index
