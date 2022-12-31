import numpy as np

from yt.funcs import get_pbar
from yt.units._numpy_wrapper_functions import uconcatenate
from yt.utilities.lib.particle_mesh_operations import CICSample_3


class ParticleGenerator:
    def __init__(self, ds, num_particles, field_list, ptype="io"):
        """
        Base class for generating particle fields which may be applied to
        streams. Normally this would not be called directly, since it doesn't
        really do anything except allocate memory. Takes a *ds* to serve as the
        basis for determining grids, the number of particles *num_particles*,
        a list of fields, *field_list*, and the particle type *ptype*, which
        has a default value of "io".
        """
        self.ds = ds
        self.num_particles = num_particles
        self.field_list = [
            (ptype, fd) if isinstance(fd, str) else fd for fd in field_list
        ]
        self.field_list.append((ptype, "particle_index"))
        self.field_units = {
            (ptype, f"particle_position_{ax}"): "code_length" for ax in "xyz"
        }
        self.field_units[ptype, "particle_index"] = ""
        self.ptype = ptype

        self._set_default_fields()

        try:
            self.posx_index = self.field_list.index(self.default_fields[0])
            self.posy_index = self.field_list.index(self.default_fields[1])
            self.posz_index = self.field_list.index(self.default_fields[2])
        except Exception as e:
            raise KeyError(
                "You must specify position fields: "
                + " ".join("particle_position_%s" % ax for ax in "xyz")
            ) from e
        self.index_index = self.field_list.index((ptype, "particle_index"))

        self.num_grids = self.ds.index.num_grids
        self.NumberOfParticles = np.zeros(self.num_grids, dtype="int64")
        self.ParticleGridIndices = np.zeros(self.num_grids + 1, dtype="int64")

        self.num_fields = len(self.field_list)

        self.particles = np.zeros(
            (self.num_particles, self.num_fields), dtype="float64"
        )

    def _set_default_fields(self):
        self.default_fields = [
            (self.ptype, "particle_position_x"),
            (self.ptype, "particle_position_y"),
            (self.ptype, "particle_position_z"),
        ]

    def has_key(self, key):
        """
        Check to see if *key* is in the particle field list.
        """
        return key in self.field_list

    def keys(self):
        """
        Return the list of particle fields.
        """
        return self.field_list

    def __getitem__(self, key):
        """
        Get the field associated with key.
        """
        return self.particles[:, self.field_list.index(key)]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value. Note that we assume
        that the particles have been sorted by grid already, so
        make sure the setting of the field is consistent with this.
        """
        self.particles[:, self.field_list.index(key)] = val[:]

    def __len__(self):
        """
        The number of particles
        """
        return self.num_particles

    def get_for_grid(self, grid):
        """
        Return a dict containing all of the particle fields in the specified grid.
        """
        ind = grid.id - grid._id_offset
        start = self.ParticleGridIndices[ind]
        end = self.ParticleGridIndices[ind + 1]
        tr = {}
        for field in self.field_list:
            fi = self.field_list.index(field)
            if field in self.field_units:
                tr[field] = self.ds.arr(
                    self.particles[start:end, fi], self.field_units[field]
                )
            else:
                tr[field] = self.particles[start:end, fi]
        return tr

    def _setup_particles(self, x, y, z, setup_fields=None):
        """
        Assigns grids to particles and sets up particle positions. *setup_fields* is
        a dict of fields other than the particle positions to set up.
        """
        particle_grids, particle_grid_inds = self.ds.index._find_points(x, y, z)
        idxs = np.argsort(particle_grid_inds)
        self.particles[:, self.posx_index] = x[idxs]
        self.particles[:, self.posy_index] = y[idxs]
        self.particles[:, self.posz_index] = z[idxs]
        self.NumberOfParticles = np.bincount(
            particle_grid_inds.astype("intp"), minlength=self.num_grids
        )
        if self.num_grids > 1:
            np.add.accumulate(
                self.NumberOfParticles.squeeze(), out=self.ParticleGridIndices[1:]
            )
        else:
            self.ParticleGridIndices[1] = self.NumberOfParticles.squeeze()
        if setup_fields is not None:
            for key, value in setup_fields.items():
                field = (self.ptype, key) if isinstance(key, str) else key
                if field not in self.default_fields:
                    self.particles[:, self.field_list.index(field)] = value[idxs]

    def assign_indices(self, function=None, **kwargs):
        """
        Assign unique indices to the particles. The default is to just use
        numpy.arange, but any function may be supplied with keyword arguments.
        """
        if function is None:
            self.particles[:, self.index_index] = np.arange(self.num_particles)
        else:
            self.particles[:, self.index_index] = function(**kwargs)

    def map_grid_fields_to_particles(self, mapping_dict):
        r"""
        For the fields in  *mapping_dict*, map grid fields to the particles
        using CIC sampling.

        Examples
        --------
        >>> field_map = {
        ...     "density": "particle_density",
        ...     "temperature": "particle_temperature",
        ... }
        >>> particles.map_grid_fields_to_particles(field_map)
        """
        pbar = get_pbar("Mapping fields to particles", self.num_grids)
        for i, grid in enumerate(self.ds.index.grids):
            pbar.update(i + 1)
            if self.NumberOfParticles[i] > 0:
                start = self.ParticleGridIndices[i]
                end = self.ParticleGridIndices[i + 1]
                # Note we add one ghost zone to the grid!
                cube = grid.retrieve_ghost_zones(1, list(mapping_dict.keys()))
                le = np.array(grid.LeftEdge).astype(np.float64)
                dims = np.array(grid.ActiveDimensions).astype(np.int32)
                for gfield, pfield in mapping_dict.items():
                    self.field_units[pfield] = cube[gfield].units
                    field_index = self.field_list.index(pfield)
                    CICSample_3(
                        self.particles[start:end, self.posx_index],
                        self.particles[start:end, self.posy_index],
                        self.particles[start:end, self.posz_index],
                        self.particles[start:end, field_index],
                        np.int64(self.NumberOfParticles[i]),
                        cube[gfield],
                        le,
                        dims,
                        grid.dds[0],
                    )
        pbar.finish()

    def apply_to_stream(self, overwrite=False, **kwargs):
        """
        Apply the particles to a grid-based stream dataset. If particles
        already exist, and overwrite=False, do not overwrite them, but add
        the new ones to them.
        """
        grid_data = []
        for i, g in enumerate(self.ds.index.grids):
            data = {}
            number_of_particles = self.NumberOfParticles[i]
            if not overwrite:
                number_of_particles += g.NumberOfParticles
            grid_particles = self.get_for_grid(g)
            for field in self.field_list:
                if number_of_particles > 0:
                    if (
                        g.NumberOfParticles > 0
                        and not overwrite
                        and field in self.ds.field_list
                    ):
                        # We have particles in this grid, we're not
                        # overwriting them, and the field is in the field
                        # list already
                        data[field] = uconcatenate([g[field], grid_particles[field]])
                    else:
                        # Otherwise, simply add the field in
                        data[field] = grid_particles[field]
                else:
                    # We don't have particles in this grid
                    data[field] = np.array([], dtype="float64")
            grid_data.append(data)
        self.ds.index.update_data(grid_data)


class FromListParticleGenerator(ParticleGenerator):
    def __init__(self, ds, num_particles, data, ptype="io"):
        r"""
        Generate particle fields from array-like lists contained in a dict.

        Parameters
        ----------
        ds : `Dataset`
            The dataset which will serve as the base for these particles.
        num_particles : int
            The number of particles in the dict.
        data : dict of NumPy arrays
            The particle fields themselves.
        ptype : string, optional
            The particle type for these particle fields. Default: "io"

        Examples
        --------
        >>> num_p = 100000
        >>> posx = np.random.random(num_p)
        >>> posy = np.random.random(num_p)
        >>> posz = np.random.random(num_p)
        >>> mass = np.ones(num_p)
        >>> data = {
        ...     "particle_position_x": posx,
        ...     "particle_position_y": posy,
        ...     "particle_position_z": posz,
        ...     "particle_mass": mass,
        ... }
        >>> particles = FromListParticleGenerator(ds, num_p, data)
        """

        field_list = list(data.keys())
        if "particle_position_x" in data:
            x = data.pop("particle_position_x")
            y = data.pop("particle_position_y")
            z = data.pop("particle_position_z")
        elif (ptype, "particle_position_x") in data:
            x = data.pop((ptype, "particle_position_x"))
            y = data.pop((ptype, "particle_position_y"))
            z = data.pop((ptype, "particle_position_z"))

        xcond = np.logical_or(x < ds.domain_left_edge[0], x >= ds.domain_right_edge[0])
        ycond = np.logical_or(y < ds.domain_left_edge[1], y >= ds.domain_right_edge[1])
        zcond = np.logical_or(z < ds.domain_left_edge[2], z >= ds.domain_right_edge[2])
        cond = np.logical_or(xcond, ycond)
        cond = np.logical_or(zcond, cond)

        if np.any(cond):
            raise ValueError("Some particles are outside of the domain!!!")

        super().__init__(ds, num_particles, field_list, ptype=ptype)
        self._setup_particles(x, y, z, setup_fields=data)


class LatticeParticleGenerator(ParticleGenerator):
    def __init__(
        self,
        ds,
        particles_dims,
        particles_left_edge,
        particles_right_edge,
        field_list,
        ptype="io",
    ):
        r"""
        Generate particles in a lattice arrangement.

        Parameters
        ----------
        ds : `Dataset`
            The dataset which will serve as the base for these particles.
        particles_dims : int, array-like
            The number of particles along each dimension
        particles_left_edge : float, array-like
            The 'left-most' starting positions of the lattice.
        particles_right_edge : float, array-like
             The 'right-most' ending positions of the lattice.
        field_list : list of strings
             A list of particle fields
        ptype : string, optional
            The particle type for these particle fields. Default: "io"

        Examples
        --------
        >>> dims = (128, 128, 128)
        >>> le = np.array([0.25, 0.25, 0.25])
        >>> re = np.array([0.75, 0.75, 0.75])
        >>> fields = [
        ...     ("all", "particle_position_x"),
        ...     ("all", "particle_position_y"),
        ...     ("all", "particle_position_z"),
        ...     ("all", "particle_density"),
        ...     ("all", "particle_temperature"),
        ... ]
        >>> particles = LatticeParticleGenerator(ds, dims, le, re, fields)
        """

        num_x = particles_dims[0]
        num_y = particles_dims[1]
        num_z = particles_dims[2]
        xmin = particles_left_edge[0]
        ymin = particles_left_edge[1]
        zmin = particles_left_edge[2]
        xmax = particles_right_edge[0]
        ymax = particles_right_edge[1]
        zmax = particles_right_edge[2]
        DLE = ds.domain_left_edge.in_units("code_length").d
        DRE = ds.domain_right_edge.in_units("code_length").d

        xcond = (xmin < DLE[0]) or (xmax >= DRE[0])
        ycond = (ymin < DLE[1]) or (ymax >= DRE[1])
        zcond = (zmin < DLE[2]) or (zmax >= DRE[2])
        cond = xcond or ycond or zcond

        if cond:
            raise ValueError("Proposed bounds for particles are outside domain!!!")

        super().__init__(ds, num_x * num_y * num_z, field_list, ptype=ptype)

        dx = (xmax - xmin) / (num_x - 1)
        dy = (ymax - ymin) / (num_y - 1)
        dz = (zmax - zmin) / (num_z - 1)
        inds = np.indices((num_x, num_y, num_z))
        xpos = inds[0] * dx + xmin
        ypos = inds[1] * dy + ymin
        zpos = inds[2] * dz + zmin

        self._setup_particles(xpos.flat[:], ypos.flat[:], zpos.flat[:])


class WithDensityParticleGenerator(ParticleGenerator):
    def __init__(
        self,
        ds,
        data_source,
        num_particles,
        field_list,
        density_field=("gas", "density"),
        ptype="io",
    ):
        r"""
        Generate particles based on a density field.

        Parameters
        ----------
        ds : `Dataset`
            The dataset which will serve as the base for these particles.
        data_source :
            `yt.data_objects.selection_objects.base_objects.YTSelectionContainer`
            The data source containing the density field.
        num_particles : int
            The number of particles to be generated
        field_list : list of strings
            A list of particle fields
        density_field : tuple, optional
            A density field which will serve as the distribution function for the
            particle positions. Theoretically, this could be any 'per-volume' field.
        ptype : string, optional
            The particle type for these particle fields. Default: "io"

        Examples
        --------
        >>> sphere = ds.sphere(ds.domain_center, 0.5)
        >>> num_p = 100000
        >>> fields = [
        ...     ("all", "particle_position_x"),
        ...     ("all", "particle_position_y"),
        ...     ("all", "particle_position_z"),
        ...     ("all", "particle_density"),
        ...     ("all", "particle_temperature"),
        ... ]
        >>> particles = WithDensityParticleGenerator(
        ...     ds, sphere, num_particles, fields, density_field="Dark_Matter_Density"
        ... )
        """

        super().__init__(ds, num_particles, field_list, ptype=ptype)

        num_cells = len(data_source[("index", "x")].flat)
        max_mass = (
            data_source[density_field] * data_source[("gas", "cell_volume")]
        ).max()
        num_particles_left = num_particles
        all_x = []
        all_y = []
        all_z = []

        pbar = get_pbar("Generating Particles", num_particles)
        tot_num_accepted = int(0)

        while num_particles_left > 0:

            m = np.random.uniform(high=1.01 * max_mass, size=num_particles_left)
            idxs = np.random.randint(low=0, high=num_cells, size=num_particles_left)
            m_true = (
                data_source[density_field] * data_source[("gas", "cell_volume")]
            ).flat[idxs]
            accept = m <= m_true
            num_accepted = accept.sum()
            accepted_idxs = idxs[accept]

            xpos = (
                data_source[("index", "x")].flat[accepted_idxs]
                + np.random.uniform(low=-0.5, high=0.5, size=num_accepted)
                * data_source[("index", "dx")].flat[accepted_idxs]
            )
            ypos = (
                data_source[("index", "y")].flat[accepted_idxs]
                + np.random.uniform(low=-0.5, high=0.5, size=num_accepted)
                * data_source[("index", "dy")].flat[accepted_idxs]
            )
            zpos = (
                data_source[("index", "z")].flat[accepted_idxs]
                + np.random.uniform(low=-0.5, high=0.5, size=num_accepted)
                * data_source[("index", "dz")].flat[accepted_idxs]
            )

            all_x.append(xpos)
            all_y.append(ypos)
            all_z.append(zpos)

            num_particles_left -= num_accepted
            tot_num_accepted += num_accepted
            pbar.update(tot_num_accepted)

        pbar.finish()

        x = uconcatenate(all_x)
        y = uconcatenate(all_y)
        z = uconcatenate(all_z)

        self._setup_particles(x, y, z)
