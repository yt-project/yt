import numpy as np

from yt.config import ytcfg
from yt.data_objects.field_data import YTFieldData
from yt.funcs import get_pbar, mylog
from yt.units.yt_array import array_like_field
from yt.utilities.exceptions import YTIllDefinedParticleData
from yt.utilities.lib.particle_mesh_operations import CICSample_3
from yt.utilities.on_demand_imports import _h5py as h5py
from yt.utilities.parallel_tools.parallel_analysis_interface import parallel_root_only


class ParticleTrajectories:
    r"""A collection of particle trajectories in time over a series of
    datasets.

    Parameters
    ----------
    outputs : ~yt.data_objects.time_series.DatasetSeries
        DatasetSeries object from which to draw the particles.
    indices : array_like
        An integer array of particle indices whose trajectories we
        want to track. If they are not sorted they will be sorted.
    fields : list of strings, optional
        A set of fields that is retrieved when the trajectory
        collection is instantiated. Default: None (will default
        to the fields 'particle_position_x', 'particle_position_y',
        'particle_position_z')
    suppress_logging : boolean
        Suppress yt's logging when iterating over the simulation time
        series. Default: False
    ptype : str, optional
        Only use this particle type. Default: None, which uses all particle type.

    Examples
    --------
    >>> my_fns = glob.glob("orbit_hdf5_chk_00[0-9][0-9]")
    >>> my_fns.sort()
    >>> fields = [
    ...     ("all", "particle_position_x"),
    ...     ("all", "particle_position_y"),
    ...     ("all", "particle_position_z"),
    ...     ("all", "particle_velocity_x"),
    ...     ("all", "particle_velocity_y"),
    ...     ("all", "particle_velocity_z"),
    ... ]
    >>> ds = load(my_fns[0])
    >>> init_sphere = ds.sphere(ds.domain_center, (0.5, "unitary"))
    >>> indices = init_sphere[("all", "particle_index")].astype("int")
    >>> ts = DatasetSeries(my_fns)
    >>> trajs = ts.particle_trajectories(indices, fields=fields)
    >>> for t in trajs:
    ...     print(
    ...         t[("all", "particle_velocity_x")].max(),
    ...         t[("all", "particle_velocity_x")].min(),
    ...     )
    """

    def __init__(
        self, outputs, indices, fields=None, suppress_logging=False, ptype=None
    ):

        indices.sort()  # Just in case the caller wasn't careful
        self.field_data = YTFieldData()
        self.data_series = outputs
        self.masks = []
        self.sorts = []
        self.array_indices = []
        self.indices = indices
        self.num_indices = len(indices)
        self.num_steps = len(outputs)
        self.times = []
        self.suppress_logging = suppress_logging
        self.ptype = ptype if ptype else "all"

        if fields is None:
            fields = []

        if self.suppress_logging:
            old_level = int(ytcfg.get("yt", "log_level"))
            mylog.setLevel(40)
        ds_first = self.data_series[0]
        dd_first = ds_first.all_data()

        fds = {}
        for field in (
            "particle_index",
            "particle_position_x",
            "particle_position_y",
            "particle_position_z",
        ):
            fds[field] = dd_first._determine_fields((self.ptype, field))[0]

        # Note: we explicitly pass dynamic=False to prevent any change in piter from
        # breaking the assumption that the same processors load the same datasets
        my_storage = {}
        pbar = get_pbar("Constructing trajectory information", len(self.data_series))
        for i, (sto, ds) in enumerate(
            self.data_series.piter(storage=my_storage, dynamic=False)
        ):
            dd = ds.all_data()
            newtags = dd[fds["particle_index"]].d.astype("int64")
            mask = np.in1d(newtags, indices, assume_unique=True)
            sort = np.argsort(newtags[mask])
            array_indices = np.where(np.in1d(indices, newtags, assume_unique=True))[0]
            self.array_indices.append(array_indices)
            self.masks.append(mask)
            self.sorts.append(sort)

            pfields = {}
            for field in (f"particle_position_{ax}" for ax in "xyz"):
                pfields[field] = dd[fds[field]].ndarray_view()[mask][sort]

            sto.result_id = ds.parameter_filename
            sto.result = (ds.current_time, array_indices, pfields)
            pbar.update(i + 1)
        pbar.finish()

        if self.suppress_logging:
            mylog.setLevel(old_level)

        sorted_storage = sorted(my_storage.items())
        _fn, (time, *_) = sorted_storage[0]
        time_units = time.units
        times = [time.to(time_units) for _fn, (time, *_) in sorted_storage]
        self.times = self.data_series[0].arr([time.value for time in times], time_units)

        self.particle_fields = []
        output_field = np.empty((self.num_indices, self.num_steps))
        output_field.fill(np.nan)
        for field in (f"particle_position_{ax}" for ax in "xyz"):
            for i, (_fn, (_time, indices, pfields)) in enumerate(sorted_storage):
                try:
                    # This will fail if particles ids are
                    # duplicate. This is due to the fact that the rhs
                    # would then have a different shape as the lhs
                    output_field[indices, i] = pfields[field]
                except ValueError as e:
                    raise YTIllDefinedParticleData(
                        "This dataset contains duplicate particle indices!"
                    ) from e
            self.field_data[field] = array_like_field(
                dd_first, output_field.copy(), fds[field]
            )
            self.particle_fields.append(field)

        # Instantiate fields the caller requested
        self._get_data(fields)

    def has_key(self, key):
        return key in self.field_data

    def keys(self):
        return self.field_data.keys()

    def __getitem__(self, key):
        """
        Get the field associated with key.
        """
        if key == "particle_time":
            return self.times
        if key not in self.field_data:
            self._get_data([key])
        return self.field_data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.field_data[key] = val

    def __delitem__(self, key):
        """
        Delete the field from the trajectory
        """
        del self.field_data[key]

    def __iter__(self):
        """
        This iterates over the trajectories for
        the different particles, returning dicts
        of fields for each trajectory
        """
        for idx in range(self.num_indices):
            traj = {}
            traj["particle_index"] = self.indices[idx]
            traj["particle_time"] = self.times
            for field in self.field_data.keys():
                traj[field] = self[field][idx, :]
            yield traj

    def __len__(self):
        """
        The number of individual trajectories
        """
        return self.num_indices

    def add_fields(self, fields):
        """
        Add a list of fields to an existing trajectory

        Parameters
        ----------
        fields : list of strings
            A list of fields to be added to the current trajectory
            collection.

        Examples
        ________
        >>> trajs = ParticleTrajectories(my_fns, indices)
        >>> trajs.add_fields([("all", "particle_mass"), ("all", "particle_gpot")])
        """
        self._get_data(fields)

    def _get_data(self, fields):
        """
        Get a list of fields to include in the trajectory collection.
        The trajectory collection itself is a dict of 2D numpy arrays,
        with shape (num_indices, num_steps)
        """

        missing_fields = [field for field in fields if field not in self.field_data]
        if not missing_fields:
            return

        if self.suppress_logging:
            old_level = int(ytcfg.get("yt", "log_level"))
            mylog.setLevel(40)
        ds_first = self.data_series[0]
        dd_first = ds_first.all_data()

        fds = {}
        new_particle_fields = []
        for field in missing_fields:
            fds[field] = dd_first._determine_fields(field)[0]
            if field not in self.particle_fields:
                ftype = fds[field][0]
                if ftype in self.data_series[0].particle_types:
                    self.particle_fields.append(field)
                    new_particle_fields.append(field)

        grid_fields = [
            field for field in missing_fields if field not in self.particle_fields
        ]
        step = int(0)
        fields_str = ", ".join(str(f) for f in missing_fields)
        pbar = get_pbar(
            f"Generating [{fields_str}] fields in trajectories",
            self.num_steps,
        )

        # Note: we explicitly pass dynamic=False to prevent any change in piter from
        # breaking the assumption that the same processors load the same datasets
        my_storage = {}
        for i, (sto, ds) in enumerate(
            self.data_series.piter(storage=my_storage, dynamic=False)
        ):
            mask = self.masks[i]
            sort = self.sorts[i]
            pfield = {}

            if new_particle_fields:  # there's at least one particle field
                dd = ds.all_data()
                for field in new_particle_fields:
                    # This is easy... just get the particle fields
                    pfield[field] = dd[fds[field]].d[mask][sort]

            if grid_fields:
                # This is hard... must loop over grids
                for field in grid_fields:
                    pfield[field] = np.zeros(self.num_indices)
                x = self["particle_position_x"][:, step].d
                y = self["particle_position_y"][:, step].d
                z = self["particle_position_z"][:, step].d
                particle_grids, particle_grid_inds = ds.index._find_points(x, y, z)

                # This will fail for non-grid index objects
                for grid in particle_grids:
                    cube = grid.retrieve_ghost_zones(1, grid_fields)
                    for field in grid_fields:
                        CICSample_3(
                            x,
                            y,
                            z,
                            pfield[field],
                            self.num_indices,
                            cube[fds[field]],
                            np.array(grid.LeftEdge).astype(np.float64),
                            np.array(grid.ActiveDimensions).astype(np.int32),
                            grid.dds[0],
                        )
            sto.result_id = ds.parameter_filename
            sto.result = (self.array_indices[i], pfield)
            pbar.update(step)
            step += 1
        pbar.finish()

        output_field = np.empty((self.num_indices, self.num_steps))
        output_field.fill(np.nan)
        for field in missing_fields:
            fd = fds[field]
            for i, (_fn, (indices, pfield)) in enumerate(sorted(my_storage.items())):
                output_field[indices, i] = pfield[field]
            self.field_data[field] = array_like_field(dd_first, output_field.copy(), fd)

        if self.suppress_logging:
            mylog.setLevel(old_level)

    def trajectory_from_index(self, index):
        """
        Retrieve a single trajectory corresponding to a specific particle
        index

        Parameters
        ----------
        index : int
            This defines which particle trajectory from the
            ParticleTrajectories object will be returned.

        Returns
        -------
        A dictionary corresponding to the particle's trajectory and the
        fields along that trajectory

        Examples
        --------
        >>> from yt.mods import *
        >>> import matplotlib.pyplot as plt
        >>> trajs = ParticleTrajectories(my_fns, indices)
        >>> traj = trajs.trajectory_from_index(indices[0])
        >>> plt.plot(
        ...     traj[("all", "particle_time")],
        ...     traj[("all", "particle_position_x")],
        ...     "-x",
        ... )
        >>> plt.savefig("orbit")
        """
        mask = np.in1d(self.indices, (index,), assume_unique=True)
        if not np.any(mask):
            print("The particle index %d is not in the list!" % (index))
            raise IndexError
        fields = [field for field in sorted(self.field_data.keys())]
        traj = {}
        traj["particle_time"] = self.times
        traj["particle_index"] = index
        for field in fields:
            traj[field] = self[field][mask, :][0]
        return traj

    @parallel_root_only
    def write_out(self, filename_base):
        """
        Write out particle trajectories to tab-separated ASCII files (one
        for each trajectory) with the field names in the file header. Each
        file is named with a basename and the index number.

        Parameters
        ----------
        filename_base : string
            The prefix for the outputted ASCII files.

        Examples
        --------
        >>> trajs = ParticleTrajectories(my_fns, indices)
        >>> trajs.write_out("orbit_trajectory")
        """
        fields = [field for field in sorted(self.field_data.keys())]
        num_fields = len(fields)
        first_str = "# particle_time\t" + "\t".join(fields) + "\n"
        template_str = "%g\t" * num_fields + "%g\n"
        for ix in range(self.num_indices):
            outlines = [first_str]
            for it in range(self.num_steps):
                outlines.append(
                    template_str
                    % tuple(
                        [self.times[it]] + [self[field][ix, it] for field in fields]
                    )
                )
            fid = open(filename_base + "_%d.dat" % self.indices[ix], "w")
            fid.writelines(outlines)
            fid.close()
            del fid

    @parallel_root_only
    def write_out_h5(self, filename):
        """
        Write out all the particle trajectories to a single HDF5 file
        that contains the indices, the times, and the 2D array for each
        field individually

        Parameters
        ----------

        filename : string
            The output filename for the HDF5 file

        Examples
        --------
        >>> trajs = ParticleTrajectories(my_fns, indices)
        >>> trajs.write_out_h5("orbit_trajectories")
        """
        fid = h5py.File(filename, mode="w")
        fid.create_dataset("particle_indices", dtype=np.int64, data=self.indices)
        fid.close()
        self.times.write_hdf5(filename, dataset_name="particle_times")
        fields = [field for field in sorted(self.field_data.keys())]
        for field in fields:
            self[field].write_hdf5(filename, dataset_name=f"{field}")
