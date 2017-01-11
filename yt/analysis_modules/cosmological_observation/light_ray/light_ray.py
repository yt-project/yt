"""
LightRay class and member functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.analysis_modules.cosmological_observation.cosmology_splice import \
    CosmologySplice
from yt.convenience import \
    load
from yt.frontends.ytdata.utilities import \
    save_as_dataset
from yt.units.yt_array import \
    YTArray
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.logger import \
    ytLogger as mylog
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, \
    parallel_root_only
from yt.utilities.physical_constants import speed_of_light_cgs
from yt.data_objects.static_output import Dataset

class LightRay(CosmologySplice):
    """
    A 1D object representing the path of a light ray passing through a 
    simulation.  LightRays can be either simple, where they pass through a 
    single dataset, or compound, where they pass through consecutive 
    datasets from the same cosmological simulation.  One can sample any of 
    the fields intersected by the LightRay object as it passed through 
    the dataset(s).
    
    For compound rays, the LightRay stacks together multiple datasets in a time
    series in order to approximate a LightRay's path through a volume
    and redshift interval larger than a single simulation data output.
    The outcome is something akin to a synthetic QSO line of sight.

    Once the LightRay object is set up, use LightRay.make_light_ray to
    begin making rays.  Different randomizations can be created with a
    single object by providing different random seeds to make_light_ray.

    Parameters
    ----------
    parameter_filename : string or :class:`yt.data_objects.static_output.Dataset`
        For simple rays, one may pass either a loaded dataset object or
        the filename of a dataset.
        For compound rays, one must pass the filename of the simulation
        parameter file.
    simulation_type : optional, string
        This refers to the simulation frontend type.  Do not use for simple 
        rays.
        Default: None
    near_redshift : optional, float
        The near (lowest) redshift for a light ray containing multiple
        datasets.  Do not use for simple rays.
        Default: None
    far_redshift : optional, float
        The far (highest) redshift for a light ray containing multiple
        datasets.  Do not use for simple rays.
        Default: None
    use_minimum_datasets : optional, bool
        If True, the minimum number of datasets is used to connect the
        initial and final redshift.  If false, the light ray solution
        will contain as many entries as possible within the redshift
        interval.  Do not use for simple rays.
        Default: True.
    max_box_fraction : optional, float
        In terms of the size of the domain, the maximum length a light
        ray segment can be in order to span the redshift interval from
        one dataset to another.  If using a zoom-in simulation, this
        parameter can be set to the length of the high resolution
        region so as to limit ray segments to that size.  If the
        high resolution region is not cubical, the smallest side
        should be used.
        Default: 1.0 (the size of the box)
    deltaz_min : optional, float
        Specifies the minimum :math:`\Delta z` between consecutive
        datasets in the returned list.  Do not use for simple rays.
        Default: 0.0.
    minimum_coherent_box_fraction : optional, float
        Use to specify the minimum length of a ray, in terms of the
        size of the domain, before the trajectory is re-randomized.
        Set to 0 to have ray trajectory randomized for every dataset.
        Set to np.inf (infinity) to use a single trajectory for the
        entire ray.
        Default: 0.
    time_data : optional, bool
        Whether or not to include time outputs when gathering
        datasets for time series.  Do not use for simple rays.
        Default: True.
    redshift_data : optional, bool
        Whether or not to include redshift outputs when gathering
        datasets for time series.  Do not use for simple rays.
        Default: True.
    find_outputs : optional, bool
        Whether or not to search for datasets in the current
        directory.  Do not use for simple rays.
        Default: False.
    load_kwargs : optional, dict
        If you are passing a filename of a dataset to LightRay rather than an 
        already loaded dataset, then you can optionally provide this dictionary 
        as keywords when the dataset is loaded by yt with the "load" function.
        Necessary for use with certain frontends.  E.g.
        Tipsy using "bounding_box"
        Gadget using "unit_base", etc.
        Default : None

    """
    def __init__(self, parameter_filename, simulation_type=None,
                 near_redshift=None, far_redshift=None,
                 use_minimum_datasets=True, max_box_fraction=1.0,
                 deltaz_min=0.0, minimum_coherent_box_fraction=0.0,
                 time_data=True, redshift_data=True,
                 find_outputs=False, load_kwargs=None):

        if near_redshift is not None and far_redshift is not None and \
          near_redshift >= far_redshift:
            raise RuntimeError(
                "near_redshift must be less than far_redshift.")

        self.near_redshift = near_redshift
        self.far_redshift = far_redshift
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        self.parameter_filename = parameter_filename
        if load_kwargs is None:
            self.load_kwargs = {}
        else:
            self.load_kwargs = load_kwargs
        self.light_ray_solution = []
        self._data = {}

        # The options here are:
        # 1) User passed us a dataset: use it to make a simple ray
        # 2) User passed us a dataset filename: use it to make a simple ray
        # 3) User passed us a simulation filename: use it to make a compound ray

        # Make a light ray from a single, given dataset: #1, #2
        if simulation_type is None:     
            self.simulation_type = simulation_type
            if isinstance(self.parameter_filename, Dataset):
                self.ds = self.parameter_filename
                self.parameter_filename = self.ds.basename
            elif isinstance(self.parameter_filename, str):
                self.ds = load(self.parameter_filename, **self.load_kwargs)
            if self.ds.cosmological_simulation:
                redshift = self.ds.current_redshift
                self.cosmology = Cosmology(
                    hubble_constant=self.ds.hubble_constant,
                    omega_matter=self.ds.omega_matter,
                    omega_lambda=self.ds.omega_lambda)
            else:
                redshift = 0.
            self.light_ray_solution.append({"filename": self.parameter_filename,
                                            "redshift": redshift})

        # Make a light ray from a simulation time-series. #3
        else:
            self.ds = None
            assert isinstance(self.parameter_filename, str)
            # Get list of datasets for light ray solution.
            CosmologySplice.__init__(self, self.parameter_filename, simulation_type,
                                     find_outputs=find_outputs)
            self.light_ray_solution = \
              self.create_cosmology_splice(
                  self.near_redshift, self.far_redshift,
                  minimal=self.use_minimum_datasets,
                  max_box_fraction=max_box_fraction,
                  deltaz_min=self.deltaz_min,
                  time_data=time_data,
                  redshift_data=redshift_data)

    def _calculate_light_ray_solution(self, seed=None,
                                      left_edge=None, right_edge=None,
                                      min_level=None, periodic=True,
                                      start_position=None, end_position=None,
                                      trajectory=None, filename=None):
        "Create list of datasets to be added together to make the light ray."

        # Calculate dataset sizes, and get random dataset axes and centers.
        my_random = np.random.RandomState(seed)

        # If using only one dataset, set start and stop manually.
        if start_position is not None:
            if self.near_redshift is not None or self.far_redshift is not None:
                raise RuntimeError("LightRay Error: cannot specify both " + \
                                   "start_position and a redshift range.")
            if not ((end_position is None) ^ (trajectory is None)):
                raise RuntimeError("LightRay Error: must specify either end_position " + \
                                   "or trajectory, but not both.")
            self.light_ray_solution[0]['start'] = start_position
            if end_position is not None:
                self.light_ray_solution[0]['end'] = end_position
            else:
                # assume trajectory given as r, theta, phi
                if len(trajectory) != 3:
                    raise RuntimeError("LightRay Error: trajectory must have length 3.")
                r, theta, phi = trajectory
                self.light_ray_solution[0]['end'] = self.light_ray_solution[0]['start'] + \
                  r * np.array([np.cos(phi) * np.sin(theta),
                                np.sin(phi) * np.sin(theta),
                                np.cos(theta)])
            self.light_ray_solution[0]['traversal_box_fraction'] = \
              vector_length(self.light_ray_solution[0]['start'],
                            self.light_ray_solution[0]['end'])

        # the normal way (random start positions and trajectories for each dataset)
        else:

            # For box coherence, keep track of effective depth travelled.
            box_fraction_used = 0.0

            for q in range(len(self.light_ray_solution)):
                if (q == len(self.light_ray_solution) - 1):
                    z_next = self.near_redshift
                else:
                    z_next = self.light_ray_solution[q+1]['redshift']

                # Calculate fraction of box required for a depth of delta z
                self.light_ray_solution[q]['traversal_box_fraction'] = \
                    self.cosmology.comoving_radial_distance(z_next, \
                        self.light_ray_solution[q]['redshift']).in_units("Mpccm / h") / \
                        self.simulation.box_size

                # Get dataset axis and center.
                # If using box coherence, only get start point and vector if
                # enough of the box has been used.
                if (q == 0) or (box_fraction_used >=
                                self.minimum_coherent_box_fraction):
                    if periodic:
                        self.light_ray_solution[q]['start'] = left_edge + \
                          (right_edge - left_edge) * my_random.random_sample(3)
                        theta = np.pi * my_random.random_sample()
                        phi = 2 * np.pi * my_random.random_sample()
                        box_fraction_used = 0.0
                    else:
                        ds = load(self.light_ray_solution[q]["filename"])
                        ray_length = \
                          ds.quan(self.light_ray_solution[q]['traversal_box_fraction'],
                                  "unitary")
                        self.light_ray_solution[q]['start'], \
                          self.light_ray_solution[q]['end'] = \
                          non_periodic_ray(ds, left_edge, right_edge, ray_length,
                                           my_random=my_random, min_level=min_level)
                        del ds
                else:
                    # Use end point of previous segment, adjusted for periodicity,
                    # and the same trajectory.
                    self.light_ray_solution[q]['start'] = \
                      periodic_adjust(self.light_ray_solution[q-1]['end'][:],
                                      left=left_edge, right=right_edge)

                if "end" not in self.light_ray_solution[q]:
                    self.light_ray_solution[q]['end'] = \
                      self.light_ray_solution[q]['start'] + \
                        self.light_ray_solution[q]['traversal_box_fraction'] * \
                        self.simulation.box_size * \
                        np.array([np.cos(phi) * np.sin(theta),
                                  np.sin(phi) * np.sin(theta),
                                  np.cos(theta)])
                box_fraction_used += \
                  self.light_ray_solution[q]['traversal_box_fraction']

        if filename is not None:
            self._write_light_ray_solution(filename,
                extra_info={'parameter_filename':self.parameter_filename,
                            'random_seed':seed,
                            'far_redshift':self.far_redshift,
                            'near_redshift':self.near_redshift})

    def make_light_ray(self, seed=None, periodic=True,
                       left_edge=None, right_edge=None, min_level=None,
                       start_position=None, end_position=None,
                       trajectory=None,
                       fields=None, setup_function=None,
                       solution_filename=None, data_filename=None,
                       get_los_velocity=None, use_peculiar_velocity=True,
                       redshift=None, field_parameters=None, njobs=-1):
        """
        make_light_ray(seed=None, periodic=True,
                       left_edge=None, right_edge=None, min_level=None,
                       start_position=None, end_position=None,
                       trajectory=None, fields=None, setup_function=None,
                       solution_filename=None, data_filename=None,
                       use_peculiar_velocity=True, redshift=None,
                       njobs=-1)

        Create a light ray and get field values for each lixel.  A light
        ray consists of a list of field values for cells intersected by
        the ray and the path length of the ray through those cells.
        Light ray data must be written out to an hdf5 file.

        Parameters
        ----------
        seed : optional, int
            Seed for the random number generator.
            Default: None.
        periodic : optional, bool
            If True, ray trajectories will make use of periodic
            boundaries.  If False, ray trajectories will not be
            periodic.
            Default : True.
        left_edge : optional, iterable of floats or YTArray
            The left corner of the region in which rays are to be
            generated.  If None, the left edge will be that of the
            domain.  If specified without units, it is assumed to
            be in code units.
            Default: None.
        right_edge : optional, iterable of floats or YTArray
            The right corner of the region in which rays are to be
            generated.  If None, the right edge will be that of the
            domain.  If specified without units, it is assumed to
            be in code units.
            Default: None.
        min_level : optional, int
            The minimum refinement level of the spatial region in which
            the ray passes.  This can be used with zoom-in simulations
            where the high resolution region does not keep a constant
            geometry.
            Default: None.
        start_position : optional, iterable of floats or YTArray.
            Used only if creating a light ray from a single dataset.
            The coordinates of the starting position of the ray.
            If specified without units, it is assumed to be in code units.
            Default: None.
        end_position : optional, iterable of floats or YTArray.
            Used only if creating a light ray from a single dataset.
            The coordinates of the ending position of the ray.
            If specified without units, it is assumed to be in code units.
            Default: None.
        trajectory : optional, list of floats
            Used only if creating a light ray from a single dataset.
            The (r, theta, phi) direction of the light ray.  Use either
            end_position or trajectory, not both.
            Default: None.
        fields : optional, list
            A list of fields for which to get data.
            Default: None.
        setup_function : optional, callable, accepts a ds
            This function will be called on each dataset that is loaded
            to create the light ray.  For, example, this can be used to
            add new derived fields.
            Default: None.
        solution_filename : optional, string
            Path to a text file where the trajectories of each
            subray is written out.
            Default: None.
        data_filename : optional, string
            Path to output file for ray data.
            Default: None.
        use_peculiar_velocity : optional, bool
            If True, the peculiar velocity along the ray will be sampled for
            calculating the effective redshift combining the cosmological
            redshift and the doppler redshift.
            Default: True.
        redshift : optional, float
            Used with light rays made from single datasets to specify a
            starting redshift for the ray.  If not used, the starting
            redshift will be 0 for a non-cosmological dataset and
            the dataset redshift for a cosmological dataset.
            Default: None.
        njobs : optional, int
            The number of parallel jobs over which the segments will
            be split.  Choose -1 for one processor per segment.
            Default: -1.

        Examples
        --------

        Make a light ray from multiple datasets:

        >>> import yt
        >>> from yt.analysis_modules.cosmological_observation.light_ray.api import \
        ...     LightRay
        >>> my_ray = LightRay("enzo_tiny_cosmology/32Mpc_32.enzo", "Enzo",
        ...                   0., 0.1, time_data=False)
        ...
        >>> my_ray.make_light_ray(seed=12345,
        ...                       solution_filename="solution.txt",
        ...                       data_filename="my_ray.h5",
        ...                       fields=["temperature", "density"],
        ...                       use_peculiar_velocity=True)

        Make a light ray from a single dataset:

        >>> import yt
        >>> from yt.analysis_modules.cosmological_observation.light_ray.api import \
        ...     LightRay
        >>> my_ray = LightRay("IsolatedGalaxy/galaxy0030/galaxy0030")
        ...
        >>> my_ray.make_light_ray(start_position=[0., 0., 0.],
        ...                       end_position=[1., 1., 1.],
        ...                       solution_filename="solution.txt",
        ...                       data_filename="my_ray.h5",
        ...                       fields=["temperature", "density"],
        ...                       use_peculiar_velocity=True)

        """
        if self.simulation_type is None:
            domain = self.ds
        else:
            domain = self.simulation

        assumed_units = "code_length"
        if left_edge is None:
            left_edge = domain.domain_left_edge
        elif not hasattr(left_edge, 'units'):
            left_edge = domain.arr(left_edge, assumed_units)
        left_edge.convert_to_units('unitary')

        if right_edge is None:
            right_edge = domain.domain_right_edge
        elif not hasattr(right_edge, 'units'):
            right_edge = domain.arr(right_edge, assumed_units)
        right_edge.convert_to_units('unitary')

        if start_position is not None:
            if hasattr(start_position, 'units'):
                start_position = start_position
            else:
                start_position = self.ds.arr(start_position, assumed_units)
            start_position.convert_to_units('unitary')

        if end_position is not None:
            if hasattr(end_position, 'units'):
                end_position = end_position
            else:
                end_position = self.ds.arr(end_position, assumed_units)
            end_position.convert_to_units('unitary')

        if get_los_velocity is not None:
            use_peculiar_velocity = get_los_velocity
            mylog.warn("'get_los_velocity' kwarg is deprecated. " + \
                       "Use 'use_peculiar_velocity' instead.")

        # Calculate solution.
        self._calculate_light_ray_solution(seed=seed,
                                           left_edge=left_edge,
                                           right_edge=right_edge,
                                           min_level=min_level, periodic=periodic,
                                           start_position=start_position,
                                           end_position=end_position,
                                           trajectory=trajectory,
                                           filename=solution_filename)

        if field_parameters is None:
            field_parameters = {}

        # Initialize data structures.
        self._data = {}
        # temperature field is automatically added to fields
        if fields is None: fields = []
        if (('gas', 'temperature') not in fields) and \
           ('temperature' not in fields):
           fields.append(('gas', 'temperature'))
        data_fields = fields[:]
        all_fields = fields[:]
        all_fields.extend(['dl', 'dredshift', 'redshift'])
        all_fields.extend(['x', 'y', 'z', 'dx', 'dy', 'dz'])
        data_fields.extend(['x', 'y', 'z', 'dx', 'dy', 'dz'])
        if use_peculiar_velocity:
            all_fields.extend(['velocity_x', 'velocity_y', 'velocity_z', 
                               'velocity_los', 'redshift_eff', 
                               'redshift_dopp'])
            data_fields.extend(['velocity_x', 'velocity_y', 'velocity_z'])

        all_ray_storage = {}
        for my_storage, my_segment in parallel_objects(self.light_ray_solution,
                                                       storage=all_ray_storage,
                                                       njobs=njobs):

            # In case of simple rays, use the already loaded dataset: self.ds, 
            # otherwise, load dataset for segment.
            if self.ds is None:
                ds = load(my_segment['filename'], **self.load_kwargs)
            else:
                ds = self.ds

            my_segment['unique_identifier'] = ds.unique_identifier
            if redshift is not None:
                if ds.cosmological_simulation and redshift != ds.current_redshift:
                    mylog.warn("Generating light ray with different redshift than " +
                               "the dataset itself.")
                my_segment["redshift"] = redshift

            if setup_function is not None:
                setup_function(ds)

            if not ds.cosmological_simulation:
                next_redshift = my_segment["redshift"]
            elif self.near_redshift == self.far_redshift:
                if isinstance(my_segment["traversal_box_fraction"], YTArray) and \
                  not my_segment["traversal_box_fraction"].units.is_dimensionless:
                    segment_length = \
                      my_segment["traversal_box_fraction"].in_units("Mpccm / h")
                else:
                    segment_length = my_segment["traversal_box_fraction"] * \
                      ds.domain_width[0].in_units("Mpccm / h")
                next_redshift = my_segment["redshift"] - \
                  self._deltaz_forward(my_segment["redshift"],
                                       segment_length)
            elif my_segment.get("next", None) is None:
                next_redshift = self.near_redshift
            else:
                next_redshift = my_segment['next']['redshift']

            # Make sure start, end, left, right
            # are using the dataset's unit system.
            my_start = ds.arr(my_segment['start'])
            my_end   = ds.arr(my_segment['end'])
            my_left  = ds.arr(left_edge)
            my_right = ds.arr(right_edge)
            mylog.info("Getting segment at z = %s: %s to %s." %
                       (my_segment['redshift'], my_start, my_end))

            # Break periodic ray into non-periodic segments.
            sub_segments = periodic_ray(my_start, my_end,
                                        left=my_left, right=my_right)

            # Prepare data structure for subsegment.
            sub_data = {}
            sub_data['segment_redshift'] = my_segment['redshift']
            for field in all_fields:
                sub_data[field] = []

            # Get data for all subsegments in segment.
            for sub_segment in sub_segments:
                mylog.info("Getting subsegment: %s to %s." %
                           (list(sub_segment[0]), list(sub_segment[1])))
                sub_ray = ds.ray(sub_segment[0], sub_segment[1])
                for key, val in field_parameters.items():
                    sub_ray.set_field_parameter(key, val)
                asort = np.argsort(sub_ray["t"])
                sub_data['dl'].extend(sub_ray['dts'][asort] *
                                      vector_length(sub_ray.start_point,
                                                    sub_ray.end_point))

                for field in data_fields:
                    sub_data[field].extend(sub_ray[field][asort])

                if use_peculiar_velocity:
                    line_of_sight = sub_segment[0] - sub_segment[1]
                    line_of_sight /= ((line_of_sight**2).sum())**0.5
                    sub_vel = ds.arr([sub_ray['velocity_x'],
                                      sub_ray['velocity_y'],
                                      sub_ray['velocity_z']])
                    # Line of sight velocity = vel_los
                    sub_vel_los = (np.rollaxis(sub_vel, 1) * \
                                   line_of_sight).sum(axis=1)
                    sub_data['velocity_los'].extend(sub_vel_los[asort])

                    # doppler redshift:
                    # See https://en.wikipedia.org/wiki/Redshift and 
                    # Peebles eqns: 5.48, 5.49

                    # 1 + redshift_dopp = (1 + v*cos(theta)/c) / 
                    # sqrt(1 - v**2/c**2)

                    # where v is the peculiar velocity (ie physical velocity
                    # without the hubble flow, but no hubble flow in sim, so
                    # just the physical velocity).

                    # the bulk of the doppler redshift is from line of sight 
                    # motion, but there is a small amount from time dilation 
                    # of transverse motion, hence the inclusion of theta (the 
                    # angle between line of sight and the velocity). 
                    # theta is the angle between the ray vector (i.e. line of 
                    # sight) and the velocity vectors: a dot b = ab cos(theta)

                    sub_vel_mag = sub_ray['velocity_magnitude']
                    cos_theta = line_of_sight.dot(sub_vel) / sub_vel_mag
                    # Protect against stituations where velocity mag is exactly
                    # zero, in which case zero / zero = NaN.
                    cos_theta = np.nan_to_num(cos_theta)
                    redshift_dopp = \
                        (1 + sub_vel_mag * cos_theta / speed_of_light_cgs) / \
                         np.sqrt(1 - sub_vel_mag**2 / speed_of_light_cgs**2) - 1
                    sub_data['redshift_dopp'].extend(redshift_dopp[asort])
                    del sub_vel, sub_vel_los, sub_vel_mag, cos_theta, \
                        redshift_dopp

                sub_ray.clear_data()
                del sub_ray, asort

            for key in sub_data:
                sub_data[key] = ds.arr(sub_data[key]).in_cgs()

            # Get redshift for each lixel.  Assume linear relation between l 
            # and z.
            sub_data['dredshift'] = (my_segment['redshift'] - next_redshift) * \
                (sub_data['dl'] / vector_length(my_start, my_end).in_cgs())
            sub_data['redshift'] = my_segment['redshift'] - \
              sub_data['dredshift'].cumsum() + sub_data['dredshift']

            # When using the peculiar velocity, create effective redshift 
            # (redshift_eff) field combining cosmological redshift and 
            # doppler redshift.
            
            # then to add cosmological redshift and doppler redshifts, follow
            # eqn 3.75 in Peacock's Cosmological Physics:
            # 1 + z_eff = (1 + z_cosmo) * (1 + z_doppler)

            if use_peculiar_velocity:
               sub_data['redshift_eff'] = ((1 + sub_data['redshift_dopp']) * \
                                            (1 + sub_data['redshift'])) - 1

            # Remove empty lixels.
            sub_dl_nonzero = sub_data['dl'].nonzero()
            for field in all_fields:
                sub_data[field] = sub_data[field][sub_dl_nonzero]
            del sub_dl_nonzero

            # Add to storage.
            my_storage.result = sub_data

            del ds

        # Reconstruct ray data from parallel_objects storage.
        all_data = [my_data for my_data in all_ray_storage.values()]
        # This is now a list of segments where each one is a dictionary
        # with all the fields.
        all_data.sort(key=lambda a:a['segment_redshift'], reverse=True)
        # Flatten the list into a single dictionary containing fields
        # for the whole ray.
        all_data = _flatten_dict_list(all_data, exceptions=['segment_redshift'])
        self._data = all_data

        if data_filename is not None:
            self._write_light_ray(data_filename, all_data)
            ray_ds = load(data_filename)
            return ray_ds
        else:
            return None

    def __getitem__(self, field):
        return self._data[field]

    @parallel_root_only
    def _write_light_ray(self, filename, data):
        """
        _write_light_ray(filename, data)

        Write light ray data to hdf5 file.
        """
        if self.simulation_type is None:
            ds = self.ds
        else:
            ds = {}
            ds["dimensionality"] = self.simulation.dimensionality
            ds["domain_left_edge"] = self.simulation.domain_left_edge
            ds["domain_right_edge"] = self.simulation.domain_right_edge
            ds["cosmological_simulation"] = self.simulation.cosmological_simulation
            ds["periodicity"] = (True, True, True)
            ds["current_redshift"] = self.near_redshift
            for attr in ["omega_lambda", "omega_matter", "hubble_constant"]:
                ds[attr] = getattr(self.cosmology, attr)
            ds["current_time"] = \
              self.cosmology.t_from_z(ds["current_redshift"])
            if isinstance(ds["hubble_constant"], YTArray):
                ds["hubble_constant"] = \
                  ds["hubble_constant"].to("100*km/(Mpc*s)").d
        extra_attrs = {"data_type": "yt_light_ray"}

        # save the light ray solution
        if len(self.light_ray_solution) > 0:
            # Convert everything to base unit system now to avoid
            # problems with different units for each ds.
            for s in self.light_ray_solution:
                for f in s:
                    if isinstance(s[f], YTArray):
                        s[f].convert_to_base()
            for key in self.light_ray_solution[0]:
                if key in ["next", "previous", "index"]:
                    continue
                lrsa = [sol[key] for sol in self.light_ray_solution]
                if isinstance(lrsa[-1], YTArray):
                    to_arr = YTArray
                else:
                    to_arr = np.array
                extra_attrs["light_ray_solution_%s" % key] = to_arr(lrsa)

        field_types = dict([(field, "grid") for field in data.keys()])

        # Only return LightRay elements with non-zero density
        if 'temperature' in data: f = 'temperature'
        if ('gas', 'temperature') in data: f = ('gas', 'temperature')
        if 'temperature' in data or ('gas', 'temperature') in data:
            mask = data[f] > 0
            if not np.any(mask):
                raise RuntimeError(
                    "No zones along light ray with nonzero %s. "
                    "Please modify your light ray trajectory." % (f,))
            for key in data.keys():
                data[key] = data[key][mask]
        save_as_dataset(ds, filename, data, field_types=field_types,
                        extra_attrs=extra_attrs)

    @parallel_root_only
    def _write_light_ray_solution(self, filename, extra_info=None):
        """
        _write_light_ray_solution(filename, extra_info=None)

        Write light ray solution to a file.
        """

        mylog.info("Writing light ray solution to %s." % filename)
        f = open(filename, 'w')
        if extra_info is not None:
            for par, val in extra_info.items():
                f.write("%s = %s\n" % (par, val))
        f.write("\nSegment Redshift dl/box    Start x       y             " + \
                "z             End x         y             z            Dataset\n")
        for q, my_segment in enumerate(self.light_ray_solution):
            f.write("%04d    %.6f %.6f % .10f % .10f % .10f % .10f % .10f % .10f %s\n" % \
                    (q, my_segment['redshift'], my_segment['traversal_box_fraction'],
                     my_segment['start'][0], my_segment['start'][1], my_segment['start'][2],
                     my_segment['end'][0], my_segment['end'][1], my_segment['end'][2],
                     my_segment['filename']))
        f.close()

def _flatten_dict_list(data, exceptions=None):
    """
    _flatten_dict_list(data, exceptions=None)

    Flatten the list of dicts into one dict.
    """

    if exceptions is None: exceptions = []
    new_data = {}
    for datum in data:
        for field in [field for field in datum.keys()
                      if field not in exceptions]:
            if field not in new_data:
                new_data[field] = []
            new_data[field].extend(datum[field])
    for field in new_data:
        new_data[field] = YTArray(new_data[field])
    return new_data

def vector_length(start, end):
    """
    vector_length(start, end)

    Calculate vector length.
    """

    return np.sqrt(np.power((end - start), 2).sum())

def periodic_adjust(p, left=None, right=None):
    """
    Return the point p adjusted for periodic boundaries.

    """
    if isinstance(p, YTArray):
        p.convert_to_units("unitary")
    if left is None:
        left = np.zeros_like(p)
    if right is None:
        right = np.ones_like(p)

    w = right - left
    p -= left
    return np.mod(p, w)

def periodic_distance(coord1, coord2):
    """
    periodic_distance(coord1, coord2)

    Calculate length of shortest vector between to points in periodic domain.
    """
    dif = coord1 - coord2

    dim = np.ones(coord1.shape,dtype=int)
    def periodic_bind(num):
        pos = np.abs(num % dim)
        neg = np.abs(num % -dim)
        return np.min([pos,neg],axis=0)

    dif = periodic_bind(dif)
    return np.sqrt((dif * dif).sum(axis=-1))

def periodic_ray(start, end, left=None, right=None):
    """
    periodic_ray(start, end, left=None, right=None)

    Break up periodic ray into non-periodic segments.
    Accepts start and end points of periodic ray as YTArrays.
    Accepts optional left and right edges of periodic volume as YTArrays.
    Returns a list of lists of coordinates, where each element of the
    top-most list is a 2-list of start coords and end coords of the
    non-periodic ray:

    [[[x0start,y0start,z0start], [x0end, y0end, z0end]],
     [[x1start,y1start,z1start], [x1end, y1end, z1end]],
     ...,]

    """

    if left is None:
        left = np.zeros(start.shape)
    if right is None:
        right = np.ones(start.shape)
    dim = right - left

    vector = end - start
    wall = np.zeros_like(start)
    close = np.zeros(start.shape, dtype=object)

    left_bound = vector < 0
    right_bound = vector > 0
    no_bound = vector == 0.0
    bound = vector != 0.0

    wall[left_bound] = left[left_bound]
    close[left_bound] = np.max
    wall[right_bound] = right[right_bound]
    close[right_bound] = np.min
    wall[no_bound] = np.inf
    close[no_bound] = np.min

    segments = []
    this_start = start.copy()
    this_end = end.copy()
    t = 0.0
    tolerance = 1e-6
    while t < 1.0 - tolerance:
        hit_left = (this_start <= left) & (vector < 0)
        if (hit_left).any():
            this_start[hit_left] += dim[hit_left]
            this_end[hit_left] += dim[hit_left]
        hit_right = (this_start >= right) & (vector > 0)
        if (hit_right).any():
            this_start[hit_right] -= dim[hit_right]
            this_end[hit_right] -= dim[hit_right]

        nearest = vector.unit_array * \
          np.array([close[q]([this_end[q], wall[q]]) \
                    for q in range(start.size)])
        dt = ((nearest - this_start) / vector)[bound].min()
        now = this_start + vector * dt
        close_enough = np.abs(now - nearest) / np.abs(vector.max()) < 1e-10
        now[close_enough] = nearest[close_enough]
        segments.append([this_start.copy(), now.copy()])
        this_start = now.copy()
        t += dt

    return segments

def non_periodic_ray(ds, left_edge, right_edge, ray_length, max_iter=5000,
                     min_level=None, my_random=None):

    max_length = vector_length(left_edge, right_edge)
    if ray_length > max_length:
        raise RuntimeError(
            ("The maximum segment length in the region %s to %s is %s, " +
             "but the ray length requested is %s.  Decrease ray length.") %
             (left_edge, right_edge, max_length, ray_length))

    if my_random is None:
        my_random = np.random.RandomState()
    i = 0
    while True:
        start = my_random.random_sample(3) * \
          (right_edge - left_edge) + left_edge
        theta = np.pi * my_random.random_sample()
        phi = 2 * np.pi * my_random.random_sample()
        end = start + ray_length * \
          np.array([np.cos(phi) * np.sin(theta),
                    np.sin(phi) * np.sin(theta),
                    np.cos(theta)])
        i += 1
        test_ray = ds.ray(start, end)
        if (end >= left_edge).all() and (end <= right_edge).all() and \
          (min_level is None or min_level <= 0 or
           (test_ray["grid_level"] >= min_level).all()):
            mylog.info("Found ray after %d attempts." % i)
            del test_ray
            return start, end
        del test_ray
        if i > max_iter:
            raise RuntimeError(
                ("Failed to create segment in %d attempts.  " +
                 "Decreasing ray length is recommended") % i)
