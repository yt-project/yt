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

import copy
import h5py
import numpy as np

from yt.analysis_modules.cosmological_observation.cosmology_splice import \
    CosmologySplice
from yt.convenience import \
    load
from yt.funcs import \
    mylog
from yt.units.yt_array import \
    YTArray
from yt.utilities.cosmology import \
    Cosmology
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, \
    parallel_root_only

class LightRay(CosmologySplice):
    """
    Create a LightRay object.  A light ray is much like a light cone,
    in that it stacks together multiple datasets in order to extend a
    redshift interval.  Unlike a light cone, which does randomly
    oriented projections for each dataset, a light ray consists of
    randomly oriented single rays.  The purpose of these is to create
    synthetic QSO lines of sight.

    Light rays can also be made from single datasets.
    
    Once the LightRay object is set up, use LightRay.make_light_ray to
    begin making rays.  Different randomizations can be created with a
    single object by providing different random seeds to make_light_ray.

    Parameters
    ----------
    parameter_filename : string
        The path to the simulation parameter file or dataset.
    simulation_type : optional, string
        The simulation type.  If None, the first argument is assumed to 
        refer to a single dataset.
        Default: None
    near_redshift : optional, float
        The near (lowest) redshift for a light ray containing multiple 
        datasets.  Do not use is making a light ray from a single 
        dataset.
        Default: None
    far_redshift : optional, float
        The far (highest) redshift for a light ray containing multiple 
        datasets.  Do not use is making a light ray from a single 
        dataset.
        Default: None
    use_minimum_datasets : optional, bool
        If True, the minimum number of datasets is used to connect the
        initial and final redshift.  If false, the light ray solution
        will contain as many entries as possible within the redshift
        interval.
        Default: True.
    deltaz_min : optional, float
        Specifies the minimum :math:`\Delta z` between consecutive
        datasets in the returned list.
        Default: 0.0.
    minimum_coherent_box_fraction : optional, float
        Used with use_minimum_datasets set to False, this parameter
        specifies the fraction of the total box size to be traversed
        before rerandomizing the projection axis and center.  This
        was invented to allow light rays with thin slices to sample
        coherent large scale structure, but in practice does not work
        so well.  Try setting this parameter to 1 and see what happens.
        Default: 0.0.
    time_data : optional, bool
        Whether or not to include time outputs when gathering
        datasets for time series.
        Default: True.
    redshift_data : optional, bool
        Whether or not to include redshift outputs when gathering
        datasets for time series.
        Default: True.
    find_outputs : optional, bool
        Whether or not to search for datasets in the current 
        directory.
        Default: False.

    """
    def __init__(self, parameter_filename, simulation_type=None,
                 near_redshift=None, far_redshift=None,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0,
                 time_data=True, redshift_data=True,
                 find_outputs=False):

        self.near_redshift = near_redshift
        self.far_redshift = far_redshift
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        self.parameter_filename = parameter_filename

        self.light_ray_solution = []
        self._data = {}

        # Make a light ray from a single, given dataset.        
        if simulation_type is None:
            ds = load(parameter_filename)
            if ds.cosmological_simulation:
                redshift = ds.current_redshift
                self.cosmology = Cosmology(
                    hubble_constant=ds.hubble_constant,
                    omega_matter=ds.omega_matter,
                    omega_lambda=ds.omega_lambda,
                    unit_registry=ds.unit_registry)
            else:
                redshift = 0.
            self.light_ray_solution.append({"filename": parameter_filename,
                                            "redshift": redshift})

        # Make a light ray from a simulation time-series.
        else:
            # Get list of datasets for light ray solution.
            CosmologySplice.__init__(self, parameter_filename, simulation_type,
                                     find_outputs=find_outputs)
            self.light_ray_solution = \
              self.create_cosmology_splice(self.near_redshift, self.far_redshift,
                                           minimal=self.use_minimum_datasets,
                                           deltaz_min=self.deltaz_min,
                                           time_data=time_data,
                                           redshift_data=redshift_data)

    def _calculate_light_ray_solution(self, seed=None, 
                                      start_position=None, end_position=None,
                                      trajectory=None, filename=None):
        "Create list of datasets to be added together to make the light ray."

        # Calculate dataset sizes, and get random dataset axes and centers.
        np.random.seed(seed)

        # If using only one dataset, set start and stop manually.
        if start_position is not None:
            if len(self.light_ray_solution) > 1:
                raise RuntimeError("LightRay Error: cannot specify start_position " + \
                                   "if light ray uses more than one dataset.")
            if not ((end_position is None) ^ (trajectory is None)):
                raise RuntimeError("LightRay Error: must specify either end_position " + \
                                   "or trajectory, but not both.")
            self.light_ray_solution[0]['start'] = np.array(start_position)
            if end_position is not None:
                self.light_ray_solution[0]['end'] = np.array(end_position)
            else:
                # assume trajectory given as r, theta, phi
                if len(trajectory) != 3:
                    raise RuntimeError("LightRay Error: trajectory must have lenght 3.")
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

                # Simple error check to make sure more than 100% of box depth
                # is never required.
                if (self.light_ray_solution[q]['traversal_box_fraction'] > 1.0):
                    mylog.error("Warning: box fraction required to go from z = %f to %f is %f" %
                                (self.light_ray_solution[q]['redshift'], z_next,
                                 self.light_ray_solution[q]['traversal_box_fraction']))
                    mylog.error("Full box delta z is %f, but it is %f to the next data dump." %
                                (self.light_ray_solution[q]['dz_max'],
                                 self.light_ray_solution[q]['redshift']-z_next))

                # Get dataset axis and center.
                # If using box coherence, only get start point and vector if
                # enough of the box has been used,
                # or if box_fraction_used will be greater than 1 after this slice.
                if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
                        (box_fraction_used >
                         self.minimum_coherent_box_fraction) or \
                        (box_fraction_used +
                         self.light_ray_solution[q]['traversal_box_fraction'] > 1.0):
                    # Random start point
                    self.light_ray_solution[q]['start'] = np.random.random(3)
                    theta = np.pi * np.random.random()
                    phi = 2 * np.pi * np.random.random()
                    box_fraction_used = 0.0
                else:
                    # Use end point of previous segment and same theta and phi.
                    self.light_ray_solution[q]['start'] = \
                      self.light_ray_solution[q-1]['end'][:]

                self.light_ray_solution[q]['end'] = \
                  self.light_ray_solution[q]['start'] + \
                    self.light_ray_solution[q]['traversal_box_fraction'] * \
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

    def make_light_ray(self, seed=None,
                       start_position=None, end_position=None,
                       trajectory=None,
                       fields=None, setup_function=None,
                       solution_filename=None, data_filename=None,
                       get_los_velocity=True, redshift=None,
                       njobs=-1):
        """
        Create a light ray and get field values for each lixel.  A light
        ray consists of a list of field values for cells intersected by
        the ray and the path length of the ray through those cells.
        Light ray data can be written out to an hdf5 file.

        Parameters
        ----------
        seed : optional, int
            Seed for the random number generator.
            Default: None.
        start_position : optional, list of floats
            Used only if creating a light ray from a single dataset.
            The coordinates of the starting position of the ray.
            Default: None.
        end_position : optional, list of floats
            Used only if creating a light ray from a single dataset.
            The coordinates of the ending position of the ray.
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
        get_los_velocity : optional, bool
            If True, the line of sight velocity is calculated for
            each point in the ray.
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
        ...                       get_los_velocity=True)

        Make a light ray from a single dataset:

        >>> import yt
        >>> from yt.analysis_modules.cosmological_observation.light_ray.api import \
        ...     LightRay
        >>> my_ray = LightRay("enzo_tiny_cosmology/DD0040/DD0040")
        ...
        >>> my_ray.make_light_ray(start_position=[0., 0., 0.],
        ...                       end_position=[1., 1., 1.],
        ...                       solution_filename="solution.txt",
        ...                       data_filename="my_ray.h5",
        ...                       fields=["temperature", "density"],
        ...                       get_los_velocity=True)
        
        """

        # Calculate solution.
        self._calculate_light_ray_solution(seed=seed, 
                                           start_position=start_position, 
                                           end_position=end_position,
                                           trajectory=trajectory,
                                           filename=solution_filename)

        # Initialize data structures.
        self._data = {}
        if fields is None: fields = []
        data_fields = fields[:]
        all_fields = fields[:]
        all_fields.extend(['dl', 'dredshift', 'redshift'])
        if get_los_velocity:
            all_fields.extend(['x-velocity', 'y-velocity',
                               'z-velocity', 'los_velocity'])
            data_fields.extend(['x-velocity', 'y-velocity', 'z-velocity'])

        all_ray_storage = {}
        for my_storage, my_segment in parallel_objects(self.light_ray_solution,
                                                       storage=all_ray_storage,
                                                       njobs=njobs):

            # Load dataset for segment.
            ds = load(my_segment['filename'])

            if redshift is not None:
                if ds.cosmological_simulation and redshift != ds.current_redshift:
                    mylog.warn("Generating light ray with different redshift than " +
                               "the dataset itself.")
                my_segment["redshift"] = redshift

            if setup_function is not None:
                setup_function(ds)
            
            my_segment["start"] = ds.domain_width * my_segment["start"] + \
                ds.domain_left_edge
            my_segment["end"] = ds.domain_width * my_segment["end"] + \
                ds.domain_left_edge

            if not ds.cosmological_simulation:
                next_redshift = my_segment["redshift"]
            elif self.near_redshift == self.far_redshift:
                next_redshift = my_segment["redshift"] - \
                  self._deltaz_forward(my_segment["redshift"], 
                                       ds.domain_width[0].in_units("Mpccm / h") *
                                       my_segment["traversal_box_fraction"])
            elif my_segment.get("next") is None:
                next_redshift = self.near_redshift
            else:
                next_redshift = my_segment['next']['redshift']

            mylog.info("Getting segment at z = %s: %s to %s." %
                       (my_segment['redshift'], my_segment['start'],
                        my_segment['end']))

            # Break periodic ray into non-periodic segments.
            sub_segments = periodic_ray(my_segment['start'], my_segment['end'],
                                        left=ds.domain_left_edge,
                                        right=ds.domain_right_edge)

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
                asort = np.argsort(sub_ray["t"])
                sub_data['dl'].extend(sub_ray['dts'][asort] *
                                      vector_length(sub_ray.start_point,
                                                    sub_ray.end_point))
                for field in data_fields:
                    sub_data[field].extend(sub_ray[field][asort])

                if get_los_velocity:
                    line_of_sight = sub_segment[1] - sub_segment[0]
                    line_of_sight /= ((line_of_sight**2).sum())**0.5
                    sub_vel = ds.arr([sub_ray['x-velocity'],
                                      sub_ray['y-velocity'],
                                      sub_ray['z-velocity']])
                    sub_data['los_velocity'].extend((np.rollaxis(sub_vel, 1) *
                                                     line_of_sight).sum(axis=1)[asort])
                    del sub_vel

                sub_ray.clear_data()
                del sub_ray, asort

            for key in sub_data:
                if key in "xyz": continue
                sub_data[key] = ds.arr(sub_data[key]).in_cgs()

            # Get redshift for each lixel.  Assume linear relation between l and z.
            sub_data['dredshift'] = (my_segment['redshift'] - next_redshift) * \
                (sub_data['dl'] / vector_length(my_segment['start'], 
                                                my_segment['end']).in_cgs())
            sub_data['redshift'] = my_segment['redshift'] - \
              sub_data['dredshift'].cumsum() + sub_data['dredshift']

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

        if data_filename is not None:
            self._write_light_ray(data_filename, all_data)

        self._data = all_data
        return all_data

    @parallel_root_only
    def _write_light_ray(self, filename, data):
        "Write light ray data to hdf5 file."

        mylog.info("Saving light ray data to %s." % filename)
        output = h5py.File(filename, 'w')
        for field in data.keys():
            output.create_dataset(field, data=data[field])
            output[field].attrs["units"] = str(data[field].units)
        output.close()

    @parallel_root_only
    def _write_light_ray_solution(self, filename, extra_info=None):
        "Write light ray solution to a file."

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
    "Flatten the list of dicts into one dict."

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
    "Calculate vector length."

    return np.sqrt(np.power((end - start), 2).sum())

def periodic_distance(coord1, coord2):
    "Calculate length of shortest vector between to points in periodic domain."
    dif = coord1 - coord2

    dim = np.ones(coord1.shape,dtype=int)
    def periodic_bind(num):
        pos = np.abs(num % dim)
        neg = np.abs(num % -dim)
        return np.min([pos,neg],axis=0)

    dif = periodic_bind(dif)
    return np.sqrt((dif * dif).sum(axis=-1))

def periodic_ray(start, end, left=None, right=None):
    "Break up periodic ray into non-periodic segments."

    if left is None:
        left = np.zeros(start.shape)
    if right is None:
        right = np.ones(start.shape)
    dim = right - left

    vector = end - start
    wall = np.zeros(start.shape)
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
        segments.append([np.copy(this_start), np.copy(now)])
        this_start = now.copy()
        t += dt

    return segments
