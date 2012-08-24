"""
LightRay class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2012 Britton Smith.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import copy
import h5py
import numpy as na

from yt.funcs import *

from yt.analysis_modules.cosmological_observation.cosmology_splice import \
     CosmologySplice
from yt.analysis_modules.halo_profiler.multi_halo_profiler import \
     HaloProfiler
from yt.convenience import load
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    only_on_root, \
    parallel_objects, \
    parallel_root_only

class LightRay(CosmologySplice):
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0,
                 time_data=True, redshift_data=True,
                 find_outputs=False):
        """
        Create a LightRay object.  A light ray is much like a light cone,
        in that it stacks together multiple datasets in order to extend a
        redshift interval.  Unlike a light cone, which does randomly
        oriented projections for each dataset, a light ray consists of
        randomly oriented single rays.  The purpose of these is to create
        synthetic QSO lines of sight.

        Once the LightRay object is set up, use LightRay.make_light_ray to
        begin making rays.  Different randomizations can be created with a
        single object by providing different random seeds to make_light_ray.

        Parameters
        ----------
        parameter_filename : string
            The simulation parameter file.
        simulation_type : string
            The simulation type.
        near_redshift : float
            The near (lowest) redshift for the light ray.
        far_redshift : float
            The far (highest) redshift for the light ray.
        use_minimum_datasets : bool
            If True, the minimum number of datasets is used to connect the
            initial and final redshift.  If false, the light ray solution
            will contain as many entries as possible within the redshift
            interval.
            Default: True.
        deltaz_min : float
            Specifies the minimum :math:`\Delta z` between consecutive
            datasets in the returned list.
            Default: 0.0.
        minimum_coherent_box_fraction : float
            Used with use_minimum_datasets set to False, this parameter
            specifies the fraction of the total box size to be traversed
            before rerandomizing the projection axis and center.  This
            was invented to allow light rays with thin slices to sample
            coherent large scale structure, but in practice does not work
            so well.  Try setting this parameter to 1 and see what happens.
            Default: 0.0.
        time_data : bool
            Whether or not to include time outputs when gathering
            datasets for time series.
            Default: True.
        redshift_data : bool
            Whether or not to include redshift outputs when gathering
            datasets for time series.
            Default: True.
        find_outputs : bool
            Whether or not to search for parameter files in the current 
            directory.
            Default: False.

        """

        self.near_redshift = near_redshift
        self.far_redshift = far_redshift
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction

        self.light_ray_solution = []
        self._data = {}

        # Get list of datasets for light ray solution.
        CosmologySplice.__init__(self, parameter_filename, simulation_type,
                                 find_outputs=find_outputs)
        self.light_ray_solution = \
          self.create_cosmology_splice(self.near_redshift, self.far_redshift,
                                       minimal=self.use_minimum_datasets,
                                       deltaz_min=self.deltaz_min,
                                       time_data=time_data,
                                       redshift_data=redshift_data)

    def _calculate_light_ray_solution(self, seed=None, filename=None):
        "Create list of datasets to be added together to make the light ray."

        # Calculate dataset sizes, and get random dataset axes and centers.
        na.random.seed(seed)

        # For box coherence, keep track of effective depth travelled.
        box_fraction_used = 0.0

        for q in range(len(self.light_ray_solution)):
            if (q == len(self.light_ray_solution) - 1):
                z_next = self.near_redshift
            else:
                z_next = self.light_ray_solution[q+1]['redshift']

            # Calculate fraction of box required for a depth of delta z
            self.light_ray_solution[q]['traversal_box_fraction'] = \
                self.cosmology.ComovingRadialDistance(\
                z_next, self.light_ray_solution[q]['redshift']) * \
                self.simulation.hubble_constant / \
                self.simulation.box_size

            # Simple error check to make sure more than 100% of box depth
            # is never required.
            if (self.light_ray_solution[q]['traversal_box_fraction'] > 1.0):
                mylog.error("Warning: box fraction required to go from z = %f to %f is %f" %
                            (self.light_ray_solution[q]['redshift'], z_next,
                             self.light_ray_solution[q]['traversal_box_fraction']))
                mylog.error("Full box delta z is %f, but it is %f to the next data dump." %
                            (self.light_ray_solution[q]['deltazMax'],
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
                self.light_ray_solution[q]['start'] = na.random.random(3)
                theta = na.pi * na.random.random()
                phi = 2 * na.pi * na.random.random()
                box_fraction_used = 0.0
            else:
                # Use end point of previous segment and same theta and phi.
                self.light_ray_solution[q]['start'] = \
                  self.light_ray_solution[q-1]['end'][:]

            self.light_ray_solution[q]['end'] = \
              self.light_ray_solution[q]['start'] + \
                self.light_ray_solution[q]['traversal_box_fraction'] * \
                na.array([na.cos(phi) * na.sin(theta),
                          na.sin(phi) * na.sin(theta),
                          na.cos(theta)])
            box_fraction_used += \
              self.light_ray_solution[q]['traversal_box_fraction']

        if filename is not None:
            self._write_light_ray_solution(filename,
                extra_info={'parameter_filename':self.parameter_filename,
                            'random_seed':seed,
                            'far_redshift':self.far_redshift,
                            'near_redshift':self.near_redshift})

    def make_light_ray(self, seed=None, fields=None,
                       solution_filename=None, data_filename=None,
                       get_los_velocity=False,
                       get_nearest_halo=False,
                       nearest_halo_fields=None,
                       halo_profiler_parameters=None,
                       njobs=1, dynamic=False):
        """
        Create a light ray and get field values for each lixel.  A light
        ray consists of a list of field values for cells intersected by
        the ray and the path length of the ray through those cells.
        Light ray data can be written out to an hdf5 file.

        Parameters
        ----------
        seed : int
            Seed for the random number generator.
            Default: None.
        fields : list
            A list of fields for which to get data.
            Default: None.
        solution_filename : string
            Path to a text file where the trajectories of each
            subray is written out.
            Default: None.
        data_filename : string
            Path to output file for ray data.
            Default: None.
        get_los_velocity : bool
            If True, the line of sight velocity is calculated for
            each point in the ray.
            Default: False.
        get_nearest_halo : bool
            If True, the HaloProfiler will be used to calculate the
            distance and mass of the nearest halo for each point in the
            ray.  This option requires additional information to be
            included.  See below for an example.
            Default: False.
        nearest_halo_fields : list
            A list of fields to be calculated for the halos nearest to
            every lixel in the ray.
            Default: None.
        halo_profiler_parameters: dict
            A dictionary of parameters to be passed to the HaloProfiler
            to create the appropriate data used to get properties for
            the nearest halos.
            Default: None.
        njobs : int
            The number of parallel jobs over which the slices for the
            halo mask will be split.  Choose -1 for one processor per
            individual slice and 1 to have all processors work together
            on each projection.
            Default: 1
        dynamic : bool
            If True, use dynamic load balancing to create the projections.
            Default: False.

        Getting the Nearest Galaxies
        ----------------------------
        The light ray tool will use the HaloProfiler to calculate the
        distance and mass of the nearest halo to that pixel.  In order
        to do this, a dictionary called halo_profiler_parameters is used
        to pass instructions to the HaloProfiler.  This dictionary has
        three additional keywords.

        halo_profiler_kwargs : dict
            A dictionary of standard HaloProfiler keyword
            arguments and values to be given to the HaloProfiler.

        halo_profiler_actions : list
            A list of actions to be performed by the HaloProfiler.
            Each item in the list should be a dictionary with the following
            entries: "function", "args", and "kwargs", for the function to
            be performed, the arguments supplied to that function, and the
            keyword arguments.

        halo_list : string
            'all' to use the full halo list, or 'filtered' to use the
            filtered halo list created after calling make_profiles.

        Examples
        --------

        from yt.mods import *
        from yt.analysis_modules.halo_profiler.api import *
        from yt.analysis_modules.cosmological_analysis.light_ray.api import LightRay

        halo_profiler_kwargs = {'halo_list_file': 'HopAnalysis.out'}

        halo_profiler_actions = []
        # Add a virial filter.
        halo_profiler_actions.append({'function': add_halo_filter,
                                      'args': VirialFilter,
                                      'kwargs': {'overdensity_field': 'ActualOverdensity',
                                                 'virial_overdensity': 200,
                                                 'virial_filters': \
                                                     [['TotalMassMsun','>=','1e14']],
                                                 'virial_quantities': \
                                                     ['TotalMassMsun','RadiusMpc']}})
        # Make the profiles.
        halo_profiler_actions.append({'function': make_profiles,
                                      'args': None,
                                      'kwargs': {'filename': 'VirializedHalos.out'}})

        halo_list = 'filtered'

        halo_profiler_parameters = dict(halo_profiler_kwargs=halo_profiler_kwargs,
                                        halo_profiler_actions=halo_profiler_actions,
                                        halo_list=halo_list)

        my_ray = LightRay('simulation.par', 'Enzo', 0., 0.1,
                          use_minimum_datasets=True,
                          time_data=False)

        my_ray.make_light_ray(seed=12345,
                              solution_filename='solution.txt',
                              data_filename='my_ray.h5',
                              fields=['Temperature', 'Density'],
                              get_nearest_halo=True,
                              nearest_halo_fields=['TotalMassMsun_100',
                                                   'RadiusMpc_100'],
                              halo_profiler_parameters=halo_profiler_parameters,
                              get_los_velocity=True)

        """

        if halo_profiler_parameters is None:
            halo_profiler_parameters = {}
        if nearest_halo_fields is None:
            nearest_halo_fields = []

        # Calculate solution.
        self._calculate_light_ray_solution(seed=seed, filename=solution_filename)

        # Initialize data structures.
        self._data = {}
        if fields is None: fields = []
        all_fields = [field for field in fields]
        all_fields.extend(['dl', 'dredshift', 'redshift'])
        if get_nearest_halo:
            all_fields.extend(['x', 'y', 'z', 'nearest_halo'])
            all_fields.extend(['nearest_halo_%s' % field \
                               for field in nearest_halo_fields])
            fields.extend(['x', 'y', 'z'])
        if get_los_velocity:
            all_fields.extend(['x-velocity', 'y-velocity',
                               'z-velocity', 'los_velocity'])
            fields.extend(['x-velocity', 'y-velocity', 'z-velocity'])

        all_ray_storage = {}
        for my_storage, my_segment in parallel_objects(self.light_ray_solution,
                                                       storage=all_ray_storage,
                                                       njobs=njobs, dynamic=dynamic):
            mylog.info("Creating ray segment at z = %f." %
                       my_segment['redshift'])
            if my_segment['next'] is None:
                next_redshift = self.near_redshift
            else:
                next_redshift = my_segment['next']['redshift']

            mylog.info("Getting segment at z = %s: %s to %s." %
                       (my_segment['redshift'], my_segment['start'],
                        my_segment['end']))

            if get_nearest_halo:
                halo_list = self._get_halo_list(my_segment['filename'],
                                                **halo_profiler_parameters)

            # Load dataset for segment.
            pf = load(my_segment['filename'])

            # Break periodic ray into non-periodic segments.
            sub_segments = periodic_ray(my_segment['start'], my_segment['end'])

            # Prepare data structure for subsegment.
            sub_data = {}
            sub_data['segment_redshift'] = my_segment['redshift']
            for field in all_fields:
                sub_data[field] = na.array([])

            # Get data for all subsegments in segment.
            for sub_segment in sub_segments:
                mylog.info("Getting subsegment: %s to %s." %
                           (list(sub_segment[0]), list(sub_segment[1])))
                sub_ray = pf.h.ray(sub_segment[0], sub_segment[1])
                sub_data['dl'] = na.concatenate([sub_data['dl'],
                                                 (sub_ray['dts'] *
                                                  vector_length(sub_segment[0],
                                                                sub_segment[1]))])
                for field in fields:
                    sub_data[field] = na.concatenate([sub_data[field],
                                                      (sub_ray[field])])

                if get_los_velocity:
                    line_of_sight = sub_segment[1] - sub_segment[0]
                    line_of_sight /= ((line_of_sight**2).sum())**0.5
                    sub_vel = na.array([sub_ray['x-velocity'],
                                        sub_ray['y-velocity'],
                                        sub_ray['z-velocity']])
                    sub_data['los_velocity'] = \
                      na.concatenate([sub_data['los_velocity'],
                                      (na.rollaxis(sub_vel, 1) *
                                       line_of_sight).sum(axis=1)])
                    del sub_vel

                sub_ray.clear_data()
                del sub_ray

            # Get redshift for each lixel.  Assume linear relation between l and z.
            sub_data['dredshift'] = (my_segment['redshift'] - next_redshift) * \
                (sub_data['dl'] / vector_length(my_segment['start'], my_segment['end']))
            sub_data['redshift'] = my_segment['redshift'] \
                - sub_data['dredshift'].cumsum() + sub_data['dredshift']

            # Calculate distance to nearest object on halo list for each lixel.
            if get_nearest_halo:
                sub_data.update(self._get_nearest_halo_properties(sub_data, halo_list,
                                fields=nearest_halo_fields))
                sub_data['nearest_halo'] *= pf.units['mpccm']

            # Remove empty lixels.
            sub_dl_nonzero = sub_data['dl'].nonzero()
            for field in all_fields:
                sub_data[field] = sub_data[field][sub_dl_nonzero]
            del sub_dl_nonzero

            # Convert dl to cm.
            sub_data['dl'] *= pf.units['cm']

            # Add to storage.
            my_storage.result = sub_data

            pf.h.clear_all_data()
            del pf

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

    def _get_halo_list(self, dataset, halo_profiler_kwargs=None,
                       halo_profiler_actions=None, halo_list='all'):
        "Load a list of halos for the dataset."

        if halo_profiler_kwargs is None: halo_profiler_kwargs = {}
        if halo_profiler_actions is None: halo_profiler_actions = []

        hp = HaloProfiler(dataset, **halo_profiler_kwargs)
        for action in halo_profiler_actions:
            if not action.has_key('args'): action['args'] = ()
            if not action.has_key('kwargs'): action['kwargs'] = {}
            action['function'](hp, *action['args'], **action['kwargs'])

        if halo_list == 'all':
            return_list = copy.deepcopy(hp.all_halos)
        elif halo_list == 'filtered':
            return_list = copy.deepcopy(hp.filtered_halos)
        else:
            mylog.error("Keyword, halo_list, must be either 'all' or 'filtered'.")
            return_list = None

        del hp
        return return_list

    def _get_nearest_halo_properties(self, data, halo_list, fields=None):
        """
        Calculate distance to nearest object in halo list for each lixel in data.
        Return list of distances and masses of nearest objects.
        """

        if fields is None: fields = []

        # Create position array from halo list.
        halo_centers = na.array(map(lambda halo: halo['center'], halo_list))
        halo_field_values = dict([(field, na.array(map(lambda halo: halo[field],
                                                       halo_list))) \
                                  for field in fields])

        nearest_distance = na.zeros(data['x'].shape)
        field_data = dict([(field, na.zeros(data['x'].shape)) \
                           for field in fields])
        for index in xrange(nearest_distance.size):
            nearest = na.argmin(periodic_distance(na.array([data['x'][index],
                                                            data['y'][index],
                                                            data['z'][index]]),
                                                  halo_centers))
            nearest_distance[index] = periodic_distance(na.array([data['x'][index],
                                                                  data['y'][index],
                                                                  data['z'][index]]),
                                                        halo_centers[nearest])
            for field in fields:
                field_data[field][index] = halo_field_values[field][nearest]

        return_data = {'nearest_halo': nearest_distance}
        for field in fields:
            return_data['nearest_halo_%s' % field] = field_data[field]
        return return_data

    @parallel_root_only
    def _write_light_ray(self, filename, data):
        "Write light ray data to hdf5 file."

        mylog.info("Saving light ray data to %s." % filename)
        output = h5py.File(filename, 'w')
        for field in data.keys():
            output.create_dataset(field, data=data[field])
        output.close()

    @parallel_root_only
    def _write_light_ray_solution(self, filename, extra_info=None):
        "Write light ray solution to a file."

        mylog.info("Writing light ray solution to %s." % filename)
        f = open(filename, 'w')
        if extra_info is not None:
            for par, val in extra_info.items():
                f.write("%s = %s\n" % (par, val))
        f.write("\nSegment Redshift dl/box    Start x       y             z             End x         y             z            Dataset\n")
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
            if field in new_data:
                new_data[field] = na.concatenate([new_data[field], datum[field]])
            else:
                new_data[field] = na.copy(datum[field])
    return new_data

def vector_length(start, end):
    "Calculate vector length."

    return na.sqrt(na.power((end - start), 2).sum())

def periodic_distance(coord1, coord2):
    "Calculate length of shortest vector between to points in periodic domain."
    dif = coord1 - coord2

    dim = na.ones(coord1.shape,dtype=int)
    def periodic_bind(num):
        pos = na.abs(num % dim)
        neg = na.abs(num % -dim)
        return na.min([pos,neg],axis=0)

    dif = periodic_bind(dif)
    return na.sqrt((dif * dif).sum(axis=-1))

def periodic_ray(start, end, left=None, right=None):
    "Break up periodic ray into non-periodic segments."

    if left is None:
        left = na.zeros(start.shape)
    if right is None:
        right = na.ones(start.shape)
    dim = right - left

    vector = end - start
    wall = na.zeros(start.shape)
    close = na.zeros(start.shape, dtype=object)

    left_bound = vector < 0
    right_bound = vector > 0
    no_bound = vector == 0.0
    bound = vector != 0.0

    wall[left_bound] = left[left_bound]
    close[left_bound] = na.max
    wall[right_bound] = right[right_bound]
    close[right_bound] = na.min
    wall[no_bound] = na.inf
    close[no_bound] = na.min

    segments = []
    this_start = na.copy(start)
    this_end = na.copy(end)
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

        nearest = na.array([close[q]([this_end[q], wall[q]]) \
                                for q in range(start.size)])
        dt = ((nearest - this_start) / vector)[bound].min()
        now = this_start + vector * dt
        close_enough = na.abs(now - nearest) < 1e-10
        now[close_enough] = nearest[close_enough]
        segments.append([na.copy(this_start), na.copy(now)])
        this_start = na.copy(now)
        t += dt

    return segments
