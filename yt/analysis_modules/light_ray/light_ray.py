"""
LightRay class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

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

try:
    from mpi4py import MPI
    parallel_light_ray = True
    my_rank = MPI.COMM_WORLD.rank
    my_size = MPI.COMM_WORLD.size
except ImportError:
    parallel_light_ray = False
    my_rank = 0
    my_size = 1

import copy
import h5py
import numpy as na

from yt.funcs import *

from yt.analysis_modules.simulation_handler.enzo_simulation import \
    EnzoSimulation
from yt.analysis_modules.halo_profiler.multi_halo_profiler import \
    HaloProfiler
from yt.convenience import load

class LightRay(EnzoSimulation):
    def __init__(self, enzo_parameter_file, final_redshift, initial_redshift, 
                 deltaz_min=0.0, use_minimum_datasets=True, 
                 minimum_coherent_box_fraction=0.0, **kwargs):
        """
        Create a LightRay object.  A light ray is much like a light cone, in that 
        it stacks together multiple datasets in order to extend a redshift interval.  
        Unlike a light cone, which does randomly oriented projections for each dataset, 
        a light ray consists of randomly oriented single rays.  The purpose of these 
        is to create synthetic QSO lines of sight.

        Once the LightRay object is set up, use LightRay.make_light_ray to begin making 
        rays.  Different randomizations can be created with a single object by providing 
        different random seeds to make_light_ray.

        :param enzo_parameter_file (string): path to simulation parameter file.
        :param final_redshift (float): lower bound of the ray redshift interval.
        :param initial_redshift (float): upper bound of the ray redshift interval.
        :param deltaz_min (float): minimum delta z between consecutive datasets.
        Default: 0.0.
        :param use_minimum_datasets (bool): if True, the minimum number of datasets is 
        used to connect the initial and final redshift.  If false, the light ray 
        solution will contain as many entries as possible within the redshift interval.  
        Default: True.
        minimum_coherent_box_fraction (float): used with use_minimum_datasets set to False, 
        this parameter specifies the fraction of the total box size to be traversed before 
        rerandomizing the projection axis and center.  This was invented to allow light cones 
        with thin slices to sample coherent large scale structure, but in practice does not 
        work so well.  It is not very clear what this will do to a light ray.  Default: 0.0.
        """

        EnzoSimulation.__init__(self, enzo_parameter_file, 
                                initial_redshift=initial_redshift,
                                final_redshift=final_redshift, links=True,
                                enzo_parameters={'CosmologyComovingBoxSize':float}, 
                                **kwargs)

        self.deltaz_min = deltaz_min
        self.use_minimum_datasets = use_minimum_datasets
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        self.dtype = '>f8'

        self.light_ray_solution = []
        self._data = {}

        # Don't use box coherence with maximum dataset depths.
        if self.use_minimum_datasets and self.minimum_coherent_box_fraction > 0:
            mylog.info("Setting minimum_coherent_box_fraction to 0 with minimal light ray.")
            self.minimum_coherent_box_fraction = 0

        # Get list of datasets for light ray solution.
        self.light_ray_solution = \
            self.create_cosmology_splice(minimal=self.use_minimum_datasets,
                                         deltaz_min=self.deltaz_min)

    def _calculate_light_ray_solution(self, seed=None, filename=None):
        "Create list of datasets to be added together to make the light ray."

        # Calculate dataset sizes, and get random dataset axes and centers.
        na.random.seed(seed)

        # For box coherence, keep track of effective depth travelled.
        boxFractionUsed = 0.0

        for q in range(len(self.light_ray_solution)):
            if (q == len(self.light_ray_solution) - 1):
                z_next = self.final_redshift
            else:
                z_next = self.light_ray_solution[q+1]['redshift']

            # Calculate fraction of box required for a depth of delta z
            self.light_ray_solution[q]['TraversalBoxFraction'] = \
                self.cosmology.ComovingRadialDistance(\
                z_next, self.light_ray_solution[q]['redshift']) * \
                self.enzoParameters['CosmologyHubbleConstantNow'] / \
                self.enzoParameters['CosmologyComovingBoxSize']

            # Simple error check to make sure more than 100% of box depth 
            # is never required.
            if (self.light_ray_solution[q]['TraversalBoxFraction'] > 1.0):
                mylog.error("Warning: box fraction required to go from z = %f to %f is %f" % 
                            (self.light_ray_solution[q]['redshift'], z_next,
                             self.light_ray_solution[q]['TraversalBoxFraction']))
                mylog.error("Full box delta z is %f, but it is %f to the next data dump." % 
                            (self.light_ray_solution[q]['deltazMax'],
                             self.light_ray_solution[q]['redshift']-z_next))

            # Get dataset axis and center.
            # If using box coherence, only get start point and vector if 
            # enough of the box has been used, 
            # or if boxFractionUsed will be greater than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
                    (boxFractionUsed > self.minimum_coherent_box_fraction) or \
                    (boxFractionUsed + self.light_ray_solution[q]['TraversalBoxFraction'] > 1.0):
                # Random star point
                self.light_ray_solution[q]['start'] = na.random.random(3)
                theta = na.pi * na.random.random()
                phi = 2 * na.pi * na.random.random()
                boxFractionUsed = 0.0
            else:
                # Use end point of previous segment and same theta and phi.
                self.light_ray_solution[q]['start'] = self.light_ray_solution[q-1]['end'][:]

            self.light_ray_solution[q]['end'] = self.light_ray_solution[q]['start'] + \
                self.light_ray_solution[q]['TraversalBoxFraction'] * \
                na.array([na.cos(phi) * na.sin(theta),
                          na.sin(phi) * na.sin(theta),
                          na.cos(theta)])
            boxFractionUsed += self.light_ray_solution[q]['TraversalBoxFraction']

        if filename is not None:
            self._write_light_ray_solution(filename, 
                                           extra_info={'enzo_parameter_file':self.enzo_parameter_file, 
                                                       'RandomSeed':seed,
                                                       'initial_redshift':self.initial_redshift, 
                                                       'final_redshift':self.final_redshift})

    def make_light_ray(self, seed=None, fields=None, 
                       solution_filename=None, data_filename=None,
                       get_nearest_galaxy=False, get_los_velocity=False, 
                       halo_mass_field='TotalMassMsun_200', **kwargs):
        """
        Create a light ray and get field values for each lixel.  A light ray consists of 
        a list of field values for cells intersected by the ray and the path length of 
        the ray through those cells.  Light ray data can be written out to an hdf5 file.

        :param seed (int): seed for the random number generator.  Default: None.
        :param fields (list): a list of fields for which to get data.  Default: None.
        :param solution_filename (string): path to a text file where the trajectories of each 
        subray is written out.  Default: None.
        :param data_filename (string): path to output file for ray data.  Default: None.
        :param get_nearest_galaxy (bool): if True, the HaloProfiler will be used to calculate 
        the distance and mass of the nearest halo for each point in the ray.  This option 
        requires additional information to be included.  See below for an example.  
        Default: False.
        :param get_los_velocity (bool): if True, the line of sight velocity is calculated for 
        each point in the ray.  Default: False.

        GETTING THE NEAREST GALAXIES
        The light ray tool will use the HaloProfiler to calculate the distance and mass 
        of the nearest halo to that pixel.  In order to do this, four additional keyword 
        arguments must be supplied to tell the HaloProfiler what to do.

        :param halo_profiler_kwargs (dict): a dictionary of standard HaloProfiler keyword 
        arguments and values to be given to the HaloProfiler.
               EXAMPLE: halo_profiler_kwargs = {'halo_list_format': {'id':0, 
                                                                     'center':[4, 5, 6]},
                                                                     'TotalMassMsun':1},
                                                'halo_list_file': 'HopAnalysis.out'}

        :param halo_profiler_actions (list): a list of actions to be performed by the 
        HaloProfiler.  Each item in the list should be a dictionary with the following 
        entries: "function", "args", and "kwargs", for the function to be performed, 
        the arguments supplied to that function, and the keyword arguments.
               EXAMPLE: halo_profiler_actions = [{'function': make_profiles,
                                                  'args': None,
                                                  'kwargs': {'filename': 'VirializedHalos.out'}},
                                                 {'function': add_halo_filter,
                                                  'args': VirialFilter,
                                                  'kwargs': {'overdensity_field': 'ActualOverdensity',
                                                             'virial_overdensity': 200,
                                                             'virial_filters': [['TotalMassMsun','>=','1e14']],
                                                             'virial_quantities': ['TotalMassMsun','RadiusMpc']}}]

        :param halo_list (string): 'all' to use the full halo list, or 'filtered' to use 
        the filtered halo list created after calling make_profiles.
               EXAMPLE: halo_list = 'filtered'

        :param halo_mass_field (string): the field from the halo list to use for mass.  
        Default: 'TotalMassMsun_200'.
        """

        # Calculate solution.
        self._calculate_light_ray_solution(seed=seed, filename=solution_filename)

        # Initialize data structures.
        self._data = {}
        data = []
        if fields is None: fields = []
        all_fields = [field for field in fields]
        all_fields.extend(['dl', 'dredshift', 'redshift'])
        if get_nearest_galaxy:
            all_fields.extend(['x', 'y', 'z', 'nearest_galaxy', 'nearest_galaxy_mass'])
            fields.extend(['x', 'y', 'z'])
        if get_los_velocity:
            all_fields.extend(['x-velocity', 'y-velocity', 'z-velocity', 'los_velocity'])
            fields.extend(['x-velocity', 'y-velocity', 'z-velocity'])

        todo = na.arange(my_rank, len(self.light_ray_solution), my_size)
        for index in todo:
            segment = self.light_ray_solution[index]
            mylog.info("Proc %04d: creating ray segment at z = %f." % 
                       (my_rank, segment['redshift']))
            if segment['next'] is None:
                next_redshift = self.final_redshift
            else:
                next_redshift = segment['next']['redshift']

            mylog.info("Getting segment at z = %s: %s to %s." % 
                       (segment['redshift'], segment['start'], segment['end']))

            if get_nearest_galaxy:
                halo_list = self._get_halo_list(segment['filename'], **kwargs)

            # Load dataset for segment.
            pf = load(segment['filename'])

            # Break periodic ray into non-periodic segments.
            sub_segments = periodic_ray(segment['start'], segment['end'])

            # Prepare data structure for subsegment.
            sub_data = {}
            sub_data['segment_redshift'] = segment['redshift']
            for field in all_fields:
                sub_data[field] = na.array([], dtype=self.dtype)

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
                    sub_data['los_velocity'] = na.concatenate([sub_data['los_velocity'], 
                                                               (na.rollaxis(sub_vel, 1) * 
                                                                line_of_sight).sum(axis=1)])
                    del sub_vel

                sub_ray.clear_data()
                del sub_ray

            # Get redshift for each lixel.  Assume linear relation between l and z.
            sub_data['dredshift'] = (segment['redshift'] - next_redshift) * \
                (sub_data['dl'] / vector_length(segment['start'], segment['end']))
            sub_data['redshift'] = segment['redshift'] \
                - sub_data['dredshift'].cumsum() + sub_data['dredshift']

            # Calculate distance to nearest object on halo list for each lixel.
            if get_nearest_galaxy:
                sub_data['nearest_galaxy'], sub_data['nearest_galaxy_mass'] = \
                    self._get_nearest_galaxy_distance(sub_data, halo_list,
                                                      halo_mass_field=halo_mass_field)
                sub_data['nearest_galaxy'] *= pf.units['mpccm']

            # Remove empty lixels.
            sub_dl_nonzero = sub_data['dl'].nonzero()
            for field in all_fields:
                sub_data[field] = sub_data[field][sub_dl_nonzero]
            del sub_dl_nonzero

            # Convert dl to Mpc comoving.
            sub_data['dl'] *= pf.units['cm']

            # Add segment to list.
            data.append(sub_data)

            pf.h.clear_all_data()
            del pf

        if parallel_light_ray:
            MPI.COMM_WORLD.Barrier()
            if my_rank == 0:
                for proc in range(1, my_size):
                    buf = MPI.COMM_WORLD.recv(source=proc, tag=0)
                    data += buf
            else:
                MPI.COMM_WORLD.send(data, dest=0, tag=0)
                del data
            MPI.COMM_WORLD.Barrier()

        if my_rank == 0:
            data.sort(key=lambda a:a['segment_redshift'], reverse=True)
            data = self._flatten_ray_data(data, exceptions=['segment_redshift'])

            if data_filename is not None:
                self._write_light_ray(data_filename, data)

            self._data = data
            return data

    def _flatten_ray_data(self, data, exceptions=None):
        "Flatten the list of dicts into one dict."

        if exceptions is None: exceptions = []
        new_data = {}
        for datum in data:
            for field in [field for field in datum.keys() 
                          if field not in exceptions]:
                if new_data.has_key(field):
                    new_data[field] = na.concatenate([new_data[field], datum[field]])
                else:
                    new_data[field] = na.copy(datum[field])

        return new_data                

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

    def _get_nearest_galaxy_distance(self, data, halo_list, 
                                     halo_mass_field='TotalMassMsun_200'):
        """
        Calculate distance to nearest object in halo list for each lixel in data.
        Return list of distances and masses of nearest objects.
        """

        # Create position array from halo list.
        halo_centers = na.array(map(lambda halo: halo['center'], halo_list))
        halo_mass = na.array(map(lambda halo: halo[halo_mass_field], halo_list))

        nearest_distance = na.zeros(data['x'].shape)
        nearest_mass = na.zeros(data['x'].shape)
        for index in xrange(nearest_distance.size):
            nearest = na.argmin(periodic_distance(na.array([data['x'][index], 
                                                            data['y'][index],
                                                            data['z'][index]]), 
                                                  halo_centers))
            nearest_distance[index] = periodic_distance(na.array([data['x'][index], 
                                                                  data['y'][index],
                                                                  data['z'][index]]), 
                                                        halo_centers[nearest])
            nearest_mass[index] = halo_mass[nearest]

        return (nearest_distance, nearest_mass)

    def _write_light_ray(self, filename, data):
        "Write light ray data to hdf5 file."

        mylog.info("Saving light ray data to %s." % filename)
        output = h5py.File(filename, 'w')
        for field in data.keys():
            output.create_dataset(field, data=data[field], dtype=self.dtype)
        output.close()

    def _write_light_ray_solution(self, filename, extra_info=None):
        "Write light ray solution to a file."

        if my_rank != 0: return

        mylog.info("Writing light ray solution to %s." % filename)
        f = open(filename, 'w')
        if extra_info is not None:
            for par, val in extra_info.items():
                f.write("%s = %s\n" % (par, val))
        f.write("\nSegment Redshift dl/box    Start x       y             z             End x         y             z            Dataset\n")
        for q, segment in enumerate(self.light_ray_solution):
            f.write("%04d    %.6f %.6f % .10f % .10f % .10f % .10f % .10f % .10f %s\n" % \
                        (q, segment['redshift'], segment['TraversalBoxFraction'],
                         segment['start'][0], segment['start'][1], segment['start'][2],
                         segment['end'][0], segment['end'][1], segment['end'][2], 
                         segment['filename']))
        f.close()

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
        nearest = na.array([close[q]([this_end[q], wall[q]]) \
                                for q in range(start.size)])
        dt = ((nearest - this_start) / vector)[bound].min()
        now = this_start + vector * dt
        segments.append([na.copy(this_start), na.copy(now)])
        this_start = na.copy(now)
        hit_left = (this_start <= left) & (vector < 0)
        if (hit_left).any():
            this_start[hit_left] += dim[hit_left]
            this_end[hit_left] += dim[hit_left]
        hit_right = (this_start >= right) & (vector > 0)
        if (hit_right).any():
            this_start[hit_right] -= dim[hit_right]
            this_end[hit_right] -= dim[hit_right]
        t += dt

    return segments
