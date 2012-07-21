"""
LightCone class and member functions.

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
import os

from yt.funcs import *
from yt.analysis_modules.cosmology.cosmology_splice import \
     CosmologySplice
from yt.convenience import load
from yt.utilities.cosmology import Cosmology
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    only_on_root, \
    parallel_objects, \
    parallel_root_only
from yt.visualization.image_writer import write_image

from .common_n_volume import common_volume
from .halo_mask import \
    _light_cone_halo_mask
from .light_cone_projection import \
     _light_cone_projection

class LightCone(CosmologySplice):
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift,
                 observer_redshift=0.0,
                 field_of_view_in_arcminutes=600.0,
                 image_resolution_in_arcseconds=60.0,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0,
                 set_parameters=None,
                 output_dir='LC', output_prefix='LightCone'):
        """
        Initialize a LightCone object.

        Parameters
        ----------
        near_redshift : float
            The near (lowest) redshift for the light cone.
        far_redshift : float
            The far (highest) redshift for the light cone.
        observer_redshift : float
            The redshift of the observer.
            Default: 0.0.
        field_of_view_in_arcminutes : float
            The field of view of the image in units of arcminutes.
            Default: 600.0.
        image_resolution_in_arcseconds : float
            The size of each image pixel in units of arcseconds.
            Default: 60.0.
        use_minimum_datasets : bool
            If True, the minimum number of datasets is used to connect the initial
            and final redshift.  If false, the light cone solution will contain
            as many entries as possible within the redshift interval.
            Default: True.
        deltaz_min : float
            Specifies the minimum :math:`\Delta z` between consecutive datasets in
            the returned list.
            Default: 0.0.
        minimum_coherent_box_fraction : float
            Used with use_minimum_datasets set to False, this parameter specifies
            the fraction of the total box size to be traversed before rerandomizing
            the projection axis and center.  This was invented to allow light cones
            with thin slices to sample coherent large scale structure, but in
            practice does not work so well.  Try setting this parameter to 1 and
            see what happens.
            Default: 0.0.
        set_parameters : dict
            Dictionary of parameters to attach to pf.parameters.
            Default: None.
        output_dir : string
            The directory in which images and data files will be written.
            Default: 'LC'.
        output_prefix : string
            The prefix of all images and data files.
            Default: 'LightCone'.

        """

        self.near_redshift = near_redshift
        self.far_redshift = far_redshift
        self.observer_redshift = observer_redshift
        self.field_of_view_in_arcminutes = field_of_view_in_arcminutes
        self.image_resolution_in_arcseconds = image_resolution_in_arcseconds
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        if set_parameters is None:
            self.set_parameters = {}
        else:
            self.set_parameters = set_parameters
        self.output_dir = output_dir
        self.output_prefix = output_prefix

        self.master_solution = [] # kept to compare with recycled solutions
        self.projection_stack = []
        self.projection_weight_field_stack = []
        self.halo_mask = []

        # Original random seed of the first solution.
        self.original_random_seed = 0

        # Parameters for recycling light cone solutions.
        self.recycle_solution = False
        self.recycle_random_seed = 0

        # Calculate number of pixels.
        self.pixels = int(self.field_of_view_in_arcminutes * 60.0 / \
                          self.image_resolution_in_arcseconds)

        # Create output directory.
        if not os.path.exists(self.output_dir):
            only_on_root(os.mkdir, self.output_dir)

        # Calculate light cone solution.
        CosmologySplice.__init__(self, parameter_filename, simulation_type)
        self.light_cone_solution = \
          self.create_cosmology_splice(self.near_redshift, self.far_redshift,
                                       minimal=self.use_minimum_datasets,
                                       deltaz_min=self.deltaz_min)

    def calculate_light_cone_solution(self, seed=None, filename=None):
        r"""Create list of projections to be added together to make the light cone.

        Several sentences providing an extended description. Refer to
        variables using back-ticks, e.g. `var`.

        Parameters
        ----------
        seed : int
            The seed for the random number generator.  Any light cone solution
            can be reproduced by giving the same random seed.  Default: None
            (each solution will be distinct).
        filename : string
            If given, a text file detailing the solution will be written out.
            Default: None.
        """

        # Don't use box coherence with maximum projection depths.
        if self.use_minimum_datasets and \
                self.minimum_coherent_box_fraction > 0:
            mylog.info("Setting minimum_coherent_box_fraction to 0 with minimal light cone.")
            self.minimum_coherent_box_fraction = 0

        # Make sure recycling flag is off.
        self.recycle_solution = False

        # Get rid of old halo mask, if one was there.
        self.halo_mask = []

        if seed is not None:
            self.original_random_seed = int(seed)

        # Calculate projection sizes, and get
        # random projection axes and centers.
        na.random.seed(self.original_random_seed)

        # For box coherence, keep track of effective depth travelled.
        box_fraction_used = 0.0

        for q in range(len(self.light_cone_solution)):
            if self.light_cone_solution[q].has_key('previous'):
                del self.light_cone_solution[q]['previous']
            if self.light_cone_solution[q].has_key('next'):
                del self.light_cone_solution[q]['next']
            if (q == len(self.light_cone_solution) - 1):
                z_next = self.near_redshift
            else:
                z_next = self.light_cone_solution[q+1]['redshift']

            # Calculate fraction of box required for a depth of delta z
            self.light_cone_solution[q]['box_depth_fraction'] = \
                self.cosmology.ComovingRadialDistance(z_next, \
                        self.light_cone_solution[q]['redshift']) * \
                        self.simulation.hubble_constant / \
                        self.simulation.box_size

            # Simple error check to make sure more than 100% of box depth
            # is never required.
            if (self.light_cone_solution[q]['box_depth_fraction'] > 1.0):
                mylog.debug("Warning: box fraction required to go from z = %f to %f is %f" %
                            (self.light_cone_solution[q]['redshift'], z_next,
                             self.light_cone_solution[q]['box_depth_fraction']))
                mylog.debug("Full box delta z is %f, but it is %f to the next data dump." %
                            (self.light_cone_solution[q]['deltazMax'],
                             self.light_cone_solution[q]['redshift']-z_next))

            # Calculate fraction of box required for width corresponding to
            # requested image size.
            scale = self.cosmology.AngularScale_1arcsec_kpc(self.observer_redshift,
                self.light_cone_solution[q]['redshift'])
            size = self.field_of_view_in_arcminutes * 60.0 * scale / 1000.0
            boxSizeProper = self.simulation.box_size / \
              (self.simulation.hubble_constant *
               (1.0 + self.light_cone_solution[q]['redshift']))
            self.light_cone_solution[q]['box_width_fraction'] = size / boxSizeProper

            # Get projection axis and center.
            # If using box coherence, only get random axis and center if enough
            # of the box has been used, or if box_fraction_used will be greater
            # than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
              (box_fraction_used > self.minimum_coherent_box_fraction) or \
              (box_fraction_used +
               self.light_cone_solution[q]['box_depth_fraction'] > 1.0):
                # Random axis and center.
                self.light_cone_solution[q]['projection_axis'] = \
                  na.random.randint(0, 3)
                self.light_cone_solution[q]['projection_center'] = \
                  [na.random.random() for i in range(3)]
                box_fraction_used = 0.0
            else:
                # Same axis and center as previous slice,
                # but with depth center shifted.
                self.light_cone_solution[q]['projection_axis'] = \
                  self.light_cone_solution[q-1]['projection_axis']
                self.light_cone_solution[q]['projection_center'] = \
                  copy.deepcopy(self.light_cone_solution[q-1]['projection_center'])
                self.light_cone_solution[q]['projection_center']\
                  [self.light_cone_solution[q]['projection_axis']] += \
                    0.5 * (self.light_cone_solution[q]['box_depth_fraction'] +
                           self.light_cone_solution[q-1]['box_depth_fraction'])
                if self.light_cone_solution[q]['projection_center']\
                  [self.light_cone_solution[q]['projection_axis']] >= 1.0:
                    self.light_cone_solution[q]['projection_center']\
                      [self.light_cone_solution[q]['projection_axis']] -= 1.0

            box_fraction_used += self.light_cone_solution[q]['box_depth_fraction']

        # Store this as the master solution.
        self.master_solution = [copy.deepcopy(q) \
                                for q in self.light_cone_solution]

        # Write solution to a file.
        if filename is not None:
            self._save_light_cone_solution(filename=filename)

    def get_halo_mask(self, mask_file=None,
                      cube_file=None, map_file=None,
                      virial_overdensity=200,
                      halo_profiler_parameters=None,
                      njobs=1, dynamic=False):
        r"""Gets a halo mask from a file or makes a new one.

        Parameters
        ----------
        mask_file : string
            An hdf5 file to output the halo mask.
            Default: None
        cub_file : string
            An hdf5 file to output a halo mask for each slice
            of the light cone.
            Default: None
        map_file : string
            A text file to which to output the halo map (locations
            in the images of the halos.
            Default: None
        virial_overdensity : float
            The overdensity at which the virial radius is calculated
            and used as the radial for the halo mask.
            Default: 200
        halo_profiler_parameters: dict
            A dictionary of parameters to be passed to the HaloProfiler
            for each slice of the light cone.
            Default: None
        njobs : int
            The number of parallel jobs over which the slices for the
            halo mask will be split.  Choose -1 for one processor per
            individual slice and 1 to have all processors work together
            on each projection.
            Default: 1
        dynamic : bool
            If True, use dynamic load balancing to create the projections.
            Default: False.

        """

        if halo_profiler_parameters is None:
            halo_profiler_parameters = {}

        # Check if file already exists.
        if mask_file is not None and os.path.exists(mask_file):
            mylog.info('Reading halo mask from %s.' % mask_file)
            in_file = h5py.File(mask_file, 'r')
            self.halo_mask = in_file['HaloMask'].value
            in_file.close()

        # Otherwise, make a halo mask.
        else:
            halo_mask_cube = _light_cone_halo_mask(self, mask_file=mask_file,
                                                   cube_file=cube_file,
                                                   map_file=map_file,
                                                   virial_overdensity=\
                                                   virial_overdensity,
                                                   halo_profiler_parameters=\
                                                   halo_profiler_parameters,
                                                   njobs=njobs,
                                                   dynamic=dynamic)
            # Collapse cube into final mask.
            self.halo_mask = na.ones(shape=(self.pixels, self.pixels),
                                     dtype=bool)
            for mask in halo_mask_cube:
                self.halo_mask *= mask
            del halo_mask_cube

    def project_light_cone(self, field, weight_field=None, apply_halo_mask=False,
                           node=None, save_stack=True, save_slice_images=False,
                           cmap_name='algae', photon_field=False,
                           njobs=1, dynamic=False):
        r"""Create projections for light cone, then add them together.

        Parameters
        ----------
        field : string
            The projected field.
        weight_field : string
            the weight field of the projection.  This has the same meaning as
            in standard projections.
            Default: None.
        apply_halo_mask : bool
            if True, a boolean mask is apply to the light cone projection.  See
            below for a description of halo masks.
            Default: False.
        node : string
            a prefix to be prepended to the node name under which the
            projection data is serialized.
            Default: None.
        save_stack : bool
            if True, the light cone data including each individual
            slice is written to an hdf5 file.
            Default: True.
        save_slice_images : bool
            save images for each individual projection slice.
            Default: False.
        cmap_name : string
            color map for images.
            Default: 'algae'.
        photon_field : bool
            if True, the projection data for each slice is decremented by 4 Pi
            R^2`, where R is the luminosity distance between the observer and
            the slice redshift.
            Default: False.
        njobs : int
            The number of parallel jobs over which the light cone projection
            will be split.  Choose -1 for one processor per individual
            projection and 1 to have all processors work together on each
            projection.
            Default: 1.
        dynamic : bool
            If True, use dynamic load balancing to create the projections.
            Default: False.
        """

        # Clear projection stack.
        self.projection_stack = []
        self.projection_weight_field_stack = []
        if (self.light_cone_solution[-1].has_key('object')):
            del self.light_cone_solution[-1]['object']

        if not(self.output_dir.endswith("/")):
                 self.output_dir += "/"

        # for q, output in enumerate(self.light_cone_solution):
        all_storage = {}
        for my_storage, output in parallel_objects(self.light_cone_solution,
                                                   storage=all_storage,
                                                   dynamic=dynamic):
            output['object'] = load(output['filename'])
            output['object'].parameters.update(self.set_parameters)
            frb = _light_cone_projection(output, field, self.pixels,
                                         weight_field=weight_field, node=node)

            if photon_field:
                # Decrement the flux by the luminosity distance.
                # Assume field in frb is in erg/s/cm^2/Hz
                dL = self.cosmology.LuminosityDistance(self.observer_redshift,
                                                       output['redshift']) #in Mpc
                boxSizeProper = self.simulation.box_size / \
                  (self.simulation.hubble_constant * (1.0 + output['redshift']))
                pixelarea = (boxSizeProper/self.pixels)**2 #in proper cm^2
                factor = pixelarea/(4.0*na.pi*dL**2)
                mylog.info("Distance to slice = %e" % dL)
                frb[field] *= factor #in erg/s/cm^2/Hz on observer's image plane.

            if weight_field is None:
                my_storage.result = {'field': frb[field]}
            else:
                my_storage.result = {'field': (frb[field] *
                                               frb['weight_field']),
                                     'weight_field': frb['weight_field']}

            del output['object']

        # Combine results from each slice.
        all_slices = all_storage.keys()
        all_slices.sort()
        for my_slice in all_slices:
            if save_slice_images:
                if node is None:
                    name = "%s%s_%04d_%04d" % (self.output_dir,
                                               self.output_prefix,
                                               my_slice,
                                               len(self.light_cone_solution))
                else:
                    name = "%s%s_%s_%04d_%04d" % (self.output_dir,
                                                  self.output_prefix,
                                                  node, my_slice,
                                                  len(self.light_cone_solution))
                if weight_field is None:
                    my_image = all_storage[my_slice]['field']
                else:
                    my_image = all_storage[my_slice]['field'] / \
                      all_storage[my_slice]['weight_field']
                only_on_root(write_image, na.log10(my_image),
                             "%s_%s.png" % (name, field), cmap_name=cmap_name)

            self.projection_stack.append(all_storage[my_slice]['field'])
            if weight_field is not None:
                self.projection_weight_field_stack.append(all_storage[my_slice]['field'])

        # Add up slices to make light cone projection.
        if (weight_field is None):
            light_cone_projection = sum(self.projection_stack)
        else:
            light_cone_projection = sum(self.projection_stack) / \
              sum(self.projection_weight_field_stack)

        if node is None:
            filename = "%s%s" % (self.output_dir, self.output_prefix)
        else:
            filename = "%s%s_%s" % (self.output_dir, self.output_prefix, node)

        # Apply halo mask.
        if apply_halo_mask:
            if len(self.halo_mask) > 0:
                mylog.info("Applying halo mask.")
                light_cone_projection *= self.halo_mask
            else:
                mylog.error("No halo mask loaded, call get_halo_mask.")


        # Write image.
        if save_slice_images:
            only_on_root(write_image, na.log10(light_cone_projection),
                         "%s_%s.png" % (filename, field), cmap_name=cmap_name)

        # Write stack to hdf5 file.
        if save_stack:
            self._save_light_cone_stack(field=field, weight_field=weight_field,
                                        filename=filename)

    def rerandomize_light_cone_solution(self, new_seed, recycle=True, filename=None):
        """
        When making a projection for a light cone, only randomizations along the
        line of sight make any given projection unique, since the lateral shifting
        and tiling is done after the projection is made.  Therefore, multiple light
        cones can be made from a single set of projections by introducing different
        lateral random shifts and keeping all the original shifts along the line of
        sight.  This routine will take in a new random seed and rerandomize the
        parts of the light cone that do not contribute to creating a unique
        projection object.  Additionally, this routine is built such that if the
        same random seed is given for the rerandomizing, the solution will be
        identical to the original.

        This routine has now been updated to be a general solution rescrambler.
        If the keyword recycle is set to True, then it will recycle.  Otherwise, it
        will create a completely new solution.

        new_sed : float
            The new random seed.
        recycle : bool
            If True, the new solution will have the same shift in the line of
            sight as the original solution.  Since the projections of each
            slice are serialized and stored for the entire width of the box
            (even if the width used is left than the total box), the projection
            data can be deserialized instead of being remade from scratch.
            This can greatly speed up the creation of a large number of light
            cone projections.
            Default: True.
        filename : string
            If given, a text file detailing the solution will be written out.
            Default: None.

        """

        # Get rid of old halo mask, if one was there.
        self.halo_mask = []

        # Clean pf objects out of light cone solution.
        for my_slice in self.light_cone_solution:
            if my_slice.has_key('object'):
                del my_slice['object']

        if recycle:
            mylog.debug("Recycling solution made with %s with new seed %s." %
                        (self.original_random_seed, new_seed))
            self.recycle_random_seed = int(new_seed)
        else:
            mylog.debug("Creating new solution with random seed %s." % new_seed)
            self.original_random_seed = int(new_seed)
            self.recycle_random_seed = 0

        self.recycle_solution = recycle

        # Keep track of fraction of volume in common between the original and
        # recycled solution.
        common_volume = 0.0
        total_volume = 0.0

        # For box coherence, keep track of effective depth travelled.
        box_fraction_used = 0.0

        # Seed random number generator with new seed.
        na.random.seed(int(new_seed))

        for q, output in enumerate(self.light_cone_solution):
            # It is necessary to make the same number of calls to the random
            # number generator so the original solution willbe produced if the
            # same seed is given.

            # Get projection axis and center.
            # If using box coherence, only get random axis and center if enough
            # of the box has been used, or if box_fraction_used will be greater
            # than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
                    (box_fraction_used > self.minimum_coherent_box_fraction) or \
                    (box_fraction_used + self.light_cone_solution[q]['box_depth_fraction'] > 1.0):
                # Get random projection axis and center.
                # If recycling, axis will get thrown away since it is used in
                # creating a unique projection object.
                newAxis = na.random.randint(0, 3)

                newCenter = [na.random.random() for i in range(3)]
                box_fraction_used = 0.0
            else:
                # Same axis and center as previous slice, but with depth center shifted.
                newAxis = self.light_cone_solution[q-1]['projection_axis']
                newCenter = copy.deepcopy(self.light_cone_solution[q-1]['projection_center'])
                newCenter[newAxis] += \
                    0.5 * (self.light_cone_solution[q]['box_depth_fraction'] +
                           self.light_cone_solution[q-1]['box_depth_fraction'])
                if newCenter[newAxis] >= 1.0:
                    newCenter[newAxis] -= 1.0

            if recycle:
                output['projection_axis'] = self.master_solution[q]['projection_axis']
            else:
                output['projection_axis'] = newAxis

            box_fraction_used += self.light_cone_solution[q]['box_depth_fraction']

            # Make list of rectangle corners to calculate common volume.
            newCube = na.zeros(shape=(len(newCenter), 2))
            oldCube = na.zeros(shape=(len(newCenter), 2))
            for w in range(len(newCenter)):
                if (w == self.master_solution[q]['projection_axis']):
                    oldCube[w] = [self.master_solution[q]['projection_center'][w] -
                                  0.5 * self.master_solution[q]['box_depth_fraction'],
                                  self.master_solution[q]['projection_center'][w] +
                                  0.5 * self.master_solution[q]['box_depth_fraction']]
                else:
                    oldCube[w] = [self.master_solution[q]['projection_center'][w] -
                                  0.5 * self.master_solution[q]['box_width_fraction'],
                                  self.master_solution[q]['projection_center'][w] +
                                  0.5 * self.master_solution[q]['box_width_fraction']]

                if (w == output['projection_axis']):
                    if recycle:
                        newCube[w] = oldCube[w]
                    else:
                        newCube[w] = \
                          [newCenter[w] -
                           0.5 * self.master_solution[q]['box_depth_fraction'],
                           newCenter[w] +
                           0.5 * self.master_solution[q]['box_depth_fraction']]
                else:
                    newCube[w] = [newCenter[w] -
                                  0.5 * self.master_solution[q]['box_width_fraction'],
                                  newCenter[w] +
                                  0.5 * self.master_solution[q]['box_width_fraction']]

            common_volume += common_volume(oldCube, newCube,
                                           periodic=na.array([[0, 1],
                                                              [0, 1],
                                                              [0, 1]]))
            total_volume += output['box_depth_fraction'] * \
              output['box_width_fraction']**2

            # Replace centers for every axis except the line of sight axis.
            for w in range(len(newCenter)):
                if not(recycle and
                       (w == self.light_cone_solution[q]['projection_axis'])):
                    self.light_cone_solution[q]['projection_center'][w] = \
                      newCenter[w]

        if recycle:
            mylog.debug("Fractional common volume between master and recycled solution is %.2e" % \
                        (common_volume / total_volume))
        else:
            mylog.debug("Fraction of total volume in common with old solution is %.2e." % \
                        (common_volume / total_volume))
            self.master_solution = [copy.deepcopy(q) \
                                    for q in self.light_cone_solution]

        # Write solution to a file.
        if filename is not None:
            self._save_light_cone_solution(filename=filename)

    def restore_master_solution(self):
        "Reset the active light cone solution to the master solution."
        self.light_cone_solution = [copy.deepcopy(q) \
                                    for q in self.master_solution]

    @parallel_root_only
    def _save_light_cone_solution(self, filename="light_cone.dat"):
        "Write out a text file with information on light cone solution."

        mylog.info("Saving light cone solution to %s." % filename)

        f = open(filename, 'w')
        if self.recycle_solution:
            f.write("Recycled Solution\n")
            f.write("OriginalRandomSeed = %s\n" % self.original_random_seed)
            f.write("RecycleRandomSeed = %s\n" % self.recycle_random_seed)
        else:
            f.write("Original Solution\n")
            f.write("OriginalRandomSeed = %s\n" % self.original_random_seed)
        f.write("parameter_filename = %s\n" % self.parameter_filename)
        f.write("\n")
        for q, output in enumerate(self.light_cone_solution):
            f.write("Proj %04d, %s, z = %f, depth/box = %f, width/box = %f, axis = %d, center = %f, %f, %f\n" %
                    (q, output['filename'], output['redshift'],
                     output['box_depth_fraction'], output['box_width_fraction'],
                     output['projection_axis'], output['projection_center'][0],
                     output['projection_center'][1], output['projection_center'][2]))
        f.close()

    @parallel_root_only
    def _save_light_cone_stack(self, field=None, weight_field=None,
                               filename=None, over_write=True):
        "Save the light cone projection stack as a 3d array in and hdf5 file."

        # Make list of redshifts to include as a dataset attribute.
        redshiftList = na.array([my_slice['redshift'] \
                                 for my_slice in self.light_cone_solution])

        field_node = "%s_%s" % (field, weight_field)
        weight_field_node = "weight_field_%s" % weight_field

        if (filename is None):
            filename = "%s/%s_data" % (self.output_dir, self.output_prefix)
        if not(filename.endswith('.h5')):
               filename += ".h5"

        if (len(self.projection_stack) == 0):
            mylog.debug("save_light_cone_stack: no projection data loaded.")
            return

        mylog.info("Writing light cone data to %s." % filename)

        output = h5py.File(filename, "a")

        node_exists = field_node in output

        if node_exists:
            if over_write:
                mylog.info("Dataset, %s, already exists, overwriting." %
                           field_node)
                write_data = True
                del output[field_node]
            else:
                mylog.info("Dataset, %s, already exists in %s, not saving." %
                           (field_node, filename))
                write_data = False
        else:
            write_data = True

        if write_data:
            mylog.info("Saving %s to %s." % (field_node, filename))
            self.projection_stack = na.array(self.projection_stack)
            field_dataset = output.create_dataset(field_node,
                                                  data=self.projection_stack)
            field_dataset.attrs['redshifts'] = redshiftList
            field_dataset.attrs['observer_redshift'] = \
              na.float(self.observer_redshift)
            field_dataset.attrs['field_of_view_in_arcminutes'] = \
              na.float(self.field_of_view_in_arcminutes)
            field_dataset.attrs['image_resolution_in_arcseconds'] = \
              na.float(self.image_resolution_in_arcseconds)

        if (len(self.projection_weight_field_stack) > 0):
            if node_exists:
                if over_write:
                    mylog.info("Dataset, %s, already exists, overwriting." %
                               weight_field_node)
                    del output[field_node]
                else:
                    mylog.info("Dataset, %s, already exists in %s, not saving." %
                               (weight_field_node, filename))
                    write_data = False
            else:
                write_data = True

            if write_data:
                mylog.info("Saving %s to %s." % (weight_field_node, filename))
                self.projection_weight_field_stack = \
                  na.array(self.projection_weight_field_stack)
                weight_field_dataset = \
                  output.create_dataset(weight_field_node,
                                        data=self.projection_weight_field_stack)
                weight_field_dataset.attrs['redshifts'] = redshiftList
                weight_field_dataset.attrs['observer_redshift'] = \
                  na.float(self.observer_redshift)
                weight_field_dataset.attrs['field_of_view_in_arcminutes'] = \
                  na.float(self.field_of_view_in_arcminutes)
                weight_field_dataset.attrs['image_resolution_in_arcseconds'] = \
                  na.float(self.image_resolution_in_arcseconds)

        output.close()
