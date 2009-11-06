"""
LightCone class and member functions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Britton Smith.  All Rights Reserved.

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

from yt.extensions.lightcone import *
from yt.extensions.EnzoSimulation import *
from yt.logger import lagosLogger as mylog
from yt.config import ytcfg
from yt.funcs import *
from Common_nVolume import *
from HaloMask import *
import copy
import os
import numpy as na

class LightCone(EnzoSimulation):
    def __init__(self, EnzoParameterFile, initial_redshift=1.0, final_redshift=0.0, observer_redshift=0.0,
                 field_of_view_in_arcminutes=600.0, image_resolution_in_arcseconds=60.0, 
                 use_minimum_datasets=True, deltaz_min=0.0, minimum_coherent_box_fraction=0.0,
                 output_dir='LC', output_prefix='LightCone'):
        """
        Initialize a LightCone object.
        :param initial_redshift (float): the initial (highest) redshift for the light cone.  Default: 1.0.
        :param final_redshift (float): the final (lowest) redshift for the light cone.  Default: 0.0.
        :param observer_redshift (float): the redshift of the observer.  Default: 0.0.
        :param field_of_view_in_arcminutes (float): the field of view of the image in units of arcminutes.  
               Default: 600.0.
        :param image_resolution_in_arcseconds (float): the size of each image pixel in units of arcseconds.  
               Default: 60.0.
        :param use_minimum_datasets (bool): if True, the minimum number of datasets is used to connect the 
               initial and final redshift.  If false, the light cone solution will contain as many entries 
               as possible within the redshift interval.  Default: True.
        :param deltaz_min (float): specifies the minimum :math:`\Delta z` between consecutive datasets in 
               the returned list.  Default: 0.0.
        :param minimum_coherent_box_fraction (float): used with use_minimum_datasets set to False, this 
               parameter specifies the fraction of the total box size to be traversed before rerandomizing 
               the projection axis and center.  This was invented to allow light cones with thin slices to 
               sample coherent large scale structure, but in practice does not work so well.  Try setting 
               this parameter to 1 and see what happens.  Default: 0.0.
        :param output_dir (str): the directory in which images and data files will be written.  Default: 'LC'.
        :param output_prefix (str): the prefix of all images and data files.  Default: 'LightCone'.
        """

        self.initial_redshift = initial_redshift
        self.final_redshift = final_redshift
        self.observer_redshift = observer_redshift
        self.field_of_view_in_arcminutes = field_of_view_in_arcminutes
        self.image_resolution_in_arcseconds = image_resolution_in_arcseconds
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        self.output_dir = output_dir
        self.output_prefix = output_prefix

        self.master_solution = [] # kept to compare with recycled solutions
        self.projection_stack = []
        self.projection_weight_field_stack = []
        self.halo_mask = []

        # Original random seed of the first solution.
        self.originalRandomSeed = 0

        # Parameters for recycling light cone solutions.
        self.recycleSolution = False
        self.recycleRandomSeed = 0

        # Initialize EnzoSimulation machinery for getting dataset list.
        EnzoSimulation.__init__(self, EnzoParameterFile, initial_redshift=self.initial_redshift,
                                final_redshift=self.final_redshift, links=True,
                                enzo_parameters={'CosmologyComovingBoxSize':float})

        # Calculate number of pixels.
        self.pixels = int(self.field_of_view_in_arcminutes * 60.0 / \
                          self.image_resolution_in_arcseconds)

        if ytcfg.getint("yt", "__parallel_rank") == 0:
            # Create output directory.
            if (os.path.exists(self.output_dir)):
                if not(os.path.isdir(self.output_dir)):
                    mylog.error("Output directory exists, but is not a directory: %s." % self.output_dir)
                    self.output_dir = './'
            else:
                os.mkdir(self.output_dir)

        # Get list of datasets for light cone solution.
        self.light_cone_solution = self.create_cosmology_splice(minimal=self.use_minimum_datasets,
                                                                deltaz_min=self.deltaz_min)

    def calculate_light_cone_solution(self, seed=None, filename=None):
        """
        Create list of projections to be added together to make the light cone.
        :param seed (int): the seed for the random number generator.  Any light cone solution 
               can be reproduced by giving the same random seed.  Default: None (each solution 
               will be distinct).
        :param filename (str): if given, a text file detailing the solution will be written out.  Default: None.
        """

        # Don't use box coherence with maximum projection depths.
        if self.use_minimum_datasets and \
                self.minimum_coherent_box_fraction > 0:
            mylog.info("Setting minimum_coherent_box_fraction to 0 with minimal light cone.")
            self.minimum_coherent_box_fraction = 0

        # Make sure recycling flag is off.
        self.recycleSolution = False

        # Get rid of old halo mask, if one was there.
        self.halo_mask = []

        if seed is not None:
            self.originalRandomSeed = int(seed)

        # Calculate projection sizes, and get random projection axes and centers.
        na.random.seed(self.originalRandomSeed)

        # For box coherence, keep track of effective depth travelled.
        boxFractionUsed = 0.0

        for q in range(len(self.light_cone_solution)):
            del self.light_cone_solution[q]['previous']
            del self.light_cone_solution[q]['next']
            if (q == len(self.light_cone_solution) - 1):
                z_next = self.final_redshift
            else:
                z_next = self.light_cone_solution[q+1]['redshift']

            # Calculate fraction of box required for a depth of delta z
            self.light_cone_solution[q]['DepthBoxFraction'] = self.cosmology.ComovingRadialDistance(z_next, self.light_cone_solution[q]['redshift']) * \
                self.enzoParameters['CosmologyHubbleConstantNow'] / self.enzoParameters['CosmologyComovingBoxSize']

            # Simple error check to make sure more than 100% of box depth is never required.
            if (self.light_cone_solution[q]['DepthBoxFraction'] > 1.0):
                mylog.debug("Warning: box fraction required to go from z = %f to %f is %f" % 
                            (self.light_cone_solution[q]['redshift'], z_next,
                             self.light_cone_solution[q]['DepthBoxFraction']))
                mylog.debug("Full box delta z is %f, but it is %f to the next data dump." % 
                            (self.light_cone_solution[q]['deltazMax'],
                             self.light_cone_solution[q]['redshift']-z_next))

            # Calculate fraction of box required for width corresponding to requested image size.
            scale = self.cosmology.AngularScale_1arcsec_kpc(self.observer_redshift, self.light_cone_solution[q]['redshift'])
            size = self.field_of_view_in_arcminutes * 60.0 * scale / 1000.0
            boxSizeProper = self.enzoParameters['CosmologyComovingBoxSize'] / (self.enzoParameters['CosmologyHubbleConstantNow'] * 
                                                                               (1.0 + self.light_cone_solution[q]['redshift']))
            self.light_cone_solution[q]['WidthBoxFraction'] = size / boxSizeProper

            # Get projection axis and center.
            # If using box coherence, only get random axis and center if enough of the box has been used, 
            # or if boxFractionUsed will be greater than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
                    (boxFractionUsed > self.minimum_coherent_box_fraction) or \
                    (boxFractionUsed + self.light_cone_solution[q]['DepthBoxFraction'] > 1.0):
                # Random axis and center.
                self.light_cone_solution[q]['ProjectionAxis'] = na.random.randint(0, 3)
                self.light_cone_solution[q]['ProjectionCenter'] = [na.random.random(), na.random.random(), na.random.random()]
                boxFractionUsed = 0.0
            else:
                # Same axis and center as previous slice, but with depth center shifted.
                self.light_cone_solution[q]['ProjectionAxis'] = self.light_cone_solution[q-1]['ProjectionAxis']
                self.light_cone_solution[q]['ProjectionCenter'] = copy.deepcopy(self.light_cone_solution[q-1]['ProjectionCenter'])
                self.light_cone_solution[q]['ProjectionCenter'][self.light_cone_solution[q]['ProjectionAxis']] += \
                    0.5 * (self.light_cone_solution[q]['DepthBoxFraction'] + self.light_cone_solution[q-1]['DepthBoxFraction'])
                if self.light_cone_solution[q]['ProjectionCenter'][self.light_cone_solution[q]['ProjectionAxis']] >= 1.0:
                    self.light_cone_solution[q]['ProjectionCenter'][self.light_cone_solution[q]['ProjectionAxis']] -= 1.0

            boxFractionUsed += self.light_cone_solution[q]['DepthBoxFraction']

        # Store this as the master solution.
        self.master_solution = [copy.deepcopy(q) for q in self.light_cone_solution]

        # Write solution to a file.
        if filename is not None:
            self._save_light_cone_solution(filename=filename)

    def get_halo_mask(self, mask_file=None, map_file=None, **kwargs):
        """
        Gets a halo mask from a file or makes a new one.
        :param mask_file (str): specify an hdf5 file to output the halo mask.
        :param map_file (str): specify a text file to output the halo map 
               (locations in image of halos).
        """

        # Get halo map if map_file given.
        if map_file is not None and not os.path.exists(map_file):
            light_cone_halo_map(self, map_file=map_file, **kwargs)

        # Check if file already exists.
        if mask_file is not None and os.path.exists(mask_file):
            mylog.info('Reading halo mask from %s.' % mask_file)
            input = h5py.File(mask_file, 'r')
            self.halo_mask = input['HaloMask'].value
            input.close()

        # Otherwise, make a halo mask.
        else:
            halo_mask_cube = light_cone_halo_mask(self, mask_file=mask_file, **kwargs)
            # Collapse cube into final mask.
            if ytcfg.getint("yt", "__parallel_rank") == 0:
                self.halo_mask = na.ones(shape=(self.pixels, self.pixels), dtype=bool)
                for mask in halo_mask_cube:
                    self.halo_mask *= mask
            del halo_mask_cube

    def project_light_cone(self, field, weight_field=None, apply_halo_mask=False, node=None,
                           save_stack=True, save_slice_images=False, flatten_stack=False, photon_field=False, **kwargs):
        """
        Create projections for light cone, then add them together.
        :param weight_field (str): the weight field of the projection.  This has the same meaning as in standard 
               projections.  Default: None.
        :param apply_halo_mask (bool): if True, a boolean mask is apply to the light cone projection.  See below for a 
               description of halo masks.  Default: False.
        :param node (str): a prefix to be prepended to the node name under which the projection data is serialized.  
               Default: None.
        :param save_stack (bool): if True, the unflatted light cone data including each individual slice is written to 
               an hdf5 file.  Default: True.
        :param save_slice_images (bool): save images for each individual projection slice.  Default: False.
        :param flatten_stack (bool): if True, the light cone stack is continually flattened each time a slice is added 
               in order to save memory.  This is generally not necessary.  Default: False.
        :param photon_field (bool): if True, the projection data for each slice is decremented by 4 Pi R^2`, where R 
               is the luminosity distance between the observer and the slice redshift.  Default: False.
        """

        # Clear projection stack.
        self.projection_stack = []
        self.projection_weight_field_stack = []
        if (self.light_cone_solution[-1].has_key('object')):
            del self.light_cone_solution[-1]['object']

        if not(self.output_dir.endswith("/")):
                 self.output_dir += "/"

        for q, output in enumerate(self.light_cone_solution):
            if node is None:
                name = "%s%s_%04d_%04d" % (self.output_dir, self.output_prefix,
                                           q, len(self.light_cone_solution))
            else:
                name = "%s%s_%s_%04d_%04d" % (self.output_dir, self.output_prefix,
                                              node, q, len(self.light_cone_solution))
            output['object'] = lagos.EnzoStaticOutput(output['filename'])
            frb = LightConeProjection(output, field, self.pixels, weight_field=weight_field,
                                      save_image=save_slice_images,
                                      name=name, node=node, **kwargs)
            if ytcfg.getint("yt", "__parallel_rank") == 0:
                if photon_field:
                    # Decrement the flux by the luminosity distance. Assume field in frb is in erg/s/cm^2/Hz
                    co = lagos.Cosmology(HubbleConstantNow = (100.0 * self.enzoParameters['CosmologyHubbleConstantNow']),
                                         OmegaMatterNow = self.enzoParameters['CosmologyOmegaMatterNow'],
                                         OmegaLambdaNow = self.enzoParameters['CosmologyOmegaLambdaNow'])
                    dL = self.cosmology.LuminosityDistance(self.observer_redshift, output['redshift']) #in Mpc
                    boxSizeProper = self.enzoParameters['CosmologyComovingBoxSize'] / (self.enzoParameters['CosmologyHubbleConstantNow'] * 
                                                                                       (1.0 + output['redshift']))
                    pixelarea = (boxSizeProper/self.pixels)**2 #in proper cm^2
                    factor = pixelarea/(4.0*na.pi*dL**2)
                    mylog.info("Distance to slice = %e" % dL)
                    frb[field] *= factor #in erg/s/cm^2/Hz on observer's image plane.

            if ytcfg.getint("yt", "__parallel_rank") == 0:
                if weight_field is not None:
                    # Data come back normalized by the weight field.
                    # Undo that so it can be added up for the light cone.
                    self.projection_stack.append(frb[field]*frb['weight_field'])
                    self.projection_weight_field_stack.append(frb['weight_field'])
                else:
                    self.projection_stack.append(frb[field])

                # Delete the frb.  This saves a decent amount of ram.
                if (q < len(self.light_cone_solution) - 1):
                    del frb

                # Flatten stack to save memory.
                if flatten_stack and (len(self.projection_stack) > 1):
                    self.projection_stack = [sum(self.projection_stack)]
                    if weight_field is not None:
                        self.projection_weight_field_stack = [sum(self.projection_weight_field_stack)]

            # Delete the plot collection now that the frb is deleted.
            del output['pc']

            # Unless this is the last slice, delete the dataset object.
            # The last one will be saved to make the plot collection.
            if (q < len(self.light_cone_solution) - 1):
                del output['object']

        if ytcfg.getint("yt", "__parallel_rank") == 0:
            # Add up slices to make light cone projection.
            if (weight_field is None):
                lightConeProjection = sum(self.projection_stack)
            else:
                lightConeProjection = sum(self.projection_stack) / sum(self.projection_weight_field_stack)

            if node is None:
                filename = "%s%s" % (self.output_dir, self.output_prefix)
            else:
                filename = "%s%s_%s" % (self.output_dir, self.output_prefix, node)

            # Save the last fixed resolution buffer for the plot collection, 
            # but replace the data with the full light cone projection data.
            frb.data[field] = lightConeProjection

            # Write stack to hdf5 file.
            if save_stack:
                self._save_light_cone_stack(field=field, weight_field=weight_field, filename=filename)

            # Apply halo mask.
            if apply_halo_mask:
                if len(self.halo_mask) > 0:
                    mylog.info("Applying halo mask.")
                    frb.data[field] *= self.halo_mask
                else:
                    mylog.error("No halo mask loaded, call get_halo_mask.")

            # Make a plot collection for the light cone projection.
            center = [0.5 * (self.light_cone_solution[-1]['object'].parameters['DomainLeftEdge'][w] + 
                             self.light_cone_solution[-1]['object'].parameters['DomainRightEdge'][w])
                      for w in range(self.light_cone_solution[-1]['object'].parameters['TopGridRank'])]
            pc = raven.PlotCollection(self.light_cone_solution[-1]['object'], center=center)
            pc.add_fixed_resolution_plot(frb, field, **kwargs)
            pc.save(filename)

            # Return the plot collection so the user can remake the plot if they want.
            return pc

    def rerandomize_light_cone_solution(self, newSeed, recycle=True, filename=None):
        """
        When making a projection for a light cone, only randomizations along the line of sight make any 
        given projection unique, since the lateral shifting and tiling is done after the projection is made.
        Therefore, multiple light cones can be made from a single set of projections by introducing different 
        lateral random shifts and keeping all the original shifts along the line of sight.
        This routine will take in a new random seed and rerandomize the parts of the light cone that do not contribute 
        to creating a unique projection object.  Additionally, this routine is built such that if the same random 
        seed is given for the rerandomizing, the solution will be identical to the original.

        This routine has now been updated to be a general solution rescrambler.  If the keyword recycle is set to 
        True, then it will recycle.  Otherwise, it will create a completely new solution.

        :param recycle (bool): if True, the new solution will have the same shift in the line of sight as the original 
               solution.  Since the projections of each slice are serialized and stored for the entire width of the 
               box (even if the width used is left than the total box), the projection data can be deserialized 
               instead of being remade from scratch.  This can greatly speed up the creation of a large number of light 
               cone projections.  Default: True.
        :param filename (str): if given, a text file detailing the solution will be written out.  Default: None.
        """

        # Get rid of old halo mask, if one was there.
        self.halo_mask = []

        # Clean pf objects out of light cone solution.
        for slice in self.light_cone_solution:
            if slice.has_key('object'):
                del slice['object']

        if recycle:
            mylog.debug("Recycling solution made with %s with new seed %s." % 
                        (self.originalRandomSeed, newSeed))
            self.recycleRandomSeed = int(newSeed)
        else:
            mylog.debug("Creating new solution with random seed %s." % newSeed)
            self.originalRandomSeed = int(newSeed)
            self.recycleRandomSeed = 0

        self.recycleSolution = recycle

        # Keep track of fraction of volume in common between the original and recycled solution.
        commonVolume = 0.0
        totalVolume = 0.0

        # For box coherence, keep track of effective depth travelled.
        boxFractionUsed = 0.0

        # Seed random number generator with new seed.
        na.random.seed(int(newSeed))

        for q, output in enumerate(self.light_cone_solution):
            # It is necessary to make the same number of calls to the random number generator
            # so the original solution willbe produced if the same seed is given.

            # Get projection axis and center.
            # If using box coherence, only get random axis and center if enough of the box has been used, 
            # or if boxFractionUsed will be greater than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
                    (boxFractionUsed > self.minimum_coherent_box_fraction) or \
                    (boxFractionUsed + self.light_cone_solution[q]['DepthBoxFraction'] > 1.0):
                # Get random projection axis and center.
                # If recycling, axis will get thrown away since it is used in creating a unique projection object.
                newAxis = na.random.randint(0, 3)

                newCenter = [na.random.random(), na.random.random(), na.random.random()]
                boxFractionUsed = 0.0
            else:
                # Same axis and center as previous slice, but with depth center shifted.
                newAxis = self.light_cone_solution[q-1]['ProjectionAxis']
                newCenter = copy.deepcopy(self.light_cone_solution[q-1]['ProjectionCenter'])
                newCenter[newAxis] += \
                    0.5 * (self.light_cone_solution[q]['DepthBoxFraction'] + self.light_cone_solution[q-1]['DepthBoxFraction'])
                if newCenter[newAxis] >= 1.0:
                    newCenter[newAxis] -= 1.0

            if recycle:
                output['ProjectionAxis'] = self.master_solution[q]['ProjectionAxis']
            else:
                output['ProjectionAxis'] = newAxis

            boxFractionUsed += self.light_cone_solution[q]['DepthBoxFraction']

            # Make list of rectangle corners to calculate common volume.
            newCube = na.zeros(shape=(len(newCenter), 2))
            oldCube = na.zeros(shape=(len(newCenter), 2))
            for w in range(len(newCenter)):
                if (w == self.master_solution[q]['ProjectionAxis']):
                    oldCube[w] = [self.master_solution[q]['ProjectionCenter'][w] - 0.5 * self.master_solution[q]['DepthBoxFraction'],
                                  self.master_solution[q]['ProjectionCenter'][w] + 0.5 * self.master_solution[q]['DepthBoxFraction']]
                else:
                    oldCube[w] = [self.master_solution[q]['ProjectionCenter'][w] - 0.5 * self.master_solution[q]['WidthBoxFraction'],
                                  self.master_solution[q]['ProjectionCenter'][w] + 0.5 * self.master_solution[q]['WidthBoxFraction']]

                if (w == output['ProjectionAxis']):
                    if recycle:
                        newCube[w] = oldCube[w]
                    else:
                        newCube[w] = [newCenter[w] - 0.5 * self.master_solution[q]['DepthBoxFraction'],
                                      newCenter[w] + 0.5 * self.master_solution[q]['DepthBoxFraction']]
                else:
                    newCube[w] = [newCenter[w] - 0.5 * self.master_solution[q]['WidthBoxFraction'],
                                  newCenter[w] + 0.5 * self.master_solution[q]['WidthBoxFraction']]

            commonVolume += commonNVolume(oldCube, newCube, periodic=na.array([[0, 1], [0, 1], [0, 1]]))
            totalVolume += output['DepthBoxFraction'] * output['WidthBoxFraction']**2

            # Replace centers for every axis except the line of sight axis.
            for w in range(len(newCenter)):
                if not(recycle and (w == self.light_cone_solution[q]['ProjectionAxis'])):
                    self.light_cone_solution[q]['ProjectionCenter'][w] = newCenter[w]

        if recycle:
            mylog.debug("Fractional common volume between master and recycled solution is %.2e" % (commonVolume/totalVolume))
        else:
            mylog.debug("Fraction of total volume in common with old solution is %.2e." % (commonVolume/totalVolume))
            self.master_solution = [copy.deepcopy(q) for q in self.light_cone_solution]

        # Write solution to a file.
        if filename is not None:
            self._save_light_cone_solution(filename=filename)

    def restore_master_solution(self):
        "Reset the active light cone solution to the master solution."
        self.light_cone_solution = [copy.deepcopy(q) for q in self.master_solution]

    @rootonly
    def _save_light_cone_solution(self, filename="light_cone.dat"):
        "Write out a text file with information on light cone solution."

        mylog.info("Saving light cone solution to %s." % filename)

        f = open(filename, 'w')
        if self.recycleSolution:
            f.write("Recycled Solution\n")
            f.write("OriginalRandomSeed = %s\n" % self.originalRandomSeed)
            f.write("RecycleRandomSeed = %s\n" % self.recycleRandomSeed)
        else:
            f.write("Original Solution\n")
            f.write("OriginalRandomSeed = %s\n" % self.originalRandomSeed)
        f.write("EnzoParameterFile = %s\n" % self.EnzoParameterFile)
        f.write("\n")
        for q, output in enumerate(self.light_cone_solution):
            f.write("Proj %04d, %s, z = %f, depth/box = %f, width/box = %f, axis = %d, center = %f, %f, %f\n" %
                    (q, output['filename'], output['redshift'], output['DepthBoxFraction'], output['WidthBoxFraction'],
                    output['ProjectionAxis'], output['ProjectionCenter'][0], output['ProjectionCenter'][1], output['ProjectionCenter'][2]))
        f.close()

    def _save_light_cone_stack(self, field=None, weight_field=None, filename=None, over_write=True):
        "Save the light cone projection stack as a 3d array in and hdf5 file."

        # Make list of redshifts to include as a dataset attribute.
        redshiftList = na.array([slice['redshift'] for slice in self.light_cone_solution])

        field_node = "%s_%s" % (field, weight_field)
        weight_field_node = "weight_field_%s" % weight_field

        import h5py
        if (filename is None):
            filename = "%s/%s_data" % (self.output_dir, self.output_prefix)
        if not(filename.endswith('.h5')):
               filename += ".h5"

        if (len(self.projection_stack) == 0):
            mylog.debug("save_light_cone_stack: no projection data loaded.")
            return

        mylog.info("Writing light cone data to %s." % filename)

        output = h5py.File(filename, "a")

        node_exists = field_node in output.listnames()

        if node_exists:
            if over_write:
                mylog.info("Dataset, %s, already exists, overwriting." % field_node)
                write_data = True
                del output[field_node]
            else:
                mylog.info("Dataset, %s, already exists in %s, not saving." % (field_node, filename))
                write_data = False
        else:
            write_data = True

        if write_data:
            mylog.info("Saving %s to %s." % (field_node, filename))
            self.projection_stack = na.array(self.projection_stack)
            field_dataset = output.create_dataset(field_node, data=self.projection_stack)
            field_dataset.attrs['redshifts'] = redshiftList
            field_dataset.attrs['observer_redshift'] = na.float(self.observer_redshift)
            field_dataset.attrs['field_of_view_in_arcminutes'] = na.float(self.field_of_view_in_arcminutes)
            field_dataset.attrs['image_resolution_in_arcseconds'] = na.float(self.image_resolution_in_arcseconds)

        if (len(self.projection_weight_field_stack) > 0):
            if node_exists:
                if over_write:
                    mylog.info("Dataset, %s, already exists, overwriting." % weight_field_node)
                    del output[field_node]
                else:
                    mylog.info("Dataset, %s, already exists in %s, not saving." % (weight_field_node, filename))
                    write_data = False
            else:
                write_data = True

            if write_data:
                mylog.info("Saving %s to %s." % (weight_field_node, filename))
                self.projection_weight_field_stack = na.array(self.projection_weight_field_stack)
                weight_field_dataset = output.create_dataset(weight_field_node, data=self.projection_weight_field_stack)
                weight_field_dataset.attrs['redshifts'] = redshiftList
                weight_field_dataset.attrs['observer_redshift'] = na.float(self.observer_redshift)
                weight_field_dataset.attrs['field_of_view_in_arcminutes'] = na.float(self.field_of_view_in_arcminutes)
                weight_field_dataset.attrs['image_resolution_in_arcseconds'] = na.float(self.image_resolution_in_arcseconds)

        output.close()
