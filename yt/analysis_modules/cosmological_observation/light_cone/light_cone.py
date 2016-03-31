"""
LightCone class and member functions.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import numpy as np
import os

from yt.config import \
    ytcfg
from yt.funcs import \
    mylog, \
    only_on_root
from yt.analysis_modules.cosmological_observation.cosmology_splice import \
    CosmologySplice
from yt.convenience import \
    load
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_objects, \
    parallel_root_only
from yt.visualization.image_writer import \
    write_image
from yt.units.yt_array import \
    YTArray
from .light_cone_projection import \
    _light_cone_projection

class LightCone(CosmologySplice):
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
    time_data : bool
        Whether or not to include time outputs when gathering
        datasets for time series.
        Default: True.
    redshift_data : bool
        Whether or not to include redshift outputs when gathering
        datasets for time series.
        Default: True.
    find_outputs : bool
        Whether or not to search for datasets in the current 
        directory.
        Default: False.
    set_parameters : dict
        Dictionary of parameters to attach to ds.parameters.
        Default: None.
    output_dir : string
        The directory in which images and data files will be written.
        Default: "LC".
    output_prefix : string
        The prefix of all images and data files.
        Default: "LightCone".

    """
    def __init__(self, parameter_filename, simulation_type,
                 near_redshift, far_redshift,
                 observer_redshift=0.0,
                 use_minimum_datasets=True, deltaz_min=0.0,
                 minimum_coherent_box_fraction=0.0,
                 time_data=True, redshift_data=True,
                 find_outputs=False, set_parameters=None,
                 output_dir="LC", output_prefix="LightCone"):

        self.near_redshift = near_redshift
        self.far_redshift = far_redshift
        self.observer_redshift = observer_redshift
        self.use_minimum_datasets = use_minimum_datasets
        self.deltaz_min = deltaz_min
        self.minimum_coherent_box_fraction = minimum_coherent_box_fraction
        if set_parameters is None:
            self.set_parameters = {}
        else:
            self.set_parameters = set_parameters
        self.output_dir = output_dir
        self.output_prefix = output_prefix

        # Create output directory.
        if not os.path.exists(self.output_dir):
            only_on_root(os.mkdir, self.output_dir)

        # Calculate light cone solution.
        CosmologySplice.__init__(self, parameter_filename, simulation_type,
                                 find_outputs=find_outputs)
        self.light_cone_solution = \
          self.create_cosmology_splice(self.near_redshift, self.far_redshift,
                                       minimal=self.use_minimum_datasets,
                                       deltaz_min=self.deltaz_min,
                                       time_data=time_data,
                                       redshift_data=redshift_data)

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

        # Don"t use box coherence with maximum projection depths.
        if self.use_minimum_datasets and \
                self.minimum_coherent_box_fraction > 0:
            mylog.info("Setting minimum_coherent_box_fraction to 0 with " +
                       "minimal light cone.")
            self.minimum_coherent_box_fraction = 0

        # Calculate projection sizes, and get
        # random projection axes and centers.
        seed = int(seed)
        np.random.seed(seed)

        # For box coherence, keep track of effective depth travelled.
        box_fraction_used = 0.0

        for q in range(len(self.light_cone_solution)):
            if "previous" in self.light_cone_solution[q]:
                del self.light_cone_solution[q]["previous"]
            if "next" in self.light_cone_solution[q]:
                del self.light_cone_solution[q]["next"]
            if q == len(self.light_cone_solution) - 1:
                z_next = self.near_redshift
            else:
                z_next = self.light_cone_solution[q+1]["redshift"]

            # Calculate fraction of box required for a depth of delta z
            self.light_cone_solution[q]["box_depth_fraction"] = \
                (self.cosmology.comoving_radial_distance(z_next, \
                        self.light_cone_solution[q]["redshift"]) / \
                        self.simulation.box_size).in_units("")

            # Calculate fraction of box required for width corresponding to
            # requested image size.
            proper_box_size = self.simulation.box_size / \
              (1.0 + self.light_cone_solution[q]["redshift"])
            self.light_cone_solution[q]["box_width_per_angle"] = \
              (self.cosmology.angular_scale(self.observer_redshift,
               self.light_cone_solution[q]["redshift"]) /
               proper_box_size).in_units("1 / degree")

            # Simple error check to make sure more than 100% of box depth
            # is never required.
            if self.light_cone_solution[q]["box_depth_fraction"] > 1.0:
                mylog.error(("Warning: box fraction required to go from " +
                             "z = %f to %f is %f") %
                            (self.light_cone_solution[q]["redshift"], z_next,
                             self.light_cone_solution[q]["box_depth_fraction"]))
                mylog.error(("Full box delta z is %f, but it is %f to the " +
                             "next data dump.") %
                            (self.light_cone_solution[q]["dz_max"],
                             self.light_cone_solution[q]["redshift"]-z_next))

            # Get projection axis and center.
            # If using box coherence, only get random axis and center if enough
            # of the box has been used, or if box_fraction_used will be greater
            # than 1 after this slice.
            if (q == 0) or (self.minimum_coherent_box_fraction == 0) or \
              (box_fraction_used > self.minimum_coherent_box_fraction) or \
              (box_fraction_used +
               self.light_cone_solution[q]["box_depth_fraction"] > 1.0):
                # Random axis and center.
                self.light_cone_solution[q]["projection_axis"] = \
                  np.random.randint(0, 3)
                self.light_cone_solution[q]["projection_center"] = \
                  np.random.random(3)
                box_fraction_used = 0.0
            else:
                # Same axis and center as previous slice,
                # but with depth center shifted.
                self.light_cone_solution[q]["projection_axis"] = \
                  self.light_cone_solution[q-1]["projection_axis"]
                self.light_cone_solution[q]["projection_center"] = \
                  self.light_cone_solution[q-1]["projection_center"].copy()
                self.light_cone_solution[q]["projection_center"]\
                  [self.light_cone_solution[q]["projection_axis"]] += \
                    0.5 * (self.light_cone_solution[q]["box_depth_fraction"] +
                           self.light_cone_solution[q-1]["box_depth_fraction"])
                if self.light_cone_solution[q]["projection_center"]\
                  [self.light_cone_solution[q]["projection_axis"]] >= 1.0:
                    self.light_cone_solution[q]["projection_center"]\
                      [self.light_cone_solution[q]["projection_axis"]] -= 1.0

            box_fraction_used += self.light_cone_solution[q]["box_depth_fraction"]

        # Write solution to a file.
        if filename is not None:
            self._save_light_cone_solution(filename=filename)

    def project_light_cone(self, field_of_view, image_resolution, field,
                           weight_field=None, photon_field=False,
                           save_stack=True, save_final_image=True,
                           save_slice_images=False,
                           cmap_name=None,
                           njobs=1, dynamic=False):
        r"""Create projections for light cone, then add them together.

        Parameters
        ----------
        field_of_view : YTQuantity or tuple of (float, str)
            The field of view of the image and the units.
        image_resolution : YTQuantity or tuple of (float, str)
            The size of each image pixel and the units.
        field : string
            The projected field.
        weight_field : string
            the weight field of the projection.  This has the same meaning as
            in standard projections.
            Default: None.
        photon_field : bool
            if True, the projection data for each slice is decremented by 4 Pi
            R^2`, where R is the luminosity distance between the observer and
            the slice redshift.
            Default: False.
        save_stack : bool
            if True, the light cone data including each individual
            slice is written to an hdf5 file.
            Default: True.
        save_final_image : bool
            if True, save an image of the final light cone projection.
            Default: True.
        save_slice_images : bool
            save images for each individual projection slice.
            Default: False.
        cmap_name : string
            color map for images.
            Default: your default colormap.
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

        if cmap_name is None:
            cmap_name = ytcfg.get("yt", "default_colormap")

        if isinstance(field_of_view, tuple) and len(field_of_view) == 2:
            field_of_view = self.simulation.quan(field_of_view[0],
                                                 field_of_view[1])
        elif not isinstance(field_of_view, YTArray):
          raise RuntimeError("field_of_view argument must be either a YTQauntity " +
                             "or a tuple of type (float, str).")
        if isinstance(image_resolution, tuple) and len(image_resolution) == 2:
            image_resolution = self.simulation.quan(image_resolution[0],
                                                    image_resolution[1])
        elif not isinstance(image_resolution, YTArray):
          raise RuntimeError("image_resolution argument must be either a YTQauntity " +
                             "or a tuple of type (float, str).")
        
        # Calculate number of pixels on a side.
        pixels = (field_of_view / image_resolution).in_units("")

        # Clear projection stack.
        projection_stack = []
        projection_weight_stack = []
        if "object" in self.light_cone_solution[-1]:
            del self.light_cone_solution[-1]["object"]

        # for q, output in enumerate(self.light_cone_solution):
        all_storage = {}
        for my_storage, output in parallel_objects(self.light_cone_solution,
                                                   storage=all_storage,
                                                   dynamic=dynamic):
            output["object"] = load(output["filename"])
            output["object"].parameters.update(self.set_parameters)

            # Calculate fraction of box required for width corresponding to
            # requested image size.
            proper_box_size = self.simulation.box_size / \
              (1.0 + output["redshift"])
            output["box_width_fraction"] = (output["box_width_per_angle"] *
                                            field_of_view).in_units("")
            
            frb = _light_cone_projection(output, field, pixels,
                                         weight_field=weight_field)

            if photon_field:
                # Decrement the flux by the luminosity distance.
                # Assume field in frb is in erg/s/cm^2/Hz
                dL = self.cosmology.luminosity_distance(self.observer_redshift,
                                                        output["redshift"])
                proper_box_size = self.simulation.box_size / \
                  (1.0 + output["redshift"])
                pixel_area = (proper_box_size.in_cgs() / pixels)**2 #in proper cm^2
                factor = pixel_area / (4.0 * np.pi * dL.in_cgs()**2)
                mylog.info("Distance to slice = %s" % dL)
                frb[field] *= factor #in erg/s/cm^2/Hz on observer"s image plane.

            if weight_field is None:
                my_storage.result = {"field": frb[field]}
            else:
                my_storage.result = {"field": (frb[field] *
                                               frb["weight_field"]),
                                     "weight_field": frb["weight_field"]}

            del output["object"]

        # Combine results from each slice.
        all_slices = list(all_storage.keys())
        all_slices.sort()
        for my_slice in all_slices:
            if save_slice_images:
                name = os.path.join(self.output_dir,
                                    "%s_%04d_%04d" %
                                    (self.output_prefix,
                                     my_slice, len(self.light_cone_solution)))
                if weight_field is None:
                    my_image = all_storage[my_slice]["field"]
                else:
                    my_image = all_storage[my_slice]["field"] / \
                      all_storage[my_slice]["weight_field"]
                only_on_root(write_image, np.log10(my_image),
                             "%s_%s.png" % (name, field), cmap_name=cmap_name)

            projection_stack.append(all_storage[my_slice]["field"])
            if weight_field is not None:
                projection_weight_stack.append(all_storage[my_slice]["field"])

        projection_stack = self.simulation.arr(projection_stack)
        projection_weight_stack = self.simulation.arr(projection_weight_stack)
                
        # Add up slices to make light cone projection.
        if (weight_field is None):
            light_cone_projection = projection_stack.sum(axis=0)
        else:
            light_cone_projection = \
              projection_stack.sum(axis=0) / \
              self.simulation.arr(projection_weight_stack).sum(axis=0)

        filename = os.path.join(self.output_dir, self.output_prefix)

        # Write image.
        if save_final_image:
            only_on_root(write_image, np.log10(light_cone_projection),
                         "%s_%s.png" % (filename, field), cmap_name=cmap_name)

        # Write stack to hdf5 file.
        if save_stack:
            self._save_light_cone_stack(field, weight_field,
                projection_stack, projection_weight_stack,
                filename=filename,
                attrs={"field_of_view": str(field_of_view),
                       "image_resolution": str(image_resolution)})

    @parallel_root_only
    def _save_light_cone_solution(self, filename="light_cone.dat"):
        "Write out a text file with information on light cone solution."

        mylog.info("Saving light cone solution to %s." % filename)

        f = open(filename, "w")
        f.write("# parameter_filename = %s\n" % self.parameter_filename)
        f.write("\n")
        f.write("# Slice    Dataset    Redshift    depth/box    " + \
                "width/degree    axis    center\n")
        for q, output in enumerate(self.light_cone_solution):
            f.write(("%04d %s %f %f %f %d %f %f %f\n") %
                    (q, output["filename"], output["redshift"],
                     output["box_depth_fraction"], output["box_width_per_angle"],
                     output["projection_axis"], output["projection_center"][0],
                     output["projection_center"][1], output["projection_center"][2]))
        f.close()

    @parallel_root_only
    def _save_light_cone_stack(self, field, weight_field,
                               pstack, wstack,
                               filename=None, attrs=None):
        "Save the light cone projection stack as a 3d array in and hdf5 file."

        if attrs is None:
            attrs = {}
        
        # Make list of redshifts to include as a dataset attribute.
        redshift_list = np.array([my_slice["redshift"] \
                                 for my_slice in self.light_cone_solution])

        field_node = "%s_%s" % (field, weight_field)
        weight_field_node = "weight_field_%s" % weight_field

        if (filename is None):
            filename = os.path.join(self.output_dir, "%s_data" % self.output_prefix)
        if not(filename.endswith(".h5")):
               filename += ".h5"

        if pstack.size == 0:
            mylog.info("save_light_cone_stack: light cone projection is empty.")
            return

        mylog.info("Writing light cone data to %s." % filename)

        fh = h5py.File(filename, "a")

        if field_node in fh:
            del fh[field_node]

        mylog.info("Saving %s to %s." % (field_node, filename))
        dataset = fh.create_dataset(field_node,
                                          data=pstack)
        dataset.attrs["units"] = str(pstack.units)
        dataset.attrs["redshifts"] = redshift_list
        dataset.attrs["observer_redshift"] = np.float(self.observer_redshift)
        for key, value in attrs.items():
            dataset.attrs[key] = value

        if wstack.size > 0:
            if weight_field_node in fh:
                del fh[weight_field_node]

            mylog.info("Saving %s to %s." % (weight_field_node, filename))
            dataset = fh.create_dataset(weight_field_node,
                                        data=wstack)
            dataset.attrs["units"] = str(wstack.units)
            dataset.attrs["redshifts"] = redshift_list
            dataset.attrs["observer_redshift"] = np.float(self.observer_redshift)
            for key, value in attrs.items():
                dataset.attrs[key] = value

        fh.close()
