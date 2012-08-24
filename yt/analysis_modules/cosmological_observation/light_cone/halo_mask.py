"""
Light cone halo mask functions.

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

import copy
import h5py
import numpy as np

from yt.funcs import *
from yt.analysis_modules.halo_profiler.api import \
     HaloProfiler
from yt.convenience import load
from yt.utilities.parallel_tools.parallel_analysis_interface import \
     parallel_objects, \
     parallel_root_only

def _light_cone_halo_mask(lightCone, cube_file=None,
                          mask_file=None, map_file=None,
                          halo_profiler_parameters=None,
                          virial_overdensity=200,
                          njobs=1, dynamic=False):
    "Make a boolean mask to cut clusters out of light cone projections."

    if halo_profiler_parameters is None:
        halo_profiler_parameters = {}

    pixels = int(lightCone.field_of_view_in_arcminutes * 60.0 /
                 lightCone.image_resolution_in_arcseconds)

    # Loop through files in light cone solution and get virial quantities.
    halo_map_storage = {}
    for my_storage, my_slice in \
      parallel_objects(lightCone.light_cone_solution,
                       njobs=njobs, dynamic=dynamic,
                       storage=halo_map_storage):
        halo_list = _get_halo_list(my_slice['filename'],
                                   **halo_profiler_parameters)
        my_storage.result = \
          {'mask': _make_slice_mask(my_slice, halo_list, pixels,
                                    virial_overdensity)}
        if map_file is not None:
            my_storage.result['map'] = \
              _make_slice_halo_map(my_slice, halo_list,
                                   virial_overdensity)

    # Reassemble halo mask and map lists.
    light_cone_mask = []
    halo_map = []
    all_slices = halo_map_storage.keys()
    all_slices.sort()
    for i in all_slices:
        light_cone_mask.append(halo_map_storage[i]['mask'])
        if map_file is not None:
            halo_map.extend(halo_map_storage[i]['map'])
    del halo_map_storage

    # Write out cube of masks from each slice.
    if cube_file is not None:
        _write_halo_mask(cube_file, np.array(light_cone_mask))

    # Write out a text list of all halos in the image.
    if map_file is not None:
        _write_halo_map(map_file, halo_map)

    # Write out final mask.
    if mask_file is not None:
        # Final mask is simply the product of the mask from each slice.
        final_mask = np.ones(shape=(pixels, pixels))
        for mask in light_cone_mask:
            final_mask *= mask
        _write_halo_mask(mask_file, final_mask)

    return light_cone_mask

@parallel_root_only
def _write_halo_mask(filename, halo_mask):
    r"""Write out an hdf5 file with the halo mask that
    can be applied to an image.
    """

    mylog.info("Saving halo mask to %s." % filename)
    output = h5py.File(filename, 'a')
    if 'HaloMask' in output.keys():
        del output['HaloMask']
    output.create_dataset('HaloMask', data=np.array(halo_mask))
    output.close()

@parallel_root_only
def _write_halo_map(filename, halo_map):
    "Write a text list of halos in a light cone image."

    mylog.info("Saving halo map to %s." % filename)
    f = open(filename, 'w')
    f.write("#z       x         y        r_image   r_mpc     m_Msun\n")
    for halo in halo_map:
        f.write("%7.4f %9.6f %9.6f %9.3e %9.3e %9.3e\n" % \
                    (halo['redshift'], halo['x'], halo['y'],
                     halo['radius_image'], halo['radius_mpc'],
                     halo['mass']))
    f.close()

def _get_halo_list(dataset, halo_profiler_kwargs=None,
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

def _make_slice_mask(slice, halo_list, pixels, virial_overdensity):
    "Make halo mask for one slice in light cone solution."

    # Get shifted, tiled halo list.
    all_halo_x, all_halo_y, \
      all_halo_radius, all_halo_mass = \
      _make_slice_halo_list(slice, halo_list, virial_overdensity)

    # Make boolean mask and cut out halos.
    dx = slice['box_width_fraction'] / pixels
    x = [(q + 0.5) * dx for q in range(pixels)]
    haloMask = np.ones(shape=(pixels, pixels), dtype=bool)

    # Cut out any pixel that has any part at all in the circle.
    for q in range(len(all_halo_radius)):
        dif_xIndex = np.array(int(all_halo_x[q]/dx) -
                              np.array(range(pixels))) != 0
        dif_yIndex = np.array(int(all_halo_y[q]/dx) -
                              np.array(range(pixels))) != 0

        xDistance = (np.abs(x - all_halo_x[q]) -
                     (0.5 * dx)) * dif_xIndex
        yDistance = (np.abs(x - all_halo_y[q]) -
                     (0.5 * dx)) * dif_yIndex

        distance = np.array([np.sqrt(w**2 + xDistance**2)
                             for w in yDistance])
        haloMask *= (distance >= all_halo_radius[q])

    return haloMask

def _make_slice_halo_map(slice, halo_list, virial_overdensity):
    "Make list of halos for one slice in light cone solution."

    # Get units to convert virial radii back to physical units.
    dataset_object = load(slice['filename'])
    Mpc_units = dataset_object.units['mpc']
    del dataset_object

    # Get shifted, tiled halo list.
    all_halo_x, all_halo_y, \
      all_halo_radius, all_halo_mass = \
      _make_slice_halo_list(slice, halo_list, virial_overdensity)

    # Construct list of halos
    halo_map = []

    for q in range(len(all_halo_x)):
        # Give radius in both physics units and
        # units of the image (0 to 1).
        radius_mpc = all_halo_radius[q] * Mpc_units
        radius_image = all_halo_radius[q] / slice['box_width_fraction']

        halo_map.append({'x': all_halo_x[q] / slice['box_width_fraction'],
                         'y': all_halo_y[q] / slice['box_width_fraction'],
                         'redshift': slice['redshift'],
                         'radius_mpc': radius_mpc,
                         'radius_image': radius_image,
                         'mass': all_halo_mass[q]})

    return halo_map

def _make_slice_halo_list(slice, halo_list, virial_overdensity):
    "Make shifted, tiled list of halos for halo mask and halo map."

   # Make numpy arrays for halo centers and virial radii.
    halo_x = []
    halo_y = []
    halo_depth = []
    halo_radius = []
    halo_mass = []

    # Get units to convert virial radii to code units.
    dataset_object = load(slice['filename'])
    Mpc_units = dataset_object.units['mpc']
    del dataset_object

    for halo in halo_list:
        if halo is not None:
            center = copy.deepcopy(halo['center'])
            halo_depth.append(center.pop(slice['projection_axis']))
            halo_x.append(center[0])
            halo_y.append(center[1])
            halo_radius.append(halo['RadiusMpc_%d' % virial_overdensity] /
                               Mpc_units)
            halo_mass.append(halo['TotalMassMsun_%d' % virial_overdensity])

    halo_x = np.array(halo_x)
    halo_y = np.array(halo_y)
    halo_depth = np.array(halo_depth)
    halo_radius = np.array(halo_radius)
    halo_mass = np.array(halo_mass)

    # Adjust halo centers along line of sight.
    depth_center = slice['projection_center'][slice['projection_axis']]
    depth_left = depth_center - 0.5 * slice['box_depth_fraction']
    depth_right = depth_center + 0.5 * slice['box_depth_fraction']

    # Make boolean mask to pick out centers in region along line of sight.
    # Halos near edges may wrap around to other side.
    add_left = (halo_depth + halo_radius) > 1 # should be box width
    add_right = (halo_depth - halo_radius) < 0

    halo_depth = np.concatenate([halo_depth,
                                 (halo_depth[add_left]-1),
                                 (halo_depth[add_right]+1)])
    halo_x = np.concatenate([halo_x, halo_x[add_left], halo_x[add_right]])
    halo_y = np.concatenate([halo_y, halo_y[add_left], halo_y[add_right]])
    halo_radius = np.concatenate([halo_radius,
                                  halo_radius[add_left],
                                  halo_radius[add_right]])
    halo_mass = np.concatenate([halo_mass,
                                halo_mass[add_left],
                                halo_mass[add_right]])

    del add_left, add_right

    # Cut out the halos outside the region of interest.
    if (slice['box_depth_fraction'] < 1):
        if (depth_left < 0):
            mask = ((halo_depth + halo_radius >= 0) &
                    (halo_depth - halo_radius <= depth_right)) | \
                ((halo_depth + halo_radius >= depth_left + 1) &
                 (halo_depth - halo_radius <= 1))
        elif (depth_right > 1):
            mask = ((halo_depth + halo_radius >= 0) &
                    (halo_depth - halo_radius <= depth_right - 1)) | \
                ((halo_depth + halo_radius >= depth_left) &
                 (halo_depth - halo_radius <= 1))
        else:
            mask = (halo_depth + halo_radius >= depth_left) & \
              (halo_depth - halo_radius <= depth_right)

        halo_x = halo_x[mask]
        halo_y = halo_y[mask]
        halo_radius = halo_radius[mask]
        halo_mass = halo_mass[mask]
        del mask
    del halo_depth

    all_halo_x = np.array([])
    all_halo_y = np.array([])
    all_halo_radius = np.array([])
    all_halo_mass = np.array([])

    # Tile halos of width box fraction is greater than one.
    # Copy original into offset positions to make tiles.
    for x in range(int(np.ceil(slice['box_width_fraction']))):
        for y in range(int(np.ceil(slice['box_width_fraction']))):
            all_halo_x = np.concatenate([all_halo_x, halo_x+x])
            all_halo_y = np.concatenate([all_halo_y, halo_y+y])
            all_halo_radius = np.concatenate([all_halo_radius, halo_radius])
            all_halo_mass = np.concatenate([all_halo_mass, halo_mass])

    del halo_x, halo_y, halo_radius, halo_mass

    # Shift centers laterally.
    offset = copy.deepcopy(slice['projection_center'])
    del offset[slice['projection_axis']]

    # Shift x and y positions.
    all_halo_x -= offset[0]
    all_halo_y -= offset[1]

    # Wrap off-edge centers back around to
    # other side (periodic boundary conditions).
    all_halo_x[all_halo_x < 0] += np.ceil(slice['box_width_fraction'])
    all_halo_y[all_halo_y < 0] += np.ceil(slice['box_width_fraction'])

    # After shifting, some centers have fractional coverage
    # on both sides of the box.
    # Find those centers and make copies to be placed on the other side.

    # Centers hanging off the right edge.
    add_x_right = all_halo_x + all_halo_radius > \
      np.ceil(slice['box_width_fraction'])
    add_x_halo_x = all_halo_x[add_x_right]
    add_x_halo_x -= np.ceil(slice['box_width_fraction'])
    add_x_halo_y = all_halo_y[add_x_right]
    add_x_halo_radius = all_halo_radius[add_x_right]
    add_x_halo_mass = all_halo_mass[add_x_right]
    del add_x_right

    # Centers hanging off the left edge.
    add_x_left = all_halo_x - all_halo_radius < 0
    add2_x_halo_x = all_halo_x[add_x_left]
    add2_x_halo_x += np.ceil(slice['box_width_fraction'])
    add2_x_halo_y = all_halo_y[add_x_left]
    add2_x_halo_radius = all_halo_radius[add_x_left]
    add2_x_halo_mass = all_halo_mass[add_x_left]
    del add_x_left

    # Centers hanging off the top edge.
    add_y_right = all_halo_y + all_halo_radius > \
      np.ceil(slice['box_width_fraction'])
    add_y_halo_x = all_halo_x[add_y_right]
    add_y_halo_y = all_halo_y[add_y_right]
    add_y_halo_y -= np.ceil(slice['box_width_fraction'])
    add_y_halo_radius = all_halo_radius[add_y_right]
    add_y_halo_mass = all_halo_mass[add_y_right]
    del add_y_right

    # Centers hanging off the bottom edge.
    add_y_left = all_halo_y - all_halo_radius < 0
    add2_y_halo_x = all_halo_x[add_y_left]
    add2_y_halo_y = all_halo_y[add_y_left]
    add2_y_halo_y += np.ceil(slice['box_width_fraction'])
    add2_y_halo_radius = all_halo_radius[add_y_left]
    add2_y_halo_mass = all_halo_mass[add_y_left]
    del add_y_left

    # Add the hanging centers back to the projection data.
    all_halo_x = np.concatenate([all_halo_x,
                                 add_x_halo_x, add2_x_halo_x,
                                 add_y_halo_x, add2_y_halo_x])
    all_halo_y = np.concatenate([all_halo_y,
                                 add_x_halo_y, add2_x_halo_y,
                                 add_y_halo_y, add2_y_halo_y])
    all_halo_radius = np.concatenate([all_halo_radius,
                                      add_x_halo_radius,
                                      add2_x_halo_radius,
                                      add_y_halo_radius,
                                      add2_y_halo_radius])
    all_halo_mass = np.concatenate([all_halo_mass,
                                    add_x_halo_mass,
                                    add2_x_halo_mass,
                                    add_y_halo_mass,
                                    add2_y_halo_mass])

    del add_x_halo_x, add_x_halo_y, add_x_halo_radius
    del add2_x_halo_x, add2_x_halo_y, add2_x_halo_radius
    del add_y_halo_x, add_y_halo_y, add_y_halo_radius
    del add2_y_halo_x, add2_y_halo_y, add2_y_halo_radius

    # Cut edges to proper width.
    cut_mask = (all_halo_x - all_halo_radius <
                slice['box_width_fraction']) & \
        (all_halo_y - all_halo_radius <
         slice['box_width_fraction'])
    all_halo_x = all_halo_x[cut_mask]
    all_halo_y = all_halo_y[cut_mask]
    all_halo_radius = all_halo_radius[cut_mask]
    all_halo_mass = all_halo_mass[cut_mask]
    del cut_mask

    return (all_halo_x, all_halo_y,
            all_halo_radius, all_halo_mass)
