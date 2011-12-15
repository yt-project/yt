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

from yt.analysis_modules.halo_profiler.api import HaloProfiler
from yt.convenience import load
from yt.utilities.logger import ytLogger as mylog
from yt.config import ytcfg
import copy
import numpy as na
import h5py

def light_cone_halo_mask(lightCone, cube_file=None, mask_file=None, **kwargs):
    "Make a boolean mask to cut clusters out of light cone projections."

    pixels = int(lightCone.field_of_view_in_arcminutes * 60.0 /
                 lightCone.image_resolution_in_arcseconds)

    light_cone_mask = []

    # Loop through files in light cone solution and get virial quantities.
    for slice in lightCone.light_cone_solution:
        halo_list = _get_halo_list(slice['filename'], **kwargs)
        light_cone_mask.append(_make_slice_mask(slice, halo_list, pixels))

    # Write out cube of masks from each slice.
    if cube_file is not None and ytcfg.getint("yt", "__topcomm_parallel_rank") == 0:
        mylog.info("Saving halo mask cube to %s." % cube_file)
        output = h5py.File(cube_file, 'a')
        output.create_dataset('haloMaskCube', data=na.array(light_cone_mask))
        output.close()

    # Write out final mask.
    if mask_file is not None and ytcfg.getint("yt", "__topcomm_parallel_rank") == 0:
        # Final mask is simply the product of the mask from each slice.
        mylog.info("Saving halo mask to %s." % mask_file)
        finalMask = na.ones(shape=(pixels, pixels))
        for mask in light_cone_mask:
            finalMask *= mask

        output = h5py.File(mask_file, 'a')
        output.create_dataset('HaloMask', data=na.array(finalMask))
        output.close()

    return light_cone_mask

def light_cone_halo_map(lightCone, map_file='halo_map.out', **kwargs):
    "Make a text list of location of halos in a light cone image with virial quantities."

    haloMap = []

    # Loop through files in light cone solution and get virial quantities.
    for slice in lightCone.light_cone_solution:
        halo_list = _get_halo_list(slice['filename'], **kwargs)
        haloMap.extend(_make_slice_halo_map(slice, halo_list))

    # Write out file.
    if ytcfg.getint("yt", "__topcomm_parallel_rank") == 0:
        mylog.info("Saving halo map to %s." % map_file)
        f = open(map_file, 'w')
        f.write("#z       x         y        M [Msun]  R [Mpc]   R [image]\n")
        for halo in haloMap:
            f.write("%7.4f %9.6f %9.6f %9.3e %9.3e %9.3e\n" % \
                        (halo['redshift'], halo['x'], halo['y'],
                         halo['mass'], halo['radiusMpc'], halo['radiusImage']))
        f.close()

def _get_halo_list(dataset, halo_profiler_kwargs=None, halo_profiler_actions=None, halo_list='all'):
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

def _make_slice_mask(slice, halo_list, pixels):
    "Make halo mask for one slice in light cone solution."

    # Get shifted, tiled halo list.
    all_halo_x, all_halo_y, all_halo_radius, all_halo_mass = _make_slice_halo_list(slice, halo_list)
 
    # Make boolean mask and cut out halos.
    dx = slice['WidthBoxFraction'] / pixels
    x = [(q + 0.5) * dx for q in range(pixels)]
    haloMask = na.ones(shape=(pixels, pixels), dtype=bool)

    # Cut out any pixel that has any part at all in the circle.
    for q in range(len(all_halo_radius)):
        dif_xIndex = na.array(int(all_halo_x[q]/dx) - na.array(range(pixels))) != 0
        dif_yIndex = na.array(int(all_halo_y[q]/dx) - na.array(range(pixels))) != 0

        xDistance = (na.abs(x - all_halo_x[q]) - (0.5 * dx)) * dif_xIndex
        yDistance = (na.abs(x - all_halo_y[q]) - (0.5 * dx)) * dif_yIndex

        distance = na.array([na.sqrt(w**2 + xDistance**2) for w in yDistance])
        haloMask *= (distance >= all_halo_radius[q])

    return haloMask

def _make_slice_halo_map(slice, halo_list):
    "Make list of halos for one slice in light cone solution."

    # Get units to convert virial radii back to physical units.
    dataset_object = load(slice['filename'])
    Mpc_units = dataset_object.units['mpc']
    del dataset_object

    # Get shifted, tiled halo list.
    all_halo_x, all_halo_y, all_halo_radius, all_halo_mass = _make_slice_halo_list(slice, halo_list)

    # Construct list of halos
    haloMap = []

    for q in range(len(all_halo_x)):
        # Give radius in both physics units and units of the image (0 to 1).
        radiusMpc = all_halo_radius[q] * Mpc_units
        radiusImage = all_halo_radius[q] / slice['WidthBoxFraction']

        haloMap.append({'x':all_halo_x[q]/slice['WidthBoxFraction'],
                        'y':all_halo_y[q]/slice['WidthBoxFraction'],
                        'redshift':slice['redshift'],
                        'radiusMpc':radiusMpc,
                        'radiusImage':radiusImage,
                        'mass':all_halo_mass[q]})

    return haloMap

def _make_slice_halo_list(slice, halo_list):
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
            halo_depth.append(center.pop(slice['ProjectionAxis']))
            halo_x.append(center[0])
            halo_y.append(center[1])
            halo_radius.append(halo['RadiusMpc']/Mpc_units)
            halo_mass.append(halo['TotalMassMsun'])

    halo_x = na.array(halo_x)
    halo_y = na.array(halo_y)
    halo_depth = na.array(halo_depth)
    halo_radius = na.array(halo_radius)
    halo_mass = na.array(halo_mass)

    # Adjust halo centers along line of sight.
    depthCenter = slice['ProjectionCenter'][slice['ProjectionAxis']]
    depthLeft = depthCenter - 0.5 * slice['DepthBoxFraction']
    depthRight = depthCenter + 0.5 * slice['DepthBoxFraction']

    # Make boolean mask to pick out centers in region along line of sight.
    # Halos near edges may wrap around to other side.
    add_left = (halo_depth + halo_radius) > 1 # should be box width
    add_right = (halo_depth - halo_radius) < 0

    halo_depth = na.concatenate([halo_depth, (halo_depth[add_left]-1), (halo_depth[add_right]+1)])
    halo_x = na.concatenate([halo_x, halo_x[add_left], halo_x[add_right]])
    halo_y = na.concatenate([halo_y, halo_y[add_left], halo_y[add_right]])
    halo_radius = na.concatenate([halo_radius, halo_radius[add_left], halo_radius[add_right]])
    halo_mass = na.concatenate([halo_mass, halo_mass[add_left], halo_mass[add_right]])

    del add_left, add_right

    # Cut out the halos outside the region of interest.
    if (slice['DepthBoxFraction'] < 1):
        if (depthLeft < 0):
            mask = ((halo_depth + halo_radius >= 0) & (halo_depth - halo_radius <= depthRight)) | \
                ((halo_depth + halo_radius >= depthLeft + 1) & (halo_depth - halo_radius <= 1))
        elif (depthRight > 1):
            mask = ((halo_depth + halo_radius >= 0) & (halo_depth - halo_radius <= depthRight - 1)) | \
                ((halo_depth + halo_radius >= depthLeft) & (halo_depth - halo_radius <= 1))
        else:
            mask = (halo_depth + halo_radius >= depthLeft) & (halo_depth - halo_radius <= depthRight)

        halo_x = halo_x[mask]
        halo_y = halo_y[mask]
        halo_radius = halo_radius[mask]
        halo_mass = halo_mass[mask]
        del mask
    del halo_depth

    all_halo_x = na.array([])
    all_halo_y = na.array([])
    all_halo_radius = na.array([])
    all_halo_mass = na.array([])

    # Tile halos of width box fraction is greater than one.
    # Copy original into offset positions to make tiles.
    for x in range(int(na.ceil(slice['WidthBoxFraction']))):
        for y in range(int(na.ceil(slice['WidthBoxFraction']))):
            all_halo_x = na.concatenate([all_halo_x, halo_x+x])
            all_halo_y = na.concatenate([all_halo_y, halo_y+y])
            all_halo_radius = na.concatenate([all_halo_radius, halo_radius])
            all_halo_mass = na.concatenate([all_halo_mass, halo_mass])

    del halo_x, halo_y, halo_radius, halo_mass

    # Shift centers laterally.
    offset = copy.deepcopy(slice['ProjectionCenter'])
    del offset[slice['ProjectionAxis']]

    # Shift x and y positions.
    all_halo_x -= offset[0]
    all_halo_y -= offset[1]

    # Wrap off-edge centers back around to other side (periodic boundary conditions).
    all_halo_x[all_halo_x < 0] += na.ceil(slice['WidthBoxFraction'])
    all_halo_y[all_halo_y < 0] += na.ceil(slice['WidthBoxFraction'])

    # After shifting, some centers have fractional coverage on both sides of the box.
    # Find those centers and make copies to be placed on the other side.

    # Centers hanging off the right edge.
    add_x_right = all_halo_x + all_halo_radius > na.ceil(slice['WidthBoxFraction'])
    add_x_halo_x = all_halo_x[add_x_right]
    add_x_halo_x -= na.ceil(slice['WidthBoxFraction'])
    add_x_halo_y = all_halo_y[add_x_right]
    add_x_halo_radius = all_halo_radius[add_x_right]
    add_x_halo_mass = all_halo_mass[add_x_right]
    del add_x_right

    # Centers hanging off the left edge.
    add_x_left = all_halo_x - all_halo_radius < 0
    add2_x_halo_x = all_halo_x[add_x_left]
    add2_x_halo_x += na.ceil(slice['WidthBoxFraction'])
    add2_x_halo_y = all_halo_y[add_x_left]
    add2_x_halo_radius = all_halo_radius[add_x_left]
    add2_x_halo_mass = all_halo_mass[add_x_left]
    del add_x_left

    # Centers hanging off the top edge.
    add_y_right = all_halo_y + all_halo_radius > na.ceil(slice['WidthBoxFraction'])
    add_y_halo_x = all_halo_x[add_y_right]
    add_y_halo_y = all_halo_y[add_y_right]
    add_y_halo_y -= na.ceil(slice['WidthBoxFraction'])
    add_y_halo_radius = all_halo_radius[add_y_right]
    add_y_halo_mass = all_halo_mass[add_y_right]
    del add_y_right

    # Centers hanging off the bottom edge.
    add_y_left = all_halo_y - all_halo_radius < 0
    add2_y_halo_x = all_halo_x[add_y_left]
    add2_y_halo_y = all_halo_y[add_y_left]
    add2_y_halo_y += na.ceil(slice['WidthBoxFraction'])
    add2_y_halo_radius = all_halo_radius[add_y_left]
    add2_y_halo_mass = all_halo_mass[add_y_left]
    del add_y_left

    # Add the hanging centers back to the projection data.
    all_halo_x = na.concatenate([all_halo_x, add_x_halo_x, add2_x_halo_x, add_y_halo_x, add2_y_halo_x])
    all_halo_y = na.concatenate([all_halo_y, add_x_halo_y, add2_x_halo_y, add_y_halo_y, add2_y_halo_y])
    all_halo_radius = na.concatenate([all_halo_radius, add_x_halo_radius, add2_x_halo_radius,
                                      add_y_halo_radius, add2_y_halo_radius])
    all_halo_mass = na.concatenate([all_halo_mass, add_x_halo_mass, add2_x_halo_mass,
                                    add_y_halo_mass, add2_y_halo_mass])

    del add_x_halo_x, add_x_halo_y, add_x_halo_radius
    del add2_x_halo_x, add2_x_halo_y, add2_x_halo_radius
    del add_y_halo_x, add_y_halo_y, add_y_halo_radius
    del add2_y_halo_x, add2_y_halo_y, add2_y_halo_radius

    # Cut edges to proper width.
    cut_mask = (all_halo_x - all_halo_radius < slice['WidthBoxFraction']) & \
        (all_halo_y - all_halo_radius < slice['WidthBoxFraction'])
    all_halo_x = all_halo_x[cut_mask]
    all_halo_y = all_halo_y[cut_mask]
    all_halo_radius = all_halo_radius[cut_mask]
    all_halo_mass = all_halo_mass[cut_mask]
    del cut_mask

    return (all_halo_x, all_halo_y, all_halo_radius, all_halo_mass)
