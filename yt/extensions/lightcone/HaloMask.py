"""
Light cone halo mask functions.

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

from yt.extensions.HaloProfiler import *
import yt.lagos as lagos
import copy
import numpy as na
import tables as h5

#### Note: assumption of box width 1.  I'll fix it someday.

def MakeLightConeHaloMask(lightCone,HaloMaskParameterFile,cube_file=None,mask_file=None):
    "Make a boolean mask to cut clusters out of light cone projections."

    pixels = int(lightCone.lightConeParameters['FieldOfViewInArcMinutes'] * 60.0 / \
            lightCone.lightConeParameters['ImageResolutionInArcSeconds'])

    lightConeMask = []

    # Loop through files in light cone solution and get virial quantities.
    for slice in lightCone.lightConeSolution:
        hp = HaloProfiler(slice['filename'],HaloMaskParameterFile)
        hp._LoadVirialData()

        lightConeMask.append(_MakeHaloMask(slice,hp.virialQuantities,pixels))

    # Write out cube of masks from each slice.
    if cube_file is not None:
        output = h5.openFile(cube_file,'a')
        output.createArray("/",'haloMaskCube',na.array(lightConeMask))
        output.close()

    # Write out final mask.
    if mask_file is not None:
        # Final mask is simply the product of the mask from each slice.
        finalMask = na.ones(shape=(pixels,pixels))
        for mask in lightConeMask:
            finalMask *= mask

        output = h5.openFile(mask_file,'a')
        output.createArray("/",'haloMask',na.array(finalMask))
        output.close()

    return lightConeMask

def _MakeHaloMask(slice,virialQuantities,pixels):
    "Make halo mask for one slice in light cone solution."

    # Make numpy arrays for halo centers and virial radii.
    halo_x = []
    halo_y = []
    halo_depth = []
    halo_radius = []

    # Get units to convert virial radii to code units.
    dataset_object = lagos.EnzoStaticOutput(slice['filename'])
    Mpc_units = dataset_object.units['mpc']
    del dataset_object

    for halo in virialQuantities:
        if halo is not None:
            center = copy.deepcopy(halo['center'])
            halo_depth.append(center.pop(slice['ProjectionAxis']))
            halo_x.append(center[0])
            halo_y.append(center[1])
            halo_radius.append(halo['RadiusMpc']/Mpc_units)

    halo_x = na.array(halo_x)
    halo_y = na.array(halo_y)
    halo_depth = na.array(halo_depth)
    halo_radius = na.array(halo_radius)

    # Adjust halo centers along line of sight.
    depthCenter = slice['ProjectionCenter'][slice['ProjectionAxis']]
    depthLeft = depthCenter - 0.5 * slice['DepthBoxFraction']
    depthRight = depthCenter + 0.5 * slice['DepthBoxFraction']

    # Make boolean mask to pick out centers in region along line of sight.
    # Halos near edges may wrap around to other side.
    add_left = (halo_depth + halo_radius) > 1 # should be box width
    add_right = (halo_depth - halo_radius) < 0

    halo_depth = na.concatenate([halo_depth,(halo_depth[add_left]-1),(halo_depth[add_right]+1)])
    halo_x = na.concatenate([halo_x,halo_x[add_left],halo_x[add_right]])
    halo_y = na.concatenate([halo_y,halo_y[add_left],halo_y[add_right]])
    halo_radius = na.concatenate([halo_radius,halo_radius[add_left],halo_radius[add_right]])

    del add_left,add_right

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
        del mask
    del halo_depth

    all_halo_x = na.array([])
    all_halo_y = na.array([])
    all_halo_radius = na.array([])

    # Tile halos of width box fraction is greater than one.
    # Copy original into offset positions to make tiles.
    for x in range(int(na.ceil(slice['WidthBoxFraction']))):
        for y in range(int(na.ceil(slice['WidthBoxFraction']))):
            all_halo_x = na.concatenate([all_halo_x,halo_x+x])
            all_halo_y = na.concatenate([all_halo_y,halo_y+y])
            all_halo_radius = na.concatenate([all_halo_radius,halo_radius])

    del halo_x,halo_y,halo_radius

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
    del add_x_right

    # Centers hanging off the left edge.
    add_x_left = all_halo_x - all_halo_radius < 0
    add2_x_halo_x = all_halo_x[add_x_left]
    add2_x_halo_x += na.ceil(slice['WidthBoxFraction'])
    add2_x_halo_y = all_halo_y[add_x_left]
    add2_x_halo_radius = all_halo_radius[add_x_left]
    del add_x_left

    # Centers hanging off the top edge.
    add_y_right = all_halo_y + all_halo_radius > na.ceil(slice['WidthBoxFraction'])
    add_y_halo_x = all_halo_x[add_y_right]
    add_y_halo_y = all_halo_y[add_y_right]
    add_y_halo_y -= na.ceil(slice['WidthBoxFraction'])
    add_y_halo_radius = all_halo_radius[add_y_right]
    del add_y_right

    # Centers hanging off the bottom edge.
    add_y_left = all_halo_y - all_halo_radius < 0
    add2_y_halo_x = all_halo_x[add_y_left]
    add2_y_halo_y = all_halo_y[add_y_left]
    add2_y_halo_y += na.ceil(slice['WidthBoxFraction'])
    add2_y_halo_radius = all_halo_radius[add_y_left]
    del add_y_left

    # Add the hanging centers back to the projection data.
    all_halo_x = na.concatenate([all_halo_x,add_x_halo_x,add2_x_halo_x,add_y_halo_x,add2_y_halo_x])
    all_halo_y = na.concatenate([all_halo_y,add_x_halo_y,add2_x_halo_y,add_y_halo_y,add2_y_halo_y])
    all_halo_radius = na.concatenate([all_halo_radius,add_x_halo_radius,add2_x_halo_radius,
                                      add_y_halo_radius,add2_y_halo_radius])

    del add_x_halo_x,add_x_halo_y,add_x_halo_radius
    del add2_x_halo_x,add2_x_halo_y,add2_x_halo_radius
    del add_y_halo_x,add_y_halo_y,add_y_halo_radius
    del add2_y_halo_x,add2_y_halo_y,add2_y_halo_radius

    # Cut edges to proper width.
    cut_mask = (all_halo_x - all_halo_radius < slice['WidthBoxFraction']) & \
        (all_halo_y - all_halo_radius < slice['WidthBoxFraction'])
    all_halo_x = all_halo_x[cut_mask]
    all_halo_y = all_halo_y[cut_mask]
    all_halo_radius = all_halo_radius[cut_mask]
    del cut_mask

    # Make boolean mask and cut out halos.
    dx = slice['WidthBoxFraction'] / pixels
    x = [(q + 0.5) * dx for q in range(pixels)]
    haloMask = na.ones(shape=(pixels,pixels),dtype=bool)

    # Cut out any pixel that has any part at all in the circle.
    for q in range(len(all_halo_radius)):
        dif_xIndex = na.array(int(all_halo_x[q]/dx) - na.array(range(pixels))) != 0
        dif_yIndex = na.array(int(all_halo_y[q]/dx) - na.array(range(pixels))) != 0

        xDistance = (na.abs(x - all_halo_x[q]) - (0.5 * dx)) * dif_xIndex
        yDistance = (na.abs(x - all_halo_y[q]) - (0.5 * dx)) * dif_yIndex

        distance = na.array([na.sqrt(w**2 + xDistance**2) for w in yDistance])
        haloMask *= (distance >= all_halo_radius[q])

    return haloMask
