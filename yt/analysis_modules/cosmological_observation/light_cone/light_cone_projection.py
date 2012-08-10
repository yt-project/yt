"""
Create randomly centered, tiled projections to be used in light cones.

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

import numpy as na
import copy

from yt.funcs import *
from yt.visualization.fixed_resolution import \
    FixedResolutionBuffer
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_blocking_call

@parallel_blocking_call
def _light_cone_projection(lightConeSlice, field, pixels, weight_field=None,
                           save_image=False, node=None, field_cuts=None):
    "Create a single projection to be added into the light cone stack."

    # Use some projection parameters to seed random number generator to make unique node name.
    # We are just saving the projection object, so only the projection axis needs to be considered
    # since the lateral shifting and tiling occurs after the projection object is made.
    # Likewise, only the box_depth_fraction needs to be considered.

    # projection_axis
    # projection_center[projection_axis]
    # box_depth_fraction

    # Name node with user specified keyword if given with 'node' keyword.
    node_name = "LightCone_%s_%d_%f_%f" % (node, lightConeSlice['projection_axis'],
                                           lightConeSlice['projection_center'][lightConeSlice['projection_axis']],
                                           lightConeSlice['box_depth_fraction'])

    mylog.info("Making projection at z = %f from %s." % (lightConeSlice['redshift'], lightConeSlice['filename']))

    region_center = [0.5 * (lightConeSlice['object'].domain_right_edge[q] +
                            lightConeSlice['object'].domain_left_edge[q]) \
                         for q in range(lightConeSlice['object'].dimensionality)]

    # 1. The Depth Problem
    # Use coordinate field cut in line of sight to cut projection to proper depth.
    if field_cuts is None:
        these_field_cuts = []
    else:
        these_field_cuts = copy.deepcopy(field_cuts)

    if (lightConeSlice['box_depth_fraction'] < 1):
        axis = ('x', 'y', 'z')[lightConeSlice['projection_axis']]
        depthLeft = lightConeSlice['projection_center'][lightConeSlice['projection_axis']] \
            - 0.5 * lightConeSlice['box_depth_fraction']
        depthRight = lightConeSlice['projection_center'][lightConeSlice['projection_axis']] \
            + 0.5 * lightConeSlice['box_depth_fraction']
        if (depthLeft < 0):
            cut_mask = "((grid[\"%s\"] + 0.5*grid[\"d%s\"] >= 0) & (grid[\"%s\"] - 0.5*grid[\"d%s\"] <= %f)) | ((grid[\"%s\"] + 0.5*grid[\"d%s\"] >= %f) & (grid[\"%s\"] - 0.5*grid[\"d%s\"] <= 1))" % \
                (axis, axis, axis, axis, depthRight, axis, axis, (depthLeft+1), axis, axis)
        elif (depthRight > 1):
            cut_mask = "((grid[\"%s\"] + 0.5*grid[\"d%s\"] >= 0) & (grid[\"%s\"] - 0.5*grid[\"d%s\"] <= %f)) | ((grid[\"%s\"] + 0.5*grid[\"d%s\"] >= %f) & (grid[\"%s\"] - 0.5*grid[\"d%s\"] <= 1))" % \
                (axis, axis, axis, axis, (depthRight-1), axis, axis, depthLeft, axis, axis)
        else:
            cut_mask = "(grid[\"%s\"] + 0.5*grid[\"d%s\"] >= %f) & (grid[\"%s\"] - 0.5*grid[\"%s\"] <= %f)" % (axis, axis, depthLeft, axis, axis, depthRight)

        these_field_cuts.append(cut_mask)

    # Make projection.
    proj = lightConeSlice['object'].h.proj(field, lightConeSlice['projection_axis'], 
                                           weight_field, center=region_center,
                                           field_cuts=these_field_cuts, node_name=node_name)

    # 2. The Tile Problem
    # Tile projection to specified width.

    # Original projection data.
    original_px = copy.deepcopy(proj['px'])
    original_py = copy.deepcopy(proj['py'])
    original_pdx = copy.deepcopy(proj['pdx'])
    original_pdy = copy.deepcopy(proj['pdy'])
    original_field = copy.deepcopy(proj[field])
    original_weight_field = copy.deepcopy(proj['weight_field'])

    # Copy original into offset positions to make tiles.
    for x in range(int(na.ceil(lightConeSlice['box_width_fraction']))):
        for y in range(int(na.ceil(lightConeSlice['box_width_fraction']))):
            if ((x + y) > 0):
                proj['px'] = na.concatenate([proj['px'], original_px+x])
                proj['py'] = na.concatenate([proj['py'], original_py+y])
                proj['pdx'] = na.concatenate([proj['pdx'], original_pdx])
                proj['pdy'] = na.concatenate([proj['pdy'], original_pdy])
                proj[field] = na.concatenate([proj[field], original_field])
                proj['weight_field'] = na.concatenate([proj['weight_field'],
                                                       original_weight_field])

    # Delete originals.
    del original_px
    del original_py
    del original_pdx
    del original_pdy
    del original_field
    del original_weight_field

    # 3. The Shift Problem
    # Shift projection by random x and y offsets.

    offset = copy.deepcopy(lightConeSlice['projection_center'])
    # Delete depth coordinate.
    del offset[lightConeSlice['projection_axis']]

    # Shift x and y positions.
    proj['px'] -= offset[0]
    proj['py'] -= offset[1]

    # Wrap off-edge cells back around to other side (periodic boundary conditions).
    proj['px'][proj['px'] < 0] += na.ceil(lightConeSlice['box_width_fraction'])
    proj['py'][proj['py'] < 0] += na.ceil(lightConeSlice['box_width_fraction'])

    # After shifting, some cells have fractional coverage on both sides of the box.
    # Find those cells and make copies to be placed on the other side.

    # Cells hanging off the right edge.
    add_x_right = proj['px'] + 0.5 * proj['pdx'] > \
      na.ceil(lightConeSlice['box_width_fraction'])
    add_x_px = proj['px'][add_x_right]
    add_x_px -= na.ceil(lightConeSlice['box_width_fraction'])
    add_x_py = proj['py'][add_x_right]
    add_x_pdx = proj['pdx'][add_x_right]
    add_x_pdy = proj['pdy'][add_x_right]
    add_x_field = proj[field][add_x_right]
    add_x_weight_field = proj['weight_field'][add_x_right]
    del add_x_right

    # Cells hanging off the left edge.
    add_x_left = proj['px'] - 0.5 * proj['pdx'] < 0
    add2_x_px = proj['px'][add_x_left]
    add2_x_px += na.ceil(lightConeSlice['box_width_fraction'])
    add2_x_py = proj['py'][add_x_left]
    add2_x_pdx = proj['pdx'][add_x_left]
    add2_x_pdy = proj['pdy'][add_x_left]
    add2_x_field = proj[field][add_x_left]
    add2_x_weight_field = proj['weight_field'][add_x_left]
    del add_x_left

    # Cells hanging off the top edge.
    add_y_right = proj['py'] + 0.5 * proj['pdy'] > \
      na.ceil(lightConeSlice['box_width_fraction'])
    add_y_px = proj['px'][add_y_right]
    add_y_py = proj['py'][add_y_right]
    add_y_py -= na.ceil(lightConeSlice['box_width_fraction'])
    add_y_pdx = proj['pdx'][add_y_right]
    add_y_pdy = proj['pdy'][add_y_right]
    add_y_field = proj[field][add_y_right]
    add_y_weight_field = proj['weight_field'][add_y_right]
    del add_y_right

    # Cells hanging off the bottom edge.
    add_y_left = proj['py'] - 0.5 * proj['pdy'] < 0
    add2_y_px = proj['px'][add_y_left]
    add2_y_py = proj['py'][add_y_left]
    add2_y_py += na.ceil(lightConeSlice['box_width_fraction'])
    add2_y_pdx = proj['pdx'][add_y_left]
    add2_y_pdy = proj['pdy'][add_y_left]
    add2_y_field = proj[field][add_y_left]
    add2_y_weight_field = proj['weight_field'][add_y_left]
    del add_y_left

    # Add the hanging cells back to the projection data.
    proj['px'] = na.concatenate([proj['px'], add_x_px, add_y_px,
                                 add2_x_px, add2_y_px])
    proj['py'] = na.concatenate([proj['py'], add_x_py, add_y_py,
                                 add2_x_py, add2_y_py])
    proj['pdx'] = na.concatenate([proj['pdx'], add_x_pdx, add_y_pdx,
                                  add2_x_pdx, add2_y_pdx])
    proj['pdy'] = na.concatenate([proj['pdy'], add_x_pdy, add_y_pdy,
                                  add2_x_pdy, add2_y_pdy])
    proj[field] = na.concatenate([proj[field], add_x_field, add_y_field,
                                  add2_x_field, add2_y_field])
    proj['weight_field'] = na.concatenate([proj['weight_field'],
                                           add_x_weight_field, add_y_weight_field,
                                           add2_x_weight_field, add2_y_weight_field])

    # Delete original copies of hanging cells.
    del add_x_px, add_y_px, add2_x_px, add2_y_px
    del add_x_py, add_y_py, add2_x_py, add2_y_py
    del add_x_pdx, add_y_pdx, add2_x_pdx, add2_y_pdx
    del add_x_pdy, add_y_pdy, add2_x_pdy, add2_y_pdy
    del add_x_field, add_y_field, add2_x_field, add2_y_field
    del add_x_weight_field, add_y_weight_field, add2_x_weight_field, add2_y_weight_field

    # Tiles were made rounding up the width to the nearest integer.
    # Cut off the edges to get the specified width.
    # Cut in the x direction.
    cut_x = proj['px'] - 0.5 * proj['pdx'] < lightConeSlice['box_width_fraction']
    proj['px'] = proj['px'][cut_x]
    proj['py'] = proj['py'][cut_x]
    proj['pdx'] = proj['pdx'][cut_x]
    proj['pdy'] = proj['pdy'][cut_x]
    proj[field] = proj[field][cut_x]
    proj['weight_field'] = proj['weight_field'][cut_x]
    del cut_x

    # Cut in the y direction.
    cut_y = proj['py'] - 0.5 * proj['pdy'] < lightConeSlice['box_width_fraction']
    proj['px'] = proj['px'][cut_y]
    proj['py'] = proj['py'][cut_y]
    proj['pdx'] = proj['pdx'][cut_y]
    proj['pdy'] = proj['pdy'][cut_y]
    proj[field] = proj[field][cut_y]
    proj['weight_field'] = proj['weight_field'][cut_y]
    del cut_y

    # Create fixed resolution buffer to return back to the light cone object.
    # These buffers will be stacked together to make the light cone.
    frb = FixedResolutionBuffer(proj, (0, lightConeSlice['box_width_fraction'],
                                       0, lightConeSlice['box_width_fraction']),
                                (pixels, pixels), antialias=False)

    return frb
