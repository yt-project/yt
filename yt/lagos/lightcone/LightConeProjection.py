"""
Create randomly centered, tiled projections to be used in light cones.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Britton Smith.  All Rights Reserved.

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

from yt.lagos.lightcone import *
import random as rand

def LightConeProjection(lightConeSlice,field,pixels,weight_field=None,save_image=False,name=""):
    "Create a single projection to be added into the light cone stack."

    # Use some projection parameters to seed random number generator to make unique node name.
    rand.seed(lightConeSlice['ProjectionCenter'][0] + lightConeSlice['ProjectionCenter'][1] + lightConeSlice['ProjectionCenter'][2] + 
              lightConeSlice['ProjectionAxis'] + lightConeSlice['DepthBoxFraction'] + lightConeSlice['WidthBoxFraction'])
    node_name = "%s_%s_%09d" % (field,weight_field,rand.randint(1,999999999))

    print "Making projection at z = %f from %s." % (lightConeSlice['redshift'],lightConeSlice['filename'])

    # Make an Enzo data object.
    dataset_object = lightConeSlice['object']

    # 1. The Depth Problem
    # Excise a region from the box with the specified depth and make projection with that.
    region_center = [0.5,0.5,0.5]

    region_left = [0.0,0.0,0.0]
    region_left[lightConeSlice['ProjectionAxis']] = lightConeSlice['ProjectionCenter'][lightConeSlice['ProjectionAxis']] - \
        0.5 * lightConeSlice['DepthBoxFraction']

    region_right = [1.0,1.0,1.0]
    region_right[lightConeSlice['ProjectionAxis']] = lightConeSlice['ProjectionCenter'][lightConeSlice['ProjectionAxis']] + \
        0.5 * lightConeSlice['DepthBoxFraction']

    periodic_region = dataset_object.h.periodic_region(region_center,region_left,region_right)

    # Make plot collection.
    pc = raven.PlotCollection(dataset_object,center=region_center)

    # Make projection.
    pc.add_projection(field,lightConeSlice['ProjectionAxis'],weight_field=weight_field,source=periodic_region,use_colorbar=True,
                      node_name=node_name)

    # Serialize projection data.
    pc.plots[0].data._serialize(node_name,force=True)

    # Delete region, maybe save some ram.
    del periodic_region

    # 2. The Tile Problem
    # Tile projection to specified width.

    # Original projection data.
    original_px = pc.plots[0].data['px']
    original_py = pc.plots[0].data['py']
    original_pdx = pc.plots[0].data['pdx']
    original_pdy = pc.plots[0].data['pdy']
    original_field = pc.plots[0].data[field]
    original_weight_field = pc.plots[0].data['weight_field']

    # Copy original into offset positions to make tiles.
    for x in range(int(na.ceil(lightConeSlice['WidthBoxFraction']))):
        for y in range(int(na.ceil(lightConeSlice['WidthBoxFraction']))):
            if ((x + y) > 0):
                pc.plots[0].data['px'] = na.concatenate([pc.plots[0].data['px'],original_px+x])
                pc.plots[0].data['py'] = na.concatenate([pc.plots[0].data['py'],original_py+y])
                pc.plots[0].data['pdx'] = na.concatenate([pc.plots[0].data['pdx'],original_pdx])
                pc.plots[0].data['pdy'] = na.concatenate([pc.plots[0].data['pdy'],original_pdy])
                pc.plots[0].data[field] = na.concatenate([pc.plots[0].data[field],original_field])
                pc.plots[0].data['weight_field'] = na.concatenate([pc.plots[0].data['weight_field'],original_weight_field])

    # Delete originals.
    del original_px
    del original_py
    del original_pdx
    del original_pdy
    del original_field
    del original_weight_field

    # 3. The Shift Problem
    # Shift projection by random x and y offsets.

    offset = lightConeSlice['ProjectionCenter']
    # Delete depth coordinate.
    del offset[lightConeSlice['ProjectionAxis']]

    # Shift x and y positions.
    pc.plots[0]['px'] -= offset[0]
    pc.plots[0]['py'] -= offset[1]

    # Wrap off-edge cells back around to other side (periodic boundary conditions).
    pc.plots[0]['px'][pc.plots[0]['px'] < 0] += na.ceil(lightConeSlice['WidthBoxFraction'])
    pc.plots[0]['py'][pc.plots[0]['py'] < 0] += na.ceil(lightConeSlice['WidthBoxFraction'])

    # After shifting, some cells have fractional coverage on both sides of the box.
    # Find those cells and make copies to be placed on the other side.

    # Cells hanging off the right edge.
    add_x_px = pc.plots[0]['px'][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_x_px -= na.ceil(lightConeSlice['WidthBoxFraction'])
    add_x_py = pc.plots[0]['py'][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_x_pdx = pc.plots[0]['pdx'][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_x_pdy = pc.plots[0]['pdy'][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_x_field = pc.plots[0][field][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_x_weight_field = pc.plots[0]['weight_field'][pc.plots[0]['px'] + 0.5 * pc.plots[0]['pdx'] > na.ceil(lightConeSlice['WidthBoxFraction'])]

    # Cells hanging off the left edge.
    add2_x_px = pc.plots[0]['px'][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]
    add2_x_px += na.ceil(lightConeSlice['WidthBoxFraction'])
    add2_x_py = pc.plots[0]['py'][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]
    add2_x_pdx = pc.plots[0]['pdx'][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]
    add2_x_pdy = pc.plots[0]['pdy'][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]
    add2_x_field = pc.plots[0][field][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]
    add2_x_weight_field = pc.plots[0]['weight_field'][pc.plots[0]['px'] - 0.5 * pc.plots[0]['pdx'] < 0]

    # Cells hanging off the top edge.
    add_y_px = pc.plots[0]['px'][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_y_py = pc.plots[0]['py'][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_y_py -= na.ceil(lightConeSlice['WidthBoxFraction'])
    add_y_pdx = pc.plots[0]['pdx'][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_y_pdy = pc.plots[0]['pdy'][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_y_field = pc.plots[0][field][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]
    add_y_weight_field = pc.plots[0]['weight_field'][pc.plots[0]['py'] + 0.5 * pc.plots[0]['pdy'] > na.ceil(lightConeSlice['WidthBoxFraction'])]

    # Cells hanging off the bottom edge.
    add2_y_px = pc.plots[0]['px'][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]
    add2_y_py = pc.plots[0]['py'][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]
    add2_y_py += na.ceil(lightConeSlice['WidthBoxFraction'])
    add2_y_pdx = pc.plots[0]['pdx'][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]
    add2_y_pdy = pc.plots[0]['pdy'][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]
    add2_y_field = pc.plots[0][field][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]
    add2_y_weight_field = pc.plots[0]['weight_field'][pc.plots[0]['py'] - 0.5 * pc.plots[0]['pdy'] < 0]

    # Add the hanging cells back to the projection data.
    pc.plots[0].data['px'] = na.concatenate([pc.plots[0]['px'],add_x_px,add_y_px,add2_x_px,add2_y_px])
    pc.plots[0].data['py'] = na.concatenate([pc.plots[0]['py'],add_x_py,add_y_py,add2_x_py,add2_y_py])
    pc.plots[0].data['pdx'] = na.concatenate([pc.plots[0]['pdx'],add_x_pdx,add_y_pdx,add2_x_pdx,add2_y_pdx])
    pc.plots[0].data['pdy'] = na.concatenate([pc.plots[0]['pdy'],add_x_pdy,add_y_pdy,add2_x_pdy,add2_y_pdy])
    pc.plots[0].data[field] = na.concatenate([pc.plots[0][field],add_x_field,add_y_field,add2_x_field,add2_y_field])
    pc.plots[0].data['weight_field'] = na.concatenate([pc.plots[0]['weight_field'],
                                                               add_x_weight_field,add_y_weight_field,add2_x_weight_field,add2_y_weight_field])

    # Delete original copies of hanging cells.
    del add_x_px,add_y_px,add2_x_px,add2_y_px
    del add_x_py,add_y_py,add2_x_py,add2_y_py
    del add_x_pdx,add_y_pdx,add2_x_pdx,add2_y_pdx
    del add_x_pdy,add_y_pdy,add2_x_pdy,add2_y_pdy
    del add_x_field,add_y_field,add2_x_field,add2_y_field
    del add_x_weight_field,add_y_weight_field,add2_x_weight_field,add2_y_weight_field

    # Tiles were made rounding up the width to the nearest integer.
    # Cut off the edges to get the specified width.
    pc.set_xlim(0,lightConeSlice['WidthBoxFraction'])
    pc.set_ylim(0,lightConeSlice['WidthBoxFraction'])

    # Save an image if requested.
    if (save_image):
        pc.save(name)

    # Create fixed resolution buffer to return back to the light cone object.
    # These buffers will be stacked together to make the light cone.
    frb = raven.FixedResolutionBuffer(pc.plots[0].data,(0,lightConeSlice['WidthBoxFraction'],0,lightConeSlice['WidthBoxFraction']),
                                      (pixels,pixels),antialias=True)

    return frb
