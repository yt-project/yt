"""
Create randomly centered, tiled projections to be used in light cones.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import \
     mylog
from yt.visualization.fixed_resolution import \
    FixedResolutionBuffer
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    parallel_blocking_call

@parallel_blocking_call
def _light_cone_projection(my_slice, field, pixels, weight_field=None,
                           save_image=False, field_cuts=None):
    "Create a single projection to be added into the light cone stack."

    # We are just saving the projection object, so only the projection axis 
    # needs to be considered since the lateral shifting and tiling occurs after 
    # the projection object is made.
    # Likewise, only the box_depth_fraction needs to be considered.

    mylog.info("Making projection at z = %f from %s." % \
               (my_slice["redshift"], my_slice["filename"]))

    region_center = [0.5 * (my_slice["object"].domain_right_edge[q] +
                            my_slice["object"].domain_left_edge[q]) \
                         for q in range(my_slice["object"].dimensionality)]

    # 1. The Depth Problem
    # Use coordinate field cut in line of sight to cut projection to proper depth.
    if field_cuts is None:
        these_field_cuts = []
    else:
        these_field_cuts = field_cuts.copy()

    if (my_slice["box_depth_fraction"] < 1):
        axis = ("x", "y", "z")[my_slice["projection_axis"]]
        depthLeft = \
          my_slice["projection_center"][my_slice["projection_axis"]] \
            - 0.5 * my_slice["box_depth_fraction"]
        depthRight = \
          my_slice["projection_center"][my_slice["projection_axis"]] \
            + 0.5 * my_slice["box_depth_fraction"]
        if (depthLeft < 0):
            cut_mask = ("((obj[\"%s\"] + 0.5*obj[\"d%s\"] >= 0) & " + \
              "(obj[\"%s\"] - 0.5*obj[\"d%s\"] <= %f)) | " + \
              "((obj[\"%s\"] + 0.5*obj[\"d%s\"] >= %f) & " + \
              "(obj[\"%s\"] - 0.5*obj[\"d%s\"] <= 1))") % \
                (axis, axis, axis, axis, depthRight, 
                 axis, axis, (depthLeft+1), axis, axis)
        elif (depthRight > 1):
            cut_mask = ("((obj[\"%s\"] + 0.5*obj[\"d%s\"] >= 0) & " + \
              "(obj[\"%s\"] - 0.5*obj[\"d%s\"] <= %f)) | " + \
              "((obj[\"%s\"] + 0.5*obj[\"d%s\"] >= %f) & " + \
              "(obj[\"%s\"] - 0.5*obj[\"d%s\"] <= 1))") % \
                (axis, axis, axis, axis, (depthRight-1),
                 axis, axis, depthLeft, axis, axis)
        else:
            cut_mask = ("(obj[\"%s\"] + 0.5*obj[\"d%s\"] >= %f) & " + \
              "(obj[\"%s\"] - 0.5*obj[\"%s\"] <= %f)") % \
              (axis, axis, depthLeft, axis, axis, depthRight)

        these_field_cuts.append(cut_mask)

    data_source = my_slice["object"].all_data()
    cut_region = data_source.cut_region(these_field_cuts)
        
    # Make projection.
    proj = my_slice["object"].proj(field, my_slice["projection_axis"], 
        weight_field, center=region_center,
        data_source=cut_region)
    proj_field = proj.field[0]

    del data_source, cut_region
    
    # 2. The Tile Problem
    # Tile projection to specified width.

    # Original projection data.
    original_px = proj.field_data["px"].in_units("code_length").copy()
    original_py = proj.field_data["py"].in_units("code_length").copy()
    original_pdx = proj.field_data["pdx"].in_units("code_length").copy()
    original_pdy = proj.field_data["pdy"].in_units("code_length").copy()
    original_field = proj.field_data[proj_field].copy()
    original_weight_field = proj.field_data["weight_field"].copy()

    for my_field in ["px", "py", "pdx", "pdy", proj_field, "weight_field"]:
        proj.field_data[my_field] = [proj.field_data[my_field]]

    # Copy original into offset positions to make tiles.
    for x in range(int(np.ceil(my_slice["box_width_fraction"]))):
        x = my_slice["object"].quan(x, "code_length")
        for y in range(int(np.ceil(my_slice["box_width_fraction"]))):
            y = my_slice["object"].quan(y, "code_length")
            if ((x + y) > 0):
                proj.field_data["px"] += [original_px+x]
                proj.field_data["py"] += [original_py+y]
                proj.field_data["pdx"] += [original_pdx]
                proj.field_data["pdy"] += [original_pdy]
                proj.field_data["weight_field"] += [original_weight_field]
                proj.field_data[proj_field] += [original_field]

    for my_field in ["px", "py", "pdx", "pdy", proj_field, "weight_field"]:
        proj.field_data[my_field] = \
          my_slice["object"].arr(proj.field_data[my_field]).flatten()

    # Delete originals.
    del original_px
    del original_py
    del original_pdx
    del original_pdy
    del original_field
    del original_weight_field

    # 3. The Shift Problem
    # Shift projection by random x and y offsets.

    image_axes = np.roll(np.arange(3), -my_slice["projection_axis"])[1:]
    di_left_x  = my_slice["object"].domain_left_edge[image_axes[0]]
    di_right_x = my_slice["object"].domain_right_edge[image_axes[0]]
    di_left_y  = my_slice["object"].domain_left_edge[image_axes[1]]
    di_right_y = my_slice["object"].domain_right_edge[image_axes[1]]
    
    offset = my_slice["projection_center"].copy() * \
      my_slice["object"].domain_width
    offset = np.roll(offset, -my_slice["projection_axis"])[1:]

    # Shift x and y positions.
    proj.field_data["px"] -= offset[0]
    proj.field_data["py"] -= offset[1]

    # Wrap off-edge cells back around to other side (periodic boundary conditions).
    proj.field_data["px"][proj.field_data["px"] < di_left_x] += \
      np.ceil(my_slice["box_width_fraction"]) * di_right_x
    proj.field_data["py"][proj.field_data["py"] < di_left_y] += \
      np.ceil(my_slice["box_width_fraction"]) * di_right_y

    # After shifting, some cells have fractional coverage on both sides of the box.
    # Find those cells and make copies to be placed on the other side.

    # Cells hanging off the right edge.
    add_x_right = proj.field_data["px"] + 0.5 * proj.field_data["pdx"] > \
      np.ceil(my_slice["box_width_fraction"]) * di_right_x
    add_x_px = proj.field_data["px"][add_x_right]
    add_x_px -= np.ceil(my_slice["box_width_fraction"]) * di_right_x
    add_x_py = proj.field_data["py"][add_x_right]
    add_x_pdx = proj.field_data["pdx"][add_x_right]
    add_x_pdy = proj.field_data["pdy"][add_x_right]
    add_x_field = proj.field_data[proj_field][add_x_right]
    add_x_weight_field = proj.field_data["weight_field"][add_x_right]
    del add_x_right

    # Cells hanging off the left edge.
    add_x_left = proj.field_data["px"] - 0.5 * proj.field_data["pdx"] < di_left_x
    add2_x_px = proj.field_data["px"][add_x_left]
    add2_x_px += np.ceil(my_slice["box_width_fraction"]) * di_right_x
    add2_x_py = proj.field_data["py"][add_x_left]
    add2_x_pdx = proj.field_data["pdx"][add_x_left]
    add2_x_pdy = proj.field_data["pdy"][add_x_left]
    add2_x_field = proj.field_data[proj_field][add_x_left]
    add2_x_weight_field = proj.field_data["weight_field"][add_x_left]
    del add_x_left

    # Cells hanging off the top edge.
    add_y_right = proj.field_data["py"] + 0.5 * proj.field_data["pdy"] > \
      np.ceil(my_slice["box_width_fraction"]) * di_right_y
    add_y_px = proj.field_data["px"][add_y_right]
    add_y_py = proj.field_data["py"][add_y_right]
    add_y_py -= np.ceil(my_slice["box_width_fraction"]) * di_right_y
    add_y_pdx = proj.field_data["pdx"][add_y_right]
    add_y_pdy = proj.field_data["pdy"][add_y_right]
    add_y_field = proj.field_data[proj_field][add_y_right]
    add_y_weight_field = proj.field_data["weight_field"][add_y_right]
    del add_y_right

    # Cells hanging off the bottom edge.
    add_y_left = proj.field_data["py"] - 0.5 * proj.field_data["pdy"] < di_left_y
    add2_y_px = proj.field_data["px"][add_y_left]
    add2_y_py = proj.field_data["py"][add_y_left]
    add2_y_py += np.ceil(my_slice["box_width_fraction"]) * di_right_y
    add2_y_pdx = proj.field_data["pdx"][add_y_left]
    add2_y_pdy = proj.field_data["pdy"][add_y_left]
    add2_y_field = proj.field_data[proj_field][add_y_left]
    add2_y_weight_field = proj.field_data["weight_field"][add_y_left]
    del add_y_left

    # Add the hanging cells back to the projection data.
    proj.field_data["px"] = proj.field_data["px"].unit_quantity * \
      np.concatenate([proj.field_data["px"].d, add_x_px.d, add_y_px.d,
                      add2_x_px.d, add2_y_px.d])
    proj.field_data["py"] = proj.field_data["py"].unit_quantity * \
        np.concatenate([proj.field_data["py"].d, add_x_py.d, add_y_py.d,
                        add2_x_py.d, add2_y_py.d])
    proj.field_data["pdx"] = proj.field_data["pdx"].unit_quantity * \
        np.concatenate([proj.field_data["pdx"].d, add_x_pdx.d, add_y_pdx.d,
                        add2_x_pdx.d, add2_y_pdx.d])
    proj.field_data["pdy"] = proj.field_data["pdy"].unit_quantity * \
        np.concatenate([proj.field_data["pdy"].d, add_x_pdy.d, add_y_pdy.d,
                        add2_x_pdy.d, add2_y_pdy.d])
    proj.field_data[proj_field] = proj.field_data["pdy"].unit_quantity * \
        np.concatenate([proj.field_data[proj_field].d, add_x_field.d, add_y_field.d,
                        add2_x_field.d, add2_y_field.d])
    proj.field_data["weight_field"] = proj.field_data["weight_field"].unit_quantity * \
        np.concatenate([proj.field_data["weight_field"].d,
                        add_x_weight_field.d, add_y_weight_field.d,
                        add2_x_weight_field.d, add2_y_weight_field.d])

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
    cut_x = proj.field_data["px"] - 0.5 * proj.field_data["pdx"] < \
      di_right_x * my_slice["box_width_fraction"]
    proj.field_data["px"] = proj.field_data["px"][cut_x]
    proj.field_data["py"] = proj.field_data["py"][cut_x]
    proj.field_data["pdx"] = proj.field_data["pdx"][cut_x]
    proj.field_data["pdy"] = proj.field_data["pdy"][cut_x]
    proj.field_data[proj_field] = proj.field_data[proj_field][cut_x]
    proj.field_data["weight_field"] = proj.field_data["weight_field"][cut_x]
    del cut_x

    # Cut in the y direction.
    cut_y = proj.field_data["py"] - 0.5 * proj.field_data["pdy"] < \
      di_right_y * my_slice["box_width_fraction"]
    proj.field_data["px"] = proj.field_data["px"][cut_y]
    proj.field_data["py"] = proj.field_data["py"][cut_y]
    proj.field_data["pdx"] = proj.field_data["pdx"][cut_y]
    proj.field_data["pdy"] = proj.field_data["pdy"][cut_y]
    proj.field_data[proj_field] = proj.field_data[proj_field][cut_y]
    proj.field_data["weight_field"] = proj.field_data["weight_field"][cut_y]
    del cut_y

    # Create fixed resolution buffer to return back to the light cone object.
    # These buffers will be stacked together to make the light cone.
    frb = FixedResolutionBuffer(proj, 
        (di_left_x, di_right_x * my_slice["box_width_fraction"],
         di_left_y, di_right_y * my_slice["box_width_fraction"]),
        (pixels, pixels), antialias=False)

    return frb
