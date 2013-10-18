"""
This module contains a routine to search for topologically connected sets



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from itertools import chain
import numpy as np

from yt.funcs import *
import yt.utilities.data_point_utilities as data_point_utilities
import yt.utilities.lib as amr_utils

def identify_contours(data_source, field, min_val, max_val,
                          cached_fields=None):
    tree = amr_utils.ContourTree()
    gct = amr_utils.TileContourTree(min_val, max_val)
    total_contours = 0
    contours = {}
    empty_mask = np.ones((1,1,1), dtype="uint8")
    for (grid, node, (sl, dims, gi)) in data_source.tiles.slice_traverse():
        nid = node.node_id
        values = grid[field][sl].astype("float64")
        contour_ids = np.zeros(dims, "int64") - 1
        gct.identify_contours(values, contour_ids, total_contours)
        new_contours = tree.cull_candidates(contour_ids)
        total_contours += new_contours.shape[0]
        tree.add_contours(new_contours)
        # Now we can create a partitioned grid with the contours.
        pg = amr_utils.PartitionedGrid(
            [contour_ids], empty_mask, g.dds * gi, g.dds * (gi + dims), dims)
        contours[nid] = (g.Level, pg)
    trunk = data_source.tiles.tree.trunk
    amr_utils.link_node_contours(trunk, contours, tree)
    #joins = tree.cull_joins(bt)
    #tree.add_joins(joins)
    joins = tree.export()
    ff = data_source["tempContours"].astype("int64")
    amr_utils.update_joins(joins, ff)
    data_source["tempContours"] = ff.astype("float64")
    data_source._flush_data_to_grids("tempContours", -1, dtype='int64')
    del data_source.field_data["tempContours"] # Force a reload from the grids
    data_source.get_data("tempContours")
    contour_ind = {}
    i = 0
    handled = set()
    for contour_id in data_source["tempContours"]:
        if contour_id == -1 or contour_id in handled: continue
        handled.add(contour_id)
        contour_ind[i] = np.where(data_source["tempContours"] == contour_id)
        mylog.debug("Contour id %s has %s cells", i, contour_ind[i][0].size)
        i += 1
    print "TREE ENTRIES", tree.count()
    mylog.info("Identified %s contours between %0.5e and %0.5e",
               len(contour_ind.keys()),min_val,max_val)
    for grid in chain(grid_set):
        grid.field_data.pop("tempContours", None)
    del data_source.field_data["tempContours"]
    return contour_ind
