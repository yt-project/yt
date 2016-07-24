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

import numpy as np

from collections import defaultdict

from yt.funcs import mylog, get_pbar
from yt.utilities.lib.contour_finding import \
    ContourTree, TileContourTree, link_node_contours, \
    update_joins
from yt.utilities.lib.partitioned_grid import \
    PartitionedGrid

def identify_contours(data_source, field, min_val, max_val,
                          cached_fields=None):
    tree = ContourTree()
    gct = TileContourTree(min_val, max_val)
    total_contours = 0
    contours = {}
    node_ids = []
    DLE = data_source.ds.domain_left_edge
    masks = dict((g.id, m) for g, m in data_source.blocks)
    for (g, node, (sl, dims, gi)) in data_source.tiles.slice_traverse():
        node.node_ind = len(node_ids)
        nid = node.node_id
        node_ids.append(nid)
        values = g[field][sl].astype("float64")
        contour_ids = np.zeros(dims, "int64") - 1
        mask = masks[g.id][sl].astype("uint8")
        total_contours += gct.identify_contours(values, contour_ids,
                                                mask, total_contours)
        new_contours = tree.cull_candidates(contour_ids)
        tree.add_contours(new_contours)
        # Now we can create a partitioned grid with the contours.
        LE = (DLE + g.dds * gi).in_units("code_length").ndarray_view()
        RE = LE + (dims * g.dds).in_units("code_length").ndarray_view()
        pg = PartitionedGrid(g.id,
            [contour_ids.view("float64")], mask,
            LE, RE, dims.astype("int64"))
        contours[nid] = (g.Level, node.node_ind, pg, sl)
    node_ids = np.array(node_ids).astype("int64")
    if node_ids.size == 0:
        return 0, {}
    trunk = data_source.tiles.tree.trunk
    mylog.info("Linking node (%s) contours.", len(contours))
    link_node_contours(trunk, contours, tree, node_ids)
    mylog.info("Linked.")
    #joins = tree.cull_joins(bt)
    #tree.add_joins(joins)
    joins = tree.export()
    contour_ids = defaultdict(list)
    pbar = get_pbar("Updating joins ... ", len(contours))
    final_joins = np.unique(joins[:,1])
    for i, nid in enumerate(sorted(contours)):
        level, node_ind, pg, sl = contours[nid]
        ff = pg.my_data[0].view("int64")
        update_joins(joins, ff, final_joins)
        contour_ids[pg.parent_grid_id].append((sl, ff))
        pbar.update(i)
    pbar.finish()
    rv = dict()
    rv.update(contour_ids)
    # NOTE: Because joins can appear in both a "final join" and a subsequent
    # "join", we can't know for sure how many unique joins there are without
    # checking if no cells match or doing an expensive operation checking for
    # the unique set of final join values.
    return final_joins.size, rv
