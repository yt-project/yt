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

def coalesce_join_tree(jtree1):
    joins = defaultdict(set)
    nj = jtree1.shape[0]
    for i1 in range(nj):
        current_new = jtree1[i1, 0]
        current_old = jtree1[i1, 1]
        for i2 in range(nj):
            if jtree1[i2, 1] == current_new:
                current_new = max(current_new, jtree1[i2, 0])
        jtree1[i1, 0] = current_new
    for i1 in range(nj):
        joins[jtree1[i1, 0]].update([jtree1[i1, 1], jtree1[i1, 0]])
    updated = -1
    while updated != 0:
        keys = list(reversed(sorted(joins.keys())))
        updated = 0
        for k1 in keys + keys[::-1]:
            if k1 not in joins: continue
            s1 = joins[k1]
            for k2 in keys + keys[::-1]:
                if k2 >= k1: continue
                if k2 not in joins: continue
                s2 = joins[k2]
                if k2 in s1:
                    s1.update(joins.pop(k2))
                    updated += 1
                elif not s1.isdisjoint(s2):
                    s1.update(joins.pop(k2))
                    s1.update([k2])
                    updated += 1
    tr = []
    for k in joins.keys():
        v = joins.pop(k)
        tr.append((k, np.array(list(v), dtype="int64")))
    return tr

def identify_contours(data_source, field, min_val, max_val,
                          cached_fields=None):
    pbar = get_pbar("First pass", len(data_source._grids))
    grids = sorted(data_source._grids, key=lambda g: -g.Level)
    tree = amr_utils.ContourTree()
    gct = amr_utils.GridContourTree(min_val, max_val)
    total_contours = 0
    for gi,grid in enumerate(grids):
        pbar.update(gi+1)
        cm = data_source._get_cut_mask(grid)
        if cm is True: cm = na.ones(grid.ActiveDimensions, dtype='int32')
        old_field_parameters = grid.field_parameters
        grid.field_parameters = data_source.field_parameters
        values = grid[field]
        grid.field_parameters = old_field_parameters
        grid["tempContours"] = na.zeros(grid.ActiveDimensions, "int64") - 1
        gct.identify_contours(values, grid["tempContours"], cm, total_contours)
        new_contours = tree.cull_candidates(grid["tempContours"])
        total_contours += new_contours.shape[0]
        tree.add_contours(new_contours)
    pbar.finish()
    pbar = get_pbar("Calculating joins ", len(data_source._grids))
    grid_set = set()
    for gi,grid in enumerate(grids):
        pbar.update(gi)
        cg = grid.retrieve_ghost_zones(1, "tempContours", smoothed=False)
        grid_set.update(set(cg._grids))
        fd = cg["tempContours"].astype('int64')
        bt = amr_utils.construct_boundary_relationships(fd)
        # This recipe is from josef.pktd on the SciPy mailing list:
        # http://mail.scipy.org/pipermail/numpy-discussion/2009-August/044664.html
        joins = tree.cull_joins(bt)
        tree.add_joins(joins)
    pbar.finish()
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
        contour_ind[i] = na.where(data_source["tempContours"] == contour_id)
        mylog.debug("Contour id %s has %s cells", i, contour_ind[i][0].size)
        i += 1
    print "TREE ENTRIES", tree.count()
    mylog.info("Identified %s contours between %0.5e and %0.5e",
               len(contour_ind.keys()),min_val,max_val)
    for grid in chain(grid_set):
        grid.field_data.pop("tempContours", None)
    del data_source.field_data["tempContours"]
    return contour_ind
