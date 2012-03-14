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
    cur_max_id = np.sum([g.ActiveDimensions.prod() for g in data_source._grids])
    pbar = get_pbar("First pass", len(data_source._grids))
    grids = sorted(data_source._grids, key=lambda g: -g.Level)
    tree = amr_utils.ContourTree()
    for gi,grid in enumerate(grids):
        pbar.update(gi+1)
        cm = data_source._get_cut_mask(grid)
        if cm is True: cm = np.ones(grid.ActiveDimensions, dtype='bool')
        old_field_parameters = grid.field_parameters
        grid.field_parameters = data_source.field_parameters
        local_ind = np.where( (grid[field] > min_val)
                            & (grid[field] < max_val) & cm )
        grid.field_parameters = old_field_parameters
        if local_ind[0].size == 0: continue
        kk = np.arange(cur_max_id, cur_max_id-local_ind[0].size, -1)
        grid["tempContours"] = np.ones(grid.ActiveDimensions, dtype='int64') * -1
        grid["tempContours"][local_ind] = kk[:]
        cur_max_id -= local_ind[0].size
        xi_u,yi_u,zi_u = np.where(grid["tempContours"] > -1)
        cor_order = np.argsort(-1*grid["tempContours"][(xi_u,yi_u,zi_u)])
        fd_orig = grid["tempContours"].copy()
        xi = xi_u[cor_order]
        yi = yi_u[cor_order]
        zi = zi_u[cor_order]
        while data_point_utilities.FindContours(grid["tempContours"], xi, yi, zi) < 0:
            pass
        new_contours = tree.cull_candidates(grid["tempContours"])
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
