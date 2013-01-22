"""
This module contains a routine to search for topologically connected sets

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

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
    total_contours = 0
    tree = []
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
        total_contours += np.unique(grid["tempContours"][grid["tempContours"] > -1]).size
        new_contours = np.unique(grid["tempContours"][grid["tempContours"] > -1]).tolist()
        tree += zip(new_contours, new_contours)
    tree = set(tree)
    pbar.finish()
    pbar = get_pbar("Calculating joins ", len(data_source._grids))
    grid_set = set()
    for gi,grid in enumerate(grids):
        pbar.update(gi)
        cg = grid.retrieve_ghost_zones(1, "tempContours", smoothed=False)
        grid_set.update(set(cg._grids))
        fd = cg["tempContours"].astype('int64')
        boundary_tree = amr_utils.construct_boundary_relationships(fd)
        tree.update(((a, b) for a, b in boundary_tree))
    pbar.finish()
    sort_new = np.array(list(tree), dtype='int64')
    mylog.info("Coalescing %s joins", sort_new.shape[0])
    joins = coalesce_join_tree(sort_new)
    #joins = [(i, np.array(list(j), dtype="int64")) for i, j in sorted(joins.items())]
    pbar = get_pbar("Joining ", len(joins))
    # This process could and should be done faster
    print "Joining..."
    t1 = time.time()
    ff = data_source["tempContours"].astype("int64")
    amr_utils.update_joins(joins, ff)
    data_source["tempContours"] = ff.astype("float64")
    #for i, new in enumerate(sorted(joins.keys())):
    #    pbar.update(i)
    #    old_set = joins[new]
    #    for old in old_set:
    #        if old == new: continue
    #        i1 = (data_source["tempContours"] == old)
    #        data_source["tempContours"][i1] = new
    t2 = time.time()
    print "Finished joining in %0.2e seconds" % (t2-t1)
    pbar.finish()
    data_source._flush_data_to_grids("tempContours", -1, dtype='int64')
    del data_source.field_data["tempContours"] # Force a reload from the grids
    data_source.get_data("tempContours", in_grids=True)
    contour_ind = {}
    i = 0
    for contour_id in np.unique(data_source["tempContours"]):
        if contour_id == -1: continue
        contour_ind[i] = np.where(data_source["tempContours"] == contour_id)
        mylog.debug("Contour id %s has %s cells", i, contour_ind[i][0].size)
        i += 1
    mylog.info("Identified %s contours between %0.5e and %0.5e",
               len(contour_ind.keys()),min_val,max_val)
    for grid in chain(grid_set):
        grid.field_data.pop("tempContours", None)
    del data_source.field_data["tempContours"]
    return contour_ind
