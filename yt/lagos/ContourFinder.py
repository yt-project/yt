"""
This module contains a routine to search for topologically connected sets

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *

class GridConsiderationQueue:
    def __init__(self, white_list, priority_func=None):
        """
        This class exists to serve the contour finder.  It ensures that
        we can create a cascading set of queue dependencies, and if
        a grid is touched again ahead of time we can bump it to the top
        of the queue again.  It like has few uses.
        """
        self.to_consider = []
        self.considered = set()
        self.n = 0
        self.white_list = set(white_list)
        self.priority_func = priority_func

    def add(self, grids, force=False):
        if not hasattr(grids,'size'):
            grids = ensure_list(grids)
        i = self.n
        to_check = self.white_list.intersection(grids)
        if not force: to_check.difference_update(self.considered)
        for g in sorted(to_check, key=self.priority_func):
            try:
                # We only delete from subsequent checks
                del self.to_consider[self.to_consider.index(g, i)]
                self.to_consider.insert(i,g)
                i += 1
            except ValueError:
                self.to_consider.append(g)

    def __iter__(self):
        return self

    def next(self):
        if self.n >= len(self.to_consider):
            raise StopIteration
        tr = self.to_consider[self.n]
        self.considered.add(tr)
        self.n += 1
        return tr

    def progress(self):
        return self.n, len(self.to_consider)

# We want an algorithm that deals with growing a given contour to *all* the
# cells in a grid.

def identify_contours(data_source, field, min_val, max_val, cached_fields=None):
    """
    Given a *data_source*, we will search for topologically connected sets
    in *field* between *min_val* and *max_val*.
    """
    if cached_fields is None: cached_fields = defaultdict(lambda: dict())
    maxn_cells = na.sum([g.ActiveDimensions.prod() for g in data_source._grids])
    contour_ind = na.where( (data_source[field] > min_val)
                          & (data_source[field] < max_val))[0]
    np = contour_ind.size
    if np == 0:
        return {}
    cur_max_id = maxn_cells - np
    contour_ids = na.arange(maxn_cells, cur_max_id, -1) + 1 # Minimum of 1
    data_source["tempContours"] = na.ones(data_source[field].shape, dtype='int32')*-1
    mylog.info("Contouring over %s cells with %s candidates", contour_ids[0],np)
    data_source["tempContours"][contour_ind] = contour_ids[:]
    data_source._flush_data_to_grids("tempContours", -1, dtype='int32')
    my_queue = GridConsiderationQueue(white_list = data_source._grids,
                    priority_func = lambda g: -1*g["tempContours"].max())
    my_queue.add(data_source._grids)
    for i,grid in enumerate(my_queue):
        mylog.info("Examining %s of %s", *my_queue.progress())
        max_before = grid["tempContours"].max()
        to_get = ["tempContours"]
        if field in cached_fields[grid.id] and \
            not na.any( (cached_fields[grid.id][field] > min_val)
                      & (cached_fields[grid.id][field] < max_val)):
            continue
        for f in [field, "GridIndices"]:
            if f not in cached_fields[grid.id]: to_get.append(f)
        cg = grid.retrieve_ghost_zones(1,to_get)
        for f in [field, "GridIndices"]:
            if f in cached_fields[grid.id]:
                cg.data[f] = cached_fields[grid.id][f]
            else:
                cached_fields[grid.id][f] = cg[f] 
        local_ind = na.where( (cg[field] > min_val)
                            & (cg[field] < max_val)
                            & (cg["tempContours"] == -1) )
        if local_ind[0].size > 0:
            kk = na.arange(cur_max_id, cur_max_id-local_ind[0].size, -1)
            cg["tempContours"][local_ind] = kk[:]
            cur_max_id -= local_ind[0].size
        fd = cg["tempContours"].astype('int64')
        fd_original = fd.copy()
        if na.all(fd > -1):
            fd[:] = fd.max()
        else:
            xi_u,yi_u,zi_u = na.where(fd > -1)
            cor_order = na.argsort(-1*fd[(xi_u,yi_u,zi_u)])
            xi = xi_u[cor_order]
            yi = yi_u[cor_order]
            zi = zi_u[cor_order]
            while PointCombine.FindContours(fd, xi, yi, zi) < 0: pass
        cg["tempContours"] = fd.copy().astype('float64')
        cg.flush_data("tempContours")
        my_queue.add(cg._grids)
        force_ind = na.unique(cg["GridIndices"][na.where(
            (cg["tempContours"] > fd_original)
          & (cg["GridIndices"] != grid.id-grid._id_offset) )])
        if len(force_ind) > 0:
            my_queue.add(data_source.hierarchy.grids[force_ind.astype('int32')], force=True)
        for ax in 'xyz':
            if not iterable(grid['d%s'%ax]):
                grid['d%s'%ax] = grid['d%s'%ax]*na.ones(grid.ActiveDimensions)
    del data_source.data["tempContours"] # Force a reload from the grids
    data_source.get_data("tempContours", in_grids=True)
    i = 0
    contour_ind = {}
    for contour_id in na.unique(data_source["tempContours"]):
        if contour_id == -1: continue
        contour_ind[i] = na.where(data_source["tempContours"] == contour_id)
        mylog.debug("Contour id %s has %s cells", i, contour_ind[i][0].size)
        i += 1
    mylog.info("Identified %s contours between %0.5e and %0.5e",
               len(contour_ind.keys()),min_val,max_val)
    for grid in data_source._grids:
        grid.data.pop("tempContours", None)
    del data_source.data["tempContours"]
    return contour_ind

def check_neighbors(data_object, field="Contours"):
    """
    This method is a means of error checking in the contour finder.
    """
    n_bad = na.zeros(1, dtype='int32')
    for cid in na.unique(data_object[field]):
        if cid == -1: continue
        ids = na.where(data_object[field] == cid)
        mx = data_object['x'][ids].copy()
        my = data_object['y'][ids].copy()
        mz = data_object['z'][ids].copy()
        mdx = data_object['dx'][ids].copy()
        mdy = data_object['dy'][ids].copy()
        mdz = data_object['dz'][ids].copy()
        grid_ids_m = data_object['GridIndices'][ids].copy()
        grid_levels_m = data_object['GridLevel'][ids].copy()
        mp = mx.size
        ids = na.where( (data_object[field] != cid)
                      & (data_object[field] >=  0 ))
        nx = data_object['x'][ids].copy()
        ny = data_object['y'][ids].copy()
        nz = data_object['z'][ids].copy()
        ndx = data_object['dx'][ids].copy()
        ndy = data_object['dy'][ids].copy()
        ndz = data_object['dz'][ids].copy()
        grid_ids_n = data_object['GridIndices'][ids].copy()
        grid_levels_n = data_object['GridLevel'][ids].copy()
        np = nx.size
        weave.inline(check_cell_distance,
                   ['mx','my','mz','mdx','mdy','mdz','mp',
                    'nx','ny','nz','ndx','ndy','ndz','np','n_bad',
                    'grid_ids_m', 'grid_levels_m', 'grid_ids_n', 'grid_levels_n'],
                    compiler='gcc', type_converters=converters.blitz,
                    auto_downcast=0, verbose=2)
    return n_bad[0]

check_cell_distance = \
r"""
using namespace std;
int i, j, k;
long double cell_dist, rad_m, rad_n;
k=0;
for(i=0;i<mp;i++){
  for(j=0;j<np;j++){
    /*
   cell_dist = sqrtl(pow(mx(i)-nx(j),2) +
                     pow(my(i)-ny(j),2) +
                     pow(mz(i)-nz(j),2));
   rad_m = sqrtl(pow(mdx(i),2) +
                 pow(mdy(i),2) +
                 pow(mdz(i),2));
   rad_n = sqrtl(pow(ndx(j),2) +
                 pow(ndy(j),2) +
                 pow(ndz(j),2));
    */
   //if(cell_dist > 1.01 * (rad_n/2.0+rad_m/2.0)) continue;
   if(fabsl(mx(i)-nx(j))>(mdx(i)+ndx(j))/2.0) continue;
   if(fabsl(my(i)-ny(j))>(mdy(i)+ndy(j))/2.0) continue;
   if(fabsl(mz(i)-nz(j))>(mdz(i)+ndz(j))/2.0) continue;
   k++;
   break;
   cout << cell_dist << "\t" << 1.01*(rad_n/2.0+rad_m/2.0) << "\t";
   cout << grid_ids_m(i) << "\t" << grid_ids_n(j) << endl;
  }
}
n_bad(0) += k;
"""
