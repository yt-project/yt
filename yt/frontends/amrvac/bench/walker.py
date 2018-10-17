"""This is a prototype subproject

GOAL: write a main() function that can fill in the gaps in the input
data in terms of grids: amrvac .dat file doesn't mention a level L
grid if it's been fully refined into (L+1) level grids, so we insert them as we go


issues : this routine need an initial guess of the final number of grids (ngrids)
this is feasible by hand on a small data set but I don't have a solution ready made for that.

"""

from math import ceil
import numpy as np


def bundle_grids(levels, ixs):
    grids = [{'lvl': l, 'ix': np.array(ix)} for l,ix in zip(levels, ixs)]
    return grids

def get_englobante(ixs):
    """guess the ix coordinates of the grid one level above the input one"""
    return np.array([ceil(i/2) for i in ixs])

def main(patches:list, nblocks:np.ndarray, ngrids:int):
    """return a completes list of grids: patches + missing (fully refined) grids"""
    assert len(patches) <= ngrids
    def xmin(grid:dict):
        """get the lower boundaries of a grid as a fraction of the domain width"""
        block_width = 1/nblocks
        return (grid['ix']-1) / (nblocks * 2**(grid['lvl']-1))

    def xmax(grid:dict):
        """get the upper boundaries of a grid as a fraction of the domain width"""
        return grid['ix'] / (nblocks * 2**(grid['lvl']-1))

    levels_out = [0]     * ngrids
    ixs_out    = [(0,0)] * ngrids
    grids = bundle_grids(levels_out, ixs_out)
    igrid = ipatch = 0
    current_level = 1
    while igrid < ngrids:
        patch = patches[ipatch]
        current_level = min(current_level, patch['lvl'])
        while patch['lvl'] > current_level+1:
            grids[igrid] = {'lvl': current_level, 'ix': get_englobante(patch['ix'])}
            current_level += 1
            igrid += 1
        if patch['lvl'] == current_level+1:
            out_of_bounds_inf = np.any(xmin(patch) < xmin(grids[igrid-1]))
            out_of_bounds_sup = np.any(xmax(patch) > xmax(grids[igrid-1]))
            print(out_of_bounds_inf, out_of_bounds_sup)#; return
            print(patch, grids[igrid-1])
            print(xmin(patch),  xmin(grids[igrid-1]))
            print(xmax(patch),  xmax(grids[igrid-1]))
            #return
            if out_of_bounds_inf or out_of_bounds_sup :
                grids[igrid] = {'lvl': current_level, 'ix': get_englobante(patch['ix'])}
                current_level += 1
                igrid += 1
                
        grids[igrid] = patches[ipatch]
        igrid  += 1
        ipatch += 1
    return grids
