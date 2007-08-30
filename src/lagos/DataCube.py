# Let's take this a single grid at a time.
#
# Grid1 -> 16x18x20 (for instance)
# Let's pre-suppose it's completely within the region we want.
# We have a refinement factor (let's call it nref) which is really just
# currentdx/finaldx
#
# for i in range(nref):
#   for j in range(nref):
#     for k in range(nref):
#       outGrid[nref*indices[0]+i, nref*indices[1]+j, nref*indices[2]+k] =
#                  ingrid[indices]
#
# Now, what do we do in case we need to clip the refined grid's edges?
# (Note that we will, actually, have to run this process for every single grid,
# thus making the previous algorithm completely stupid!)
#
# i,j,k = indices(na.array(outGrid.shape)*nref)
# toGrab = na.where( \
#   na.all([ na.logical_and(i >= beginIndices[0], i <= endIndices[0]), \
#            na.logical_and(j >= beginIndices[1], j <= endIndices[1]), \
#            na.logical_and(k >= beginIndices[2], k <= endIndices[2])]) )
# That gives the indices that we want to transform.
# Now we work in reverse!  We have the refined cell locations we want, so we
# just have to convert that to the coarse cell.
# 
# coarseCell = (floor(float(i)/nref),
#               floor(float(j)/nref),
#               floor(float(k)/nref))
# gridOut[i,j,k] = gridIn[coarseCell]
#

from yt.lagos import *
from math import log10
import types

def dataCube(pf, LE, RE, cubeDim, fields):
    if not isinstance(fields, types.ListType): fields = [fields]
    dx = ((RE-LE)/cubeDim).max()
    grids, gridI = pf.hierarchy.getBoxGrids(LE, RE)
    cubesOut = [na.zeros([cubeDim]*3, dtype=nT.Float32)] * len(fields)
    gridIbyLevel = na.argsort(pf.hierarchy.gridLevels[gridI], axis=0)
    for grid in grids[gridIbyLevel].ravel():
        if grid.dx < dx:
            print "Grid dx too small", grid
            continue
        # Let's get our start and end indices
        overlapShape = (na.minimum(RE, grid.RightEdge) 
                      - na.maximum(LE, grid.LeftEdge))/dx
        overlapShape = na.floor(overlapShape).astype(nT.Int32)
        if na.any(overlapShape < 0):
            continue
        nref = grid.dx / dx
        print "Working on", grid, nref, grid.dx, dx, grid.Level
        cubeOffset = na.maximum(0, grid.LeftEdge - LE)
        gridOffset = ((na.maximum(grid.LeftEdge, LE)) - grid.LeftEdge)/grid.dx
        print overlapShape, cubeOffset
        i,j,k = na.indices(overlapShape)
        i += (cubeOffset[0]/dx).astype(nT.Int32)
        j += (cubeOffset[1]/dx).astype(nT.Int32)
        k += (cubeOffset[2]/dx).astype(nT.Int32)
        gi = (na.floor(i.astype("float32")/nref)+gridOffset[0]).astype(nT.Int32)
        gj = (na.floor(j.astype("float32")/nref)+gridOffset[1]).astype(nT.Int32)
        gk = (na.floor(k.astype("float32")/nref)+gridOffset[2]).astype(nT.Int32)
        for i, field in enumerate(fields):
            print "Copying", field
            cubesOut[i][i,j,k] = grid[field][gi,gj,gk]
    return cubesOut
