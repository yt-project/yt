from yt.lagos import *

from nws.client import NetWorkSpace

currentESO = {}

def DispatchedProjection(fn, Workspace = "testing", axis=0, field = "Density"):
    ws = NetWorkSpace(Workspace)
    if not currentESO.has_key(fn):
        a = EnzoStaticOutput(fn)
        currentESO[fn] = a
    else: a = currentESO[fn]
    myGrids = a.h.grids[ws.fetch('myGrids')]
    RE = a.h.gridRightEdge
    LE = a.h.gridLeftEdge
    grids = a.h.grids
    for grid in myGrids:
        grid.generateOverlapMasks(axis, LE, RE)
        grid.myOverlapGrids[axis] = grids[na.where(grid.myOverlapMasks[0] == 1)]
        grid.retVal = grid.getProjection(axis, field, False)
    for grid1 in myGrids:
        print "Combining %s / %s" % (grid1, len(myGrids))
        rv1 = grid1.retVal
        if rv1[0].shape[0] == 0: continue
        for grid2 in grid.myOverlapGrids[axis]:
            try:
                rv2 = grid2.retVal
            except AttributeError:
                continue
            if rv2[0].shape[0] == 0 or \
               grid1.id == grid2.id:
                continue
            PointCombine.CombineData( \
                rv1[0], rv1[1], rv1[2], rv1[3], rv1[4],
                rv2[0], rv2[1], rv2[2], rv2[3], rv2[4],
                0)
            goodI = na.where(rv2[0] > -1)
            grid2.retVal = [rv2[i][goodI] for i in range(5)]

    all_data = [ [ grid.retVal[j] for grid in myGrids] for j in range(5)]
    ws.store('myData',all_data)
