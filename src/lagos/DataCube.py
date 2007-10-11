"""
Method for extracting a datacube.  When this goes out of testing phase it will
become a fully-featured class, deriving from EnzoData.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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
from yt.progressbar import *
from math import log10
import types

def dataCube(pf, LE, RE, cubeDim, fields):
    """
    Return a datacube at a fixed resolution of a given field or fields, with
    given left edge and right edge.

    @param pf: ParameterFile to extract from
    @param LE: Left edge (list or array)
    @param RE: Right edge (list or array)
    @param cubeDim: The number of cells on a side in the cube
    @param fields: The field(s) to extract.
    """
    if not isinstance(fields, types.ListType): fields = [fields]
    LE = na.array(LE)
    RE = na.array(RE)
    dx = ((RE-LE)/cubeDim).max()
    grids, gridI = pf.hierarchy.getBoxGrids(LE, RE)
    #grids, gridI = pf.hierarchy.getBoxGrids([0,0,0],[1,1,1])
    #grids = pf.hierarchy.grids
    #gridI = na.indices(grids.shape)
    cubesOut = [na.zeros([cubeDim]*3, dtype=nT.Float32) - 999] * len(fields)
    gridIbyLevel = na.argsort(pf.hierarchy.gridLevels[gridI], axis=0)
    widgets = [ 'Getting grid data  ',
                Percentage(), ' ',
                Bar(marker=RotatingMarker()),
                ' ', ETA(), ' ']
    pbar = ProgressBar(widgets=widgets,
                             maxval=len(grids)).start()
    for gi,grid in enumerate(grids[gridIbyLevel].ravel()):
        pbar.update(gi)
        #print grid, dx, grid.dx#, grid.myChildMask
        if grid.dx < dx:
            #print "Grid dx too small", grid, dx, grid.dx
            continue
        lastLevel = False
        if grid.dx / 2.0 < dx:
            lastLevel = True
            #print "LASTLEVEL", grid.dx / 2.0, dx
        for i,field in enumerate(fields):
            #print "Getting field", field
            weaveDict = {
                'nxc': int(grid.ActiveDimensions[0]),
                'nyc': int(grid.ActiveDimensions[1]),
                'nzc': int(grid.ActiveDimensions[2]),
                'leftEdgeCoarse': na.array(grid.LeftEdge),
                'rf': int(grid.dx / dx),
                'dx_c': float(grid.dx),
                'dy_c': float(grid.dx),
                'dz_c': float(grid.dx),
                'dx_f': float(dx),
                'dy_f': float(dx),
                'dz_f': float(dx),
                'cubeLeftEdge': na.array(LE),
                'cubeRightEdge': na.array(RE),
                'nxf': int(cubeDim),
                'nyf': int(cubeDim),
                'nzf': int(cubeDim),
                'fieldData': grid[field],
                'cubeData': cubesOut[i],
                'childMask' : grid.myChildMask,
                'lastLevel': int(lastLevel)
            }
            #print "Inlining", weaveDict['rf']
            weave.inline(DataCubeRefineCoarseData,
                         weaveDict.keys(), local_dict=weaveDict,
                         compiler='gcc',
                         type_converters=converters.blitz,
                         auto_downcast=0, verbose=2)
            #print "Done inlining"
    return cubesOut
