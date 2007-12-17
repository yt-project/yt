"""
Various non-grid data containers.

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

class EnzoData:
    """
    Generic EnzoData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    """
    _grids = None

    def __init__(self, pf, fields):
        """
        Sets up EnzoData instance

        @param hierarchy: hierarchy we're associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param fields: Fields represented in the data
        @type fields: list of strings
        """
        if pf != None:
            self.pf = pf
            self.hierarchy = pf.hierarchy
        if not isinstance(fields, types.ListType):
            fields = [fields]
        self.fields = fields
        self.data = {}

    def clear_data(self):
        """
        Clears out all data from the EnzoData instance, freeing memory.
        """
        del self.data
        self.data = {}

    def has_key(self, key):
        return self.data.has_key(key)

    def _refresh_data(self):
        """
        Wipes data and rereads/regenerates it from the self.fields.
        """
        self.clear_data()
        self.get_data()

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        if not self.data.has_key(key):
            if key not in self.fields:
                self.fields.append(key)
            self.get_data(key)
        return self.data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.data[key] = val

    def _generate_field(self, fieldName):
        """
        Generates, or attempts to generate, a field not found in the data file

        @param fieldName: field to generate
        @type fieldName: string
        """
        if fieldInfo.has_key(fieldName):
            # First we check the validator
            try:
                fieldInfo[fieldName].Validate(self)
            except NeedsGridType:
                # We leave this to be implementation-specific
                self._generate_field_from_grids(fieldName)
            else:
                self[fieldName] = fieldInfo[fieldName](self)
        else: # Can't find the field, try as it might
            raise exceptions.KeyError, fieldName

    def _generate_field_from_grids(self, fieldName):
        pass

class Enzo2DData(EnzoData):
    """
    Class to represent a set of EnzoData that's 2-D in nature.  Slices and
    projections, primarily.
    """
    def __init__(self, axis, fields, pf=None):
        """
        Prepares the Enzo2DData.

        @param hierarchy: hierarchy associated with this data
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to which this data is parallel
        @type axis: integer
        @param fields: fields to be processed or generated
        @type fields: list of strings
        """
        self.axis = axis
        EnzoData.__init__(self, pf, fields)

    def interpolate_discretize(self, LE, RE, field, side, logSpacing=True):
        """
        This returns a uniform grid of points, interpolated using the nearest
        neighbor method.

        @note: Requires Delaunay triangulation, which is not included in
        most/all scipy binaries.
        @param LE: Left Edge
        @type LE: array of Floats
        @param RE: Right Edge
        @type RE: array of Floats
        @param field: The field to discretize
        @type field: string
        @param side: The number of points on a side
        @type side: int
        """
        import scipy.sandbox.delaunay as de
        if logSpacing:
            zz = na.log10(self[field])
        else:
            zz = self[field]
        xi, yi = na.array( \
                 na.mgrid[LE[0]:RE[0]:side*1j, \
                          LE[1]:RE[1]:side*1j], nT.Float64)
        zi = de.Triangulation(self.x,self.y).nn_interpolator(zz)\
                 [LE[0]:RE[0]:side*1j, \
                  LE[1]:RE[1]:side*1j]
        if logSpacing:
            zi = 10**(zi)
        return [xi,yi,zi]
        #return [xx,yy,zz]

class EnzoSliceBase(Enzo2DData):
    """
    A slice at a given coordinate along a given axis through the entire grid
    hierarchy.
    """

    @time_execution
    def __init__(self, axis, coord, fields = None, center=None, fRet=False,
                 pf=None):
        """
        Returns an instance of EnzoSlice.

        We are not mandating a field be passed in
        The field and coordinate we want to be able to change -- however, the
        axis we do NOT want to change.  So think of EnzoSlice as defining a slice
        operator, rather than a set piece of data.

        Arguments:
        @param hierarchy: hierarchy associated with this data
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to which this data is parallel
        @type axis: integer
        @param coord: three points defining the center
        @type coord: na.array
        @keyword fields: fields to be processed or generated
        @type fields: list of strings
        @keyword fRet: Do we want a false return?  As in, do we want all of our
        data points?
        @type fRet: bool
        """
        Enzo2DData.__init__(self, axis, fields, pf)
        self.center = center
        self.coord = coord
        self.fRet = fRet
        self._refresh_data()

    def reslice(self, coord):
        mylog.debug("Setting coordinate to %0.5e" % coord)
        self.coord = coord
        self.refresh_data()

    def shift(self, val):
        """
        Moves the slice coordinate up by either a floating point value, or an
        integer number of indices of the finest grid.

        @param val: shift amount
        @type val: integer (number of cells) or float (distance)
        """
        if isinstance(val, types.FloatType):
            # We add the dx
            self.coord += val
        elif isinstance(val, types.IntType):
            # Here we assume that the grid is the max level
            level = self.hierarchy.max_level
            self.coord
            dx = self.hierarchy.gridDxs[self.hierarchy.levelIndices[level][0]]
            self.coord += dx * val
        self.refreshData()

    def _generate_coords(self):
        points = []
        for grid in self.grids:
            points.append(self._generate_grid_coords(grid))
        t = na.concatenate(points)
        self['x'] = t[:,0]
        self['y'] = t[:,1]
        self['z'] = t[:,2]
        self['dx'] = t[:,3]
        self['dy'] = t[:,3]
        self['dz'] = t[:,3]
        self.ActiveDimensions = (t.shape[0], 1, 1)

    @time_execution
    def get_data(self, field = None):
        """
        Iterates over the list of fields and generates/reads them all.
        Probably shouldn't be called directly.

        @keyword field: field to add (list or single)
        @type field: string or list of strings
        """
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        self.grids, ind = self.hierarchy.findSliceGrids(self.coord, self.axis)
        if not self.has_key('dx'):
            self._generate_coords()
        if isinstance(field, types.StringType):
            fieldsToGet = [field]
        elif isinstance(field, types.ListType):
            fieldsToGet = field
        elif field == None:
            fieldsToGet = self.fields
        for field in fieldsToGet:
            if self.data.has_key(field):
                continue
            rvs=[]
            if field in self.hierarchy.fieldList:
                self[field] = na.concatenate(
                    [self._get_data_from_grid(grid, field)
                     for grid in self.grids])
            else:
                self._generate_field(field)

    def _generate_grid_coords(self, grid):
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        #sl.reverse()
        sl = tuple(sl)
        nx = grid.myChildMask.shape[xaxis]
        ny = grid.myChildMask.shape[yaxis]
        cm = na.where(grid.myChildMask[sl].ravel() == 1)
        cmI = na.indices((nx,ny))
        xind = cmI[0,:].ravel()
        xpoints = na.ones(cm[0].shape, nT.Float64)
        xpoints *= xind[cm]*grid.dx+(grid.LeftEdge[xaxis] + 0.5*grid.dx)
        yind = cmI[1,:].ravel()
        ypoints = na.ones(cm[0].shape, nT.Float64)
        ypoints *= yind[cm]*grid.dx+(grid.LeftEdge[yaxis] + 0.5*grid.dx)
        zpoints = na.ones(xpoints.shape, nT.Float64) * self.coord
        dx = na.ones(xpoints.shape, nT.Float64) * grid.dx/2.0
        t = na.array([xpoints, ypoints, zpoints, dx]).swapaxes(0,1)
        return t

    def _get_data_from_grid(self, grid, field):
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        slHERE = tuple(sl)
        sl.reverse()
        slHDF = tuple(sl)
        if not grid.has_key(field):
            dv = self._read_data_slice(grid, field, slHDF)
        else:
            dv = grid[field][slHERE]
        cm = na.where(grid.myChildMask[slHERE].ravel() == 1)
        dv=dv.swapaxes(0,2)
        dataVals = dv.ravel()[cm]
        return dataVals

    def _generate_field_from_grids(self, field):
        data = []
        for grid in self.grids:
            temp = grid[field]
            data.append(self._get_data_from_grid(grid,field))
        self[field] = na.concatenate(data)

class EnzoProjBase(Enzo2DData):
    def __init__(self, axis, field, weightField = None,
                 max_level = None, center = None, pf = None,
                 type=0):
        """
        Returns an instance of EnzoProj.

        @param hierarchy: the hierarchy we are projecting
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to project along
        @type axis: integer
        @param field: the field to project (NOT multiple)
        @type field: string
        @keyword weightField: the field to weight by
        @keyword max_level: the maximum level to project through
        @keyword type: The type of projection: 0 for sum, 1 for MIP
        """
        Enzo2DData.__init__(self, axis, field, pf)
        if max_level == None:
            ###CHANGE###
            #max_level = self.hierarchy.max_level
            max_level = self.hierarchy.maxLevel
        self._max_level = max_level
        self._weight = weightField
        self.center = center
        if type == 1:
            self.type="MIP"
            self.func = na.max
        else:
            self.type="SUM"
            self.func = na.sum
        self.__calculate_memory()
        self._refresh_data()

    @time_execution
    def __calculate_memory(self):
        """
        Here we simply calculate how much memory is needed, which speeds up
        allocation later
        """
        level_mem = {}
        h = self.hierarchy
        i = 0
        mylog.info("Calculating memory usage")
        for level in range(self._max_level+1):
            level_mem[level] = 0
            mylog.debug("Examining level %s", level)
            grids = h.levelIndices[level]
            numGrids = len(grids)
            RE = h.gridRightEdge[grids].copy()
            LE = h.gridLeftEdge[grids].copy()
            for grid in h.grids[grids]:
                if (i%1e3) == 0:
                    mylog.debug("Reading and masking %s / %s", i, h.numGrids)
                for ax in [0,1,2]:
                    grid.generateOverlapMasks(ax, LE, RE)
                    grid.myOverlapGrids[ax] = \
                      h.grids[grids[na.where(grid.myOverlapMasks[ax] == 1)]]
                level_mem[level] += \
                          grid.ActiveDimensions.prod() / \
                          grid.ActiveDimensions[ax]
                i += 1
        for level in range(self._max_level+1):
            gI = h.selectLevel(level)
            mylog.debug("%s cells and %s grids for level %s", \
                level_mem[level], len(gI), level)
        mylog.debug("We need %s cells total",
                    na.add.reduce(level_mem.values()))
        self.__memory_per_level = level_mem

    def _serialize(self):
        mylog.info("Serializing data...")
        nodeName = "%s_%s_%s" % (self.fields[0], self.weightField, self.axis)
        mylog.info("nodeName: %s", nodeName)
        projArray = na.array([self.x, self.y, self.dx, self.dy, self[self.fields[0]]])
        self.hierarchy.saveData(projArray, "/Projections", nodeName)
        mylog.info("Done serializing...")

    def _deserialize(self):
        nodeName = "%s_%s_%s" % (self.fields[0], self.weightField, self.axis)
        array=self.hierarchy.getData("/Projections", nodeName)
        if array == None:
            return
        self['x'] = array[0,:]
        self['y'] = array[1,:]
        self['dx'] = array[2,:]
        self['dy']= array[3,:]
        self[self.fields[0]] = array[4,:]
        return True

    def __project_level(self, level, field):
        grids_to_project = self.hierarchy.selectGrids(level)
        zeroOut = (level != self._max_level)
        pbar = get_pbar('Projecting level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids_to_project))
        for pi, grid in enumerate(grids_to_project):
            grid.retVal = grid.getProjection(self.axis, field, zeroOut,
                                             weight=self._weight,
                                             func = self.func)
            pbar.update(pi)
        pbar.finish()
        self.__combine_grids(level) # In-place
        all_data = [ [grid.retVal[j] for grid in grids_to_project] for j in range(5)]
        for grid in grids_to_project:
            cI = na.where(grid.retVal[3]==0) # Where childmask = 0
            grid.coarseData = [grid.retVal[j][cI] for j in range(5)]
        levelData = [na.concatenate(all_data[i]) for i in range(5)]
        mylog.debug("All done combining and refining with a final %s points",
                    levelData[0].shape[0])
        dblI = na.where((levelData[0]>-1) & (levelData[3]==1))
        if self._weight != None:
            weightedData = levelData[2][dblI] / levelData[4][dblI]
        else:
            weightedData = levelData[2][dblI]
        mylog.debug("Level %s done: %s final of %s", \
                   level, len(dblI[0]), \
                   levelData[0].shape[0])
        dx = grids_to_project[0].dx * na.ones(len(dblI[0]), dtype=nT.Float64)
        return [levelData[0][dblI], levelData[1][dblI], weightedData, dx]

    def __combine_grids(self, level):
        grids = self.hierarchy.selectGrids(level)
        pbar = get_pbar('Combining level % 2i / % 2i ' \
                          % (level, self._max_level), len(grids))
        # We have an N^2 check, so we try to be as quick as possible
        # and to skip as many as possible
        for pi, grid1 in enumerate(grids):
            pbar.update(pi)
            if grid1.retVal[0].shape[0] == 0: continue
            for grid2 in grid1.myOverlapGrids[self.axis]:
                if grid2.retVal[0].shape[0] == 0 \
                  or grid1.id == grid2.id:
                    continue
                args = grid1.retVal + grid2.retVal + [0]
                PointCombine.CombineData(*args)
                goodI = na.where(grid2.retVal[0] > -1)
                grid2.retVal = [grid2.retVal[i][goodI] for i in range(5)]
            numRefined = 0
            if level <= self._max_level and level > 0:
                for grid2 in grid1.Parent.myOverlapGrids[self.axis]:
                    if grid2.coarseData[0].shape[0] == 0: continue # Already refined
                    args = grid1.retVal[:3] + [grid1.retVal[4]] + \
                           grid2.coarseData[:3] + [grid2.coarseData[4]] + [2]
                    numRefined += PointCombine.RefineCoarseData(*args)
        pbar.finish()

    @time_execution
    def get_data(self, field = None):
        if not field: field = self.fields[0]
        all_data = []
        h = self.hierarchy
        for level in range(0, self._max_level+1):
            all_data.append(self.__project_level(level, field))
        all_data = na.concatenate(all_data, axis=1)
        # We now convert to half-widths and center-points
        self['dx'] = all_data[3,:]
        self['x'] = (all_data[0,:]+0.5) * self['dx']
        self['y'] = (all_data[1,:]+0.5) * self['dx']
        self['dx'] *= 0.5
        self['dy'] = self['dx'].copy()
        self.data[field] = all_data[2,:]
        # Now, we should clean up after ourselves...
        [grid.clearAll for grid in h.grids]
