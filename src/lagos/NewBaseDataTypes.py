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

class EnzoData:
    """
    Generic EnzoData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    """
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
        @todo: We want d[xyz] and [xyz] to be data fields from now on
        Clears out all data from the EnzoData instance, freeing memory.
        """
        del self.data
        self.data = {}

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
        if self.data.has_key(key):
            return self.data[key]
        else:
            if field not in self.fields:
                self.fields.append(field)
            self.get_data(field)
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
        else:
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
    grids = None
#    @time_execution
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
        integer number of indices of the finest grid.  Note that hippodraw
        doesn't like it if we change the number of values in a column, as of
        right now, so if the number of points in the slice changes, it will
        kill HD.

        @param val: shift amount
        @type val: integer (number of cells) or float (distance)
        """
        if isinstance(val, types.FloatType):
            # We add the dx
            self.coord += val
        elif isinstance(val, types.IntType):
            # Here we assume that the grid is the max level
            level = self.hierarchy.maxLevel
            self.coord
            dx = self.hierarchy.gridDxs[self.hierarchy.levelIndices[level][0]]
            self.coord += dx * val
        self.refreshData()

    def _generate_coords(self):
        for grid in g:
            #print "Generating coords for grid %s" % (grid.id)
            points.append(self._generate_grid_coords(grid))
        t = na.concatenate(points)
        self['x'] = t[:,0]
        self['y'] = t[:,1]
        self['z'] = t[:,2]
        self['dx'] = t[:,3]
        self['dy'] = t[:,3]
        self['dz'] = t[:,3]
        self.ActiveDimensions = (len(self.x), 1, 1)

#    @time_execution
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
        #print g
        points = []
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
