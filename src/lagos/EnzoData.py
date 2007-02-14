#
# raven:
#   A module for dealing with Enzo data
#   Currently isolated fromall HippoDraw classes
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from yt.lagos import *

class EnzoData:
    """
    Generic EnzoData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    """
    def __init__(self, hierarchy, fields):
        """
        Sets up EnzoData instance

        Arguments:
            hierarchy -- EnzoHierarchy we're associated with
            fields -- Fields represented in the data
        """
        self.hierarchy = hierarchy
        if not isinstance(fields, types.ListType):
            fields = [fields]
        self.fields = fields
        self.data = {}
        
    def clearData(self):
        """
        Clears out all data from the EnzoData instance, freeing memory.
        """
        del self.data 
        self.data = {}
        self.dx = None
        self.dy = None
        self.dz = None
        self.x = None
        self.y = None
        self.z = None

    def refreshData(self):
        """
        Wipes data and rereads/regenerates it from the self.fields.
        """
        self.clearData()
        self.getData()

    def readAllData(self):
        """
        Reads all fields defined in self.hierarchy.fieldList except those
        containing "particle" or "Dark"
        """
        for field in self.hierarchy.fieldList:
            if field.find("particle") > -1 or field.find("Dark") > -1:
                continue
            t = self[field]

    def __getitem__(self, key):
        """
        Returns a single field.  Will add if necessary.
        """
        #print "Getting item %s" % (key)
        if self.data.has_key(key):
            return self.data[key]
        else:
            self.addField(key)
            return self.data[key]

    def __setitem__(self, key, val):
        """
        Sets a field to be some other value.
        """
        self.data[key] = val

    def addField(self, field):
        """
        Adds a field to the current fields already in memory.  Will
        read/generate as necessary.
        """
        if field not in self.fields:
            self.fields.append(field)
        if not self.data.has_key(field):
            self.getData(field)

    def generateField(self, fieldName):
        """
        Generates, or attempts to generate,  a field not found in the file

        Arguments:
            fieldName -- string, field name

        See fields.py for more information.  fieldInfo.keys() will list all of
        the available derived fields.  Note that we also make available the
        suffices _Fraction and _Squared here.  All fields prefixed with 'k'
        will force an attempt to use the chemistry tables to generate them from
        temperature.  All fields used in generation will remain resident in
        memory.

        (should be moved to a standalone?  Or should EnzoGrid subclass EnzoData?)
        """
        if fieldName.endswith("Fraction"):
            # Very simple mass fraction here.  Could be modified easily,
            # but that would require a dict lookup, which is expensive, or
            # an elif block, which is inelegant
            baryonField = "%s_Density" % (fieldName[:-9])
            self[fieldName] = self[baryonField] / self["Density"]
        elif fieldName.endswith("Squared"):
            baryonField = fieldName[:-8]
            self[fieldName] = (self[baryonField])**2.0
        elif fieldInfo.has_key(fieldName):
            # We do a fallback to checking the fieldInfo dict
            # Note that it'll throw an exception here if it's not found...
            # ...which I'm cool with
            fieldInfo[fieldName][3](self, fieldName)
        elif fieldName.startswith("k"):
            self[fieldName] = abs(self.hierarchy.rates[self["Temperature"],fieldName])
        else:
            raise exceptions.KeyError, fieldName

class Enzo2DData(EnzoData):
    """
    Class to represent a set of EnzoData that's 2-D in nature.  Slices and
    projections, primarily.
    """
    def __init__(self, hierarchy, axis, fields):
        """
        Prepares the Enzo2DData.

        Arguments:
            hierarchy -- EnzoHierarchy associated with this data
            axis -- axis to which this data is parallel
            fields -- list of fields to be processed or generated
        """
        self.axis = axis
        EnzoData.__init__(self, hierarchy, fields)

    def output(self, fields, filename):
        """
        Doesn't work yet -- allPoints doesn't quite exist.
        """
        if not isinstance(fields, types.ListType):
            fields = [fields]
        f = open(filename,"w")
        f.write("x\ty\tdx\tdy")
        for field in fields:
            f.write("\t%s" % (field))
        f.write("\n")
        for i in range(self.x.shape[0]):
            f.write("%0.20f %0.20f %0.20f %0.20f" % \
                    (self.x[i], \
                     self.y[i], \
                     self.dx[i], \
                     self.dy[i]))
            for field in fields:
                f.write("\t%0.7e" % (self[field][i]))
            f.write("\n")
        f.close()

class EnzoProj(Enzo2DData):
    """
    EnzoProj data.  Doesn't work yet.  Needs to be implemented, following the
    method of EnzoHierarchy.getProjection
    """
    def __init__(self, hierarchy, axis, fields = None):
        self.coords = None
        self.refreshData()

class EnzoSlice(Enzo2DData):
    """
    A slice at a given coordinate along a given axis through the entire grid
    hierarchy.
    """
    def __init__(self, hierarchy, axis, coord, fields = None):
        """
        Returns an instance of EnzoSlice.

        Arguments:
            hierarchy -- EnzoHierarchy we are associated with
            axis -- axis along which we are slicing (0,1,2)
            coord -- Array of *three* points defining the center
        Keyword Arguments:
            fields -- fields to grab/generated. (Defaults to None.)

        We are not mandating a field be passed in
        The field and coordinate we want to be able to change -- however, the
        axis we do NOT want to change.
        So think of EnzoSlice as defining a slice operator, rather than a
        set piece of data.
        """
        Enzo2DData.__init__(self, hierarchy, axis, fields)
        self.coord = coord
        self.coords = None
        self.refreshData()

    def shift(self, val):
        """
        Moves the slice coordinate up by either a floating point value, or an
        integer number of indices of the finest grid.  Note that hippodraw 
        doesn't like it if we change the number of values in a column, as of
        right now, so if the number of points in the slice changes, it will
        kill HD.

        Arguments:
            val -- integer or float.
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

    def generateCoords(self):
        """
        Place holder, so that generating fields doesn't die.  Note that
        coordinates are generated at initialization.
        """
        pass

    def getData(self, field = None):
        """
        Iterates over the list of fields and generates/reads them all.

        Keyword Arguments:
            field -- field to add (list or single)
        """
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        g,ind = self.hierarchy.findSliceGrids(self.coord, self.axis)
        self.grids = g
        #print g
        points = []
        if self.dx == None:
            print "Generating coordinates"
            for grid in g:
                #print "Generating coords for grid %s" % (grid.id)
                points.append(self.generateGridCoords(grid))
            t = concatenate(points)
            print t.shape
            self.x = t[:,0]
            self.y = t[:,1]
            self.z = t[:,2]
            self.dx = t[:,3]
            self.dy = t[:,3]
            self.dz = t[:,3]
            self.coords = array([self.x, self.y, self.z])
            self.ActiveDimensions = (len(self.x), 1, 1)
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
                for grid in g:
                    #print "Getting %s from grid %s" % (field, grid.id)
                    rvs.append(self.getDataFromGrid(grid, field))
                    self[field] = concatenate(rvs)
            else:
                i = 0
                self.generateField(field)

    def generateGridCoords(self, grid):
        """
        Returns coordinates (x,y,z,dx) for a single grid.

        Arguments:
            grid -- EnzoGrid object
        """
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        if grid.myChildMask == None:
            #print "Generating child mask"
            grid.generateChildMask()
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        #sl.reverse()
        sl = tuple(sl)
        nx = grid.myChildMask.shape[xaxis]
        ny = grid.myChildMask.shape[yaxis]
        cm = where(grid.myChildMask[sl].flat == 1)
        cmI = indices((nx,ny))
        xind = cmI[0,:].flat
        xpoints = ones(cm[0].shape, Float64)
        xpoints *= xind[cm]*grid.dx+(grid.LeftEdge[xaxis] + 0.5*grid.dx)
        yind = cmI[1,:].flat
        ypoints = ones(cm[0].shape, Float64)
        ypoints *= yind[cm]*grid.dx+(grid.LeftEdge[yaxis] + 0.5*grid.dx)
        zpoints = ones(xpoints.shape, Float64) * self.coord
        dx = ones(xpoints.shape, Float64) * grid.dx/2.0
        t = array([xpoints, ypoints, zpoints, dx])
        t.swapaxes(0,1)
        return t

    def getDataFromGrid(self, grid, field):
        """
        Reads a slice from a single grid

        Arguments:
            grid -- EnzoGrid object
            field -- field to read

        Could probably be more efficient, possibly by accepting a list of
        fields.
        """
        if grid.myChildMask == None:
            #print "Generating child mask"
            grid.generateChildMask()
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        # NOTE: This should be fixed.  I don't think it works properly or
        # intelligently.
        wantedIndex = int(((self.coord-grid.LeftEdge[self.axis])/grid.dx))
        #print self.coord, grid.LeftEdge, grid.dx
        sl = [slice(None), slice(None), slice(None)]
        sl[self.axis] = slice(wantedIndex, wantedIndex + 1)
        slHERE = tuple(sl)
        sl.reverse()
        slHDF = tuple(sl)
        #print sl
        dv = self.readDataSlice(grid, field, slHDF)
        #print dv.flat.shape, grid.myChildMask[slHERE].flat.shape, grid.ActiveDimensions
        cm = where(grid.myChildMask[slHERE].flat == 1)
        #print field, cm[0].shape, dv.shape
        dv.swapaxes(0,2)
        dataVals = dv.flat[cm]
        #print cm[0].shape, grid.myChildMask.shape, dataVals.shape
        # We now have a couple one dimensional arrays.  We will
        # make these one array, and return them as [x y val dx dy]
        return dataVals

class Enzo3DData(EnzoData):
    """
    Class describing a cluster of data points, not necessarily sharing a
    coordinate.
    """
    def __init__(self, hierarchy, center, fields):
        """
        Returns an instance of Enzo3DData, or prepares one.

        Arguments:
            hierarchy -- EnzoHierarchy we are associated with
            center -- array, center of the region
            fields -- fields to read/generate
        """
        # We are not mandating a field be passed in
        # The field and coordinate we want to be able to change -- however, the
        # center we do NOT want to change.
        EnzoData.__init__(self, hierarchy, fields)
        self.center = center
        self.coords = None
        self.dx = None

    def output(self, fields, filename):
        """
        Doesn't work yet -- allPoints doesn't quite exist
        """
        if not isinstance(fields, types.ListType):
            fields = [fields]
        f = open(filename,"w")
        f.write("x\ty\tz\tdx\tdy\tdz")
        for field in fields:
            f.write("\t%s" % (field))
        f.write("\n")
        for i in range(self.x.shape[0]):
            f.write("%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f" % \
                    (self.x[i], \
                     self.y[i], \
                     self.z[i], \
                     self.dx[i], \
                     self.dy[i], \
                     self.dz[i]))
            for field in fields:
                f.write("\t%0.7e" % (self[field][i]))
            f.write("\n")
        f.close()

    def generateCoords(self):
        """
        Basically a placeholder for the generated fields.
        """
        self.ActiveDimensions = (len(self.x), 1, 1)
        if self.coords == None:
            self.coords = array([self.x, self.y, self.z])
            print "Returning from generateCoords"

class EnzoRegion(Enzo3DData):
    def __init__(self, hierarchy, center, leftEdge, rightEdge, fields = None):
        Enzo3DData.__init__(self, hierarchy, center, fields)
        self.fields = ["Radius"] + self.fields
        self.leftEdge = leftEdge
        self.rightEdge = rightEdge
        self.radii = None # Still have radii
        self.refreshData()

    def getData(self, field = None):
        ind = where(sum(transpose(logical_or(greater(self.gridLeftEdge, rightEdge), \
                                               less(self.gridRightEdge, leftEdge)))) == 0)
        g = self.hierarchy.grids[ind]
        self.grids = g
        points = []
        if self.dx == None:
            for grid in g:
                #print "Generating coords for grid %s" % (grid.id)
                points.append(self.generateGridCoords(grid))
            t = concatenate(points)
            print t.shape
            self.x = t[:,0]
            self.y = t[:,1]
            self.z = t[:,2]
            self.radii = t[:,3]
            self.dx = t[:,4]
            self.dy = t[:,4]
            self.dz = t[:,4]
        if isinstance(field, types.StringType):
            fieldsToGet = [field]
        elif isinstance(field, types.ListType):
            fieldsToGet = field
        elif field == None:
            fieldsToGet = self.fields
        for field in fieldsToGet:
            if self.data.has_key(field):
                continue
            print "Getting field %s from %s" % (field, len(g))
            sf = zeros(self.x.shape[0],Float64)
            i = 0
            if field in self.hierarchy.fieldList:
                for grid in g:
                    #print "\tGetting %s from grid %s" % (field, grid.id)
                    ta = self.getDataFromGrid(grid, field)
                    sf[i:i+ta.shape[0]] = ta
                    i += ta.shape[0]
                self[field] = sf
            else:
                self.generateField(field)

    def getDataFromGrid(self, grid, field):
        #print "\tGetting data"
        if grid.myChildMask == None or grid.myChildIndices == None:
            #print "\tGenerating child mask"
            grid.generateChildMask()
        # First we find the cells that are within the sphere
        i0 = map(int, (self.leftEdge  - grid.LeftEdge) / grid.dx)
        i1 = map(int, map(ceil, (self.rightEdge - grid.LeftEdge) / grid.dx))
        i0 = choose(less(i0,0), (i0,0))
        i1 = choose(greater(i1,grid.ActiveDimensions-1), (i1,grid.ActiveDimensions-1))
        cutMask = zeros(grid.ActiveDimensions, Int64)
        cutMask[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] = 1
        pointI = where(logical_and(cutMask, grid.myChildMask==1) == 1)
        return grid[field][pointI]

    def generateGridCoords(self, grid):
        l = time.time()
        if grid.coords == None:
            grid.generateCoords()
        if grid.myChildMask == None or grid.myChildIndices == None:
            grid.generateChildMask()
        # First we find the cells that are within the sphere
        i0 = map(int, (self.leftEdge  - grid.LeftEdge) / grid.dx)
        i1 = map(int, map(ceil, (self.rightEdge - grid.LeftEdge) / grid.dx))
        i0 = choose(less(i0,0), (i0,0))
        i1 = choose(greater(i1,grid.ActiveDimensions-1), (i1,grid.ActiveDimensions-1))
        cutMask = zeros(grid.ActiveDimensions, Int64)
        cutMask[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] = 1
        pointI = where(logical_and(cutMask, grid.myChildMask==1) == 1)
        dx = ones(pointI[0].shape[0], Float64) * grid.dx
        tr = array([grid.coords[0,:][pointI], \
                grid.coords[1,:][pointI], \
                grid.coords[2,:][pointI], \
                grid["RadiusCode"][pointI],
                dx], Float64)
        tr.swapaxes(0,1)
        return tr


class EnzoSphere(Enzo3DData):
    def __init__(self, hierarchy, center, radius, fields = None):
        Enzo3DData.__init__(self, hierarchy, center, fields)
        self.fields = ["Radius"] + self.fields
        self.radius = radius
        self.radii = None
        self.refreshData()

    def getData(self, field = None):
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        g,ind = self.hierarchy.findSphereGrids(self.center, self.radius)
        self.grids = g
        points = []
        if self.dx == None:
            for grid in g:
                #print "Generating coords for grid %s" % (grid.id)
                points.append(self.generateGridCoords(grid))
            t = concatenate(points)
            print t.shape
            self.x = t[:,0]
            self.y = t[:,1]
            self.z = t[:,2]
            self.radii = t[:,3]
            self.dx = t[:,4]
            self.dy = t[:,4]
            self.dz = t[:,4]
        if isinstance(field, types.StringType):
            fieldsToGet = [field]
        elif isinstance(field, types.ListType):
            fieldsToGet = field
        elif field == None:
            fieldsToGet = self.fields
        for field in fieldsToGet:
            if self.data.has_key(field):
                continue
            print "Getting field %s from %s" % (field, len(g))
            sf = zeros(self.x.shape[0],Float64)
            i = 0
            if field in self.hierarchy.fieldList:
                for grid in g:
                    #print "\tGetting %s from grid %s" % (field, grid.id)
                    ta = self.getDataFromGrid(grid, field)
                    sf[i:i+ta.shape[0]] = ta
                    i += ta.shape[0]
                self[field] = sf
            else:
                self.generateField(field)

    def getDataFromGrid(self, grid, field):
        #print "\tGetting data"
        if grid.myChildMask == None or grid.myChildIndices == None:
            #print "\tGenerating child mask"
            grid.generateChildMask()
        # First we find the cells that are within the sphere
        pointI = where(logical_and((grid["RadiusCode"]<=self.radius),grid.myChildMask==1)==1)
        return grid[field][pointI]

    def generateGridCoords(self, grid):
        l = time.time()
        if grid.coords == None:
            grid.generateCoords()
        if grid.myChildMask == None or grid.myChildIndices == None:
            grid.generateChildMask()
        # First we find the cells that are within the sphere
        pointI = where(logical_and((grid["RadiusCode"]<=self.radius),grid.myChildMask==1)==1)
        dx = ones(pointI[0].shape[0], Float64) * grid.dx
        tr = array([grid.coords[0,:][pointI], \
                grid.coords[1,:][pointI], \
                grid.coords[2,:][pointI], \
                grid["RadiusCode"][pointI],
                dx], Float64)
        tr.swapaxes(0,1)
        return tr
