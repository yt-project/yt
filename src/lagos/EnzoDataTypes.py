"""
Various non-grid data containers.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@change: Sun Mar 11 18:10:57 PDT 2007 implemented EnzoProj
"""

from yt.lagos import *

class EnzoData:
    """
    Generic EnzoData container.  By itself, will attempt to
    generate field, read fields (method defined by derived classes)
    """
    def __init__(self, hierarchy, fields):
        """
        Sets up EnzoData instance

        @param hierarchy: hierarchy we're associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param fields: Fields represented in the data
        @type fields: list of strings
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

        @param field: the field to add
        @type field: string
        """
        if field not in self.fields:
            self.fields.append(field)
        if not self.data.has_key(field):
            self.getData(field)

    def generateField(self, fieldName):
        """
        Generates, or attempts to generate,  a field not found in the file

        See fields.py for more information.  fieldInfo.keys() will list all of
        the available derived fields.  Note that we also make available the
        suffices _Fraction and _Squared here.  All fields prefixed with 'k'
        will force an attempt to use the chemistry tables to generate them from
        temperature.  All fields used in generation will remain resident in
        memory.

        Also, it now works better with vectors.

        (should be moved to a standalone?  Or should EnzoGrid subclass EnzoData?)

        @param fieldName: field to generate
        @type fieldName: string
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
        elif fieldName.endswith("_vcomp"):
            baryonField = fieldName[:-8]
            index = int(fieldName[-7:-6])
            self[fieldName] = self[baryonField][index,:]
        elif fieldName.endswith("_tcomp"):
            baryonField = fieldName[:-9]
            ii = map(int, fieldName[-8:-6])
            self[fieldName] = self[baryonField][ii[0],ii[1],:]
        elif fieldName.endswith("Abs"):
            baryonField = fieldName[:-4]
            self[fieldName] = abs(self[baryonField])
        elif fieldInfo.has_key(fieldName):
            # We do a fallback to checking the fieldInfo dict
            # Note that it'll throw an exception here if it's not found...
            # ...which I'm cool with
            fieldInfo[fieldName][3](self, fieldName)
        elif fieldName.startswith("k"):
            self[fieldName] = abs(self.hierarchy.rates[self["Temperature"],fieldName])
        else:
            raise exceptions.KeyError, fieldName

    def writeOut(self, filename, header=False):
        """
        Writes out the generalized data to an ASCII file.  Probably quite slow.
        Note that we're doing everything in python, which is ... slow to do.
        Maybe I'll write a C-extension someday.

        @param filename: The filename to output
        @type filename: string
        @keyword header: Should we output a tab-separated header
        @type header: Boolean
        """
        # Now we first set up a new array, with all the stuff, to speed 
        # things up.
        data = []
        for field in self.fields:
            data.append(self[field])
        tw = na.array(data,nT.Float32)
        fs = "%0.6e\t" * len(data)
        #print fs, tw.shape
        f = open(filename,"w")
        if header:
            f.write("x\ty\tz\tdx\t" + "\t".join(self.fields) + "\n")
        for i in xrange(self.coords.shape[1]):
            if i % 1000 == 0:
                mylog.debug("Writing record %i / %i", i, self.coords.shape[1])
            f.write("%0.15e\t" * 3 % tuple(self.coords[:,i]))
            f.write("%0.15e\t" % (self.dx[i]))
            f.write(fs % tuple(tw[:,i]))
            f.write("\n")
        f.close()
        del data

class Enzo2DData(EnzoData):
    """
    Class to represent a set of EnzoData that's 2-D in nature.  Slices and
    projections, primarily.
    """
    def __init__(self, hierarchy, axis, fields):
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
        for i in xrange(self.x.shape[0]):
            f.write("%0.20f %0.20f %0.20f %0.20f" % \
                    (self.x[i], \
                     self.y[i], \
                     self.dx[i], \
                     self.dy[i]))
            for field in fields:
                f.write("\t%0.7e" % (self[field][i]))
            f.write("\n")
        f.close()

    def interpolateDiscretize(self, LE, RE, field, side, logSpacing=True):
        """
        This returns a uniform grid of points, interpolated using the nearest
        neighbor method.

        @note: Requires NumPy
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
        """
        points = na.where(na.logical_and( \
                            na.logical_and( \
                                self.x <= RE[0], \
                                self.x >= LE[0]),
                            na.logical_and( \
                                self.y <= RE[1], \
                                self.y >= LE[1])) == 1)
        """
        #xx = na.array(self.x[points], dtype=nT.Float64)
        #yy = na.array(self.y[points], dtype=nT.Float64)
        #zz = na.array(self[field][points], dtype=nT.Float64)
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

    def discretize(self, LE, RE, field, side, logSpacing=True):
        pass

class EnzoProj(Enzo2DData):
    """
    EnzoProj data.  Doesn't work yet.  Needs to be implemented, following the
    method of EnzoHierarchy.getProjection.

    """
    def __getitem__(self, key):
        """
        We override this because we don't want to add new fields for new keys.
        """
        return self.data[key]

    def __init__(self, hierarchy, axis, field, weightField = None, maxLevel = None):
        """
        Returns an instance of EnzoProj.  Note that this object will be fairly static.

        @todo: Implement
        @param hierarchy: the hierarchy we are projecting
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param axis: axis to project along
        @type axis: integer
        @param field: the field to project (NOT multiple)
        @type field: string
        @keyword weightField: the field to weight by
        @keyword maxLevel: the maximum level to project through
        """
        Enzo2DData.__init__(self, hierarchy, axis, field)
        if maxLevel == None:
            maxLevel = self.hierarchy.maxLevel
        self.maxLevel = maxLevel
        self.weightField = weightField
        self.coords = None
        self.calculateMemory()
        self.refreshData()

    def calculateMemory(self):
        """
        Here we simply calculate how much memory is needed, which speeds up
        allocation later
        """
        memoryPerLevel = {}
        totalProj = 0
        h = self.hierarchy
        i = 0
        for level in range(self.maxLevel+1):
            memoryPerLevel[level] = 0
            mylog.info("Working on level %s", level)
            grids = h.levelIndices[level]
            numGrids = len(grids)
            RE = h.gridRightEdge[grids].copy()
            LE = h.gridLeftEdge[grids].copy()
            for grid in h.grids[grids]:
                if (i%1e3) == 0:
                    mylog.debug("Reading and masking %s / %s", i, h.numGrids)
                # We unroll this so as to avoid instantiating a range() for
                # every grid
                grid.generateOverlapMasks(0, LE, RE)
                grid.myOverlapGrids[0] = h.grids[grids[na.where(grid.myOverlapMasks[0] == 1)]]
                grid.generateOverlapMasks(1, LE, RE)
                grid.myOverlapGrids[1] = h.grids[grids[na.where(grid.myOverlapMasks[1] == 1)]]
                grid.generateOverlapMasks(2, LE, RE)
                grid.myOverlapGrids[2] = h.grids[grids[na.where(grid.myOverlapMasks[2] == 1)]]
                myNeeds = grid.ActiveDimensions[(self.axis+1)%3]\
                         *grid.ActiveDimensions[(self.axis+2)%3]
                memoryPerLevel[level] += myNeeds
                totalProj += myNeeds
                i += 1
        for level in range(self.maxLevel+1):
            gI = na.where(h.gridLevels == level)
            mylog.debug("%s cells and %s grids for level %s", \
                memoryPerLevel[level], len(gI[0]), level)
        mylog.debug("We need %s cells total", totalProj)
        self.memoryPerLevel = memoryPerLevel

    def getData(self):
        """
        Unfortunately, projecting is more linear, less-OO than slicing, it
        seems to me.  We don't want to view this as a projection "operator" or
        anything like that, so this will be a fairly straightforward port of the old
        getProjection code.
        """
        # Note that we only look at the first field, since we know it is a list
        field = self.fields[0]
        i = 0
        zeroOut = True
        dataByLevel = {}
        totalGridsProjected = 0
        dbl_coarse = {}
        h = self.hierarchy
        weightField = self.weightField
        minLevel = 0 # We're not going to allow projections from selected
                     # levels right now
        fullLength = 0
        axis = self.axis # Speed up the access here by a miniscule amount
        dataFieldName = field
        for level in range(minLevel, self.maxLevel+1):
            zeroOut = (level != self.maxLevel)
            mylog.info("Projecting through level = %s", level)
            llen = self.memoryPerLevel[level]
            tempLevelData = [na.zeros(llen, nT.Int64), \
                             na.zeros(llen, nT.Int64), \
                             na.zeros(llen, nT.Float64), \
                             na.zeros(llen, nT.Float64)]
            myInd = na.where(h.gridLevels == level)[0]
            gridsToProject = h.grids[myInd]
            index = 0
            for grid in gridsToProject:
                i += 1
                grid.retVal = grid.getProjection(axis, field, zeroOut, weight=weightField)
            totalGridsProjected += i
            mylog.debug("Grid projecting done (%s / %s total) with %s points", \
                        totalGridsProjected, h.numGrids, index)
            mylog.debug("Combining level %s...", level)
            i = 0
            for grid1 in gridsToProject:
                i += 1
                if grid1.retVal[0].shape[0] == 0:
                    continue # We skip anything already processed
                for grid2 in grid1.myOverlapGrids[axis]:
                    if grid2.retVal[0].shape[0] == 0:
                        continue # Already processed
                    if grid1.id == grid2.id:
                        continue # Same grid!
                    index=EnzoCombine.CombineData( \
                            grid1.retVal[0], grid1.retVal[1], \
                            grid1.retVal[2], grid1.retVal[3], grid1.retVal[4], \
                            grid2.retVal[0], grid2.retVal[1], \
                            grid2.retVal[2], grid2.retVal[3], grid2.retVal[4], 0)
                    goodI = na.where(grid2.retVal[0] > -1)
                    grid2.retVal[0] = grid2.retVal[0][goodI].copy()
                    grid2.retVal[1] = grid2.retVal[1][goodI].copy()
                    grid2.retVal[2] = grid2.retVal[2][goodI].copy()
                    grid2.retVal[3] = grid2.retVal[3][goodI].copy()
                    grid2.retVal[4] = grid2.retVal[4][goodI].copy()
                numRefined = 0
                if (level > minLevel) and (level <= self.maxLevel):
                    for grid2 in grid1.Parent.myOverlapGrids[axis]:
                        if grid2.coarseData[0].shape[0] == 0:
                            continue # Already refined
                        numRefined += EnzoCombine.RefineCoarseData( \
                            grid1.retVal[0], grid1.retVal[1], \
                            grid1.retVal[2], grid1.retVal[4], \
                            grid2.coarseData[0], grid2.coarseData[1], \
                            grid2.coarseData[2], grid2.coarseData[4], 2)
            all_data = [[],[],[],[],[]]
            for grid in gridsToProject:
                all_data[0].append(grid.retVal[0])
                all_data[1].append(grid.retVal[1])
                all_data[2].append(grid.retVal[2])
                all_data[3].append(grid.retVal[3])
                all_data[4].append(grid.retVal[4])
                cI = na.where(grid.retVal[3]==0)
                grid.coarseData = [grid.retVal[0][cI], \
                                   grid.retVal[1][cI], \
                                   grid.retVal[2][cI], \
                                   grid.retVal[3][cI], \
                                   grid.retVal[4][cI]]
            levelData = []
            levelData.append(na.concatenate(all_data[0]))
            levelData.append(na.concatenate(all_data[1]))
            levelData.append(na.concatenate(all_data[2]))
            levelData.append(na.concatenate(all_data[3]))
            levelData.append(na.concatenate(all_data[4]))
            mylog.debug("All done combining and refining with a final %s points", levelData[0].shape[0])
            dx=gridsToProject[0].dx # Assume uniform dx
            dblI = na.where(na.logical_and((levelData[0]>-1), (levelData[3] == 1))==1)
            if weightField != None:
                weightedData = levelData[2][dblI] / levelData[4][dblI]
            else:
                weightedData = levelData[2][dblI]
            dataByLevel[level] = [levelData[0][dblI], levelData[1][dblI], \
                                  weightedData, dx]
            time4 = time.time()
            if dataByLevel[level][0].shape[0] > 0:
                exampleValue = dataByLevel[level][2][0]
            else:
                exampleValue = 0.0
            mylog.debug("Level %s done: %s final of %s (%s)", \
                       level, dataByLevel[level][0].shape[0], \
                       levelData[0].shape[0], exampleValue)
            fullLength += dataByLevel[level][0].shape[0]
            del levelData
        # Now we construct a simple dictionary
        all_x = []
        all_y = []
        all_z = []
        all_dx = []
        for level in range(minLevel, self.maxLevel+1):
            all_x.append(dataByLevel[level][0])
            all_y.append(dataByLevel[level][1])
            all_z.append(dataByLevel[level][2])
            all_dx.append(na.ones(dataByLevel[level][2].shape)*dataByLevel[level][3])
        # We now convert to half-widths and center-points
        self.dx = na.concatenate(all_dx) / 2.0
        self.dy = na.concatenate(all_dx) / 2.0
        self.x = (na.concatenate(all_x)+0.5) * self.dx*2.0
        self.y = (na.concatenate(all_y)+0.5) * self.dy*2.0
        self.data[field] = na.concatenate(all_z)
        

class EnzoSlice(Enzo2DData):
    """
    A slice at a given coordinate along a given axis through the entire grid
    hierarchy.
    """
    def __init__(self, hierarchy, axis, coord, fields = None, center=None, fRet=False):
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
        Enzo2DData.__init__(self, hierarchy, axis, fields)
        self.center = center
        self.coord = coord
        self.coords = None
        self.fRet = fRet
        self.refreshData()

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

    def generateCoords(self):
        """
        Place holder, so that generating fields doesn't die.  Note that
        coordinates are generated at initialization.
        """
        pass

    def getData(self, field = None):
        """
        Iterates over the list of fields and generates/reads them all.
        Probably shouldn't be called directly.

        @keyword field: field to add (list or single)
        @type field: string or list of strings
        """
        # We get it for the values in fields and coords
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        g,ind = self.hierarchy.findSliceGrids(self.coord, self.axis)
        self.grids = g
        #print g
        points = []
        if self.dx == None:
            for grid in g:
                #print "Generating coords for grid %s" % (grid.id)
                points.append(self.generateGridCoords(grid))
            t = na.concatenate(points)
            self.x = t[:,0]
            self.y = t[:,1]
            self.z = t[:,2]
            self.dx = t[:,3]
            self.dy = t[:,3]
            self.dz = t[:,3]
            #self.coords = na.array([self.x, self.y, self.z])
            self.coords = na.zeros((3,len(self.x)),nT.Float64)
            #print self.coords.shape
            self.coords[x_dict[self.axis],:] = t[:,0]
            self.coords[y_dict[self.axis],:] = t[:,1]
            self.coords[self.axis,:] = t[:,2]
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
                    self[field] = na.concatenate(rvs)
            else:
                i = 0
                self.generateField(field)

    def generateGridCoords(self, grid):
        """
        Returns coordinates (x,y,z,dx) for a single grid.

        @param grid: grid to generate coords for
        @type grid: L{EnzoGrid<EnzoGrid>}
        """
        xaxis = x_dict[self.axis]
        yaxis = y_dict[self.axis]
        if grid.myChildMask == None:
            #print "Generating child mask"
            grid.generateChildMask(fRet=self.fRet)
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

    def getDataFromGrid(self, grid, field):
        """
        Reads a slice from a single grid

        Could probably be more efficient, possibly by accepting a list of
        fields.

        @param grid: grid to slice
        @type grid: L{EnzoGrid<EnzoGrid>}
        @param field: field to read
        @type field: string
        """
        if grid.myChildMask == None:
            #print "Generating child mask"
            grid.generateChildMask(fRet=self.fRet)
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
        # The following bit breaks the finest grid.  Hopefully we won't suffer
        # too much of a performance hit by always reading in the data.
        #if not(grid.data.has_key(field)):
            #print "Reading slice", grid
            #dv = self.readDataSlice(grid, field, slHDF)
        #else:
            #print "We seem to have it already...", grid
            #dv = grid[field][slHERE]
        dv = self.readDataSlice(grid, field, slHDF)
        #print dv.flat.shape, grid.myChildMask[slHERE].flat.shape, grid.ActiveDimensions
        cm = na.where(grid.myChildMask[slHERE].ravel() == 1)
        #print field, cm[0].shape, dv.shape
        dv=dv.swapaxes(0,2)
        dataVals = dv.ravel()[cm]
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

        @param hierarchy: hierarchy we are associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param center: center of the region
        @type center: array of floats
        @param fields: fields to read/generate
        @type fields: list of strings
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
            self.coords = na.array([self.x, self.y, self.z])

    def makeProfile(self, fields, nBins, rInner, rOuter):
        """
        Here we make a profile.  Note that, for now, we do log-spacing in the
        bins, and we also assume we are given the radii in code units.
        field_weights defines the weighting for a given field.  If not given,
        assumed to be mass-weighted.
        
        Units might be screwy, but I don't think so.  Unless otherwise stated,
        returned in code units.  CellMass, by the way, is Msun.  And
        accumulated.

        Returns the bin outer radii and a dict of the profiles.

        @param fields: fields to get profiles of
        @type fields: list of strings
        @param nBins: number of bins
        @type nBins: int
        @param rInner: inner radius, code units
        @param rOuter: outer radius, code units
        """
        if not isinstance(fields, types.ListType):
            fields = [fields]
        rOuter = min(rOuter, self.findMaxRadius)
        # Let's make the bins
        lrI = log10(rInner)
        lrO = log10(rOuter)
        st = (lrO-lrI)/(nBins-1)
        bins = 10**(na.arange(lrI, lrO+st, st))
        radiiOrder = na.argsort(self["RadiusCode"])
        fieldCopies = {} # We double up our memory usage here for sorting
        radii = self["RadiusCode"][radiiOrder]
        #print radii.max()
        #print radii.min()
        #print radii.shape
        #print "BINS!", bins
        binIndices = na.searchsorted(bins, radii)
        nE = self["RadiusCode"].shape[0]
        #defaultWeight = na.ones(nE, nT.Float32)
        defaultWeight = self["CellMass"][radiiOrder]
        fieldProfiles = {}
        if "CellMass" not in fields:
            fields.append("CellMass")
        for field in fields:
            fc = self[field][radiiOrder]
            fp = na.zeros(nBins,nT.Float32)
            if field_weights.has_key(field):
                if field_weights[field] == -999:
                    ss = "Accumulation weighting"
                    weight = na.ones(nE, nT.Float32) * -999
                elif field_weights[field] != None:
                    ww = field_weights[field]
                    ss="Weighting with %s" % (ww)
                    weight = self[ww][radiiOrder]
                elif field_weights[field] == None:
                    ss="Not weighted"
                    weight = na.ones(nE, nT.Float32)
                else:
                    mylog.warning("UNDEFINED weighting for %s", field)
                    ss="Undefined weighting"
            else:
                ss="Weighting with default"
                weight = defaultWeight
            mylog.info("Getting profile of %s (%s)", field, ss)
            EnzoCombine.BinProfile(fc, binIndices, \
                                   fp, weight)
            fieldProfiles[field] = fp
        #return bins
        #print rOuter, rInner, st, nBins, bins, bins[:nBins]
        fieldProfiles["OuterRadius"] = bins[:nBins]
        co = AnalyzeClusterOutput(fieldProfiles)
        return co

class EnzoRegion(Enzo3DData):
    def __init__(self, hierarchy, center, leftEdge, rightEdge, fields = None):
        """
        Initializer for a rectangular prism of data.

        @note: Center does not have to be (rightEdge - leftEdge) / 2.0
        """
        Enzo3DData.__init__(self, hierarchy, center, fields)
        self.fields = ["Radius"] + self.fields
        self.leftEdge = leftEdge
        self.rightEdge = rightEdge
        self.radii = None # Still have radii
        self.refreshData()

    def findMaxRadius(self):
        RE = na.array(self.rightEdge)
        LE = na.array(self.leftEdge)
        return (RE-LE).min()

    def getData(self, field = None):
        ind = na.where(na.sum(na.transpose(na.logical_or(na.greater(self.hierarchy.gridLeftEdge, self.rightEdge), \
                                               na.less(self.hierarchy.gridRightEdge, self.leftEdge)))) == 0)
        g = self.hierarchy.grids[ind]
        self.grids = g
        points = []
        if self.dx == None:
            for grid in g:
                #print "Generating coords for grid %s" % (grid.id)
                c = grid.center
                grid.center = self.center
                points.append(self.generateGridCoords(grid))
                t = na.concatenate(points)
                grid.center = c
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
            mylog.info("Getting field %s from %s", field, len(g))
            sf = na.zeros(self.x.shape[0],nT.Float64)
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
        i0 = na.choose(na.less(i0,0), (i0,0))
        i1 = na.choose(na.greater(i1,grid.ActiveDimensions-1), (i1,grid.ActiveDimensions-1))
        cutMask = na.zeros(grid.ActiveDimensions, nT.Int64)
        cutMask[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] = 1
        pointI = na.where(na.logical_and(cutMask, grid.myChildMask==1) == 1)
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
        i0 = na.choose(na.less(i0,0), (i0,0))
        i1 = na.choose(na.greater(i1,grid.ActiveDimensions-1), (i1,grid.ActiveDimensions-1))
        cutMask = na.zeros(grid.ActiveDimensions, nT.Int64)
        cutMask[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] = 1
        pointI = na.where(na.logical_and(cutMask, grid.myChildMask==1) == 1)
        dx = na.ones(pointI[0].shape[0], nT.Float64) * grid.dx
        tr = na.array([grid.coords[0,:][pointI], \
                grid.coords[1,:][pointI], \
                grid.coords[2,:][pointI], \
                grid["RadiusCode"][pointI],
                dx], nT.Float64).swapaxes(0,1)
        return tr

class EnzoDataCube(Enzo3DData):
    def __init__(self, hierarchy, level, center, dim, fields = None):
        """
        This returns a cube of equal-resolution data, with 
        dimensions dim x dim x dim .

        Note that we typically generate one, then write it out.  We do not
        expect to be modifying them once created -- so I am not really
        implementing the normal set of niceties.

        @param hierarchy: the EnzoHierarchy to associate with
        @type hierarchy: L{EnzoHierarchy}
        @param level: the level we are going to set our data to be at
        @param center: the center point
        @type center: sequence of floats
        @param dim: the number of cells on a side
        @type dim: int
        @param fields: fields to snagify
        @type fields: list of strings
        """

        Enzo3DData.__init__(self, hierarchy, center, fields)
        self.fields = ["Radius"] + self.fields
        dx = hierarchy.gridDxs[hierarchy.levelIndices[level][0]]
        self.leftEdge  = center - dx*(dim/2.0)
        self.rightEdge = center + dx*(dim/2.0)
        self.level = level
        self.dx = dx
        self.refreshData()

    def getData(self, field = None):
        ind = na.where(\
            na.sum(na.transpose(na.logical_or(na.greater(self.hierarchy.gridLeftEdge, self.rightEdge), \
                                              na.less(self.hierarchy.gridRightEdge, self.leftEdge)))) == 0)
        g = self.hierarchy.grids[ind]
        # I think it's faster if we snag all the points, then cut at the end.
        # This may not be true for the root grid.  But I don't care.  Again,
        # this is not supposed to be a super-fast process, or one that's done
        # interactively.  (Very often, at least.)
        vals = {}
        for field in self.fields + ['x','y','z']:
            vals[field] = []
        for grid in g:
            if grid.Level > self.level:
                continue
            print grid
            tr = grid.atResolution(self.level, self.fields)
            for field in self.fields + ['x', 'y', 'z']:
                vals[field].append(tr[field])
        for field in self.fields + ['x','y','z']:
            vals[field] = na.concatenate(vals[field])
        self.data = vals

    def writeOut(self, filename, header=False):
        self.coords = na.array([self['x'], self['y'], self['z']])
        self.dx = na.ones(len(self['x']), dtype=self['x'].dtype)*self.dx
        EnzoData.writeOut(self, filename, header)

class EnzoSphere(Enzo3DData):
    """
    A sphere of points
    """
    def __init__(self, hierarchy, center, radius, fields = None):
        """
        @param hierarchy: hierarchy we are associated with
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param center: center of the region
        @type center: array of floats
        @param radius: radius of the sphere
        @type radius: float
        @keyword fields: fields to read/generate
        @type fields: list of strings
        """
        Enzo3DData.__init__(self, hierarchy, center, fields)
        self.fields = ["Radius"] + self.fields
        self.radius = radius
        self.radii = None
        self.refreshData()

    def findMaxRadius(self):
        return self.radius

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
                c = grid.center
                grid.center = self.center
                points.append(self.generateGridCoords(grid))
                grid.center = c
            t = na.concatenate(points)
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
            mylog.info("Getting field %s from %s", field, len(g))
            sf = na.zeros(self.x.shape[0],nT.Float64)
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
        pointI = na.where(na.logical_and((grid["RadiusCode"]<=self.radius),grid.myChildMask==1)==1)
        return grid[field][pointI]

    def generateGridCoords(self, grid):
        l = time.time()
        if grid.coords == None:
            grid.generateCoords()
        if grid.myChildMask == None or grid.myChildIndices == None:
            grid.generateChildMask()
        # First we find the cells that are within the sphere
        pointI = na.where(na.logical_and((grid["RadiusCode"]<=self.radius),grid.myChildMask==1)==1)
        dx = na.ones(pointI[0].shape[0], nT.Float64) * grid.dx
        tr = na.array([grid.coords[0,:][pointI], \
                grid.coords[1,:][pointI], \
                grid.coords[2,:][pointI], \
                grid["RadiusCode"][pointI],
                dx], nT.Float64).swapaxes(0,1)
        return tr
