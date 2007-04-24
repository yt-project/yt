"""
Enzo hierarchy container class

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Either subclass EnzoHierarchy from EnzoParameterFile, or get rid of
overlap and make Hierarchy an option property of EnzoParameterFile
"""

from yt.lagos import *
import string, re, gc
#import EnzoCombine
import yt.enki

class EnzoParameterFile:
    """
    This class is a stripped down class that simply reads and parses, without
    looking at the hierarchy.
    """
    def __init__(self, filename, hdf_version = None):
        """
        @note: We disregard hdf_version here
        """
        if filename.endswith(".hierarchy"):
            filename = filename[:-10]
        self.parameterFilename = "%s" % (filename)
        self.basename = os.path.basename(filename)
        self.directory = os.path.dirname(filename)
        self.fullpath = os.path.abspath(self.directory)
        if len(self.directory) == 0:
            self.directory = "."
        self.conversionFactors = {}
        self.parameters = {}
        self.parameters["CurrentTimeIdentifier"] = \
            int(os.stat(self.parameterFilename)[ST_CTIME])
        self.parseParameterFile()
        self.setUnits()
        
    def __getitem__(self, key):
        """
        Returns units, parameters, or conversionFactors in that order
        """
        if self.units.has_key(key):
            return self.units[key]
        elif self.parameters.has_key(key):
            return self.parameters[key]
        return self.conversionFactors[key]

    def __repr__(self):
        return self.basename

    def keys(self):
        """
        Returns a list of possible keys, from units, parameters and
        conversionFactors
        """
        return self.units.keys() \
             + self.parameters.keys() \
             + self.conversionFactors.keys()

    def has_key(self, key):
        """
        Returns true or false
        """
        return (self.units.has_key(key) or \
                self.parameters.has_key(key) or \
                self.conversionFactors.has_key(key))

    def parseParameterFile(self):
        """
        Parses the parameter file and establishes the various
        dictionaries.
        """
        # Let's read the file
        lines = open(self.parameterFilename).readlines()
        for lineI in xrange(len(lines)):
            line = lines[lineI]
            if len(line) < 2:
                continue
            param, vals = map(strip,map(rstrip,line.split("=")))
            if parameterDict.has_key(param):
                t = map(parameterDict[param], vals.split())
                if len(t) == 1:
                    self.parameters[param] = t[0]
                else:
                    self.parameters[param] = t
                if param.endswith("Units"):
                    dataType = param[:-5]
                    self.conversionFactors[dataType] = self.parameters[param]
            elif param.startswith("#DataCGS"):
                # Assume of the form: #DataCGSConversionFactor[7] = 2.38599e-26 g/cm^3
                if lines[lineI-1].find("Label") >= 0:
                    kk = lineI-1
                elif lines[lineI-2].find("Label") >= 0:
                    kk = lineI-2
                dataType = lines[kk].split("=")[-1].rstrip().strip()
                convFactor = float(line.split("=")[-1].split()[0])
                self.conversionFactors[dataType] = convFactor
            elif param.startswith("#CGSConversionFactor"):
                dataType = param[20:].rstrip()
                convFactor = float(line.split("=")[-1])
                self.conversionFactors[dataType] = convFactor

    def setUnits(self):
        """
        Generates the conversion to various physical units based on the parameter file
        """
        self.units = {}
        if len(self.parameters) == 0:
            self.parseParameterFile()
        if self["ComovingCoordinates"]:
            z = self["CosmologyCurrentRedshift"]
            boxh = self["CosmologyComovingBoxSize"]
            self.units['aye']  = (1.0 + self["CosmologyInitialRedshift"])/(z - 1.0)
            if not self.has_key("Time"):
                LengthUnit = 3.086e24 * boxh / self["CosmologyHubbleConstantNow"] \
                             / (1+self["CosmologyInitialRedshift"])
                self.conversionFactors["Time"] = LengthUnit / self["x-velocity"]
        elif self.has_key("LengthUnits"):
            # We are given LengthUnits, which is number of cm per box length
            # So we convert that to box-size in Mpc
            z = 0
            boxh = 3.24077e-25 * self["LengthUnits"]
            self.units['aye']  = 1.0
        else:
            z = 0
            boxh = 1.0
            self.units['aye'] = 1.0
        seconds = self["Time"]
        box = boxh/(1+z)
        # Units are all defined in ravenDefs, thus making it easy to add new
        # na.ones!  (Not that there are really that many more to add...
        for unit in unitList.keys():
            self.units[unit] = unitList[unit] * box
            self.units[unit+'h'] = unitList[unit] * boxh
        self.units['1']     = 1
        self.units['years'] = seconds / (365*3600*24.0)
        self.units['days']  = seconds / (3600*24.0)

class EnzoHierarchy(EnzoParameterFile):
    """
    Class for handling Enzo timestep outputs
    """
    def __init__(self, filename, hdf_version=4):
        """
        Returns a new instance of EnzoHierarchy

        
        @param filename:  the filename of the parameter file
        @type filename: string
        @keyword hdf_version: either 4 or 5, depending
        """
        EnzoParameterFile.__init__(self, filename)
        # For now, we default to HDF4, but allow specifying HDF5
        if hdf_version == 5:
            EnzoGrid.readDataFast = readDataHDF5
            EnzoGrid.readAllData = readAllDataHDF5
            EnzoSlice.readDataSlice = readDataSliceHDF5
            EnzoGrid.getFields = getFieldsHDF5
            warnings.filterwarnings("ignore",".*PyTables format.*")
        else:
            EnzoGrid.readDataFast = readDataHDF4
            EnzoGrid.readAllData = readAllDataHDF4
            EnzoSlice.readDataSlice = readDataSliceHDF4
            EnzoGrid.getFields = getFieldsHDF4
        time1=time.time()
        # Expect filename to be the name of the parameter file, not the
        # hierarchy
        self.hierarchyFilename = "%s.hierarchy" % (filename)
        self.boundaryFilename = "%s.boundary" % (filename)
        # Now we do a bit of a hack to figure out how many grids there are,
        # so we can pre-create our arrays of dimensions, indices, etc
        self.hierarchyLines = open(self.hierarchyFilename).readlines()
        self.hierarchyString = open(self.hierarchyFilename).read()
        #re_gridID = re.compile("Grid\s=\s\d*", re.M)
        for i in xrange(len(self.hierarchyLines)-1,0,-1):
            line = self.hierarchyLines[i]
            if line.startswith("Grid ="):
                self.numGrids = int(line.split("=")[-1])
                break
        if self.parameters.has_key("CompilerPrecision") \
           and (self.parameters["CompilerPrecision"] == "r8" \
                or self.parameters["CompilerPrecision"] == "r4"):
            EnzoFloatType = nT.Float32
        else:
            EnzoFloatType = nT.Float64

        self.gridDimensions = na.zeros((self.numGrids,3), nT.Int32)
        self.gridStartIndices = na.zeros((self.numGrids,3), nT.Int32)
        self.gridEndIndices = na.zeros((self.numGrids,3), nT.Int32)
        self.gridLeftEdge = na.zeros((self.numGrids,3), EnzoFloatType)
        self.gridRightEdge = na.zeros((self.numGrids,3), EnzoFloatType)
        self.gridLevels = na.zeros((self.numGrids,1), nT.Int32)
        self.gridDxs = na.zeros((self.numGrids,1), EnzoFloatType)
        self.gridTimes = na.zeros((self.numGrids,1), nT.Float64)
        self.gridNumberOfParticles = na.zeros((self.numGrids,1))
        gg = []
        for fnI in xrange(self.numGrids):
            gg.append(EnzoGrid(self, fnI+1))
        self.grids = obj.array(gg)
        del gg
        self.gridReverseTree = [None] * self.numGrids
        self.gridTree = []
        for i in xrange(self.numGrids):
            self.gridTree.append([])

        # Now some statistics:
        #   0 = number of grids
        #   1 = number of cells
        #   2 = blank
        self.levelsStats = na.zeros((MAXLEVEL,3), nT.Int32)
        for i in xrange(MAXLEVEL):
            self.levelsStats[i,2] = i

        # For use with radial plots
        self.center = None
        self.bulkVelocity = None

        # For SWIG
        self.eiTopGrid = None
        
#        #l=raw_input("Hey dumbface, stopping here to check the memory")
        self.populateHierarchy()
        time2=time.time()
        mylog.info("Took %s seconds to read the hierarchy.", time2-time1)

    def __del__(self):
        """
        Let's see if we can delete some stuff here!
        """
        #print weakref.getweakrefcount(self)
        del self.eiTopGrid
        del self.gridReverseTree
        #del self.grids
        del self.gridLeftEdge, self.gridRightEdge
        del self.gridLevels, self.gridStartIndices, self.gridEndIndices
        del self.gridTimes, self.hierarchyString, self.hierarchyLines
        #del self.grids
        for gridI in xrange(self.numGrids):
            for g in self.gridTree[gridI]:
                #print "Deleting", g
                del g
        del self.gridTree

    def clearAllGridReferences(self):
        return
        del self.grids
        #for gridI in xrange(self.numGrids):
            #for g in self.gridTree[gridI]:
                #print "Deleting", g
                #g.clearAllGridReferences()
                #del g.hierarchy
                #del g
        #print "Grid 20:", gc.get_referrers(aa)
        #print "Self:", gc.get_referrers(self)

    def populateHierarchyRE(self):
        """
        This is so much slower than the other method.  Ugh.  How unfair!
        """
        for i in xrange(self.numGrids):
            self.grids[i] = EnzoGrid(self, i+1)

        re_Grid = constructRegularExpressions("Grid",('d'))

        re_GridDimension = constructRegularExpressions("GridDimension", ('d','d','d'))
        t = re.findall(re_GridDimension, self.hierarchyString)
        self.gridDimensions[:] = \
             na.array([[int(di[0]), int(di[1]), int(di[2])] for di in t], nT.Int32)

        re_GridStartIndex = constructRegularExpressions("GridStartIndex", ('d','d','d'))
        t = re.findall(re_GridStartIndex, self.hierarchyString)
        self.gridStartIndices[:] =\
             na.array([[int(di[0]), int(di[1]), int(di[2])] for di in t], nT.Int32)

        re_GridEndIndex = constructRegularExpressions("GridEndIndex", ('d','d','d'))
        t = re.findall(re_GridEndIndex, self.hierarchyString)
        self.gridEndIndices[:] =\
             na.array([[int(di[0]), int(di[1]), int(di[2])] for di in t], nT.Int32)
            
        re_GridLeftEdge = constructRegularExpressions("GridLeftEdge", ('g','g','g'))
        t = re.findall(re_GridLeftEdge, self.hierarchyString)
        self.gridLeftEdge[:] =\
             na.array([[float(di[0]), float(di[1]), float(di[2])] for di in t], nT.Float64)
            
        re_GridRightEdge = constructRegularExpressions("GridRightEdge", ('g','g','g'))
        t = re.findall(re_GridRightEdge, self.hierarchyString)
        self.gridRightEdge[:] =\
             na.array([[float(di[0]), float(di[1]), float(di[2])] for di in t], nT.Float64)
            
        re_GridLevels = constructRegularExpressions("^Level", ('d'))
        t = re.findall(re_GridLevels, self.hierarchyString)
        self.gridLevels[:] =\
             na.array([int(di) for di in t], nT.Int32)
            
#        re_GridTimes = constructRegularExpressions("Time", ('g'))
#        t = re.findall(re_GridTimes, self.hierarchyString)
#        self.gridTimes[:] =\
#             na.array([[float(di[0]), float(di[1]), float(di[2])] for di in t], nT.Float64)

        re_BaryonFileName = constructRegularExpressions("BaryonFileName",('s'))
        t = re.findall(re_BaryonFileName, self.hierarchyString)
        for fnI in xrange(len(t)):
            self.grids[fnI].setFilename(t[fnI])

        re_GridNumberOfParticles = constructRegularExpressions("NumberOfParticles", ('d'))
        t = re.findall(re_GridNumberOfParticles, self.hierarchyString)
        self.gridNumberOfParticles[:] =\
             na.array([int(di) for di in t], nT.Int32)
            
    def populateHierarchy(self):
        """
        Instantiates all of the grid objects, with their appropriate parameters
        """
        # Now, can we do this cleverly?
        # Let's do it the unclever way, I suppose...
        time12=time.time()
        # First, we look to see if the pseudo-pickled file is available
        arrayFilename = os.path.join(self.directory,self.basename+".arrayfile")
        if os.path.exists(arrayFilename) and na.__name__ == "numarray":
#            loki = raw_input("Pausing")
            allArrays = na.fromfile(arrayFilename, type=nT.Float64, shape=(self.numGrids,18))
#            loki = raw_input("Loaded")
            self.gridDimensions[:] = allArrays[:,0:3]
            self.gridStartIndices[:] = allArrays[:,3:6]
            self.gridEndIndices[:] = allArrays[:,6:9]
            self.gridLeftEdge[:] = allArrays[:,9:12]
            self.gridRightEdge[:] = allArrays[:,12:15]
            self.gridLevels[:] = allArrays[:,15:16]
            self.gridTimes[:] = allArrays[:,16:17]
            self.gridNumberOfParticles[:] = allArrays[:,17:18]
            del allArrays
            # Now get the baryon filenames
            re_BaryonFileName = constructRegularExpressions("BaryonFileName",('s'))
            t = re.findall(re_BaryonFileName, self.hierarchyString)
#            loki = raw_input("Instantiating")
            for fnI in xrange(len(t)):
                #self.grids[fnI] = EnzoGrid(self, fnI+1)
                self.grids[fnI].setFilename(t[fnI])
#            loki = raw_input("Done")
        else:
            for line in self.hierarchyLines:
                # We can do this the slow, 'reliable' way by stripping
                # or we can manually pad all our strings, which speeds it up by a
                # factor of about ten
                #param, vals = map(strip,line.split("="))
                if len(line) < 2:
                    continue
                i = line.index("=")
                param, vals = line.split("=")
                #param, vals = line[:i], line[i+1:]
                if param == "Grid ":
                    curGrid = int(vals)
                    self.grids[curGrid-1] = EnzoGrid(self, curGrid)
                elif param == "GridDimension     ":
                    splitConvertGridParameter(vals, float, self.gridDimensions, curGrid)
                elif param == "GridStartIndex    ":
                    splitConvertGridParameter(vals, int, self.gridStartIndices, curGrid)
                elif param == "GridEndIndex      ":
                    splitConvertGridParameter(vals, int, self.gridEndIndices, curGrid)
                elif param == "GridLeftEdge      ":
                    splitConvertGridParameter(vals, float, self.gridLeftEdge, curGrid)
                elif param == "GridRightEdge     ":
                    splitConvertGridParameter(vals, float, self.gridRightEdge, curGrid)
                elif param == "Level             ":
                    splitConvertGridParameter(vals, int, self.gridLevels, curGrid)
                elif param == "Time              ":
                    splitConvertGridParameter(vals, float, self.gridTimes, curGrid)
                elif param == "NumberOfParticles   ":
                    splitConvertGridParameter(vals, float, self.gridNumberOfParticles, curGrid)
                elif param == "BaryonFileName ":
                    self.grids[curGrid-1].setFilename(vals[1:-1])
            #self.populateHierarchyRE()
            mylog.info("Creating allArray to dump to file")
            allArrays = na.zeros((self.numGrids,18),nT.Float64)
            allArrays[:,0:3] = self.gridDimensions[:]
            allArrays[:,3:6] = self.gridStartIndices[:]
            allArrays[:,6:9] = self.gridEndIndices[:]
            allArrays[:,9:12] = self.gridLeftEdge[:]
            allArrays[:,12:15] = self.gridRightEdge[:]
            allArrays[:,15:16] = self.gridLevels[:]
            allArrays[:,16:17] = self.gridTimes[:]
            allArrays[:,17:18] = self.gridNumberOfParticles[:]
            mylog.info("Outputting allArray to %s", arrayFilename)
            try:
                allArrays.tofile(arrayFilename)
            except:
                mylog.error("There was an error writing to the file.  Skipping.")
            del allArrays
        time13=time.time()
        #print "Took %s" % (time13-time12)
        self.grids[0].Level = 0
        self.gridLevels[0] = 0
        time10=time.time()
        p = re.compile(r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)$", re.M)
        for m in p.finditer(self.hierarchyString):
            secondGrid = int(m.group(3))-1
            if secondGrid == -1:
                continue
            firstGrid = int(m.group(1))-1
            if m.group(2) == "Next":
                self.gridTree[firstGrid].append(weakref.proxy(self.grids[secondGrid]))
                self.gridReverseTree[secondGrid] = firstGrid + 1
                self.grids[secondGrid].Level = self.grids[firstGrid].Level + 1
                self.gridLevels[secondGrid] = self.gridLevels[firstGrid] + 1
            elif m.group(2) == "This":
                parent = self.gridReverseTree[firstGrid]
                if parent:
                    self.gridTree[parent-1].append(weakref.proxy(self.grids[secondGrid]))
                self.gridReverseTree[secondGrid] = parent
                self.grids[secondGrid].Level = self.grids[firstGrid].Level
                self.gridLevels[secondGrid] = self.gridLevels[firstGrid]
#        loki = raw_input("All done with that shit!")
        time11=time.time()
        #print "Took %s seconds" % (time11-time10)
        self.maxLevel = self.gridLevels.max()
        # Now we do things that we need all the grids to do
        self.fieldList = self.grids[0].getFields()
        time14=time.time()
        for i in xrange(self.numGrids):
            self.levelsStats[self.gridLevels[i,0],0] += 1
            self.grids[i].prepareGrid()
#            self.levelsStats[self.gridLevels[i,0],1] += na.product(self.grids[i].ActiveDimensions)
#        loki = raw_input("All done with that shit (2)!")
        self.levelIndices = {}
        self.levelNum = {}
        time15=time.time()
        #print "Took %s" % (time15-time14)
        for level in xrange(self.maxLevel+1):
            self.levelIndices[level] = self.selectLevel(level)
            self.levelNum = len(self.levelIndices[level])

    def selectLevel(self, level):
        """
        Returns a list of indices of EnzoHierarchy.grids at the specified level

        @param level: the level
        @type level: integer
        """
        # We return a numarray of the indices of all the grids on a given level
        indices = na.where(self.gridLevels[:,0] == level)[0]
        return indices

    def getSmallestDx(self):
        for i in xrange(MAXLEVEL):
            if (self.levelsStats[i,0]) == 0:
                break
            dx = self.gridDxs[self.levelIndices[i][0]]
        return dx[0]

    def printStats(self):
        """
        Prints out relevant information about the simulation
        """
        for i in xrange(MAXLEVEL):
            if (self.levelsStats[i,0]) == 0:
                break
            print "% 3i\t% 6i\t% 11i" % \
                  (i, self.levelsStats[i,0], self.levelsStats[i,1])
            dx = self.gridDxs[self.levelIndices[i][0]]
        print "-" * 28
        print "   \t% 6i\t% 11i" % (self.levelsStats[:,0].sum(), self.levelsStats[:,1].sum())
        print "\n"
        try:
            print "z = %0.8f" % (self["CosmologyCurrentRedshift"])
        except:
            pass
        t_s = self["InitialTime"] * self["Time"]
        print "t = %0.8e = %0.8e s = %0.8e years" % \
            (self["InitialTime"], \
             t_s, t_s / (365*24*3600.0) )
        print "\nSmallest Cell:"
        u=[]
        for item in self.units.items():
            u.append((item[1],item[0]))
        u.sort()
        for unit in u:
            print "\tWidth: %0.3e %s" % (dx*unit[0], unit[1])

    def findPoint(self, coord):
        """
        Returns the objects, indices of grids containing a point

        @param coord: three floats
        @type coord: tuple of floats
        """
        # Take a floating point 3-tuple, find the grids that contain it
        # We do this the stupid way, looking along all three axes
        # Choose works like this:
        #   na.choose(condition, (false_result, true_result))
        #   so here, if gLE > coord, we get a zero, and if it's true, we get
        #   the existing mask.  That way, a single false turns the mask to
        #   zero.
        # We could do this with a 'logical and,' but it's clearer and almost as
        # fast this way.
        mask=na.ones(self.numGrids)
        for i in xrange(len(coord)):
            na.choose(na.greater(self.gridLeftEdge[:,i],coord[i]), (mask,0), mask)
            na.choose(na.greater(self.gridRightEdge[:,i],coord[i]), (0,mask), mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def findRayGrids(self, coord, axis):
        """
        Returns the objects, indices of grids that a ray intersects

        @param coord: the ray endpoint
        @type coord: tuple of floats
        @param axis: the axis the ray travels parallel to
        @type axis: integer
        """
        # Let's figure out which grids are on the slice
        mask=na.ones(self.numGrids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the two edges, we win!
        na.choose(na.greater(self.gridRightEdge[:,x_dict[axis]],coord[0]),(0,mask),mask)
        na.choose(na.greater(self.gridLeftEdge[:,x_dict[axis]],coord[0]),(mask,0),mask)
        na.choose(na.greater(self.gridRightEdge[:,y_dict[axis]],coord[1]),(0,mask),mask)
        na.choose(na.greater(self.gridLeftEdge[:,y_dict[axis]],coord[1]),(mask,0),mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def findSliceGrids(self, coord, axis):
        """
        Returns the objects, indices of grids that a slice intersects

        @param coord: three floats
        @type coord: tuple of floats
        @param axis: the axis the slice is through
        @type axis: integer
        """
        # Let's figure out which grids are on the slice
        mask=na.ones(self.numGrids)
        # So if gRE > coord, we get a mask, if not, we get a zero
        #    if gLE > coord, we get a zero, if not, mask
        # Thus, if the coordinate is between the edges, we win!
        na.choose(na.greater(self.gridRightEdge[:,axis],coord),(0,mask),mask)
        na.choose(na.greater(self.gridLeftEdge[:,axis],coord),(mask,0),mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

    def getSlice(self, center, axis, field, fileName=None, outline=False):
        """
        Returns array of vals.
        @deprecated: We now use EnzoSlice
        """
        # We take a 3-tuple of the coordinate we want to slice through, as well
        # as the axis we're slicing along
        rvs=[]
        g,ind = self.findSliceGrids(center[axis],axis)
        time1=time.time()
        for grid in g:
            mylog.info("Getting from grid %s", grid.id)
            rvs.append(grid.getSlice(center[axis],axis,field,outline))
        allPoints = na.concatenate(rvs)
        if fileName:
            time2=time.time()
            mylog.info("It took %s seconds to generate a slice through all levels", time2-time1)
            mylog.info("This means %s points in %s grids!", allPoints.shape[0], len(g))
            f=open(fileName, "w")
            f.write("x\ty\tz\tdx\tdy\n")
            for i in xrange(allPoints.shape[0]):
                f.write("%0.20f %0.20f %0.5e %0.20f %0.20f\n" % \
                        (allPoints[i,0], \
                        allPoints[i,1], \
                        allPoints[i,2], \
                        allPoints[i,3], \
                        allPoints[i,4] ) )
            f.close()
        else:
            return allPoints

    def findSphereGrids(self, center, radius):
        """
        Returns objects, indices of grids within a sphere

        @param center: coordinate of center
        @type center: tuple of floats
        @param radius: the radius of the sphere in code units!
        """
        centers = (self.gridRightEdge + self.gridLeftEdge)/2.0
        long_axis = na.maximum.reduce(self.gridRightEdge - self.gridLeftEdge, 1)
        t = centers - center
        dist = na.sqrt(t[:,0]**2+t[:,1]**2+t[:,2]**2)
        gridI = na.where(na.logical_and((self.gridDxs<=radius)[:,0],(dist < (radius + long_axis))) == 1)
        return self.grids[gridI], gridI

    def getSphere(self, center, radius, fields):
        """
        Deprecated.  Returns array of vals.  Don't use.

        @deprecated: Use EnzoSphere
        """
        # We first figure out which grids are within distance radius of center
        # If the center of the box is within the distance R+(long axis of box)
        # then we will examine them.
        # Additionally, don't consider boxes whose dx is greater than the
        # radius
        time1 = time.time()
        if not isinstance(fields, types.ListType):
            fields = [fields]
        if self.center == None:
            self.center = center
        g, ind = self.findSphereGrids(center, radius)
        xs = []
        ys = []
        zs = []
        values = []
        i = 0
        for grid in g:
            i+=1
            # Get the points now
            x,y,z, v = grid.getSphere(center, radius, fields)
            mylog.debug("Took %s / %s points from %s at level %s ( %s / %s )",  \
                len(x), na.product(grid.ActiveDimensions), grid.id, grid.Level, i, ind[0].shape[0])
            xs.append(x)
            ys.append(y)
            zs.append(z)
            values.append(v)
            grid.clearAll()
        xs = na.concatenate(xs)
        ys = na.concatenate(ys)
        zs = na.concatenate(zs)
        values = na.concatenate(values)
        time2 = time.time()
        mylog.info("Total of %s points in %0.3e seconds!", xs.shape[0], time2-time1)
        return [xs, ys, zs, values]

    def findMax(self, field, finestLevels = 1):
        """
        Returns value, center of location of maximum for a given field

        Arguments:
        @param field: field (derived or otherwise) of which to look for maximum
        @keyword finestLevels: whether or not to search NUMTOCHECK finest levels
        @type finestLevels: boolean
        """
        if finestLevels:
            gI = na.where(self.gridLevels >= self.maxLevel - NUMTOCHECK)
        else:
            gI = na.where(self.gridLevels >= 0) # Slow but pedantic
        maxVal = -1e100
        for grid in self.grids[gI[0]]:
            mylog.debug("Checking %s (level %s)", grid.id, grid.Level)
            val, coord = grid.findMax(field)
            if val > maxVal:
                maxCoord = coord
                maxVal = val
                maxGrid = grid
        mc = na.array(maxCoord)
        pos=maxGrid.getPosition(mc)
        pos[0] += 0.5*maxGrid.dx
        pos[1] += 0.5*maxGrid.dx
        pos[2] += 0.5*maxGrid.dx
        mylog.info("Max Value is %0.5e at %0.16f %0.16f %0.16f in grid %s at level %s", \
              maxVal, pos[0], pos[1], pos[2], maxGrid, maxGrid.Level)
        self.center = pos
        # This probably won't work for anyone else
        self.bulkVelocity = (maxGrid["x-velocity"][maxCoord], \
                             maxGrid["y-velocity"][maxCoord], \
                             maxGrid["z-velocity"][maxCoord])
        self.parameters["Max%sValue" % (field)] = maxVal
        self.parameters["Max%sPos" % (field)] = "%s" % (pos)
        return maxVal, pos

    def findMin(self, field):
        """
        Returns value, center of location of minimum for a given field

        Arguments:
        @param field: field (derived or otherwise) of which to look for maximum
        """
        gI = na.where(self.gridLevels >= 0) # Slow but pedantic
        minVal = 1e100
        for grid in self.grids[gI[0]]:
            mylog.debug("Checking %s (level %s)", grid.id, grid.Level)
            val, coord = grid.findMin(field)
            if val < minVal:
                minCoord = coord
                minVal = val
                minGrid = grid
        mc = na.array(minCoord)
        pos=minGrid.getPosition(mc)
        pos[0] += 0.5*minGrid.dx
        pos[1] += 0.5*minGrid.dx
        pos[2] += 0.5*minGrid.dx
        mylog.info("Min Value is %0.5e at %0.16f %0.16f %0.16f in grid %s at level %s", \
              minVal, pos[0], pos[1], pos[2], minGrid, minGrid.Level)
        self.center = pos
        # This probably won't work for anyone else
        self.binkVelocity = (minGrid["x-velocity"][minCoord], \
                             minGrid["y-velocity"][minCoord], \
                             minGrid["z-velocity"][minCoord])
        self.parameters["Min%sValue" % (field)] = minVal
        self.parameters["Min%sPos" % (field)] = "%s" % (pos)
        return minVal, pos

    def getProjection(self, axis, field, fileName=None, minLevel=0, maxLevel=None, \
                      weightField=None):
        """
        Returns projection.  Currently in flux.  Don't use directly.
        """
        # Currently weightField does nothing.
        if maxLevel == None:
            maxLevel = self.maxLevel
        # First we precalculate how much memory we will need
        totalProj = 0
        memoryPerLevel = {}
        #gridDataIndices = na.zeros((self.numGrids,4))
        i = 0
        for level in xrange(self.maxLevel+1):
            memoryPerLevel[level] = 0
            mylog.info( "Working on level %s", level)
            grids = self.levelIndices[level]
            numGrids = len(grids)
            RE = self.gridRightEdge[grids].copy()
            LE = self.gridLeftEdge[grids].copy()
            for grid in self.grids[grids]:
                if (i%1e3) == 0:
                    mylog.info("Reading and masking %s / %s", i, self.numGrids)
                for ax in xrange(3):
                    grid.generateOverlapMasks(ax, LE, RE)
                    grid.myOverlapGrids[ax] = self.grids[grids[na.where(grid.myOverlapMasks[ax] == 1)]]
                i += 1
        for grid in self.grids:
            myNeeds = grid.ActiveDimensions[(axis+1)%3]*grid.ActiveDimensions[(axis+2)%3]
            totalProj += myNeeds
            memoryPerLevel[grid.Level] += myNeeds
        for level in xrange(maxLevel+1):
            gI = na.where(self.gridLevels==level)
            mylog.info("%s cells and %s grids for level %s", \
             memoryPerLevel[level], len(gI[0]), level)
        mylog.debug("We need %s cells total", totalProj)
        # We start at the coarsest resolution levels
        i = 0
        dataByLevel = {}
        time1 = time.time()
        #aRStart = 0
        totalGridsProjected = 0
        zeroOut = True
        self.dbl_coarse = {}
        for level in xrange(minLevel,maxLevel+1):
            if level == maxLevel:
                zeroOut = False
            time3 = time.time()
            #levelData = {}
            mylog.info("Projecting through level = %s", level)
            #levelData = [na.zeros(memoryPerLevel[level], nT.Int64), \
                            #na.zeros(memoryPerLevel[level], nT.Int64), \
                            #na.zeros(memoryPerLevel[level], nT.Float64), \
                            #na.zeros(memoryPerLevel[level], nT.Float64) ]    # vals
            tempLevelData = [na.zeros(memoryPerLevel[level], nT.Int64), \
                            na.zeros(memoryPerLevel[level], nT.Int64), \
                            na.zeros(memoryPerLevel[level], nT.Float64), \
                            na.zeros(memoryPerLevel[level], nT.Float64) ]    # vals
            myLevelInd=na.where(self.gridLevels == level)[0]
            #print myLevelInd
            gridsToProject = self.grids[myLevelInd]
            ng = len(gridsToProject)
            index = 0
            i=0
            index=0
            time5=time.time()
            #print "Allocating %s for level %s" % (tempLevelData[0].shape, level)
            for grid in gridsToProject:
                i+=1
                grid.retVal=grid.getProjection(axis, field, zeroOut, weight=weightField)
            time6=time.time()
            totalGridsProjected += i
            mylog.info("Grid projecting done in %s seconds (%s / %s total) with %s points", \
                    time6-time5, totalGridsProjected, self.numGrids, index)
            time5=time.time()
            #print "\tCombining with a maximum of %s operations" % (index*index)
            mylog.info("Combining ...")
            i=0
            for grid1 in gridsToProject:
                i += 1
                #if (level > minLevel) and (level <= maxLevel):
                    #print "\tCombining grid (%s / %s) (%s, %s, %s)" % \
                        #(i,ng, len(grid1.myOverlapGrids[axis]), \
                        #len(grid1.Parent.myOverlapGrids[axis]), grid1.retVal[0].shape[0])
                if grid1.retVal[0].shape[0] == 0:
                    #print "\tSkipping grid (%s / %s) (%s, %s)" %\
                        #(i,ng, len(grid1.myOverlapGrids[axis]), grid1.retVal[0].shape[0])
                    continue
                #print grid1.myOverlapGrids[0]
                for grid2 in grid1.myOverlapGrids[axis]:
                    if grid2.retVal[0].shape[0] == 0:
                        continue
                    if grid1.id == grid2.id:
                        continue
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
                if (level > minLevel) and (level <= maxLevel):
                    for grid2 in grid1.Parent.myOverlapGrids[axis]:
                        if grid2.coarseData[0].shape[0] == 0:
                            continue
                        numRefined += EnzoCombine.RefineCoarseData( \
                            grid1.retVal[0], grid1.retVal[1], grid1.retVal[2], \
                            grid1.retVal[4], \
                            grid2.coarseData[0], grid2.coarseData[1], grid2.coarseData[2], \
                            grid2.coarseData[4], 2)
            all_data = [[],[],[],[],[]]
            for grid in gridsToProject:
                #print grid.retVal[0]
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
            # Now we concatenate our lists into an array
            #print all_data[0]
            levelData = []
            levelData.append(na.concatenate(all_data[0]))
            levelData.append(na.concatenate(all_data[1]))
            levelData.append(na.concatenate(all_data[2]))
            levelData.append(na.concatenate(all_data[3]))
            levelData.append(na.concatenate(all_data[4]))
            time6=time.time()
            mylog.info("Took %s seconds with a final %s points", time6-time5, levelData[0].shape[0])
            dx=gridsToProject[0].dx
            #print "\tdx = %s" % (dx)
            # Make new memory-aligned arrays from the old, and disregarding
            # unneeded points
            # Now we use the coarser data to add to our finer data
            # Here we also divide the weighting factor out (if used)
            time5 = time.time()
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
            mylog.info("Level %s done in %s seconds: %s final of %s (%s)", \
                       level, time4-time3, dataByLevel[level][0].shape[0], \
                       levelData[0].shape[0], exampleValue)
            del levelData
        time2 = time.time()
        mylog.info("Got all the projected points in %s seconds", time2-time1)
        self.dataByLevel = dataByLevel
        if fileName == None:
            return dataByLevel
        else:
            outputProjectionASCII(dataByLevel, fileName, minLevel, maxLevel)

    def exportParticlesPB(self, filename, filter = 1, indexboundary = 0, fields = None, scale=1.0):
        """
        Exports all the star particles, or a subset, to a pb file for viewing in
        partiview

        @param filename: filename of the .pb file to create
        @type filename: string
        @keyword filter: the particle type you want to get (assumes 2)
        @type filter: integer
        @keyword fields: the fields you want to snag.  If not supplied, it just
                      grabs the position and index.
        @keyword indexboundary: for those who want to discriminate the
                    particles with particle index
        @type indexboundary: integer
        @keyword scale: the factor to multiply the position by (defaults to 1.0)
        @type scale: float
        """
        import struct
        pbf_magic = 0xffffff98
        header_fmt = 'Iii'
        fmt = 'ifff'
        f = open(filename,"w")
        if fields:
            fmt += len(fields)*'f'
            padded_fields = string.join(fields,"\0") + "\0"
            header_fmt += "%ss" % len(padded_fields)
            args = [pbf_magic, struct.calcsize(header_fmt), len(fields), padded_fields]
            fields = ["particle_index","particle_position_x","particle_position_y","particle_position_z"] \
                   + fields
            format = 'Int32,Float32,Float32,Float32' + ',Float32'*(len(fields)-4)
        else:
            args = [pbf_magic, struct.calcsize(header_fmt), 0]
            fields = ["particle_index","particle_position_x","particle_position_y","particle_position_z"]
            format = 'Int32,Float32,Float32,Float32'
        f.write(struct.pack(header_fmt, *args))
        tot = 0
        sc = na.array([1.0] + [scale] * 3 + [1.0]*(len(fields)-4))
        gI = na.where(self.gridNumberOfParticles.ravel() > 0)
        for g in self.grids[gI]:
            pI = na.where(na.logical_and((g["particle_type"] == filter),(g["particle_index"] >= indexboundary)) == 1)
            tot += pI[0].shape[0]
            toRec = []
            for field, scale in zip(fields, sc):
                toRec.append(scale*g[field][pI])
            particle_info = rec.array(toRec,formats=format)
            particle_info.tofile(f)
        f.close()
        mylog.info("Wrote %s particles to %s", tot, filename)

    def exportBoxesPV(self, filename):
        f=open(filename,"w")
        for l in xrange(self.maxLevel):
            f.write("add object g%s = l%s\n" % (l,l))
            ind = self.selectLevel(l)
            for i in ind:
                f.write("add box -n %s -l %s %s,%s %s,%s %s,%s\n" % \
                    (i+1, self.gridLevels.ravel()[i],
                     self.gridLeftEdge[i,0], self.gridRightEdge[i,0],
                     self.gridLeftEdge[i,1], self.gridRightEdge[i,1],
                     self.gridLeftEdge[i,2], self.gridRightEdge[i,2]))

    def exportBoxes(self):
        pass

    def initializeEnzoInterface(self, idt_val = 0.0):
        """
        Here we start up the SWIG interface, grabbing what we need from it.

        @keyword idt_val: the initialdt fed to ReadParameterFile (doesn't need
                          to be set)
        @type idt_val: float
        """
        ei = yt.enki.EnzoInterface
        f = open(self.parameterFilename, "r")
        self.eiTopGridData = ei.TopGridData()
        idt = ei.new_Float()
        ei.Float_assign(idt,idt_val)
        ei.cvar.debug = 1 # Set debugging to on, for extra output!
                          # Hm, we *should* have some kind of redirection here
        ei.SetDefaultGlobalValues(self.eiTopGridData)
        # Set up an initial dt
        ei.ReadParameterFile(f,self.eiTopGridData,idt)
        ei.InitializeRateData(self.eiTopGridData.Time)
        ei.InitializeRadiationFieldData(self.eiTopGridData.Time)

def outputProjectionASCII(dataByLevel, fileName, minLevel, maxLevel):
        """
        Don't use directly.  Currently in flux.

        @deprecated: The new EnzoProj datatype already deals with outputting,
        thanks to EnzoData
        """
        mylog.debug("Now outputting to %s in five-column ascii.  How efficient!", fileName)
        k=0
        f=open(fileName,"w")
        f.write("x\ty\tz\tdx\tdy\n")
        i=0
        for level in xrange(minLevel, maxLevel+1):
            #f.write("# Level %s\n" % level)
            j=i
            numBad=0
            data = dataByLevel[level]
            dx = data[3]
            x = (data[0])*dx
            y = (data[1])*dx
            vals = data[2]
            for i in xrange(x.shape[0]):
                k+=1
                if vals[i]==0 or vals[i]==-1:
                    numBad+=1
                    print "%0.20e %0.20e %0.20e %0.20e %0.20e" % \
                            (x[i], y[i], vals[i], dx, dx)
                    continue
                f.write("%0.20e %0.20e %0.20e %0.20e %0.20e\n" % \
                        (x[i]+0.5*dx, y[i]+0.5*dx, (vals[i]), dx/2.0, dx/2.0))
        f.close()
        mylog.debug( "Wrote %s lines", k)
            
def splitConvertGridParameter(vals, func, toAdd, curGrid):
    """
    Quick function to split up a parameter and convert it and toss onto a grid
    """
    j = 0
    for v in vals.split():
        toAdd[curGrid-1,j] = func(v)
        j+=1
    #toAdd[curGrid-1,:] = map(func,vals.split())

scanf_regex = {}
scanf_regex['e'] = r"[-+]?\d+\.?\d*?|\.\d+[eE][-+]?\d+?"
scanf_regex['g'] = scanf_regex['e']
scanf_regex['f'] = scanf_regex['e']
scanf_regex['F'] = scanf_regex['e']
#scanf_regex['g'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
#scanf_regex['f'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
#scanf_regex['F'] = r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
scanf_regex['i'] = r"[-+]?(0[xX][\dA-Fa-f]+|0[0-7]*|\d+)"
scanf_regex['d'] = r"[-+]?\d+"
scanf_regex['s'] = r"\S+"

def constructRegularExpressions(param, toReadTypes):
    re_e=r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    re_i=r"[-+]?(0[xX][\dA-Fa-f]+|0[0-7]*|\d+)"
    rs = "^%s\s*=\s*" % (param)
    for t in toReadTypes:
        rs += "(%s)\s*" % (scanf_regex[t])
    rs +="$"
    return re.compile(rs,re.M)
