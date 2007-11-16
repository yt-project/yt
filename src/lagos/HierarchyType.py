"""
Enzo hierarchy container class
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
from yt.funcs import *
from collections import defaultdict
import string, re, gc, time
import cPickle
import yt.enki

dataStyleFuncs = \
   { 4: (readDataHDF4, readAllDataHDF4, getFieldsHDF4, readDataSliceHDF4), \
     5: (readDataHDF5, readAllDataHDF5, getFieldsHDF5, readDataSliceHDF5), \
     6: (readDataPacked, readAllDataPacked, getFieldsPacked, readDataSlicePacked) \
   }

class EnzoHierarchy:
    """
    Class for handling Enzo timestep outputs

    @param pf: The OutputFile we're instantiating from
    @type pf: L{EnzoOutput}
    @keyword data_style: The type of Enzo Output we're going to read from -- 
                         4 : hdf4, 5 : hdf5, 6 : packed HDF5
    @type data_style: int
    """
    eiTopGrid = None
    @time_execution
    def __init__(self, pf, data_style=4):
        # For now, we default to HDF4, but allow specifying HDF5
        # Expect filename to be the name of the parameter file, not the
        # hierarchy
        self.hierarchyFilename = os.path.abspath(pf.parameterFilename) \
                               + ".hierarchy"
        self.boundaryFilename = os.path.abspath(pf.parameterFilename) \
                               + ".boundary"
        self.directory = os.path.dirname(self.hierarchyFilename)
        self.parameterFile = pf
        self.dataFile = None
        # Now we search backwards from the end of the file to find out how many
        # grids we have, which allows us to preallocate memory
        self.hierarchyLines = open(self.hierarchyFilename).readlines()
        self.hierarchyString = open(self.hierarchyFilename).read()
        for i in xrange(len(self.hierarchyLines)-1,0,-1):
            line = self.hierarchyLines[i]
            if line.startswith("Grid ="):
                self.numGrids = int(line.split("=")[-1])
                break
        self.figureOutDataType()
        # For some reason, r8 seems to want Float64
        if self.parameters.has_key("CompilerPrecision") \
            and self.parameters["CompilerPrecision"] == "r4":
            EnzoFloatType = nT.Float32
        else:
            EnzoFloatType = nT.Float64

        self.setupClasses()
        self.gridDimensions = na.zeros((self.numGrids,3), nT.Int32)
        self.gridStartIndices = na.zeros((self.numGrids,3), nT.Int32)
        self.gridEndIndices = na.zeros((self.numGrids,3), nT.Int32)
        self.gridLeftEdge = na.zeros((self.numGrids,3), EnzoFloatType)
        self.gridRightEdge = na.zeros((self.numGrids,3), EnzoFloatType)
        self.gridLevels = na.zeros((self.numGrids,1), nT.Int32)
        self.gridDxs = na.zeros((self.numGrids,1), EnzoFloatType)
        self.gridTimes = na.zeros((self.numGrids,1), nT.Float64)
        self.gridNumberOfParticles = na.zeros((self.numGrids,1))

        self.grids = obj.array([self.grid(i+1) for i in xrange(self.numGrids)])
        self.gridReverseTree = [-1] * self.numGrids
        self.gridTree = [ [] for i in range(self.numGrids)]

        # Now some statistics:
        #   0 = number of grids
        #   1 = number of cells
        #   2 = blank
        desc = {'names': ['numgrids','numcells','level'],
                'formats':['Int32']*3}
        self.levelStats = blankRecordArray(desc, MAXLEVEL)
        self.levelStats['level'] = [i for i in range(MAXLEVEL)]

        # For use with derived quantities depending on centers
        # Although really, I think perhaps we should take a closer look
        # at how "center" is used.
        self.center = None
        self.bulkVelocity = None

        self.initializeDataFile()
        self.populateHierarchy()
    def figureOutDataType(self):
        for i in xrange(len(self.hierarchyLines)-1,0,-1):
            line = self.hierarchyLines[i]
            if line.startswith("BaryonFileName") or \
               line.startswith("FileName "):
                testGrid = line.split("=")[-1].strip().rstrip()
                break
        if testGrid[0] != os.path.sep:
            testGrid = os.path.join(self.directory, testGrid)
        try:
            a = SD.SD(testGrid)
            self.data_style = 4
            mylog.debug("Detected HDF4")
        except:
            a = tables.openFile(testGrid)
            for b in a.iterNodes("/"):
                c = "%s" % (b)
                break
            if c.startswith("/Grid"):
                mylog.debug("Detected packed HDF5")
                self.data_style = 6
            else:
                mylog.debug("Detected unpacked HDF5")
                self.data_style = 5
            a.close()

    def setupClasses(self):
        """
        This is our class factory.  It takes the base classes and assigns to
        them appropriate data-reading functions based on the data-style.
    
        @postcondition: .grid, .prof, .slice, .region, .datacube and .sphere
                        will be classes, instantiated with the appropriate 
                        methods of obtaining data.
        """
        dd = { 'readDataFast' : dataStyleFuncs[self.data_style][0],
               'readAllData' : dataStyleFuncs[self.data_style][1],
               'getFields' : dataStyleFuncs[self.data_style][2],
               'readDataSlice' : dataStyleFuncs[self.data_style][3],
               'pf' : self.parameterFile,
               'hierarchy': self }
        self.grid = classobj("EnzoGrid",(EnzoGridBase,), dd)
        self.proj = classobj("EnzoProj",(EnzoProjBase,), dd)
        self.slice = classobj("EnzoSlice",(EnzoSliceBase,), dd)
        self.region = classobj("EnzoRegion",(EnzoRegionBase,), dd)
        self.datacube = classobj("EnzoDataCube",(EnzoDataCubeBase,), dd)
        self.sphere = classobj("EnzoSphere",(EnzoSphereBase,), dd)

    def initializeDataFile(self):
        """
        We initialize our data-serialization file here.

        @precond: tables must be imported and we must have write access to the
                  directory the data is contained in.  (Otherwise silent failure.)
        """
        if not ytcfg.getboolean('lagos','serialize'): return
        fn = os.path.join(self.directory,"%s.yt" % self["CurrentTimeIdentifier"])
        try:
            self.dataFile = tables.openFile(fn, "a")
        except:
            pass

    def saveData(self, array, node, name):
        """
        Arbitrary numpy data will be saved to the region in the datafile
        described by node and name.
        @arg array: The data to be saved.
        @type array: NumPy array
        @arg node: The HDF5 node to save to
        @type node: String
        @arg name: Name of the array in the file
        @type name: String
        """
        if self.dataFile != None:
            self.dataFile.createArray(node, name, array, createparents=True)
            self.dataFile.flush()

    def getData(self, node, name):
        if self.dataFile == None:
            return None
        try:
            return self.dataFile.getNode(node, name)
        except tables.exceptions.NoSuchNodeError:
            return None

    def __xattrs__(self, mode="default"):
        return ("basename", "getTimeID()","maxLevel")

    def __del__(self):
        """
        Let's see if we can delete some stuff here!
        """
        if self.dataFile:
            self.dataFile.close()
            del self.dataFile
        try:
            del self.eiTopGrid
        except:
            pass
        del self.gridReverseTree
        del self.gridLeftEdge, self.gridRightEdge
        del self.gridLevels, self.gridStartIndices, self.gridEndIndices
        del self.gridTimes, self.hierarchyString, self.hierarchyLines
        for gridI in xrange(self.numGrids):
            for g in self.gridTree[gridI]:
                del g
        del self.gridTree

    @time_execution
    def populateHierarchy(self):
        """
        Instantiates all of the grid objects, with their appropriate
        parameters.  This is the work-horse.
        """
        harray = self.getData("/", "Hierarchy")
        if harray:
            self.gridDimensions[:] = harray[:,0:3]
            self.gridStartIndices[:] = harray[:,3:6]
            self.gridEndIndices[:] = harray[:,6:9]
            self.gridLeftEdge[:] = harray[:,9:12]
            self.gridRightEdge[:] = harray[:,12:15]
            self.gridLevels[:] = harray[:,15:16]
            self.gridTimes[:] = harray[:,16:17]
            self.gridNumberOfParticles[:] = harray[:,17:18]
            del harray
            # Now get the baryon filenames
            re_BaryonFileName = constructRegularExpressions("BaryonFileName",('s'))
            t = re.findall(re_BaryonFileName, self.hierarchyString)
            for fnI in xrange(len(t)):
                self.grids[fnI].setFilename(t[fnI])
            re_BaryonFileName = constructRegularExpressions("FileName",('s'))
            t = re.findall(re_BaryonFileName, self.hierarchyString)
            for fnI in xrange(len(t)):
                self.grids[fnI].setFilename(t[fnI])
        else:
            def splitConvertGridParameter(vals, func, toAdd, curGrid):
                """
                Quick function to split up a parameter and convert it and toss onto a grid
                """
                j = 0
                for v in vals.split():
                    toAdd[curGrid-1,j] = func(v)
                    j+=1
            for line in self.hierarchyLines:
                # We can do this the slow, 'reliable' way by stripping
                # or we can manually pad all our strings, which speeds it up by a
                # factor of about ten
                #param, vals = map(strip,line.split("="))
                if len(line) < 2:
                    continue
                param, vals = line.split("=")
                param = param.rstrip() # This slows things down considerably...
                                       # or so I used to think...
                if param == "Grid":
                    curGrid = int(vals)
                    self.grids[curGrid-1] = self.grid(curGrid)
                elif param == "GridDimension":
                    splitConvertGridParameter(vals, float, self.gridDimensions, curGrid)
                elif param == "GridStartIndex":
                    splitConvertGridParameter(vals, int, self.gridStartIndices, curGrid)
                elif param == "GridEndIndex":
                    splitConvertGridParameter(vals, int, self.gridEndIndices, curGrid)
                elif param == "GridLeftEdge":
                    splitConvertGridParameter(vals, float, self.gridLeftEdge, curGrid)
                elif param == "GridRightEdge":
                    splitConvertGridParameter(vals, float, self.gridRightEdge, curGrid)
                elif param == "Level":
                    splitConvertGridParameter(vals, int, self.gridLevels, curGrid)
                elif param == "Time":
                    splitConvertGridParameter(vals, float, self.gridTimes, curGrid)
                elif param == "NumberOfParticles":
                    splitConvertGridParameter(vals, float, self.gridNumberOfParticles, curGrid)
                elif param == "FileName":
                    self.grids[curGrid-1].setFilename(vals[1:-1])
                elif param == "BaryonFileName":
                    self.grids[curGrid-1].setFilename(vals[1:-1])
            mylog.info("Caching hierarchy information")
            allArrays = na.zeros((self.numGrids,18),nT.Float64)
            allArrays[:,0:3] = self.gridDimensions[:]
            allArrays[:,3:6] = self.gridStartIndices[:]
            allArrays[:,6:9] = self.gridEndIndices[:]
            allArrays[:,9:12] = self.gridLeftEdge[:]
            allArrays[:,12:15] = self.gridRightEdge[:]
            allArrays[:,15:16] = self.gridLevels[:]
            allArrays[:,16:17] = self.gridTimes[:]
            allArrays[:,17:18] = self.gridNumberOfParticles[:]
            self.saveData(allArrays, "/","Hierarchy")
            del allArrays
        treeArray = self.getData("/", "Tree")
        if treeArray == None:
            self.grids[0].Level = 0  # Bootstrap
            self.gridLevels[0] = 0   # Bootstrap
            p = re.compile(r"Pointer: Grid\[(\d*)\]->NextGrid(Next|This)Level = (\d*)$", re.M)
            # Now we assemble the grid tree
            # This is where all the time is spent.
            for m in p.finditer(self.hierarchyString):
                secondGrid = int(m.group(3))-1 # zero-index versus one-index
                if secondGrid == -1:
                    continue
                firstGrid = int(m.group(1))-1
                if m.group(2) == "Next":
                    self.gridTree[firstGrid].append(self.grids[secondGrid])
                    self.gridReverseTree[secondGrid] = firstGrid + 1
                    self.grids[secondGrid].Level = self.grids[firstGrid].Level + 1
                    self.gridLevels[secondGrid] = self.gridLevels[firstGrid] + 1
                elif m.group(2) == "This":
                    parent = self.gridReverseTree[firstGrid]
                    if parent:
                        self.gridTree[parent-1].append(self.grids[secondGrid])
                    self.gridReverseTree[secondGrid] = parent
                    self.grids[secondGrid].Level = self.grids[firstGrid].Level
                    self.gridLevels[secondGrid] = self.gridLevels[firstGrid]
            pTree = [ [ grid.id - 1 for grid in self.gridTree[i] ] for i in range(self.numGrids) ]
            self.gridReverseTree[0] = -1
            self.saveData(cPickle.dumps(pTree), "/", "Tree")
            self.saveData(na.array(self.gridReverseTree), "/", "ReverseTree")
            self.saveData(self.gridLevels, "/", "Levels")
        else:
            mylog.debug("Grabbing serialized")
            pTree = cPickle.loads(treeArray.read())
            self.gridReverseTree = list(self.getData("/","ReverseTree"))
            self.gridTree = [ [ self.grids[i] for i in pTree[j] ] for j in range(self.numGrids) ]
            self.gridLevels = self.getData("/","Levels")[:]
            mylog.debug("Grabbed")
        for i,v in enumerate(self.gridReverseTree):
            if v == -1: self.gridReverseTree[i] = None
        #self.gridReverseTree[0] = None
        self.maxLevel = self.gridLevels.max()
        # Now we do things that we need all the grids to do
        self.fieldList = self.grids[0].getFields()
        # The rest of this can probably be done with list comprehensions, but
        # I think this way is clearer.
        mylog.debug("Preparing grids")
        for i, grid in enumerate(self.grids):
            tlevel = self.gridLevels[i]
            grid.prepareGrid()
            self.levelStats['numgrids'][tlevel] += 1
            self.levelStats['numcells'][tlevel] += na.product(grid.ActiveDimensions)
        mylog.debug("Prepared")
        self.levelIndices = {}
        self.levelNum = {}
        for level in xrange(self.maxLevel+1):
            self.levelIndices[level] = self.selectLevel(level)
            self.levelNum[level] = len(self.levelIndices[level])

    def selectLevel(self, level):
        """
        Returns a list of indices of EnzoHierarchy.grids at the specified level

        @param level: the level
        @type level: integer
        @note: This would be more intuitive if it returned the *actual grids*.
        """
        # We return a numarray of the indices of all the grids on a given level
        indices = na.where(self.gridLevels[:,0] == level)[0]
        return indices

    def getSmallestDx(self):
        """
        Returns (in code units) the smallest dx in the simulation.
        """
        return self.gridDxs.min()

    def printStats(self):
        """
        Prints out relevant information about the simulation
        """
        for i in xrange(MAXLEVEL):
            if (self.levelStats['numgrids'][i]) == 0:
                break
            print "% 3i\t% 6i\t% 11i" % \
                  (i, self.levelStats['numgrids'][i], self.levelStats['numcells'][i])
            dx = self.gridDxs[self.levelIndices[i][0]]
        print "-" * 28
        print "   \t% 6i\t% 11i" % (self.levelStats['numgrids'].sum(), self.levelStats['numcells'].sum())
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
        for item in self.parameterFile.units.items():
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
        #ind = na.where( na.logical_and(self.gridRightEdge[:,axis] > coord, \
                                       #self.gridLeftEdge[:,axis] < coord))
        na.choose(na.greater(self.gridRightEdge[:,axis],coord),(0,mask),mask)
        na.choose(na.greater(self.gridLeftEdge[:,axis],coord),(mask,0),mask)
        ind = na.where(mask == 1)
        return self.grids[ind], ind

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

    def getBoxGrids(self, leftEdge, rightEdge):
        """
        Gets back all the grids between a left edge and right edge

        @param leftEdge: the left edge
        @type leftEdge: array
        @param rightEdge: the right edge
        @type rightEdge: array
        """
        gridI = na.where(na.logical_and(na.any(self.gridRightEdge > leftEdge, axis=1),
                                        na.any(self.gridLeftEdge < rightEdge, axis=1)) == True)
        return self.grids[gridI], gridI

    @time_execution
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

    @time_execution
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

    @time_execution
    def exportParticlesPB(self, filename, filter = 1, indexboundary = 0, fields = None, scale=1.0):
        """
        Exports all the star particles, or a subset, to a pb file for viewing in
        partiview

        @param filename: filename of the .pb file to create
        @type filename: string
        @keyword filter: the particle type you want to get (assumes 1)
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

    @time_execution
    def exportBoxesPV(self, filename):
        """
        Exports the grid structure in partiview text format.

        @arg filename: File to export to.
        @type filename: String
        """
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

    @time_execution
    def exportAmira(self, basename, fields, a5basename, timestep):
        """
        Exports the grid structure in partiview text format.
        """
        if (not iterable(fields)) or (isinstance(fields, types.StringType)):
            fields = [fields]
        for field in fields:
            tt=tables.openFile(basename % {'field':field},"w")
            k=tt.createGroup("/","Parameters and Global Attributes")
            k._f_setAttr("staggering",1)
            tt.close()
            a5=tables.openFile(a5basename % {'field':field},"a")
            a5.createGroup("/", "time-%i" % timestep)
            node = a5.getNode("/","time-%i" % timestep)
            node._f_setAttr("numLevels",self.maxLevel+1)
            node._f_setAttr("time",self["InitialTime"])
            a5.close()
        for level in range(self.maxLevel+1):
            mylog.info("Exporting level %s", level)
            for field in fields:
                a5=tables.openFile(a5basename % {'field':field},"a")
                a5.createGroup("/time-%i" % (timestep),"level-%i" % (level))
                node=a5.getNode("/time-%i" % (timestep),"level-%i" % (level))
                delta = na.array([self.gridDxs[self.levelIndices[level][0]]]*3,dtype=nT.Float64)
                node._f_setAttr("delta",delta)
                node._f_setAttr("numGrids",self.levelNum[level])
                # This next one is not necessarily true.  But, it is for
                # everyone I care about right now...
                node._f_setAttr("relativeRefinementFactor",na.array([2,2,2],dtype=nT.Int32))
                a5.close()
            gid = 0
            for grid in self.grids[self.levelIndices[level]]:
                grid.exportAmira(basename, fields, timestep, a5basename, gid)
                gid += 1

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


    """@todo: Fix these!"""
    def _get_parameters(self):
        return self.parameterFile.parameters
    parameters=property(_get_parameters)

    def __getitem__(self, item):
        return self.parameterFile[item]

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
