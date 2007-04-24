"""
Python-based grid handler, not to be confused with the SWIG-handler

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.lagos import *
import yt.enki, gc

class EnzoGrid:
    """
    Class representing a single Enzo Grid instance
    """
    def __init__(self, hierarchy, id, filename=None):
        """
        Returns an instance of EnzoGrid

        @param hierarchy: EnzoHierarchy, parent hierarchy
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        @param id: grid ID (NOT index, which is ID-1)
        @type id: int
        @keyword filename: filename holding grid data
        @type filename: string
        """
        self.id = id
        self.hierarchy = weakref.proxy(hierarchy)
        self.data = {}
        self.datasets = {}
        self.SDi = None
        self.SDi_datasets = None
        if filename:
            self.setFilename(filename)
        self.myChildMask = None
        self.myChildIndices = None
        self.myOverlapMasks = [None, None, None]
        self.myOverlapGrids = [None, None, None]

    def __getattr__(self, attr):
        """
        Attempts to grab an attribute from a grid, and failing that, grabs it
        from the hierarchy
        (Maybe dangerous?)
        """
        try:
            return eval("self.hierarchy.%s" % (attr))
        except:
            raise AttributeError, attr

    def __getitem__(self, key):
        """
        Returns a field or set of fields for a key or set of keys
        """
        if isinstance(key, types.StringType):
            if self.data.has_key(key):
                return self.data[key]
            else:
                self.readDataFast(key)
                return self.data[key]
        elif isinstance(key, types.ListType) or \
             isinstance(key, types.TupleType):
            tr = []
            for k in key:
                if self.data.has_key(k):
                    tr.append(self.data[k])
                else:
                    self.readDataFast(k)
                    tr.append(self.data[k])
            return tr
        else:
            return self.data[key]

    def __setitem__(self, key, data):
        """
        Sets a data field equal to some value or set of values
        """
        if isinstance(key, types.StringType):
            self.data[key] = data
        elif isinstance(key, types.ListType) or \
             isinstance(key, types.TupleType):
            tr = []
            for kI in range(len(key)):
                self.data[kI] = data[kI]
        else:
            self.data[key] = data

    def has_key(self, key):
        """
        Checks to see if this field *is already generated.*  Will not check to
        see if it *can*  be generated.
        """
        return self.data.has_key(key)

    def keys(self):
        """
        Returns all existing fields.
        """
        return self.data.keys()

    def clearAllGridReferences(self):
        #print "clearingAllGridReferences in EnzoGrid ", self.id
        #return
        self.clearDerivedQuantities()
        if hasattr(self, 'hierarchy'):
            del self.hierarchy
        if hasattr(self, 'Parent'):
            if self.Parent != None:
                self.Parent.clearAllGridReferences()
            del self.Parent
        if hasattr(self, 'Children'):
            for i in self.Children:
                if i != None:
                    #i.clearAllGridReferences()
                    del i
            del self.Children

    def __del__(self):
        #self.clearAllGridReferences()
        return
        #print weakref.getweakrefcount(self)
        #print gc.get_referents(self)
        #del self.Parent
        #del self.hierarchy
        #del self.Children
        #print "Really deleting ", self.id
        

    def prepareGrid(self):
        """
        Copies all the appropriate attributes from the hierarchy
        """
        # Now we give it pointers to all of its attributes
        h = self.hierarchy # cache it
        self.Dimensions = h.gridDimensions[self.id-1]
        #print gc.get_referrers(self)
        self.StartIndices = h.gridStartIndices[self.id-1]
        self.EndIndices = h.gridEndIndices[self.id-1]
        self.LeftEdge = h.gridLeftEdge[self.id-1]
        self.RightEdge = h.gridRightEdge[self.id-1]
        self.Level = h.gridLevels[self.id-1,0]
        self.Time = h.gridTimes[self.id-1,0]
        self.NumberOfParticles = h.gridNumberOfParticles[self.id-1,0]
        self.ActiveDimensions = (self.EndIndices - self.StartIndices + 1)
        #self.ActiveDimensions = self.Dimensions - 3
        self.Children = h.gridTree[self.id-1]
        pID = h.gridReverseTree[self.id-1]
        if pID != None:
            self.Parent = weakref.proxy(h.grids[pID - 1])
        else:
            self.Parent = None
        # So first we figure out what the index is.  We assume
        # that dx=dy=dz
        self.dx = (self.RightEdge[0] - self.LeftEdge[0]) / \
                  (self.EndIndices[0]-self.StartIndices[0]+1)
        self.dy = (self.RightEdge[1] - self.LeftEdge[1]) / \
                  (self.EndIndices[1]-self.StartIndices[1]+1)
        self.dz = (self.RightEdge[2] - self.LeftEdge[2]) / \
                  (self.EndIndices[2]-self.StartIndices[2]+1)
        h.gridDxs[self.id-1,0] = self.dx
        self.coords = None
        #self.generateCoords()
        del h

    def generateChildMask(self):
        """
        Generates self.myChildMask, which is zero where child grids exist (and
        thus, where higher resolution data is available.)
        """
        self.myChildMask = na.ones(self.ActiveDimensions, nT.Int32)
        for child in self.Children:
            # Now let's get our overlap
            si = [None]*3
            ei = [None]*3
            startIndex = (child.LeftEdge - self.LeftEdge)/self.dx
            endIndex = (child.RightEdge - self.LeftEdge)/self.dx
            for i in range(3):
                si[i] = int(startIndex[i])
                ei[i] = int(endIndex[i])
            self.myChildMask[si[0]:ei[0], si[1]:ei[1], si[2]:ei[2]] = 0
        #self.myIndices = na.where(self.myChildMask==1)
        self.myChildIndices = na.where(self.myChildMask==0)

    def generateOverlapMasks(self, axis, LE, RE):
        """
        Generate a mask that shows which cells overlap with other cells on
        different grids.  (If fed appropriate subsets, can be constrained to
        current level.
        Use algorithm described at http://www.gamedev.net/reference/articles/article735.asp

        @param axis: axis  along which line of sight is drawn
        @type axis: int
        @param LE: LeftEdge positions to check against
        @type LE: array of floats
        @param RE: RightEdge positions to check against
        @type RE: array of floats
        """
        x = x_dict[axis]
        y = y_dict[axis]
        cond1 = self.RightEdge[x] > LE[:,x]
        cond2 = self.LeftEdge[x] < RE[:,x]
        cond3 = self.RightEdge[y] > LE[:,y]
        cond4 = self.LeftEdge[y] < RE[:,y]
        self.myOverlapMasks[axis]=na.logical_and(na.logical_and(cond1, cond2), \
                                               na.logical_and(cond3, cond4))
    def __repr__(self):
        return "%s" % (self.id)
    def __int__(self):
        return self.id
    def setFilename(self, filename):
        if filename[0] == os.path.sep:
            self.filename = filename
        else:
            self.filename = self.hierarchy.directory + os.path.sep + filename
        return

    def findMax(self, field):
        """
        Returns value, coordinate of maximum value in this gird

        @param field: field to check
        @type field: string
        """
        coord=nd.maximum_position(self[field])
        val = self[field][coord]
        return val, coord

    def findMin(self, field):
        """
        Returns value, coordinate of minimum value in this gird

        @param field: field to check
        @type field: string
        """
        coord=nd.minimum_position(self[field])
        val = self[field][coord]
        return val, coord

    def getPosition(self, coord):
        """
        Returns position of a coordinate

        @param coord: position to check
        @type coord: array of floats
        """
        # We accept arrays here, people, not tuples
        pos = (coord + 0.0) * self.dx + self.LeftEdge
        # Should 0.0 be 0.5?
        return pos

    def getSlice(self, coord, axis, field, outline=False):
        """
        @deprecated: use EnzoSlice!
        """
        if self.myChildMask == None:
            self.generateChildMask()
        # So what's our index of slicing?  This is what we need to figure out
        # first, so we can deal with our data in the fastest way.
        # NOTE: This should be fixed.  I don't think it works properly or
        # intelligently.
        wantedIndex = int(((coord-self.LeftEdge[axis])/self.dy))
        # I can't think of a better way to do this: because we have lots of
        # different arrays, and because we don't want to mess that up, we
        # shouldn't do any axis-swapping, I think.  So, looks like we're just
        # going to have to do this the stupid way.  (I suspect there's a more
        # clever way, involving generating an array, and using %, but this
        # works.
        xaxis = x_dict[axis]
        yaxis = y_dict[axis]
        if axis == 0:
            cm = na.where(self.myChildMask[wantedIndex,:,:] == 1)
            cmI = na.indices(self.myChildMask[wantedIndex,:,:].shape)
            slicedData = self[field][wantedIndex,:,:]
        elif axis == 1:
            cm = na.where(self.myChildMask[:,wantedIndex,:] == 1)
            cmI = na.indices(self.myChildMask[:,wantedIndex,:].shape)
            slicedData = self[field][:,wantedIndex,:]
        elif axis == 2:
            cm = na.where(self.myChildMask[:,:,wantedIndex] == 1)
            cmI = na.indices(self.myChildMask[:,:,wantedIndex].shape)
            slicedData = self[field][:,:,wantedIndex]
        # So now we figure out which points we want, and their (x,y,z) values
        xind = cmI[0,:]
        xpoints = xind[cm]*self.dx+(self.LeftEdge[xaxis] + 0.5*self.dx)
        yind = cmI[1,:]
        ypoints = yind[cm]*self.dx+(self.LeftEdge[yaxis] + 0.5*self.dx)
        dataVals = slicedData[cm]
        # We now have a couple one dimensional arrays.  We will
        # make these one array, and return them as [x y val dx dy]
        if self.hierarchy.conversionFactors.has_key(field):
            conv = self.hierarchy.conversionFactors[field]
        else:
            conv = 1
        numVals = dataVals.shape[0]
        retVal = na.array(shape=(numVals,5), type=nT.Float64)
        retVal[:,0] = xpoints
        retVal[:,1] = ypoints
        retVal[:,2] = dataVals*conv
        retVal[:,3] = self.dx/2.0
        retVal[:,4] = self.dx/2.0
        if outline == True:
            del self[field]
        return retVal

    def getProjection(self, axis, field, zeroOut, weight=None):
        """
        Projects along an axis.  Currently in flux.  Shouldn't be called
        directly.
        """
        if weight == None:
            maskedData = self[field].copy()
            weightData = na.ones(maskedData.shape)
        else:
            maskedData = self[field] * self[weight]
            weightData = self[weight].copy()
        if self.myChildMask == None:
            self.generateChildMask()
        if len(self.myOverlapMasks) == 0:
            self.generateOverlapMasks()
        if zeroOut:
            maskedData[self.myChildIndices]=0
            weightData[self.myChildIndices]=0
            toCombineMask = na.logical_and.reduce(self.myChildMask, axis)
        # How do we do this the fastest?
        # We only want to project those values that don't have subgrids
        fullProj = na.sum(maskedData,axis)*self.dx # Gives correct shape
        weightProj = na.sum(weightData,axis)*self.dx
        #fullProj = na.maximum.reduce(maskedData,axis) # Gives correct shape
        if not zeroOut:
            toCombineMask = na.ones(fullProj.shape, dtype=nT.Bool)
        cmI = na.indices(fullProj.shape)
        # So now we figure out which points we want, and their (x,y,z) values
        # Note that this is currently wrong for anything other than x (axis = 0)
        xind = cmI[0,:]
        yind = cmI[1,:]
        #print "AUGH!", toCombineMask.shape, fullProj.shape, weightProj.shape, maskedData.shape, self.dx, \
                #weightData.shape
        xpoints = na.array(xind+(self.LeftEdge[x_dict[axis]]/self.dx),nT.Int64)
        ypoints = na.array(yind+(self.LeftEdge[y_dict[axis]]/self.dx),nT.Int64)
        return [xpoints.ravel(), ypoints.ravel(), fullProj.ravel(), toCombineMask.ravel(), weightProj.ravel()]

    def getSliceAll(self, coord, axis, field):
        tempMask = self.myChildMask
        self.myChildMask = na.ones(self.ActiveDimensions)
        points,dataVals=self.getSlice(coord, axis, field)
        self.myChildMask = tempMask
        return points,dataVals

    def clearAll(self):
        """
        Clears all datafields from memory.
        """
        for key in self.keys():
            del self.data[key]
        del self.data
        self.data = {}
        self.clearDerivedQuantities()

    def clearDerivedQuantities(self):
        """
        Clears coordinates, myChildIndices, myChildMask.
        """
        del self.coords
        self.coords = None
        del self.myChildIndices
        self.myChildIndices = None
        del self.myChildMask
        self.myChildMask = None

    def generateField(self, fieldName):
        """
        Generates, or attempts to generate,  a field not found in the file

        See fields.py for more information.  fieldInfo.keys() will list all of
        the available derived fields.  Note that we also make available the
        suffices _Fraction and _Squared here.  All fields prefixed with 'k'
        will force an attempt to use the chemistry tables to generate them from
        temperature.  All fields used in generation will remain resident in
        memory.

        I feel like there's a reason that EnzoGrid isn't subclassed from
        EnzoData, and I think it's related to this method.  But I can't remember now.

        @param fieldName: field name
        @type fieldName: string

        """
        # This is for making derived fields
        # Note that all fields used for derivation are kept resident in memory: probably a 
        # mistake, but it is expensive to do a lookup.  I will fix this later.
        #
        # Note that you can do a couple things: the suffices _Fraction and
        # _Squared will be dealt with appropriately.  Not sure what else to
        # add.
        if fieldName.endswith("Fraction"):
            # Very simple mass fraction here.  Could be modified easily,
            # but that would require a dict lookup, which is expensive, or
            # an elif block, which is inelegant
            baryonField = "%s_Density" % (fieldName[:-9])
            self[fieldName] = self[baryonField] / self["Density"]
        elif fieldName.endswith("Squared"):
            baryonField = fieldName[:-7]
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

    def generateCoords(self):
        """
        Creates self.coords, which is of dimensions (3,ActiveDimensions)
        """
        if self.coords != None:
            return
        ind = na.indices(self.ActiveDimensions)
        LE = na.reshape(self.LeftEdge,(3,1,1,1))
        self.coords = (ind+0.5)*self.dx+LE
    
    def getSphere(self, center, radius, fields, zeroOut = True):
        """
        Returns a sphere from within the grid.  Don't use.
        @deprecated: Use EnzoSphere!
        """
        if self.myChildMask == None or self.myChildIndices == None:
            self.generateChildMask()
        # First we find the cells that are within the sphere
        pointI = na.where(na.logical_and((self["RadiusCode"]<=radius),self.myChildMask==1)==1)
        # Note that we assumed here that all our data will be nT.Float32
        # Not a *terrible* assumption...
        trData = na.zeros((pointI[0].shape[0],len(fields)), nT.Float64)
        i = 0
        for field in fields:
            if self.hierarchy.conversionFactors.has_key(field):
                conv = self.hierarchy.conversionFactors[field]
            else:
                conv = 1
            trData[:,i] = self[field][pointI] * conv
            i+=1
        return [self.coords[0,:][pointI], \
                self.coords[1,:][pointI], \
                self.coords[2,:][pointI], \
                trData]

    def getEnzoGrid(self):
        """
        This attempts to get an instance of this particular grid from the SWIG
        interface.  Note that it first checks to see if the ParameterFile has
        been instantiated.
        """
        if self.hierarchy.eiTopGrid == None:
            self.hierarchy.initializeEnzoInterface()
        p=re.compile("Grid = %s\n" % (self.id))
        h=open(self.hierarchyFilename,"r").read()
        m=re.search(p,h)
        h=open(self.hierarchyFilename,"r")
        retVal = yt.enki.EnzoInterface.fseek(h, long(m.end()), 0)
        self.eiGrid=yt.enki.EnzoInterface.grid()
        cwd = os.getcwd() # Hate doing this, need to for relative pathnames
        os.chdir(self.hierarchy.directory)
        self.eiGrid.ReadGrid(h, 1)
        os.chdir(cwd)
        mylog.debug("Grid read with SWIG")
