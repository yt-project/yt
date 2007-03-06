"""
ProblemType class definition, including utility methods

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""
from yt.enki import *

class ProblemType:
    """
    The base class for new problem types, for initializaing data for use in
    Enzo.  Includes several helper functions that could be useful for new data
    generation.

    Note that we do *not* have helper classes for the SWIG grid.h wrapping.
    This is *important*, as we do not want to treat them too much like python
    objects.  (For now, the fact that data must be flushed back and forth means
    that we have to treat them as conduits, not objects.)
    """
    def __init__(self, MetaData, TopGrid, LevelArray, InitArgs):
        """
        @param MetaData: the MetaData associated with this problem
        @type MetaData: EnzoInterface.MetaData
        @param TopGrid: the TopGrid
        @type TopGrid: EnzoInterface.HierarchyEntry
        @param LevelArray: the array of grids on each level
        @type LevelArray: L{LevelArrayWrapper<yt.enki.LevelArrayWrapper>}
        @param InitArgs: the arguments to the init function
        @type InitArgs: dict
        """
        self.MetaData = MetaData
        self.TopGrid = TopGrid
        self.LevelArray = LevelArray
        self.Defaults = {}
        self.InitArgs = InitArgs
        self.count = 0
        self.GridFieldTypes = []
        self.FieldIndex = {}
        self.Fields = []
        self.oldFields = []

    def __getitem__(self, item):
        # We set precedence here of the InitArgs over the MetaData
        if self.InitArgs.has_key(item):
            return self.InitArgs[item]
        if self.Defaults.has_key(item):
            return self.Defaults[item]
        try:
            return self.MetaData.__getattr__(item)
        except AttributeError:
            exec("tr = EnzoInterface.cvar.%s" % (item))
            return tr

    def Initialize(self):
        """
        This function *must* be overridden.
        """
        mylog.warning("Oops, looks like your problem initializer wasn't defined!")

    def InitializeMetaData(self):
        pass

    def AddField(self, field):
        EnzoInterface.charSArray_setitem(EnzoInterface.cvar.DataLabel, self.count, field)
        self.GridFieldTypes.append(FieldTypes[field])
        self.FieldIndex[field] = self.count
        #print "Appending", field, len(self.GridFieldTypes)
        self.count += 1

    def InitializeFields(self):
        # These are the standard ones
        self.AddField("Density")
        self.AddField("TotalEnergy")
        if self["DualEnergyFormalism"]:
            self.AddField("GasEnergy")
        self.AddField("x-velocity")
        self.AddField("y-velocity")
        self.AddField("z-velocity")
        if self["CollapseTestUseColour"]:
            self.AddField("colour")
        if self["MultiSpecies"]:
            self.AddField("Electron_Density")
            self.AddField("HI_Density")
            self.AddField("HII_Density")
            self.AddField("HeI_Density")
            self.AddField("HeII_Density")
            self.AddField("HeIII_Density")
            if self["MultiSpecies"] > 1:
                self.AddField("HM_Density")
                self.AddField("H2I_Density")
                self.AddField("H2II_Density")
            if self["MultiSpecies"] > 2:
                self.AddField("DI_Density")
                self.AddField("DII_Density")
                self.AddField("HDI_Density")

    def InitializeFieldsInGrid(self, grid):
        # We now allocate memory and toss them to the thinige.
        # Note that we clear the field list  at the start of each call
        #self.oldFields.append(self.Fields) # For reference counting.  Ugly.  Hateful.
        self.Fields = []
        #print self.GridFieldTypes
        #grid.NumberOfBaryonFields = len(self.GridFieldTypes)+1
        #print "NUMBER:",grid.NumberOfBaryonFields, self.GridFieldTypes, len(self.GridFieldTypes)+1
        grid.SetNumberOfBaryonFields(len(self.GridFieldTypes))
        #grid.SetNumberOfBaryonFields(1)
        #print "NUMBER:",grid.NumberOfBaryonFields, self.GridFieldTypes, len(self.GridFieldTypes)
        for fieldI in range(len(self.GridFieldTypes)):
            EnzoInterface.intArray_setitem(grid.FieldType, fieldI, self.GridFieldTypes[fieldI])
            #self.Fields.append(grid.DirectManipulate(fieldI))
            #self.Fields.append(grid.AppendBaryonField(fieldI))
            self.Fields.append(grid.BaryonFieldToNumArray(fieldI))

    def FlushFieldsToGrid(self, grid):
        for fieldI in range(len(self.GridFieldTypes)):
            grid.NumArrayToBaryonField(self.Fields[fieldI], fieldI)

    def RefineUpTo(self, maxLevel):
        for level in range(maxLevel):
            mylog.info("Refining level %s (MaximumLevelOfRefinement = %s)", level, maxLevel)
            retVal = EnzoInterface.RebuildHierarchy(self.MetaData, self.LevelArray.la, level)
            mylog.info("Done refining level %s", level)
            if retVal != ENZO_SUCCESS:
                mylog.error("Something screwed up in RebuildHierarchy at level %s", level)
                return 0
            if self.LevelArray.GetLevel(level+1) == None:
                mylog.info("RefineUpTo: Level %s has nothing", level+1)
                break
            #print self.LevelArray.GetLevel(level+1)
            # Now we iterate over every grid on the next level, and call the
            # initializer function on each
            for Temp in self.LevelArray[level+1]:
                #print "Initializing grid on level", level+1
                self.InitializeGrid(Temp.GridData)
        # Now project back upwards
        for level in range(level+1, 0, -1):
            lUp = self.LevelArray.GetLevel(level-1)
            mylog.info("Projecting up from %s to %s", level, level-1)
            for Temp in self.LevelArray[level]:
                mylog.debug("Actual level: %s (%s) and %s (%s)", \
                    Temp.GridData.Level, Temp.GridData.GridSize(), \
                    lUp.GridData.Level, lUp.GridData.GridSize())
                Temp.GridData.ProjectSolutionToParentGrid(lUp.GridData)

    def GetCellPositions(self, grid):
        # This is probably quite slow, but, I think it is as fast as can be while still
        # retaining the required precision.  Fortunately, it's defined here, so
        # we can change it if and when I am shown to be wrong.
        ei = EnzoInterface
        size = grid.GridSize()
        # We will default to 64 here, which may be a mistake
        #print "SIZE",size
        CellPositions = zeros((size, 3), Float64)
        # Okay, we can do this the old way, or the new way
        ind = indices((ei.intArray_getitem(grid.GridDimension,0), \
                       ei.intArray_getitem(grid.GridDimension,1), \
                       ei.intArray_getitem(grid.GridDimension,2)))
        for dim in range(3):
            l = ei.EFloatArray_getitem(grid.GridLeftEdge, dim)
            w = ei.EFloatArray_getitem(ei.EFloatDimArray_getitem(grid.CellWidth, dim), 0)
            CellPositions[:,dim] = (ind[dim,:].flat + 0.5) * w + l
        return CellPositions
        # The old way, but non-functional:
        for dim in range(3):
            for i in range(size):
                tt = ei.EFloatArray_getitem(grid.CellLeftEdge, dim)
                t1 = ei.EFloatArray_getitem(ei.EFloatDimArray_getitem(grid.CellLeftEdge, dim), i)
                t2 = ei.EFloatArray_getitem(ei.EFloatDimArray_getitem(grid.CellWidth, dim), i)
                print tt, t1, t2
                CellPositions[i,dim] = t1 + 0.5 * t2
        return CellPositions
