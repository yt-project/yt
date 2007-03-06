"""
New problems are created with this function

Enki, the creator of Humankind
==============================

In Sumerian myth, Enki lay asleep in the depths of the primeval ocean, unable
to hear the lament of the gods as they complained about the difficulty of
cultivating wheat and making bread. Eventually the primeval sea, Nammu
brought the gods' tears to Enki. Enki, as the god of wisdom, was expected to
devise a solution, so he solicited Nammu and the birth-goddess Ninmah to use
clay to form the first men, who would toil and farm so that the gods could
relax.

In later Akkadian or Babylonian Cosmology there were six generations of Gods
that led to the creation of the Younger (Igigi) divinities of the Anunaki. In
the seventh generation (Akkadian "Shappatu" hence the Hebrew Shabbath =>
English Sabbath), the younger Gods went on strike, put down their tools and
refused to keep the creation working. In the Babylonian creation myth the
Enuma Elish, Abzu, the water lord, threatens to take back the creation with a
universal flood, but Enki averts the threat by imprisoning Abzu beneath the
Earth. Kingu, his son, informs his mother, Abzu's wife, the serpentine Tiamat
(Ti = Life, Ama = mother, Biblical tehwom = the deeps), and in anger she
threatens to take back the whole of creation. The Gods gather in terror, but
Enlil (his place in the Enuma Elish is later taken by Enki's son Marduk)
subdues and slays Tiamat with the arrows of his winds which he shoots down
her throat. Marduk, Enki's son, (earlier Enlil, Enki's half-brother), takes
from Tiamat the Tablet of Destinies.

In the Sumerian poem 'Ninurta and the Turtle' it is the god Enki, rather
than Enlil, who holds the tablet. Both this poem and the Akkadian AnzÃ»
poem concern the theft of the tablet by the bird Imdugud (Sumerian) or AnzÃ»
(Akkadian) . Supposedly, whoever possessed the tablets ruled the
universe.

The tablet can be compared with the concept of the Me, divine decrees that
are the special attribute of Enki.

But the problem created by the "strike of the Gods" remains, how is creation
to continue? Enki proposes that the Gods make humankind as their servant, and
give humans the task of keeping creation going. It is agreed, and Enki forms
humanity out of the red earth (Hebrew Adamah), mingled with the red blood of
the God Kingu, slain for his part in Tiamat's attack. Enlil fills his lungs
with air (Hebrew ruach, Greek pneuma, Latin spiritus), and humans are alive.
In this way, Humanity is given the task of maintaining the balance of nature
and keeping the created order in place.

Another myth, "Enki and Adapa", tells of how humanity loses the chance at
immortality. Adapa U-an (Berossus' Oannes), who is Abgallu (Ab = Water, Gal =
Great, Lu = Man) (Akkadian Apkallu), Enki's advisor, to the first king of
Eridu, Allulim, inadvertently breaks the wings of the South Wind, Ninlil (See
Lilith) (Nin = Lady, Lil = Air), daughter of Anu (the Heavens) and wife to
Enlil, King of the Gods. In terror at the thought of their retribution, Adapa
seeks the advice of Enki. Enki advises that Adapa make a deep and sincere
atonement, but advises Adapa to eat nothing given to him by the Gods, as he
will probably be given the food of death, out of their anger at his deeds.
Adapa takes Enki's advice, but the Gods, so impressed by the sincerity of
Adapa's sorrow and grief as to what he did, offered instead the fruit of
immortality. Adapa remembering Enki's words, refuses, and so misses out on
the chance of eternal life.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""


from yt.enki import *

def InitializeNew(ProblemType, ParameterFilename=None, InitArgs = {}, \
                  OutputFilename=None):
    """
    Here we accept up to two parameters.  One is a function, which needs to
    meet a set number of requirements.  The other is a n optional parameter
    file, to initialize the global variables with.

    @param ProblemType: a problem type
    @type ProblemType: L{ProblemType<ProblemType>}
    @keyword ParameterFilename: if this is not None, then it is passed in as a
            file object to ReadParameterFile.
    @type ParameterFilename: string
    @keyword MetaDataInitializerFunction: if supplied, this function gets called to set the
            various parameters (this is in lieu of ParameterFilename, or
            possibly in addition, for overriding)
    @keyword InitArgs: a dictionary of parameters to be passed to both the
            ProblemInitializerFunction and MetaDataInitializerFunction
    @type InitArgs: dict
    @keyword OutputFilename: a basename to be passed in to WriteAllData
    @type OutputFilename: string
    """
    EnzoInterface.cvar.debug = TRUE
    MetaData = EnzoInterface.TopGridData()
    TopGrid = EnzoInterface.HierarchyEntry()
    Exterior = EnzoInterface.ExternalBoundary()
    LevelArray = LevelArrayWrapper()
    Initialdt = EnzoInterface.new_Float()
    EnzoInterface.Float_assign(Initialdt, 0.0)
    EnzoInterface.SetDefaultGlobalValues(MetaData)
    # And that's it for initializing variables.  Note that we will need to
    # initialize various variables, but we haven't yet done so.  Some of them
    # will be dealt with by the parameterfile, but, again, this is optional.
    # If you want to have everything defined in python -- which you might, if
    # you are conducting statistical realizations, for instance -- you can.
    if ParameterFilename:
        f=open(ParameterFilename,"r")
        retVal = EnzoInterface.ReadParameterFile(f, MetaData, Initialdt)
        # We should probably check the return value.
        if retVal != 1:
            mylog.warning("Uh oh, parameter file (%ss) was not read correctly!", ParameterFilename)
            return 0
    MyProblem = ProblemType(MetaData, TopGrid, LevelArray, InitArgs)
    MyProblem.InitializeMetaData()
    TopGridInitialize(MetaData, TopGrid)
    # Now we do the standard Top Grid Initialization functions.
    MyProblem.Initialize()
    ExteriorInitialize(Exterior, TopGrid, MetaData)
    if OutputFilename:
        # Call WriteAllData
        EnzoInterface.WriteAllData(OutputFilename, 0, TopGrid, MetaData, Exterior, -1)

def TopGridInitialize(MetaData, TopGrid):
    """
    All the stuff enzo does in InitializeNew before calling the problem
    initializers.

    @param MetaData: the TopGridData() entry for the sim
    @type MetaData: MetaDataType
    @param TopGrid: the root HierarchyEntry
    @type TopGrid: HierarchyEntry
    """
    # We're going to ignore the parameter checking, because if you people can't
    # be grownups, I'm not going to do your work *for* you.  Besides, if you
    # set the TopGridRank < 0, don't you *deserve* to get an exception thrown
    # at you at some unpredictable point down the road?
    #
    # First we add the Ghost Zones (Pac-Man beware, as they are immune to Power
    # Pellets within these regions)
    for i in range(MetaData.TopGridRank):
        EnzoInterface.intArray_setitem(MetaData.TopGridDims, i, \
            EnzoInterface.intArray_getitem(MetaData.TopGridDims, i) + \
            2 * DEFAULT_GHOST_ZONES)
    # Okay, set the ghost zones, now to set up the grid
    TopGrid.GridData = EnzoInterface.grid()
    tg = TopGrid.GridData
    tg.PrepareGrid(MetaData.TopGridRank, MetaData.TopGridDims, \
            EnzoInterface.cvar.DomainLeftEdge, EnzoInterface.cvar.DomainRightEdge, \
            MetaData.NumberOfParticles, 0)
    tg.SetTime(MetaData.Time)
    tg.SetHydroParameters(MetaData.CourantSafetyNumber, \
            MetaData.PPMFlatteningParameter, \
            MetaData.PPMDiffusionParameter, \
            MetaData.PPMSteepeningParameter)
    tg.SetGravityParameters(MetaData.GravityBoundary)
    # Now we fix the topgrid back up, taking out Pinky, Inky, Blinky, and Clyde
    # ( http://www.flickr.com/photos/kenthenderson/92213471/ )
    for i in range(MetaData.TopGridRank):
        EnzoInterface.intArray_setitem(MetaData.TopGridDims, i, \
            EnzoInterface.intArray_getitem(MetaData.TopGridDims, i) - \
            2 * DEFAULT_GHOST_ZONES)
    TopGrid.NextGridThisLevel = None
    TopGrid.ParentGrid = None
    TopGrid.NextGridNextLevel = None  # Reset in init function
    # And now we're all done, and ready to call the problem initializer.
    # Take it away, Don Pardo!

def ExteriorInitialize(Exterior, TopGrid, MetaData):
    ei = EnzoInterface
    if Exterior.AmIPrepared():
        return
    Exterior.Prepare(TopGrid.GridData)
    if MetaData.BoundaryConditionName:
        f = open(MetaData.BoundaryConditionName, "r")
        Exterior.ReadExternalBoundary(f)
        ei.fclose(f)
    else:
        dummy = ei.new_floatArray(3)
        for dim in range(MetaData.TopGridRank):
            Exterior.InitializeExternalBoundaryFace(dim, \
                    ei.BoundaryTypeArray_getitem(MetaData.LeftFaceBoundaryCondition, dim), \
                    ei.BoundaryTypeArray_getitem(MetaData.RightFaceBoundaryCondition, dim), \
                    dummy, dummy)
        Exterior.InitializeExternalBoundaryParticles(MetaData.ParticleBoundaryType)


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
