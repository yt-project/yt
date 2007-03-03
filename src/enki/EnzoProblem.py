from yt.enki import *

class ProblemType:
    def __init__(self, MetaData, TopGrid, LevelArray, InitArgs):
        self.MetaData = MetaData
        self.TopGrid = TopGrid
        self.LevelArray = LevelArray
        self.Defaults = {}
        self.InitArgs = InitArgs
        self.count = 0
        self.FieldTypes = []

    def __getitem__(self, item):
        # We set precedence here of the InitArgs over the MetaData
        if self.InitArgs.has_key(item):
            return self.InitArgs[item]
        if self.Defaults.has_key(item):
            return self.Defaults[item]
        try:
            return MetaData.__getattr__(item)
        except AttributeError:
            exec("return EnzoInterface.cvar.%s" % (item))

    def Initialize(self):
        mylog.warning("Oops, look like your problem initializer wasn't defined!")

    def InitializeMetaData(self):
        pass

    def AddField(self, field, g):
        EnzoInterface.charSArray_setitem(EnzoInterface.cvar.DataLabel, count, field)
        self.FieldTypes.append(FieldTypes[field])
        self.count += 1

    def InitializeFields(self, grid):
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

    def GetCellPositions(self, grid):
        # This is probably quite slow, but, I think it is as fast as can be while still
        # retaining the required precision.  Fortunately, it's defined here, so
        # we can change it if and when I am shown to be wrong.
        ei = EnzoInterface
        size = grid.GridSize()
        # We will default to 64 here, which may be a mistake
        CellPositions = zeros((size, 3), Float64)
        for dim in range(3):
            for i in range(size):
                CellPositions[i,dim] = \
                    ei.arrayEFloat_getitem(ei.arrayEFloatDim_getitem(grid.CellLeftEdge, dim), i) + \
              0.5 * ei.arrayEFloat_getitem(ei.arrayEFloatDim_getitem(grid.CellWidth, dim), i)
        return CellPositions
