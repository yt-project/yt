"""
New problems are created with this function

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

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


from yt.enki import *

def InitializeNew(ThisProblemType, ParameterFilename=None, InitArgs = {}, \
                  OutputFilename=None):
    """
    Here we accept up to two parameters.  One is a function, which needs to
    meet a set number of requirements.  The other is a n optional parameter
    file, to initialize the global variables with.

    @param ThisProblemType: a problem type
    @type ThisProblemType: L{ProblemType<ProblemType>}
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
    MyProblem = ThisProblemType(MetaData, TopGrid, LevelArray, InitArgs)
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
