"""
Collapse Test intiializer

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/

Todo: Add remainder of sphere types
Todo: Add multispecies abundance setting
Todo: Get TemperatureUnits from reasonable source, not hardcoded
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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

from yt.enki.mes import *
#from EnzoProblem import ProblemType

class CollapseTestProblem(ProblemType):
    def InitializeMetaData(self):
        #self.Defaults["CollapseTestNumberOfSpheres"] = 1
        self.Defaults["CollapseTestRefineAtStart"]  = TRUE
        self.Defaults["CollapseTestUseParticles"]   = FALSE
        self.Defaults["CollapseTestUseColour"]      = FALSE
        self.Defaults["CollapseTestInitialTemperature"] = 1000

        self.Defaults["CollapseTestUniformDensity"]=0
        self.Defaults["CollapseTestVelocityShockMagnitude"]=0
        self.Defaults["CollapseTestVelocityShockWidth"]=0
        self.Defaults["CollapseTestVelocityShockDirection"]=0

    def Initialize(self):
        # First we set up the datalabels
        self.InitializeFields()
        self.InitializeGrid(self.TopGrid.GridData)
        EnzoInterface.AddLevel(self.LevelArray.la, self.TopGrid, 0)
        # Now, if we want to refine at start, we call RebuildHierarchy until it
        # doesn't refine any more
        if self["CollapseTestRefineAtStart"] == TRUE:
            self.RefineUpTo(self["MaximumRefinementLevel"])

    def InitializeGrid(self, grid):
        # First, set up the fields
        # Then, iterate over each sphere, giving it an array of coordinates and
        # getting back an array of attributes
        self.InitializeFieldsInGrid(grid)
        xyz = self.GetCellPositions(grid)
        # Alright, now we set some defaults
        dens = self.Fields[self.FieldIndex["Density"]].ravel()
        dens[:] = 1.0
        temp = na.ones(dens.shape, 'float32').ravel() * self["CollapseTestInitialTemperature"]
        vel = na.zeros((3,temp.shape[0]), 'float32')
        vel[0,:] = self.Fields[self.FieldIndex["x-velocity"]][:].ravel()
        vel[1,:] = self.Fields[self.FieldIndex["y-velocity"]][:].ravel()
        vel[2,:] = self.Fields[self.FieldIndex["z-velocity"]][:].ravel()
        #TE = self.Fields[self.FieldIndex["TotalEnergy"]].ravel()
        for sphere in self.InitArgs["Spheres"]:
            sphere.GetDensity(self, grid, xyz, dens)
            sphere.GetTemperature(self, grid, xyz, temp)
            sphere.GetVelocity(self, grid, xyz, vel)
        # These next two assignments used to be flattened.
        sh = self.Fields[self.FieldIndex["GasEnergy"]].shape
        self.Fields[self.FieldIndex["GasEnergy"]] = \
            na.reshape(temp / 206800 /((5./3.-1.0)*0.6), sh)
        self.Fields[self.FieldIndex["TotalEnergy"]] = \
            na.reshape(temp / 206800 /((5./3.-1.0)*0.6), sh)
            #temp / self["TemperatureUnits"]/((self["Gamma"]-1.0)*self["mu"])
        self.FlushFieldsToGrid(grid)

class CollapseTestSphere:
    def __init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Density = Density
        self.Temperature = Temperature
        self.Velocity = Velocity
        self.Position = Position
        self.Radius = Radius
        self.CoreRadius = CoreRadius
    def GetDensity(self, prob, grid, xyz, dens):
        r = ((xyz[:,0]-self.Position[0])**2.0 \
          +  (xyz[:,1]-self.Position[1])**2.0 \
          +  (xyz[:,2]-self.Position[2])**2.0)**0.5
        dI = na.where(r < self.Radius)
        if dI[0].shape == 0:
            return
        dens[dI] = self.MyGetDensity(prob, grid, xyz[dI], r[dI])

    def GetVelocity(self, prob, grid, xyz, vel):
        r = ((xyz[:,0]-self.Position[0])**2.0 \
          +  (xyz[:,1]-self.Position[1])**2.0 \
          +  (xyz[:,2]-self.Position[2])**2.0)**0.5
        dI = na.where(r < self.Radius)
        if dI[0].shape == 0:
            return
        for i in range(3):
            vi = vel[i,:]
            vi[dI] = self.GetMyVelocity(prob, grid, xyz[dI], r[dI], i)
            vel[i,:] = vi
        

    def GetTemperature(self, prob, grid, xyz, temp):
        r = ((xyz[:,0]-self.Position[0])**2.0 \
          +  (xyz[:,1]-self.Position[1])**2.0 \
          +  (xyz[:,2]-self.Position[2])**2.0)**0.5
        dI = na.where(r < self.Radius)
        if dI[0].shape == 0:
            return
        temp[dI] = self.GetMyTemperature(prob, grid, xyz[dI], r[dI])

    def GetMyDensity(self, prob, grid, xyz, r):
        return na.ones(r.shape,'float32')
    def GetMyTemperature(self, prob, grid, xyz, r):
        return na.ones(r.shape,'float32')*prob["CollapseTestInitialTemperature"]
    def GetMyVelocity(self, prob, grid, xyz, r, dim):
        return na.ones(r.shape,'float32')*10

class CollapseTestSphereUniform(CollapseTestSphere):
    def __init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Type = "Uniform"
        CollapseTestSphere.__init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius)
    def MyGetDensity(self, prob, grid, xyz, r):
        return na.ones(r.shape, 'float32') * self.Density

class CollapseTestSphereR2(CollapseTestSphere):
    def __init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Type = "R^-2"
        CollapseTestSphere.__init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius)
    def MyGetDensity(self, prob, grid, xyz, r):
        return self.Density * (r / self.Radius)**-2.0

