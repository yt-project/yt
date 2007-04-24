"""
Collapse Test intiializer

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Add remainder of sphere types
@todo: Add multispecies abundance setting
@todo: Get TemperatureUnits from reasonable source, not hardcoded

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
        temp = na.ones(dens.shape, nT.Float32).ravel() * self["CollapseTestInitialTemperature"]
        vel = na.zeros((3,temp.shape[0]), nT.Float32)
        vel[0,:] = self.Fields[self.FieldIndex["x-velocity"]][:].ravel()
        vel[1,:] = self.Fields[self.FieldIndex["y-velocity"]][:].ravel()
        vel[2,:] = self.Fields[self.FieldIndex["z-velocity"]][:].ravel()
        #TE = self.Fields[self.FieldIndex["TotalEnergy"]].ravel()
        for sphere in self.InitArgs["Spheres"]:
            sphere.GetDensity(self, grid, xyz, dens)
            sphere.GetTemperature(self, grid, xyz, temp)
            sphere.GetVelocity(self, grid, xyz, vel)
        self.Fields[self.FieldIndex["GasEnergy"]] = \
            temp / 206800 /((5./3.-1.0)*0.6)
            #temp / self["TemperatureUnits"]/((self["Gamma"]-1.0)*self["mu"])
        self.Fields[self.FieldIndex["TotalEnergy"]] = \
            temp / 206800 /((5./3.-1.0)*0.6)
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
        return na.ones(r.shape,nT.Float32)
    def GetMyTemperature(self, prob, grid, xyz, r):
        return na.ones(r.shape,nT.Float32)*prob["CollapseTestInitialTemperature"]
    def GetMyVelocity(self, prob, grid, xyz, r, dim):
        return na.ones(r.shape,nT.Float32)*10

class CollapseTestSphereUniform(CollapseTestSphere):
    def __init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Type = "Uniform"
        CollapseTestSphere.__init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius)
    def MyGetDensity(self, prob, grid, xyz, r):
        return na.ones(r.shape, nT.Float32) * self.Density

class CollapseTestSphereR2(CollapseTestSphere):
    def __init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Type = "R^-2"
        CollapseTestSphere.__init__(self, Density, Temperature, Velocity, Position, Radius, CoreRadius)
    def MyGetDensity(self, prob, grid, xyz, r):
        return self.Density * (r / self.Radius)**-2.0

