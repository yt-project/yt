from yt.enki import *

class CollapseTestSphere:
    def __init__(self, Type, Density, Temperature, Velocity, Position, Radius, CoreRadius):
        self.Type = Type
        self.Density = Density
        self.Temperature = Temperature
        self.Velocity = Velocity
        self.Position = Position
        self.Radius = Radius
        self.CoreRadius = CoreRadius
    def GetDensity(self, prob, grid, xyz, r, dens):
        dI = where(r < self.Radius)
        dens[dI] = maximum(self.MyGetDensity(prob, grid, xyz[dI], r[dI]), \
                       dens[dI])
    def GetTemperature(self, prob, grid, xyz, r, temp):
        dI = where(r < self.Radius)
        temp[dI] = maxximum(self.MyGetTemperature(prob, grid, xyz[dI], r[dI]), \
                       temp[dI])

class CollapseTestSphereUniform(CollapseTestSphere):
    def MyGetDensity(self, prob, grid, xyz, r):
        return ones(r.shape, Float32) * self.Density
    def MyGetTemperature(self, prob, grid, xyz, r):
        return zeros(r.shape)


class CollapseTestSphereR2(CollapseTestSphere):
    def MyGetDensity(self, prob, grid, xyz, r):
        return ones(r.shape, Float32) * self.Density * \
               (r / self.Radius)**-2.0
    def MyGetTemperature(self, prob, grid, xyz, r):
        return zeros(r.shape)


class CollapseTestProblem(EnzoProblem):
    def MetaDataInitialize(self):
        #self.Defaults["CollapseTestNumberOfSpheres"] = 1
        self.Defaults["CollapseTestRefineAtStart"]  = TRUE
        self.Defaults["CollapseTestUseParticles"]   = FALSE
        self.Defaults["CollapseTestUseColour"]      = FALSE
        self.Defaults["CollapseTestInitialTemperature"] = 1000

        self.Defaults["CollapseTestUniformDensity"]=0
        self.Defaults["CollapseTestVelocityShockMagnitude"]=0
        self.Defaults["CollapseTestVelocityShockWidth"]=0
        self.Defaults["CollapseTestVelocityShockDirection"]=0

        # Then the routine has some stuff we don't care about, since we're
        # using Python for  all the variable handling
        
        #CollapseTestSphereType[MAX_SPHERES];
        #CollapseTestSphereDensity[MAX_SPHERES],
        #CollapseTestSphereTemperature[MAX_SPHERES],
        #CollapseTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
        #CollapseTestUniformVelocity[MAX_DIMENSION];
        #CollapseTestUniformDensityGradient[MAX_DIMENSION]={0}; //sets all to 0
        #CollapseTestUniformTotalEnergyGradient[MAX_DIMENSION]={0};
  
    def Initialize(self):
        # First we set up the datalabels
        self.InitializeFieldTypes()
        self.InitializeGrid(self.TopGrid)

    def InitializeGrid(self, grid):
        # First, set up the fields
        # Then, iterate over each sphere, giving it an array of coordinates and
        # getting back an array of attributes
        xyz = self.GetCellPositions(grid)
        r = (xyz[:,0]**2.0 + xyx[:,1]**2.0 + xyz[:,2]**2.0)**0.5
        dens = zeros(r.shape, Float32)
        temp = zeros(r.shape, Float32)
        for sphere in sel.InitArgs["Sphere"]:
            dens = sphere.GetDensity(self, grid, xyz, r, dens)
            sphere.GetTemperature(self, grid, xyz, r, temp)
            #sphere.GetVelocity(self, grid, xyz, r, vel)
        # Alright, now we just need to copy in the baryon fields
