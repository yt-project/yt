"""
Fields derived from other data

Unless otherwise specified in the field name, all fields should be returned in
some form of natural, non-code units.  fieldInfo is of the form:

fieldInfo[fieldName] = (Units, ProjectedUnits, TakeLog, Function)

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@author: U{John Wise<http://www.slac.stanford.edu/~jwise/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{jwise@slac.stanford.edu<mailto:jwise@slac.stanford.edu>}
"""

from yt.lagos import *

try:
    import EnzoFortranWrapper
except:
    mylog.warning("EnzoFortranWrapper not properly installed -- no fortran fields available!")

# Add the info for any non-derived fields up here.  For any added derived
# fields, add it immediately after the function definition.

fieldInfo = {}

# fieldInfo has the following structure:
#   key == field
#   value == tuple ( naturalUnits, projUnits, takeLog, generatorFunction )

fieldInfo["Density"] = ("g cm^-3", "g cm^-2", True, None)
fieldInfo["Temperature"] = ("K", None, False, None)
fieldInfo["HII_Fraction"] = ("mass fraction", None, True, None)
fieldInfo["H2I_Fraction"] = ("mass fraction", None, True, None)

# These are all the fields that should be logged when plotted
# NOTE THAT THIS IS OVERRIDEN BY fieldInfo !
log_fields = [ "Density_Squared", "k23", "k22", "k13" ]

colormap_dict = {"Temperature":"Red Temperature"}

def Entropy(self, fieldName):
    self[fieldName] = self["Density"]**(-2./3.) * \
                           self["Temperature"]
fieldInfo["Entropy"] = (None, None, True, Entropy)

def DynamicalTime(self, fieldName):
    """
    The formulation for the dynamical time is:
    M{sqrt(3pi/(16*G*rho))} or M{sqrt(3pi/(16G))*rho^-(1/2)}
    Note that we return in our natural units already
    """
    G = self.hierarchy["GravitationalConstant"]
    t_dyn_coeff = sqrt(3*pi/(16*G)) * self.hierarchy["Time"]
    self[fieldName] = self["Density"]**(-1./2.) * t_dyn_coeff
fieldInfo["DynamicalTime"] = ("s", None, True, DynamicalTime)

def H2FormationTime(self, fieldName):
    """
    This calculates the formation time using k22 and k21
    """
    self[fieldName] = abs( (self["H2I_Density"]/2.0) / \
                       ( self["HI_Density"] \
                       * self["HI_Density"] \
                       * self.hierarchy.rates[self["Temperature"],"k22"] \
                       * self.hierarchy.rates.params["kunit_3bdy"] \
                       * self["HI_Density"] \
                       + self["HI_Density"] \
                       * self["HI_Density"] \
                       * self.hierarchy.rates[self["Temperature"],"k21"] \
                       * self.hierarchy.rates.params["kunit_3bdy"] \
                       * self["H2I_Density"]/2.0 ) ) 
fieldInfo["H2FormationTime"] = ("s", None, True, H2FormationTime)

def H2DissociationTime(self, fieldName):
    """
    This calculates the dissociation time using k23 and k13dd
    """
    self[fieldName] = abs( (self["H2I_Density"]/2.0) / \
                       ( self["H2I_Density"]/2.0 \
                       * self.hierarchy.rates[self["Temperature"],"k23"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self["H2I_Density"]/2.0 + \
                         self["HI_Density"] \
                       * self["k13DensityDependent"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self["H2I_Density"]/2.0) )
fieldInfo["H2DissociationTime"] = ("s", None, True, H2DissociationTime)

def k23DissociationTime(self, fieldName):
    """Just k23"""
    self[fieldName] = abs( self["H2I_Density"] / \
                       ( self["H2I_Density"]*2 \
                       * self.hierarchy.rates[self["Temperature"],"k23"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self["H2I_Density"]*2) )
fieldInfo["k23DissociationTime"] = ("s" , None, True, k23DissociationTime)

def k13DissociationTime(self, fieldName):
    """Just k13dd"""
    self[fieldName] = abs( self["H2I_Density"] / \
                       ( self["HI_Density"] \
                       * self["k13DensityDependent"] \
                       * self.hierarchy.rates.params["kunit"] \
                       * self["H2I_Density"] ) )
fieldInfo["k13DissociationTime"] = ("s" , None, True, k13DissociationTime)

def CIEOpticalDepthFudge(self, fieldName):
    """Optical depth from CIE"""
    tau = maximum((self["NumberDensity"]/2.0e16)**2.8, 1.0e-5)
    self[fieldName] = minimum(1.0, (1.0-exp(-tau))/tau)
fieldInfo["CIEOpticalDepthFudge"] = (None, None, False, CIEOpticalDepthFudge)

def compH2DissociationTime(self, fieldName):
    """k13/k23"""
    self[fieldName] = self["k13DissociationTime"] \
                         / self["k23DissociationTime"]
fieldInfo["compH2DissociationTime"] = ("t_k13/t_k23" , None, True, compH2DissociationTime)

def k13DensityDependent(self, fieldName):
    """k13dd"""
    dom = self.hierarchy["Density"] / 1.67e-24
    nh = minimum(self["HI_Density"]*dom, 1.0e9)
    k1 = self.hierarchy.rates[self["Temperature"],"k13_1"]
    k2 = self.hierarchy.rates[self["Temperature"],"k13_2"]
    k3 = self.hierarchy.rates[self["Temperature"],"k13_3"]
    k4 = self.hierarchy.rates[self["Temperature"],"k13_4"]
    k5 = self.hierarchy.rates[self["Temperature"],"k13_5"]
    k6 = self.hierarchy.rates[self["Temperature"],"k13_6"]
    k7 = self.hierarchy.rates[self["Temperature"],"k13_7"]
    self[fieldName] = maximum( 10.0**( \
            k1-k2/(1+(nh/k5)**k7) \
          + k3-k4/(1+(nh/k6)**k7) )\
          , 1e-30 )
fieldInfo["k13DensityDependent"] = ("cm^-3" , None, True, k13DensityDependent)

def H2EquilibriumBalance(self, fieldName):
    """L{H2FormationTime}/L{H2DissociationTime}"""
    self[fieldName] = abs(
                        self["H2FormationTime"] \
                      / self["H2DissociationTime"] )
fieldInfo["H2EquilibriumBalance"] = ("t_formation/t_dissociation" , None, True, H2EquilibriumBalance)

def H2FormationDynamicalBalance(self, fieldName):
    """t_dyn / L{H2FormationTime}"""
    self[fieldName] = abs( \
                        self["DynamicalTime"]
                     /  self["H2FormationTime"] )
fieldInfo["H2FormationDynamicalBalance"] = ("t_dyn/t_k22", None, True, H2FormationDynamicalBalance)

def H2DissociationDynamicalBalance(self, fieldName):
    """t_dyn / L{H2DissociationTime}"""
    self[fieldName] = abs( \
                        self["DynamicalTime"]
                     /  self["H2DissociationTime"] )
fieldInfo["H2DissociationDynamicalBalance"] = ("t_dyn/t_k23", None, True, H2DissociationDynamicalBalance)

def DCComp(self, fieldName):
    """t_dyn / courant_time"""
    self[fieldName] = self["DynamicalTime"]/self["CourantTimeStep"]
fieldInfo["DCComp"] = ("t_dyn/t_courant", None, True, DCComp)

def NumberDensityCode(self, fieldName):
    """amu/cc, in code units."""
    # We are going to *try* to use all the fields, and fallback when we run out
    self[fieldName] = zeros(self["HI_Density"].shape, self["HI_Density"].type())
    if self.hierarchy["MultiSpecies"] > 0:
        self[fieldName] += self["HI_Density"] / 1.0
        self[fieldName] += self["HII_Density"] / 1.0
        self[fieldName] += self["HeI_Density"] / 4.0
        self[fieldName] += self["HeII_Density"] / 4.0
        self[fieldName] += self["HeIII_Density"] / 4.0
        self[fieldName] += self["Electron_Density"] / 1.0
    if self.hierarchy["MultiSpecies"] > 1:
        self[fieldName] += self["HM_Density"] / 1.0
        self[fieldName] += self["H2I_Density"] / 2.0
        self[fieldName] += self["H2II_Density"] / 2.0
    if self.hierarchy["MultiSpecies"] > 2:
        self[fieldName] += self["DI_Density"] / 2.0
        self[fieldName] += self["DII_Density"] / 2.0
        self[fieldName] += self["HDI_Density"] / 3.0
fieldInfo["NumberDensityCode"] = (None, None, True, NumberDensityCode)

def NumberDensity(self, fieldName):
    """L{NumberDensityCode} M{* DensityConversion/m_h} S{->} particles/cc"""
    self[fieldName] = self["NumberDensityCode"] * \
                      ((self.hierarchy["Density"]/1.67e-24))
fieldInfo["NumberDensity"] = ("cm^-3", "cm^-2", True, NumberDensity)

def Outline(self, fieldName):
    """Just the level, for outlines of grid structure"""
    self[fieldName] = ones(self.myChildMask.shape) * self.Level + 1
fieldInfo["Outline"] = ("Level", "LevelSum", False, Outline)

def SoundSpeed(self, fieldName):
    """M{t_sound = sqrt(Gamma*Pressure/Density)}"""
    self[fieldName] = self.hierarchy["x-velocity"] * ( \
             self.hierarchy["Gamma"]*self["Pressure"] / \
             self["Density"] )**(1.0/2.0)
fieldInfo["SoundSpeed"] = ("cm/s", None, True, SoundSpeed)

def MachNumber(self, fieldName):
    """M{|v|/t_sound}"""
    self[fieldName] = (self.hierarchy["x-velocity"] * ( \
        self["x-velocity"]**2.0 + \
        self["y-velocity"]**2.0 + \
        self["z-velocity"]**2.0 )**(1.0/2.0)) / self["SoundSpeed"]
fieldInfo["MachNumber"] = (None, None, False, MachNumber)
        

def Pressure(self, fieldName):
    """M{(Gamma-1.0)*rho*E}"""
    self[fieldName] = (self.hierarchy["Gamma"] - 1.0) * \
                            self["Density"] * self["Gas_Energy"]
fieldInfo["Pressure"] = (None, None, True, Pressure)

def CourantTimeStep(self, fieldName):
    """
    Very simple, just a quick look at the courant timestep
    No simplification, just as done in calc_dt.[sr]c
    Additionally, we assume that GridVelocity == 0 in all dims
    """
    t1 = self.dx * self.hierarchy["cm"] / \
         (self["SoundSpeed"] + \
            (abs(self["x-velocity"])* \
             self.hierarchy["x-velocity"]))
    t2 = self.dy * self.hierarchy["cm"] / \
         (self["SoundSpeed"] +
            (abs(self["y-velocity"])* \
             self.hierarchy["y-velocity"]))
    t3 = self.dz * self.hierarchy["cm"] / \
         (self["SoundSpeed"] +
            (abs(self["z-velocity"])* \
             self.hierarchy["z-velocity"]))
    self[fieldName] = minimum(minimum(t1,t2),t3)
    del t1, t2, t3
fieldInfo["CourantTimeStep"] = ("s", None, True, CourantTimeStep)

def H2FormationCourant(self, fieldName):
    """L{H2FormationTime}/L{CourantTimeStep}"""
    self[fieldName] = self["H2FormationTime"] / self["CourantTimeStep"]
fieldInfo["H2FormationCourant"] = ("t_formation/t_courant", None, True, H2FormationCourant)

def RadialVelocity(self, fieldName):
    """
    Velocity toward object.center, subtracting off the object.bulkVelocity
    """
    # We assume a "center" variable is assigned to the hierarchy

    # Correct for the bulk velocity if it exists
    if self.hierarchy.bulkVelocity != None:
        xv = self["x-velocity"] - self.hierarchy.bulkVelocity[0]
        yv = self["y-velocity"] - self.hierarchy.bulkVelocity[1]
        zv = self["z-velocity"] - self.hierarchy.bulkVelocity[2]

    self.generateCoords()
    delx = self.coords[0,:] - self.hierarchy.center[0]
    dely = self.coords[1,:] - self.hierarchy.center[1]
    delz = self.coords[2,:] - self.hierarchy.center[2]
    radius = sqrt((delx**2.0) + (dely**2.0) + (delz**2.0))
    self[fieldName] = ( \
           ( delx * xv + \
             dely * yv +   \
             delz * zv ) / \
           ( radius) ) * self.hierarchy.conversionFactors["x-velocity"] / 100000.0
fieldInfo["RadialVelocity"] = ("km/s", None, False, RadialVelocity)

def RadiusCode(self, fieldName):
    """
    Distance from cell to object.center
    """
    self.generateCoords()
    delx = self.coords[0,:] - self.center[0]
    dely = self.coords[1,:] - self.center[1]
    delz = self.coords[2,:] - self.center[2]
    self[fieldName] = ((delx**2.0) + (dely**2.0) + (delz**2.0))**0.5
fieldInfo["RadiusCode"] = (None, None, True, RadiusCode)

def Radius(self, fieldName):
    """
    L{RadiusCode} in cm
    """
    self[fieldName] = self["RadiusCode"] * self.hierarchy["cm"]
fieldInfo["Radius"] = ("cm", None, True, Radius)

def Metallicity(self, fieldName):
    zSolar = 0.0204
    if self.hierarchy.has_key('Metal_Density'):
        self[fieldName] = self["Metal_Fraction"] / zSolar
    else:
        self[fieldName] = self["SN_Colour"] / self["Density"] / \
                          zSolar
        
fieldInfo["Metallicity"] = ("Z/Z_solar", None, True, Metallicity)

def CellMass(self, fieldName):
    """
    Calculates the mass (= rho * dx^3) in solar masses in each cell.
    Useful for mass-weighted plots.
    """
    msun = 1.989e33
    unitConversion = (self.hierarchy["Density"] / msun)
    self[fieldName] = ones(self["Density"].shape, Float64)
    self[fieldName] *= (self.dx) * self.hierarchy["cm"]
    self[fieldName] *= (self.dy) * self.hierarchy["cm"]
    self[fieldName] *= (self.dz) * self.hierarchy.units["cm"]
    self[fieldName] *= self["Density"]
    self[fieldName] *= unitConversion
fieldInfo["CellMass"] = ("Msun", None, True, CellMass)

def Volume(self, fieldName):
    """
    Useful for volume-weighted plots.
    """
    self[fieldName] = (self.dx * self.dy * self.dz) * \
                           ones(self.ActiveDimensions)
fieldInfo["Volume"] = (None, None, True, Volume)

def CoolingTime(self, fieldName):
    """
    This requires the hierarchy have an associated .rates and .cool, and
    additionally it calls the fortran wrapper.  Probably doesn't work right now.
    @todo: fix to work with yt.lagos.EnzoFortranWrapper
    """
    self[fieldName] = array(src_wrapper.runCoolMultiTime(self), Float64) \
                           * self.hierarchy["Time"]
fieldInfo["CoolingTime"] = ("s", None, True, CoolingTime)

def CoolingTimeDynamicalTimeComp(self, fieldName):
    """t_cool / t_dyn"""
    self[fieldName] = self["CoolingTime"]/self["DynamicalTime"]
fieldInfo["CoolingTimeDynamicalTimeComp"] = ("t_cool/t_dyn", None, True, CoolingTimeDynamicalTimeComp)

def CoolingTimeCourantTimeComp(self, fieldName):
    """t_cool / t_courant"""
    self[fieldName] = self["CoolingTime"]/self["CourantTimeStep"]
fieldInfo["CoolingTimeCourantTimeComp"] = ("t_cool/t_c", None, True, CoolingTimeCourantTimeComp)

def Gamma2(self, fieldName):
    """
    The minor H2 correction to the equation of state
    Note that we use the existing TempFromGE field, which is a little stupid
    since we're using Gamma2 to calculate temperature at a different time,
    but we will call it okay since Gamma2 should be a small correction
    """
    x = 6100.0 / self["TempFromGE"]
    x_gt = 2.5
    x_lt = 0.5*(5.0+2.0*x**2.0 * exp(x)/(exp(x)-1.0)**2)
    self[fieldName] = choose(greater(x,10.0),(x_lt,x_gt))
    nH2 = 0.5 * (self["H2I_Density"] + self["H2II_Density"])
    nOther = self["NumberDensityCode"] - nH2
    self[fieldName] = 1.0 + (nH2 + nOther) / \
        (nH2*self[fieldName] + nOther/(self.hierarchy["Gamma"]-1))
fieldInfo["Gamma2"] = (None, None, True, Gamma2)

def TempFromGE(self, fieldName):
    """
    We are going to calculate the temperature the same
    way it's calculated in s_r_c.  Helpful for rate solver checking.
    Note that this assumes DaulEnergyFormalism.  If you want it to work with
    RegularEnergyFormalism, do it yourself.  <smiley removed at request of
    etiquette board>
    """
    self[fieldName] = self["Pressure"] * self.hierarchy["Temp"] \
                    / (self["NumberDensityCode"])
    self[fieldName] *= (self["Gamma2"] - 1.0) \
                     / (self.hierarchy["Gamma"] - 1.0)
fieldInfo["TempFromGE"] = ("K", None, False, TempFromGE)

def TempComp(self, fieldName):
    """
    L{TempFromGE}/Temperature
    """
    self[fieldName] = self["TempFromGE"]/self["Temperature"]
fieldInfo["TempComp"] = ("t_ge/t_output",None,False,TempComp)

def HydrogenComp(self, fieldName):
    """
    (Hydrogen / Density) - 0.76
    """
    self[fieldName] = self["HI_Density"] + \
                      self["HII_Density"] + \
                      self["HM_Density"] + \
                      self["H2I_Density"] + \
                      self["H2II_Density"]
    self[fieldName] /= self["Density"]
    self[fieldName] -= 0.76
    self[fieldName] = abs(self[fieldName])
fieldInfo["HydrogenComp"] = ("rhoH/rho - 0.76",None,False,HydrogenComp)

def ColorSum(self, fieldName):
    """
    Sum up all the species fields
    """
    self[fieldName] = zeros(self["HI_Density"].shape, Float64)
    if self.hierarchy["MultiSpecies"] > 0:
        self[fieldName] += self["HI_Density"]
        self[fieldName] += self["HII_Density"]
        self[fieldName] += self["HeI_Density"]
        self[fieldName] += self["HeII_Density"]
        self[fieldName] += self["HeIII_Density"]
        self[fieldName] += self["Electron_Density"]
    if self.hierarchy["MultiSpecies"] > 1:
        self[fieldName] += self["HM_Density"]
        self[fieldName] += self["H2I_Density"]
        self[fieldName] += self["H2II_Density"]
    if self.hierarchy["MultiSpecies"] > 2:
        self[fieldName] += self["DI_Density"]
        self[fieldName] += self["DII_Density"]
        self[fieldName] += self["HDI_Density"]
fieldInfo["ColorSum"] = (None, None, True, ColorSum)

def DensityComp(self, fieldName):
    """
    Compare L{ColorSum} against the density
    """
    self[fieldName] = abs((self["ColorSum"] / self["Density"]) - 1.0)
fieldInfo["DensityComp"] = (None, None, True, DensityComp)
