import types
from collections import defaultdict
import numpy as na
from math import pi

mh = 1.67e-24
fieldInfo = {}

def add_field(name, function = None, **kwargs):
    if function == None:
        if kwargs.has_key("function"):
            function = kwargs.pop("function")
        else:
            # This will fail if it does not exist,
            # which is our desired behavior
            function = eval("_%s" % name)
    fieldInfo[name] = DerivedField(
        name, function, **kwargs
    )

class FieldDetector(defaultdict):
    params = defaultdict(lambda: 1)
    def __init__(self):
        self.requested = []
        defaultdict.__init__(self, lambda: na.ones(10))
    def __missing__(self, item):
        if fieldInfo.has_key(item):
            fieldInfo[item](self)
            return self[item]
        self.requested.append(item)
        return defaultdict.__missing__(self, item)
    def convert(self, item):
        return 1

class DerivedField:
    def __init__(self, name, function,
                 units = "", projected_units = "",
                 takeLog = True, validator = None):
        self.name = name
        self.function = function
        self.validator = validator
        self.takeLog = takeLog
    def get_values(self, data):
        pass
    def check_available(self, data):
        if self.validator != None:
            self.validator(self, data)
    def get_my_dependencies(self):
        pass
    def get_all_dependencies(self):
        pass
    def __call__(self, data):
        data[self.name] = self.function(self, data)

class FieldValidator(object):
    pass

class ValidateSpecies(FieldValidator):
    def __init__(self, species):
        FieldValidator.__init__(self)
        if not isinstance(species, types.ListType):
            species = [species]
        self.species = species
    def __call__(self, data):
        for s in species:
            if not data.has_field(s): return False
        return True

def ValidateProperty(FieldValidator):
    def __init__(self, prop):
        FieldValidator.__init__(self)
        if not isinstance(prop, types.ListType):
            prop = [prop]
        self.prop = prop
    def __call__(self, data):
        for p in prop:
            if not hasattr(data,p): return False
        return True

# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

_speciesList = ["HI","HII","Electron",
               "HeI","HeII","HeIII",
               "H2I","H2II","HM",
               "DI","DII","HDI"]
def _SpeciesFractiospecies_fractionn(field, data):
    sp = field.name.split("_")[0]
    return data[sp]/data["Density"]
for species in speciesList:
    add_field("%s_Fraction" % species,
             function=_species_fraction,
             validator=ValidateSpecies("%s_Density" % species))

def _SoundSpeed(field, data):
    return data.convert("x-velocity") * ( \
           data.params["Gamma"]*data["Pressure"] / \
           data["Density"] )**(1.0/2.0)
add_field("SoundSpeed")

def _MachNumber(field, data):
    """M{|v|/t_sound}"""
    return data["VelocityMagnitude"] / data["SoundSpeed"]
add_field("MachNumber")

def _VelocityMagnitude(field, data):
    """M{|v|}"""
    return (data.convert("x-velocity") * ( \
            data["x-velocity"]**2.0 + \
            data["y-velocity"]**2.0 + \
            data["z-velocity"]**2.0 )**(1.0/2.0))
add_field("VelocityMagnitude", takeLog=False)

def _Pressure(field, data):
    """M{(Gamma-1.0)*rho*E}"""
    return (data.params["Gamma"] - 1.0) * \
           data["Density"] * data["Gas_Energy"]
add_field("Pressure", validator=ValidateSpecies("Gas_Energy"))

def _Entropy(field, data):
    return data["Density"]**(-2./3.) * \
           data["Temperature"]
add_field("Entropy")

def _DynamicalTime(field, data):
    """
    The formulation for the dynamical time is:
    M{sqrt(3pi/(16*G*rho))} or M{sqrt(3pi/(16G))*rho^-(1/2)}
    Note that we return in our natural units already
    """
    G = data.params["GravitationalConstant"]
    t_dyn_coeff = (3*pi/(16*G))**0.5 * data.convert("Time")
    return data["Density"]**(-1./2.) * t_dyn_coeff
add_field("DynamicalTime")

def _NumberDensity(field, data):
    # We can assume that we at least have Density
    # We should actually be guaranteeing the presence of a .shape attribute,
    # but I am not currently implementing that
    fieldData = na.zeros(data["Density"].shape,
                         dtype = data["Density"].dtype)
    if data.params["MultiSpecies"] == 0:
        fieldData += data["Density"] * data.params["mu"]
    if data.params["MultiSpecies"] > 0:
        fieldData += data["HI_Density"] / 1.0
        fieldData += data["HII_Density"] / 1.0
        fieldData += data["HeI_Density"] / 4.0
        fieldData += data["HeII_Density"] / 4.0
        fieldData += data["HeIII_Density"] / 4.0
        fieldData += data["Electron_Density"] / 1.0
    if data.params["MultiSpecies"] > 1:
        fieldData += data["HM_Density"] / 1.0
        fieldData += data["H2I_Density"] / 2.0
        fieldData += data["H2II_Density"] / 2.0
    if data.params["MultiSpecies"] > 2:
        fieldData += data["DI_Density"] / 2.0
        fieldData += data["DII_Density"] / 2.0
        fieldData += data["HDI_Density"] / 3.0
    return fieldData * data.convert("Density") / mh
add_field("NumberDensity")



if __name__ == "__main__":
    k = fieldInfo.keys()
    k.sort()
    for f in k:
        e = FieldDetector()
        fieldInfo[f](e)
        print f + ":", ", ".join(e.requested)