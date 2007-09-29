fieldInfo = {}

def AddField(name, function = None, **kwargs):
    if function == None:
        if not kwargs.has_key(function):
            function = kwargs.pop("function")
        else:
            # This will fail if it does not exist,
            # which is our desired behavior
            function = eval("_%s" % name)
    fieldInfo[name] = DerivedField(
        name, function, **kwargs
    )

class FieldDetector(DefaultDict):
    pass

class DerivedField:
    def __init__(self, name, function,
                 units = "", projected_units = "",
                 takeLog = True, validator = None):
        self.function = function
        self.validator = validator
        pass
    def GetValues(self, object):
        pass
    def CheckAvailable(self, object):
        if self.validator != None:
            self.validator(self, object)
    def GetMyDependencies(self):
        pass
    def GetAllDependencies(self):
        pass
    def __call__(self, object):
        object[self.name] = self.function(self, object)

def ValidateSpecies(species):
    if not istype(species, types.ListType):
        species = [species]
    def CheckSpecies(object):
        for s in species:
            if not object.has_field(s): return False
        return True

def ValidatePropery(prop):
    if not istype(prop, types.ListType):
        prop = [prop]
    def CheckProperty(object):
        for p in prop:
            if not hasattr(object,p): return False
        return True

def _H2I_Fraction(field, object):
    return object["H2I_Density"]/object["Density"]
AddField("H2I_Fraction", validator=ValidateSpecies("H2I_Density"))

