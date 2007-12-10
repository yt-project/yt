import types
from collections import defaultdict
import numpy as na
import inspect

# All our math stuff here:
import scipy.signal
from math import pi

from yt.funcs import *

mh = 1.67e-24
fieldInfo = {}

class ValidationException(Exception):
    pass

class NeedsGridType(ValidationException):
    def __init__(self, ghost_zones = 0):
        self.ghost_zones = ghost_zones

class NeedsDataField(ValidationException):
    def __init__(self, missing_fields):
        self.missing_fields = missing_fields

class NeedsProperty(ValidationException):
    def __init__(self, missing_properties):
        self.missing_properties = missing_properties

class NeedsParameter(ValidationException):
    def __init__(self, missing_parameters):
        self.missing_parameters = missing_parameters

def add_field(name, function = None, **kwargs):
    if function == None:
        if kwargs.has_key("function"):
            function = kwargs.pop("function")
        else:
            # This will fail if it does not exist,
            # which is our desired behavior
            function = eval("_%s" % name)
    fieldInfo[name] = DerivedField(
        name, function, **kwargs)

class FieldDetector(defaultdict):
    pf = defaultdict(lambda: 1)
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
                 convert_function = None,
                 units = "", projected_units = "",
                 take_log = True, validators = None):
        self.name = name
        self._function = function
        if validators:
            self.validators = ensure_list(validators)
        else:
            self.validators = []
        self.take_log = take_log
        self._units = units
        self._projected_units = projected_units
        if not convert_function:
            convert_function = lambda a: 1.0
        self._convert_function = convert_function
    def check_available(self, data):
        for validator in self.validators:
            validator(data)
        # If we don't get an exception, we're good to go
        return True
    def get_dependencies(self):
        e = FieldDetector()
        self(e)
        return e.requested
    def get_units(self):
        return self._units
    def get_projected_units(self):
        return self._projected_units
    def __call__(self, data):
        original_fields = data.fields[:] # Copy
        dd = self._function(self, data)
        dd *= self._convert_function(data)
        for field_name in data.fields:
            if field_name not in original_fields:
                del data[field_name]
        return dd
    def get_source(self):
        return inspect.getsource(self._function)

class FieldValidator(object):
    pass

class ValidateParameter(FieldValidator):
    def __init__(self, parameters):
        FieldValidator.__init__(self)
        self.parameters = ensure_list(parameters)
    def __call__(self, data):
        doesnt_have = []
        for p in self.parameters:
            if not data.field_parameters.has_key(p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsParameter(doesnt_have)
        return True

class ValidateDataField(FieldValidator):
    def __init__(self, field):
        FieldValidator.__init__(self)
        self.fields = ensure_list(field)
    def __call__(self, data):
        doesnt_have = []
        for f in self.fields:
            if f not in data.hierarchy.field_list:
                doesnt_have.append(f)
        if len(doesnt_have) > 0:
            raise NeedsDataField(doesnt_have)
        return True

class ValidateProperty(FieldValidator):
    def __init__(self, prop):
        FieldValidator.__init__(self)
        self.prop = ensure_list(prop)
    def __call__(self, data):
        doesnt_have = []
        for p in prop:
            if not hasattr(data,p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsProperty(doesnt_have)
        return True

class ValidateSpatial(FieldValidator):
    def __init__(self, ghost_zones = 0):
        FieldValidator.__init__(self)
        self.ghost_zones = ghost_zones
    def __call__(self, data):
        # When we say spatial information, we really mean
        # that it has a three-dimensional data structure
        if self.ghost_zones == 0:
            if data._spatial:
                return True
        if self.ghost_zones == data._num_ghost_zones:
            return True
        raise NeedsGridType(self.ghost_zones)

# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

def _coordX(field, data):
    dim = data.ActiveDimensions[0]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[0]).reshape(dim,1,1)
            +0.5) * data['dx'] + data.LeftEdge[0]
add_field('x', function=_coordX,
          validators=[ValidateSpatial(0)])

def _coordY(field, data):
    dim = data.ActiveDimensions[1]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[1]).reshape(1,dim,1)
            +0.5) * data['dy'] + data.LeftEdge[1]
add_field('y', function=_coordY,
          validators=[ValidateSpatial(0)])

def _coordZ(field, data):
    dim = data.ActiveDimensions[2]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[2]).reshape(1,1,dim)
            +0.5) * data['dz'] + data.LeftEdge[2]
add_field('z', function=_coordZ,
          validators=[ValidateSpatial(0)])


_speciesList = ["HI","HII","Electron",
               "HeI","HeII","HeIII",
               "H2I","H2II","HM",
               "DI","DII","HDI"]
def _SpeciesFraction(field, data):
    sp = field.name.split("_")[0]
    return data[sp]/data["Density"]
for species in _speciesList:
    add_field("%s_Fraction" % species,
             function=_SpeciesFraction,
             validators=ValidateDataField("%s_Density" % species))

def _Ones(field, data):
    return na.ones(data["Density"].shape,
                   dtype=data["Density"].dtype)
add_field("Ones")
add_field("CellsPerBin", function=_Ones)

def _SoundSpeed(field, data):
    return ( data.pf["Gamma"]*data["Pressure"] / \
             data["Density"] )**(1.0/2.0)
add_field("SoundSpeed", units=r"$\rm{cm}/\rm{s}$")

def particle_func(p_field):
    def _Particles(field, data):
        try:
            particles = data["particle_%s" % p_field]
        except:
            particles = na.array([], dtype=data["Density"].dtype)
        return particles
    return _Particles
for pf in ["index","mass","type"] + \
          ["velocity_%s" % ax for ax in 'xyz'] + \
          ["position_%s" % ax for ax in 'xyz']:
    pfunc = particle_func(pf)
    add_field("particle_%s" % pf, function=pfunc,
              validators = [ValidateSpatial(0)])

def _MachNumber(field, data):
    """M{|v|/t_sound}"""
    return data["VelocityMagnitude"] / data["SoundSpeed"]
add_field("MachNumber")

def _RadialVelocity(field, data):
    pass

def _VelocityMagnitude(field, data):
    """M{|v|}"""
    return ( data["x-velocity"]**2.0 + \
             data["y-velocity"]**2.0 + \
             data["z-velocity"]**2.0 )**(1.0/2.0)
add_field("VelocityMagnitude", take_log=False)

def _Pressure(field, data):
    """M{(Gamma-1.0)*rho*E}"""
    return (data.pf["Gamma"] - 1.0) * \
           data["Density"] * data["Gas_Energy"]
add_field("Pressure", #validators=ValidateDataField("Gas_Energy"),
          units=r"$\rm{dyne}/\rm{cm}^{2}$")

def _Entropy(field, data):
    return data["Density"]**(-2./3.) * \
           data["Temperature"]
add_field("Entropy", units="WhoKnows")

def _DynamicalTime(field, data):
    """
    The formulation for the dynamical time is:
    M{sqrt(3pi/(16*G*rho))} or M{sqrt(3pi/(16G))*rho^-(1/2)}
    Note that we return in our natural units already
    """
    return data["Density"]**(-1./2.)
def _ConvertDynamicalTime(data):
    G = data.pf["GravitationalConstant"]
    t_dyn_coeff = (3*pi/(16*G))**0.5 \
                * data.convert("Time")
    return t_dyn_coeff
add_field("DynamicalTime", units=r"$\rm{s}$",
          convert_function=_ConvertDynamicalTime)

def _NumberDensity(field, data):
    # We can assume that we at least have Density
    # We should actually be guaranteeing the presence of a .shape attribute,
    # but I am not currently implementing that
    fieldData = na.zeros(data["Density"].shape,
                         dtype = data["Density"].dtype)
    if data.pf["MultiSpecies"] == 0:
        fieldData += data["Density"] * data.get_field_parameter("mu", 0.6)
    if data.pf["MultiSpecies"] > 0:
        fieldData += data["HI_Density"] / 1.0
        fieldData += data["HII_Density"] / 1.0
        fieldData += data["HeI_Density"] / 4.0
        fieldData += data["HeII_Density"] / 4.0
        fieldData += data["HeIII_Density"] / 4.0
        fieldData += data["Electron_Density"] / 1.0
    if data.pf["MultiSpecies"] > 1:
        fieldData += data["HM_Density"] / 1.0
        fieldData += data["H2I_Density"] / 2.0
        fieldData += data["H2II_Density"] / 2.0
    if data.pf["MultiSpecies"] > 2:
        fieldData += data["DI_Density"] / 2.0
        fieldData += data["DII_Density"] / 2.0
        fieldData += data["HDI_Density"] / 3.0
    return fieldData
def _ConvertNumberDensity(data):
    return 1.0/mh
add_field("NumberDensity", units=r"$\rm{cm}^{-3}$",
          convert_function=_ConvertNumberDensity)

def _CellMass(field, data):
    return data["Density"] * data["CellVolume"]
add_field("CellMass", units=r"$\rm{g}$")

def _CellVolume(field, data):
    return data["dx"]*data["dy"]*data["dz"]
def _ConvertCellVolume(data):
    return data.convert("cm")**3.0
add_field("CellVolume", units=r"$\rm{cm}^3$",
          convert_function=_ConvertCellVolume)

def __gauss_kern(size):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    x, y, z = na.mgrid[-size:size+1, -size:size+1, -size:size+1]
    g = na.exp(-(x**2/float(size)+y**2/float(size)+z**2/float(size)))
    return g / g.sum()

def __blur_image(im, n):
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = __gauss_kern(n)
    improc = scipy.signal.convolve(im,g, mode='same')
    return(improc)

def _SmoothedDensity(field, data):
    return __blur_image(data["Density"], 1)
add_field("SmoothedDensity", validators=[ValidateSpatial(2)])

def _AveragedDensity(field, data):
    nx, ny, nz = data["Density"].shape
    new_field = na.zeros((nx-2,ny-2,nz-2), dtype='float64')
    weight_field = na.zeros((nx-2,ny-2,nz-2), dtype='float64')
    i_i, j_i, k_i = na.mgrid[0:3,0:3,0:3]
    for i,j,k in zip(i_i.ravel(),j_i.ravel(),k_i.ravel()):
        sl = [slice(i,nx-(2-i)),slice(j,ny-(2-j)),slice(k,nz-(2-k))]
        new_field += data["Density"][sl] * data["CellMass"][sl]
        weight_field += data["CellMass"][sl]
    # Now some fancy footwork
    new_field2 = na.zeros((nx,ny,nz))
    new_field2[1:-1,1:-1,1:-1] = new_field/weight_field
    return new_field2
add_field("AveragedDensity", validators=[ValidateSpatial(1)])

def _Radius(field, data):
    center = data.get_field_parameter("center")
    radius = na.sqrt((data["x"] - center[0])**2.0 +
                     (data["y"] - center[1])**2.0 +
                     (data["z"] - center[2])**2.0)
    return radius
def _ConvertRadius(data):
    return data.convert("cm")
add_field("Radius", validators=[ValidateParameter("center")],
          convert_function = _ConvertRadius, units=r"$\rm{cm}$")
add_field("RadiusCode", function=_Radius,
          validators=[ValidateParameter("center")])

def _RadialVelocity(field, data):
    center = data.get_field_parameter("center")
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    new_field = ( (data['x']-center[0])*data["x-velocity"]
                + (data['y']-center[1])*data["y-velocity"]
                + (data['z']-center[2])*data["z-velocity"])/data["RadiusCode"]
    return new_field
def _ConvertRadialVelocity(data):
    return (data.convert("x-velocity") / 1e5)
def _ConvertRadialVelocityCGS(data):
    return (data.convert("x-velocity"))
add_field("RadialVelocity", function=_RadialVelocity,
          convert_function=_ConvertRadialVelocity, units=r"$\rm{km}/\rm{s}$",
          validators=[ValidateParameter("center"),
                      ValidateParameter("bulk_velocity")])
add_field("RadialVelocityCGS", function=_RadialVelocity,
          convert_function=_ConvertRadialVelocityCGS, units=r"$\rm{cm}/\rm{s}$")

# Now we add all the fields that we want to control, but we give a null function
# This is every Enzo field we can think of.  This will be installation-dependent,

_enzo_fields = ["Density","Temperature","Gas_Energy","Total_Energy",
                "x-velocity","y-velocity","z-velocity"]
_enzo_fields += [ "%s_Density" % sp for sp in _speciesList ]
def _returnCodeField(real_field):
    def _fieldFunction(field, data):
        return data[real_field]
    return _fieldFunction
for field in _enzo_fields:
    add_field(field, function=lambda a, b: None, take_log=True,
              validators=[ValidateDataField(field)], units=r"$\rm{g}/\rm{cm}^3$")
    add_field("_Code%s" % field, function=_returnCodeField(field),
              take_log=True, validators=[ValidateDataField(field)])

# Now we override

def _convertDensity(data):
    return data.convert("Density")
for field in ["Density"] + [ "%s_Density" % sp for sp in _speciesList ]:
    fieldInfo[field]._units = r"$\rm{g}/\rm{cm}^3$"
    fieldInfo[field]._projected_units = r"$\rm{g}/\rm{cm}^2$"
    fieldInfo[field]._convert_function=_convertDensity

def _convertEnergy(data):
    return data.convert("x-velocity")**2.0
fieldInfo["Gas_Energy"].units = r"$\rm{ergs}/\rm{g}$"
fieldInfo["Gas_Energy"]._convert_function = _convertEnergy
fieldInfo["Total_Energy"].units = r"$\rm{ergs}/\rm{g}$"
fieldInfo["Total_Energy"]._convert_function = _convertEnergy

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = fieldInfo["%s-velocity" % ax]
    f.units = r"$\rm{km}/\rm{s}$"
    f.convert_function = _convertVelocity
    f.take_log = False

fieldInfo["Temperature"].units = r"$K$"

if __name__ == "__main__":
    k = fieldInfo.keys()
    k.sort()
    for f in k:
        e = FieldDetector()
        fieldInfo[f](e)
        print f + ":", ", ".join(e.requested)