"""
Derived fields and support for existing fields included in here.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

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

import types
import numpy as na
import inspect
import copy

# All our math stuff here:
try:
    import scipy.signal
except ImportError:
    pass

from math import pi

from yt.funcs import *

mh = 1.67e-24 # g
me = 9.11e-28 # g
sigma_thompson = 6.65e-25 # cm^2
clight = 3.0e10 # cm/s
kboltz = 1.38e-16 # erg K^-1
G = 6.67e-8   # cm^3 g^-1 s^-2
Msun2g = 1.989e33
MJ_constant = (((5*kboltz)/(G*mh))**(1.5)) * (3/(4*pi))**(0.5) / Msun2g


class FieldInfoContainer: # We are all Borg.
    _shared_state = {}
    def __new__(cls, *args, **kwargs):
        self = object.__new__(cls, *args, **kwargs)
        self.__dict__ = cls._shared_state
        return self
    def __init__(self):
        self.__field_list = {}
    def __getitem__(self, key):
        if not self.__field_list.has_key(key): # Now we check it out, to see if we need to do anything to it
            if key.endswith("Code"):
                old_field = self.__field_list[key[:-4]]
                new_field = copy.copy(old_field)
                new_field.convert_function = \
                            lambda a: 1.0/old_field._convert_function(a)
            elif key.endswith("Abs"):
                old_field = self.__field_list[key[:-4]]
                new_field = copy.copy(old_field)
                new_field._function = \
                    lambda a,b: na.abs(old_field._function(a,b))
            else:
                raise KeyError
            self.__field_list[key] = new_field
        return self.__field_list[key]

fieldInfo = {}

class ValidationException(Exception):
    pass

class NeedsGridType(ValidationException):
    def __init__(self, ghost_zones = 0, fields=None):
        self.ghost_zones = ghost_zones
        self.fields = fields

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
    Level = 1
    NumberOfParticles = 0
    def __init__(self, nd = 16, pf = None):
        self.nd = nd
        self.ActiveDimensions = [nd,nd,nd]
        self.LeftEdge = [0.0,0.0,0.0]
        self.RightEdge = [1.0,1.0,1.0]
        self.dx = self.dy = self.dz = na.array([1.0])
        self.fields = []
        if pf is None:
            pf = defaultdict(lambda: 1)
        self.pf = pf
        self.requested = []
        self.requested_parameters = []
        defaultdict.__init__(self, lambda: na.ones((nd,nd,nd)))
    def __missing__(self, item):
        if fieldInfo.has_key(item) and \
            fieldInfo[item]._function.func_name != '<lambda>':
            try:
                vv = fieldInfo[item](self)
            except NeedsGridType, exc:
                ngz = exc.ghost_zones
                nfd = FieldDetector(self.nd+ngz*2)
                vv = fieldInfo[item](nfd)[ngz:-ngz,ngz:-ngz,ngz:-ngz]
                for i in vv.requested:
                    if i not in self.requested: self.requested.append(i)
                for i in vv.requested_parameters:
                    if i not in self.requested_parameters: self.requested_parameters.append(i)
            if vv is not None:
                self[item] = vv
                return self[item]
        self.requested.append(item)
        return defaultdict.__missing__(self, item)
    def get_field_parameter(self, param):
        self.requested_parameters.append(param)
        if param in ['bulk_velocity','center','height_vector']:
            return na.array([0,0,0])
        else:
            return 0.0
    _spatial = True
    _num_ghost_zones = 0
    id = 1
    def has_field_parameter(self, param): return True
    def convert(self, item): return 1

class DerivedField:
    def __init__(self, name, function,
                 convert_function = None,
                 units = "", projected_units = "",
                 take_log = True, validators = None,
                 particle_type = False, vector_field=False,
                 display_field = True, not_in_all=False,
                 projection_conversion = "cm"):
        """
        This is the base class used to describe a cell-by-cell derived field.

        :param name: is the name of the field.
        :param function: is a function handle that defines the field
        :param convert_function: must convert to CGS, if it needs to be done
        :param units: is a mathtext-formatted string that describes the field
        :param projected_units: if we display a projection, what should the units be?
        :param take_log: describes whether the field should be logged
        :param validators: is a list of :class:`FieldValidator` objects
        :param particle_type: is this field based on particles?
        :param vector_field: describes the dimensionality of the field
        :param display_field: governs its appearance in the dropdowns in reason
        :param not_in_all: is used for baryon fields from the data that are not in
                           all the grids
        :param projection_conversion: which unit should we multiply by in a
                                      projection?
        """
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
        self.particle_type = particle_type
        self.vector_field = vector_field
        self.projection_conversion = projection_conversion
        self.display_field = display_field
        self.not_in_all = not_in_all
    def check_available(self, data):
        for validator in self.validators:
            validator(data)
        # If we don't get an exception, we're good to go
        return True
    def get_dependencies(self, *args, **kwargs):
        e = FieldDetector(*args, **kwargs)
        if self._function.func_name == '<lambda>':
            e.requested.append(self.name)
        else:
            self(e)
        return e
    def get_units(self):
        return self._units
    def get_projected_units(self):
        return self._projected_units
    def __call__(self, data):
        ii = self.check_available(data)
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
            if not data.has_field_parameter(p):
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
        if isinstance(data, FieldDetector): return True
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
        for p in self.prop:
            if not hasattr(data,p):
                doesnt_have.append(p)
        if len(doesnt_have) > 0:
            raise NeedsProperty(doesnt_have)
        return True

class ValidateSpatial(FieldValidator):
    def __init__(self, ghost_zones = 0, fields=None):
        FieldValidator.__init__(self)
        self.ghost_zones = ghost_zones
        self.fields = fields
    def __call__(self, data):
        # When we say spatial information, we really mean
        # that it has a three-dimensional data structure
        if isinstance(data, FieldDetector): return True
        if not data._spatial:
            raise NeedsGridType(self.ghost_zones,self.fields)
        if self.ghost_zones == data._num_ghost_zones:
            return True
        raise NeedsGridType(self.ghost_zones)

# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

def _dx(field, data):
    return data.dx
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dx
add_field('dx', display_field=False, validators=[ValidateSpatial(0)])

def _dy(field, data):
    return data.dy
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dy
add_field('dy', display_field=False, validators=[ValidateSpatial(0)])

def _dz(field, data):
    return data.dz
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dz
add_field('dz', display_field=False, validators=[ValidateSpatial(0)])

def _coordX(field, data):
    dim = data.ActiveDimensions[0]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[0]).reshape(dim,1,1)
            +0.5) * data['dx'] + data.LeftEdge[0]
add_field('x', function=_coordX, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordY(field, data):
    dim = data.ActiveDimensions[1]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[1]).reshape(1,dim,1)
            +0.5) * data['dy'] + data.LeftEdge[1]
add_field('y', function=_coordY, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordZ(field, data):
    dim = data.ActiveDimensions[2]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[2]).reshape(1,1,dim)
            +0.5) * data['dz'] + data.LeftEdge[2]
add_field('z', function=_coordZ, display_field=False,
          validators=[ValidateSpatial(0)])


_speciesList = ["HI","HII","Electron",
               "HeI","HeII","HeIII",
               "H2I","H2II","HM",
               "DI","DII","HDI","Metal"]
def _SpeciesFraction(field, data):
    sp = field.name.split("_")[0] + "_Density"
    return data[sp]/data["Density"]
for species in _speciesList:
    add_field("%s_Fraction" % species,
             function=_SpeciesFraction,
             validators=ValidateDataField("%s_Density" % species))

def _Metallicity(field, data):
    return data["Metal_Fraction"] / 0.0204
add_field("Metallicity", units=r"Z_{\rm{Solar}}",
          validators=ValidateDataField("Metal_Density"),
          projection_conversion="1")


def _GridLevel(field, data):
    return na.ones(data["Density"].shape)*(data.Level)
add_field("GridLevel", validators=[#ValidateProperty('Level'),
                                   ValidateSpatial(0)])

def _GridIndices(field, data):
    return na.ones(data["Density"].shape)*(data.id-1)
add_field("GridIndices", validators=[#ValidateProperty('id'),
                                     ValidateSpatial(0)], take_log=False)

def _OnesOverDx(field, data):
    return na.ones(data["Density"].shape,
                   dtype=data["Density"].dtype)/data['dx']
add_field("OnesOverDx", display_field=False)

def _Ones(field, data):
    return na.ones(data.ActiveDimensions, dtype='float64')
add_field("Ones", validators=[ValidateSpatial(0)], projection_conversion="1",
          display_field = False)
add_field("CellsPerBin", function=_Ones, validators=[ValidateSpatial(0)],
          display_field = False)

def _SoundSpeed(field, data):
    return ( data.pf["Gamma"]*data["Pressure"] / \
             data["Density"] )**(1.0/2.0)
add_field("SoundSpeed", units=r"\rm{cm}/\rm{s}")

def particle_func(p_field):
    def _Particles(field, data):
        if not data.NumberOfParticles > 0:
            return na.array([], dtype='float64')
        try:
            return data._read_data(p_field).astype('float64')
        except data._read_exception:
            pass
        # This is bad.  But it's the best idea I have right now.
        return data._read_data(p_field.replace("_"," ")).astype('float64')
    return _Particles
for pf in ["index", "type", "mass"] + \
          ["velocity_%s" % ax for ax in 'xyz'] + \
          ["position_%s" % ax for ax in 'xyz']:
    pfunc = particle_func("particle_%s" % (pf))
    add_field("particle_%s" % pf, function=pfunc,
              validators = [ValidateSpatial(0)],
              particle_type=True)
for pf in ["creation_time", "dynamical_time", "metallicity_fraction"]:
    pfunc = particle_func(pf)
    add_field(pf, function=pfunc,
              validators = [ValidateSpatial(0),
                            ValidateDataField(pf)],
              particle_type=True)
add_field("particle mass", function=particle_func("particle_mass"),
          validators=[ValidateSpatial(0)], particle_type=True)

add_field("Dark matter density", function=lambda a,b: None,
          validators=[ValidateDataField("Dark matter density"),
                      ValidateSpatial(0)],
          not_in_all = True)

def _ParticleMass(field, data):
    particles = data["particle_mass"].astype('float64') * \
                just_one(data["CellVolumeCode"].ravel())
    # Note that we mandate grid-type here, so this is okay
    return particles
def _convertParticleMass(data):
    return data.convert("Density")*(data.convert("cm")**3.0)
def _convertParticleMassMsun(data):
    return data.convert("Density")*((data.convert("cm")**3.0)/1.989e33)
add_field("ParticleMass", validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMass)
add_field("ParticleMassMsun",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun)

def _MachNumber(field, data):
    """M{|v|/t_sound}"""
    return data["VelocityMagnitude"] / data["SoundSpeed"]
add_field("MachNumber")

def _CourantTimeStep(field, data):
    t1 = data['dx'] / (
        data["SoundSpeed"] + \
        abs(data["x-velocity"]))
    t2 = data['dy'] / (
        data["SoundSpeed"] + \
        abs(data["y-velocity"]))
    t3 = data['dz'] / (
        data["SoundSpeed"] + \
        abs(data["z-velocity"]))
    return na.minimum(na.minimum(t1,t2),t3)
def _convertCourantTimeStep(data):
    # SoundSpeed and z-velocity are in cm/s, dx is in code
    return data.convert("cm")
add_field("CourantTimeStep", convert_function=_convertCourantTimeStep,
          units=r"$\rm{s}$")

def _VelocityMagnitude(field, data):
    """M{|v|}"""
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    return ( (data["x-velocity"]-bulk_velocity[0])**2.0 + \
             (data["y-velocity"]-bulk_velocity[1])**2.0 + \
             (data["z-velocity"]-bulk_velocity[2])**2.0 )**(1.0/2.0)
add_field("VelocityMagnitude", take_log=False, units=r"\rm{cm}/\rm{s}")

def _TangentialOverVelocityMagnitude(field, data):
    return na.abs(data["TangentialVelocity"])/na.abs(data["VelocityMagnitude"])
add_field("TangentialOverVelocityMagnitude", take_log=False)

def _TangentialVelocity(field, data):
    return na.sqrt(data["VelocityMagnitude"]**2.0
                 - data["RadialVelocity"]**2.0)
add_field("TangentialVelocity", take_log=False, units=r"\rm{cm}/\rm{s}")

def _Pressure(field, data):
    """M{(Gamma-1.0)*rho*E}"""
    return (data.pf["Gamma"] - 1.0) * \
           data["Density"] * data["ThermalEnergy"]
add_field("Pressure", units=r"\rm{dyne}/\rm{cm}^{2}")

def _ThermalEnergy(field, data):
    if data.pf["HydroMethod"] == 2:
        return data["Total_Energy"]
    else:
        if data.pf["DualEnergyFormalism"]:
            return data["Gas_Energy"]
        else:
            return data["Total_Energy"] - 0.5*(
                   data["x-velocity"]**2.0
                 + data["y-velocity"]**2.0
                 + data["z-velocity"]**2.0 )
add_field("ThermalEnergy", units=r"\rm{ergs}/\rm{g}")

def _Entropy(field, data):
    return data["Density"]**(-2./3.) * \
           data["Temperature"]
add_field("Entropy", units="WhoKnows")

def _Height(field, data):
    # We take the dot product of the radius vector with the height-vector
    center = data.get_field_parameter("center")
    r_vec = na.array([data["x"] - center[0],
                      data["y"] - center[1],
                      data["z"] - center[2]])
    h_vec = na.array(data.get_field_parameter("height_vector"))
    h_vec = h_vec / na.sqrt(h_vec[0]**2.0+
                            h_vec[1]**2.0+
                            h_vec[2]**2.0)
    height = r_vec[0,:] * h_vec[0] \
           + r_vec[1,:] * h_vec[1] \
           + r_vec[2,:] * h_vec[2]
    return na.abs(height)
def _convertHeight(data):
    return data.convert("cm")
def _convertHeightAU(data):
    return data.convert("au")
add_field("Height", convert_function=_convertHeight,
          validators=[ValidateParameter("height_vector")],
          units=r"cm", display_field=False)
add_field("HeightAU", function=_Height,
          convert_function=_convertHeightAU,
          validators=[ValidateParameter("height_vector")],
          units=r"AU", display_field=False)

def _DiskAngle(field, data):
    # We make both r_vec and h_vec into unit vectors
    center = data.get_field_parameter("center")
    r_vec = na.array([data["x"] - center[0],
                      data["y"] - center[1],
                      data["z"] - center[2]])
    r_vec = r_vec/na.sqrt((r_vec**2.0).sum(axis=0))
    h_vec = na.array(data.get_field_parameter("height_vector"))
    dp = r_vec[0,:] * h_vec[0] \
       + r_vec[1,:] * h_vec[1] \
       + r_vec[2,:] * h_vec[2]
    return na.arccos(dp)
add_field("DiskAngle", take_log=False,
          validators=[ValidateParameter("height_vector"),
                      ValidateParameter("center")],
          display_field=False)

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
add_field("DynamicalTime", units=r"\rm{s}",
          convert_function=_ConvertDynamicalTime)

def _NumberDensity(field, data):
    # We can assume that we at least have Density
    # We should actually be guaranteeing the presence of a .shape attribute,
    # but I am not currently implementing that
    fieldData = na.zeros(data["Density"].shape,
                         dtype = data["Density"].dtype)
    if data.pf["MultiSpecies"] == 0:
        if data.has_field_parameter("mu"):
            mu = data.get_field_parameter("mu")
        else:
            mu = 0.6
        fieldData += data["Density"] * mu
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
add_field("NumberDensity", units=r"\rm{cm}^{-3}",
          convert_function=_ConvertNumberDensity)

def JeansMassMsun(field,data):
    return (MJ_constant * 
            ((data["Temperature"]/data["MeanMolecularWeight"])**(1.5)) *
            (data["Density"]**(-0.5)))
add_field("JeansMassMsun",function=JeansMassMsun,units=r"\rm{Msun}")


def _CellMass(field, data):
    return data["Density"] * data["CellVolume"]
def _convertCellMassMsun(data):
    return 5.027854e-34 # g^-1
add_field("CellMass", units=r"\rm{g}")
add_field("CellMassMsun", units=r"M_{\odot}",
          function=_CellMass,
          convert_function=_convertCellMassMsun)

def _CellVolume(field, data):
    if data['dx'].size == 1:
        try:
            return data['dx']*data['dy']*data['dx']*\
                na.ones(data.ActiveDimensions, dtype='float64')
        except AttributeError:
            return data['dx']*data['dy']*data['dx']
    return data["dx"]*data["dy"]*data["dz"]
def _ConvertCellVolumeMpc(data):
    return data.convert("mpc")**3.0
def _ConvertCellVolumeCGS(data):
    return data.convert("cm")**3.0
add_field("CellVolumeCode", units=r"\rm{BoxVolume}^3",
          function=_CellVolume)
add_field("CellVolumeMpc", units=r"\rm{Mpc}^3",
          function=_CellVolume,
          convert_function=_ConvertCellVolumeMpc)
add_field("CellVolume", units=r"\rm{cm}^3",
          function=_CellVolume,
          convert_function=_ConvertCellVolumeCGS)

def _XRayEmissivity(field, data):
    return ((data["Density"].astype('float64')**2.0) \
            *data["Temperature"]**0.5)
def _convertXRayEmissivity(data):
    return 2.168e60
add_field("XRayEmissivity",
          convert_function=_convertXRayEmissivity,
          projection_conversion="1")

def _SZKinetic(field, data):
    vel_axis = data.get_field_parameter('axis')
    if vel_axis > 2:
        raise NeedsParameter(['axis'])
    vel = data["%s-velocity" % ({0:'x',1:'y',2:'z'}[vel_axis])]
    return (vel*data["Density"])
def _convertSZKinetic(data):
    return 0.88*((sigma_thompson/mh)/clight)
add_field("SZKinetic",
          convert_function=_convertSZKinetic,
          validators=[ValidateParameter('axis')])

def _SZY(field, data):
    return (data["Density"]*data["Temperature"])
def _convertSZY(data):
    conv = (0.88/mh) * (kboltz)/(me * clight*clight) * sigma_thompson
    return conv
add_field("SZY", convert_function=_convertSZY)

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

def _DivV(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    div_x = (data["x-velocity"][sl_right,1:-1,1:-1] -
             data["x-velocity"][sl_left,1:-1,1:-1]) \
          / (div_fac*data["dx"].flat[0])
    div_y = (data["y-velocity"][1:-1,sl_right,1:-1] -
             data["y-velocity"][1:-1,sl_left,1:-1]) \
          / (div_fac*data["dy"].flat[0])
    div_z = (data["z-velocity"][1:-1,1:-1,sl_right] -
             data["z-velocity"][1:-1,1:-1,sl_left]) \
          / (div_fac*data["dz"].flat[0])
    new_field = na.zeros(data["x-velocity"].shape)
    new_field[1:-1,1:-1,1:-1] = div_x+div_y+div_z
    return na.abs(new_field)
def _convertDivV(data):
    return data.convert("cm")**-1.0
add_field("DivV", validators=[ValidateSpatial(1,
            ["x-velocity","y-velocity","z-velocity"])],
          units=r"\rm{s}^{-1}",
          convert_function=_convertDivV)

def _Contours(field, data):
    return na.ones(data["Density"].shape)*-1
add_field("Contours", validators=[ValidateSpatial(0)], take_log=False,
          display_field=False)
add_field("tempContours", function=_Contours, validators=[ValidateSpatial(0)],
          take_log=False, display_field=False)

def _SpecificAngularMomentum(field, data):
    """
    Calculate the angular velocity.  Returns a vector for each cell.
    """
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    xv = data["x-velocity"] - bv[0]
    yv = data["y-velocity"] - bv[1]
    zv = data["z-velocity"] - bv[2]
    center = data.get_field_parameter('center')
    coords = na.array([data['x'],data['y'],data['z']], dtype='float64')
    new_shape = tuple([3] + [1]*(len(coords.shape)-1))
    r_vec = coords - na.reshape(center,new_shape)
    v_vec = na.array([xv,yv,zv], dtype='float64')
    return na.cross(r_vec, v_vec, axis=0)
def _convertSpecificAngularMomentum(data):
    return data.convert("cm")
add_field("SpecificAngularMomentum",
          convert_function=_convertSpecificAngularMomentum, vector_field=True,
          units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])
def _convertSpecificAngularMomentumKMSMPC(data):
    return data.convert("mpc")/1e5
add_field("SpecificAngularMomentumKMSMPC",
          function=_SpecificAngularMomentum,
          convert_function=_convertSpecificAngularMomentumKMSMPC, vector_field=True,
          units=r"\rm{km}\rm{Mpc}/\rm{s}", validators=[ValidateParameter('center')])

def _ParticleRadius(field, data):
    center = data.get_field_parameter("center")
    radius = na.sqrt((data["particle_position_x"] - center[0])**2.0 +
                     (data["particle_position_y"] - center[1])**2.0 +
                     (data["particle_position_z"] - center[2])**2.0)
    return radius
def _Radius(field, data):
    center = data.get_field_parameter("center")
    radius = na.sqrt((data["x"] - center[0])**2.0 +
                     (data["y"] - center[1])**2.0 +
                     (data["z"] - center[2])**2.0)
    return radius
def _ConvertRadiusCGS(data):
    return data.convert("cm")
add_field("ParticleRadius", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units=r"\rm{cm}",
          particle_type = True)
add_field("Radius", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units=r"\rm{cm}")

def _ConvertRadiusMpc(data):
    return data.convert("mpc")
add_field("RadiusMpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units=r"\rm{Mpc}")
add_field("ParticleRadiusMpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units=r"\rm{Mpc}",
          particle_type=True)

def _ConvertRadiuskpc(data):
    return data.convert("kpc")
add_field("ParticleRadiuskpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}",
          particle_type=True)
add_field("Radiuskpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}")

def _ConvertRadiuskpch(data):
    return data.convert("kpch")
add_field("ParticleRadiuskpch", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}/\rm{h}",
          particle_type=True)
add_field("Radiuskpch", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}/\rm{h}")

def _ConvertRadiuspc(data):
    return data.convert("pc")
add_field("ParticleRadiuspc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units=r"\rm{pc}",
          particle_type=True)
add_field("Radiuspc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units=r"\rm{pc}")

def _ConvertRadiusAU(data):
    return data.convert("au")
add_field("ParticleRadiusAU", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units=r"\rm{AU}",
          particle_type=True)
add_field("RadiusAU", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units=r"\rm{AU}")

add_field("ParticleRadiusCode", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          particle_type=True)
add_field("RadiusCode", function=_Radius,
          validators=[ValidateParameter("center")])

def _RadialVelocity(field, data):
    center = data.get_field_parameter("center")
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    new_field = ( (data['x']-center[0])*(data["x-velocity"]-bulk_velocity[0])
                + (data['y']-center[1])*(data["y-velocity"]-bulk_velocity[1])
                + (data['z']-center[2])*(data["z-velocity"]-bulk_velocity[2])
                )/data["RadiusCode"]
    return new_field
def _RadialVelocityABS(field, data):
    return na.abs(_RadialVelocity(field, data))
def _ConvertRadialVelocityKMS(data):
    return 1e-5
add_field("RadialVelocity", function=_RadialVelocity,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center"),
                      ValidateParameter("bulk_velocity")])
add_field("RadialVelocityABS", function=_RadialVelocityABS,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center"),
                      ValidateParameter("bulk_velocity")])
add_field("RadialVelocityKMS", function=_RadialVelocity,
          convert_function=_ConvertRadialVelocityKMS, units=r"\rm{km}/\rm{s}",
          validators=[ValidateParameter("center"),
                      ValidateParameter("bulk_velocity")])

def _CuttingPlaneVelocityX(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    v_vec = na.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                - bulk_velocity[...,na.newaxis]
    return na.dot(x_vec, v_vec)
add_field("CuttingPlaneVelocityX", 
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units=r"\rm{km}/\rm{s}")
def _CuttingPlaneVelocityY(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    v_vec = na.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                - bulk_velocity[...,na.newaxis]
    return na.dot(y_vec, v_vec)
add_field("CuttingPlaneVelocityY", 
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units=r"\rm{km}/\rm{s}")

def _MeanMolecularWeight(field,data):
    return (data["Density"] / (mh *data["NumberDensity"]))
add_field("MeanMolecularWeight",function=_MeanMolecularWeight,units=r"")

def _JeansMassMsun(field,data):
    MJ_constant = (((5*kboltz)/(G*mh))**(1.5)) * \
    (3/(4*3.1415926535897931))**(0.5) / 1.989e33

    return (MJ_constant *
            ((data["Temperature"]/data["MeanMolecularWeight"])**(1.5)) *
            (data["Density"]**(-0.5)))
add_field("JeansMassMsun",function=_JeansMassMsun,
          units=r"\rm{M_{\odot}}")

# Now we add all the fields that we want to control, but we give a null function
# This is every Enzo field we can think of.  This will be installation-dependent,

_enzo_fields = ["Density","Temperature","Gas_Energy","Total_Energy",
                "x-velocity","y-velocity","z-velocity"]
_enzo_fields += [ "%s_Density" % sp for sp in _speciesList ]
for field in _enzo_fields:
    add_field(field, function=lambda a, b: None, take_log=True,
              validators=[ValidateDataField(field)], units=r"\rm{g}/\rm{cm}^3")
fieldInfo["x-velocity"].projection_conversion='1'
fieldInfo["y-velocity"].projection_conversion='1'
fieldInfo["z-velocity"].projection_conversion='1'

# Now we override

def _convertDensity(data):
    return data.convert("Density")
for field in ["Density"] + [ "%s_Density" % sp for sp in _speciesList ]:
    fieldInfo[field]._units = r"\rm{g}/\rm{cm}^3"
    fieldInfo[field]._projected_units = r"\rm{g}/\rm{cm}^2"
    fieldInfo[field]._convert_function=_convertDensity

add_field("Dark_Matter_Density", function=lambda a,b: None,
          convert_function=_convertDensity,
          validators=[ValidateDataField("Dark_Matter_Density"),
                      ValidateSpatial(0)],
          not_in_all = True)

def _convertEnergy(data):
    return data.convert("x-velocity")**2.0
fieldInfo["Gas_Energy"]._units = r"\rm{ergs}/\rm{g}"
fieldInfo["Gas_Energy"]._convert_function = _convertEnergy
fieldInfo["Total_Energy"]._units = r"\rm{ergs}/\rm{g}"
fieldInfo["Total_Energy"]._convert_function = _convertEnergy
fieldInfo["Temperature"]._units = r"\rm{K}"

def _convertVelocity(data):
    return data.convert("x-velocity")
for ax in ['x','y','z']:
    f = fieldInfo["%s-velocity" % ax]
    f._units = r"\rm{cm}/\rm{s}"
    f._convert_function = _convertVelocity
    f.take_log = False

fieldInfo["Temperature"].units = r"K"

if __name__ == "__main__":
    k = fieldInfo.keys()
    k.sort()
    for f in k:
        e = FieldDetector()
        fieldInfo[f](e)
        print f + ":", ", ".join(e.requested)
