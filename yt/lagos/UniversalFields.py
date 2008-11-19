"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.

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

from math import pi

from yt.funcs import *
from FieldInfoContainer import *

try:
    import cic_deposit
except ImportError:
    pass

mh = 1.67e-24 # g
me = 9.11e-28 # g
sigma_thompson = 6.65e-25 # cm^2
clight = 3.0e10 # cm/s
kboltz = 1.38e-16 # erg K^-1
G = 6.67e-8   # cm^3 g^-1 s^-2



# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

def _dx(field, data):
    return data.dx
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dx
add_field('dx', function=_dx, display_field=False,
          validators=[ValidateSpatial(0)])

def _dy(field, data):
    return data.dy
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dy
add_field('dy', function=_dy, display_field=False,
          validators=[ValidateSpatial(0)])

def _dz(field, data):
    return data.dz
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dz
add_field('dz', function=_dz,
          display_field=False, validators=[ValidateSpatial(0)])

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

def _GridLevel(field, data):
    return na.ones(data["Density"].shape)*(data.Level)
add_field("GridLevel", function=_GridLevel,
          validators=[#ValidateProperty('Level'),
                      ValidateSpatial(0)])

def _GridIndices(field, data):
    return na.ones(data["Density"].shape)*(data.id-data._id_offset)
add_field("GridIndices", function=_GridIndices,
          validators=[#ValidateProperty('id'),
                      ValidateSpatial(0)], take_log=False)

def _OnesOverDx(field, data):
    return na.ones(data["Density"].shape,
                   dtype=data["Density"].dtype)/data['dx']
add_field("OnesOverDx", function=_OnesOverDx,
          display_field=False)

def _Ones(field, data):
    return na.ones(data.ActiveDimensions, dtype='float64')
add_field("Ones", function=_Ones,
          validators=[ValidateSpatial(0)],
          projection_conversion="unitary",
          display_field = False)
add_field("CellsPerBin", function=_Ones, validators=[ValidateSpatial(0)],
          display_field = False)

def _SoundSpeed(field, data):
    return ( data.pf["Gamma"]*data["Pressure"] / \
             data["Density"] )**(1.0/2.0)
add_field("SoundSpeed", function=_SoundSpeed,
          units=r"\rm{cm}/\rm{s}")

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
add_field("ParticleMass",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMass)
add_field("ParticleMassMsun",
          function=_ParticleMass, validators=[ValidateSpatial(0)],
          particle_type=True, convert_function=_convertParticleMassMsun)

def _MachNumber(field, data):
    """M{|v|/t_sound}"""
    return data["VelocityMagnitude"] / data["SoundSpeed"]
add_field("MachNumber", function=_MachNumber)

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
add_field("CourantTimeStep", function=_CourantTimeStep,
          convert_function=_convertCourantTimeStep,
          units=r"$\rm{s}$")

def _VelocityMagnitude(field, data):
    """M{|v|}"""
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    return ( (data["x-velocity"]-bulk_velocity[0])**2.0 + \
             (data["y-velocity"]-bulk_velocity[1])**2.0 + \
             (data["z-velocity"]-bulk_velocity[2])**2.0 )**(1.0/2.0)
add_field("VelocityMagnitude", function=_VelocityMagnitude,
          take_log=False, units=r"\rm{cm}/\rm{s}")

def _TangentialOverVelocityMagnitude(field, data):
    return na.abs(data["TangentialVelocity"])/na.abs(data["VelocityMagnitude"])
add_field("TangentialOverVelocityMagnitude",
          function=_TangentialOverVelocityMagnitude,
          take_log=False)

def _TangentialVelocity(field, data):
    return na.sqrt(data["VelocityMagnitude"]**2.0
                 - data["RadialVelocity"]**2.0)
add_field("TangentialVelocity", 
          function=_TangentialVelocity,
          take_log=False, units=r"\rm{cm}/\rm{s}")

def _Pressure(field, data):
    """M{(Gamma-1.0)*rho*E}"""
    return (data.pf["Gamma"] - 1.0) * \
           data["Density"] * data["ThermalEnergy"]
add_field("Pressure", function=_Pressure, units=r"\rm{dyne}/\rm{cm}^{2}")

def _Entropy(field, data):
    return data["Density"]**(-2./3.) * \
           data["Temperature"]
add_field("Entropy", function=_Entropy, units="WhoKnows")

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
add_field("Height", function=_Height,
          convert_function=_convertHeight,
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
add_field("DiskAngle", function=_DiskAngle,
          take_log=False,
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
add_field("DynamicalTime", function=_DynamicalTime,
           units=r"\rm{s}",
          convert_function=_ConvertDynamicalTime)

def JeansMassMsun(field,data):
    return (MJ_constant * 
            ((data["Temperature"]/data["MeanMolecularWeight"])**(1.5)) *
            (data["Density"]**(-0.5)))
add_field("JeansMassMsun",function=JeansMassMsun,units=r"\rm{Msun}")

def _CellMass(field, data):
    return data["Density"] * data["CellVolume"]
def _convertCellMassMsun(data):
    return 5.027854e-34 # g^-1
add_field("CellMass", function=_CellMass, units=r"\rm{g}")
add_field("CellMassMsun", units=r"M_{\odot}",
          function=_CellMass,
          convert_function=_convertCellMassMsun)

def _CellMassCode(field, data):
    return data["Density"] * data["CellVolumeCode"]
def _convertCellMassCode(data):
    return 1.0/data.convert("Density")
add_field("CellMassCode", 
          function=_CellMassCode,
          convert_function=_convertCellMassCode)

def _TotalMass(field,data):
    return (data["Density"]+data["particle_density"]) * data["CellVolume"]
add_field("TotalMassMsun", units=r"M_{\odot}",
          function=_TotalMass,
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
add_field("XRayEmissivity", function=_XRayEmissivity,
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
add_field("SZKinetic", function=_SZKinetic,
          convert_function=_convertSZKinetic,
          validators=[ValidateParameter('axis')])

def _SZY(field, data):
    return (data["Density"]*data["Temperature"])
def _convertSZY(data):
    conv = (0.88/mh) * (kboltz)/(me * clight*clight) * sigma_thompson
    return conv
add_field("SZY", function=_SZY, convert_function=_convertSZY)

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
add_field("AveragedDensity",
          function=_AveragedDensity,
          validators=[ValidateSpatial(1)])

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
add_field("DivV", function=_DivV,
            validators=[ValidateSpatial(1,
            ["x-velocity","y-velocity","z-velocity"])],
          units=r"\rm{s}^{-1}",
          convert_function=_convertDivV)

def _Contours(field, data):
    return na.ones(data["Density"].shape)*-1
add_field("Contours", validators=[ValidateSpatial(0)], take_log=False,
          display_field=False, function=_Contours)
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
          function=_SpecificAngularMomentum,
          convert_function=_convertSpecificAngularMomentum, vector_field=True,
          units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])
def _convertSpecificAngularMomentumKMSMPC(data):
    return data.convert("mpc")/1e5
add_field("SpecificAngularMomentumKMSMPC",
          function=_SpecificAngularMomentum,
          convert_function=_convertSpecificAngularMomentumKMSMPC, vector_field=True,
          units=r"\rm{km}\rm{Mpc}/\rm{s}", validators=[ValidateParameter('center')])
def _AngularMomentum(field, data):
    return data["CellMass"] * data["SpecificAngularMomentum"]
add_field("AngularMomentum", function=_AngularMomentum,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", vector_field=True)
def _AngularMomentumMSUNKMSMPC(field, data):
    return data["CellMassMsun"] * data["SpecificAngularMomentumKMSMPC"]
add_field("AngularMomentumMSUNKMSMPC", function=_AngularMomentum,
          units=r"M_{\odot}\rm{km}\rm{Mpc}/\rm{s}", vector_field=True)

def _ParticleSpecificAngularMomentum(field, data):
    """
    Calculate the angular of a particle velocity.  Returns a vector for each
    particle.
    """
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    xv = data["particle_velocity_x"] - bv[0]
    yv = data["particle_velocity_y"] - bv[1]
    zv = data["particle_velocity_z"] - bv[2]
    center = data.get_field_parameter('center')
    coords = na.array([data['particle_position_x'],
                       data['particle_position_y'],
                       data['particle_position_z']], dtype='float64')
    new_shape = tuple([3] + [1]*(len(coords.shape)-1))
    r_vec = coords - na.reshape(center,new_shape)
    v_vec = na.array([xv,yv,zv], dtype='float64')
    return na.cross(r_vec, v_vec, axis=0)
add_field("ParticleSpecificAngularMomentum",
          function=_ParticleSpecificAngularMomentum,
          convert_function=_convertSpecificAngularMomentum, vector_field=True,
          units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])
def _convertSpecificAngularMomentumKMSMPC(data):
    return data.convert("mpc")/1e5
add_field("ParticleSpecificAngularMomentumKMSMPC",
          function=_ParticleSpecificAngularMomentum,
          convert_function=_convertSpecificAngularMomentumKMSMPC, vector_field=True,
          units=r"\rm{km}\rm{Mpc}/\rm{s}", validators=[ValidateParameter('center')])
def _ParticleAngularMomentum(field, data):
    return data["ParticleMass"] * data["ParticleSpecificAngularMomentum"]
add_field("ParticleAngularMomentum", 
          function=_ParticleAngularMomentum, units=r"\rm{g}\/\rm{cm}^2/\rm{s}",
          particle_type=True)
def _ParticleAngularMomentumMSUNKMSMPC(field, data):
    return data["ParticleMass"] * data["ParticleSpecificAngularMomentumKMSMPC"]
add_field("ParticleAngularMomentumMSUNKMSMPC",
          function=_ParticleAngularMomentumMSUNKMSMPC,
          units=r"M_{\odot}\rm{km}\rm{Mpc}/\rm{s}",
          particle_type=True)


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
          function=_CuttingPlaneVelocityX,
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
          function=_CuttingPlaneVelocityY,
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
