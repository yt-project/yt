"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

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

from yt.funcs import *

from yt.utilities.lib import CICDeposit_3, obtain_rvec
from yt.utilities.cosmology import Cosmology
from field_info_container import \
    add_field, \
    ValidateDataField, \
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter

from yt.utilities.physical_constants import \
     mh, \
     me, \
     sigma_thompson, \
     clight, \
     kboltz, \
     G, \
     rho_crit_now, \
     speed_of_light_cgs
     
# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

def _dx(field, data):
    return data.dds[0]
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dds[0]
add_field('dx', function=_dx, display_field=False,
          validators=[ValidateSpatial(0)])

def _dy(field, data):
    return data.dds[1]
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dds[1]
add_field('dy', function=_dy, display_field=False,
          validators=[ValidateSpatial(0)])

def _dz(field, data):
    return data.dds[2]
    return na.ones(data.ActiveDimensions, dtype='float64') * data.dds[2]
add_field('dz', function=_dz,
          display_field=False, validators=[ValidateSpatial(0)])

def _coordX(field, data):
    dim = data.ActiveDimensions[0]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[0])[:,None,None]
            +0.5) * data['dx'] + data.LeftEdge[0]
add_field('x', function=_coordX, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordY(field, data):
    dim = data.ActiveDimensions[1]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[1])[None,:,None]
            +0.5) * data['dy'] + data.LeftEdge[1]
add_field('y', function=_coordY, display_field=False,
          validators=[ValidateSpatial(0)])

def _coordZ(field, data):
    dim = data.ActiveDimensions[2]
    return (na.ones(data.ActiveDimensions, dtype='float64')
                   * na.arange(data.ActiveDimensions[2])[None,None,:]
            +0.5) * data['dz'] + data.LeftEdge[2]
add_field('z', function=_coordZ, display_field=False,
          validators=[ValidateSpatial(0)])

def _GridLevel(field, data):
    return na.ones(data.ActiveDimensions)*(data.Level)
add_field("GridLevel", function=_GridLevel,
          validators=[ValidateGridType(),
                      ValidateSpatial(0)])

def _GridIndices(field, data):
    return na.ones(data["Ones"].shape)*(data.id-data._id_offset)
add_field("GridIndices", function=_GridIndices,
          validators=[ValidateGridType(),
                      ValidateSpatial(0)], take_log=False)

def _OnesOverDx(field, data):
    return na.ones(data["Ones"].shape,
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
    if data.pf["EOSType"] == 1:
        return na.ones(data["Density"].shape, dtype='float64') * \
                data.pf["EOSSoundSpeed"]
    return ( data.pf["Gamma"]*data["Pressure"] / \
             data["Density"] )**(1.0/2.0)
add_field("SoundSpeed", function=_SoundSpeed,
          units=r"\rm{cm}/\rm{s}")

def _RadialMachNumber(field, data):
    """M{|v|/t_sound}"""
    return na.abs(data["RadialVelocity"]) / data["SoundSpeed"]
add_field("RadialMachNumber", function=_RadialMachNumber)

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

def _ParticleVelocityMagnitude(field, data):
    """M{|v|}"""
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    return ( (data["particle_velocity_x"]-bulk_velocity[0])**2.0 + \
             (data["particle_velocity_y"]-bulk_velocity[1])**2.0 + \
             (data["particle_velocity_z"]-bulk_velocity[2])**2.0 )**(1.0/2.0)
add_field("ParticleVelocityMagnitude", function=_ParticleVelocityMagnitude,
          particle_type=True, 
          take_log=False, units=r"\rm{cm}/\rm{s}")

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
    if data.has_field_parameter("mu"):
        mw = mh*data.get_field_parameter("mu")
    else :
        mw = mh
    return kboltz * data["Temperature"] / \
           ((data["Density"]/mw)**(data.pf["Gamma"] - 1.0))
add_field("Entropy", units=r"\rm{ergs}\ \rm{cm}^{3\gamma-3}",
          function=_Entropy)

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
    return (3.0*na.pi/(16*G*data["Density"]))**(1./2.)
add_field("DynamicalTime", function=_DynamicalTime,
           units=r"\rm{s}")

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
    return (data["Density"]+data["Dark_Matter_Density"]) * data["CellVolume"]
add_field("TotalMass", function=_TotalMass, units=r"\rm{g}")
add_field("TotalMassMsun", units=r"M_{\odot}",
          function=_TotalMass,
          convert_function=_convertCellMassMsun)

def _StarMass(field,data):
    return data["star_density"] * data["CellVolume"]
add_field("StarMassMsun", units=r"M_{\odot}",
          function=_StarMass,
          convert_function=_convertCellMassMsun)

def _Matter_Density(field,data):
    return (data['Density'] + data['Dark_Matter_Density'])
add_field("Matter_Density",function=_Matter_Density,units=r"\rm{g}/\rm{cm^3}")

def _ComovingDensity(field, data):
    ef = (1.0 + data.pf.current_redshift)**3.0
    return data["Density"]/ef
add_field("ComovingDensity", function=_ComovingDensity, units=r"\rm{g}/\rm{cm}^3")

# This is rho_total / rho_cr(z).
def _Convert_Overdensity(data):
    return 1 / (rho_crit_now * data.pf.hubble_constant**2 * 
                (1+data.pf.current_redshift)**3)
add_field("Overdensity",function=_Matter_Density,
          convert_function=_Convert_Overdensity, units=r"")

# This is (rho_total - <rho_total>) / <rho_total>.
def _DensityPerturbation(field, data):
    rho_bar = rho_crit_now * data.pf.omega_matter * \
        data.pf.hubble_constant**2 * \
        (1.0 + data.pf.current_redshift)**3
    return ((data['Matter_Density'] - rho_bar) / rho_bar)
add_field("DensityPerturbation",function=_DensityPerturbation,units=r"")

# This is rho_b / <rho_b>.
def _Baryon_Overdensity(field, data):
    if data.pf.has_key('omega_baryon_now'):
        omega_baryon_now = data.pf['omega_baryon_now']
    else:
        omega_baryon_now = 0.0441
    return data['Density'] / (omega_baryon_now * rho_crit_now * 
                              (data.pf['CosmologyHubbleConstantNow']**2) * 
                              ((1+data.pf['CosmologyCurrentRedshift'])**3))
add_field("Baryon_Overdensity", function=_Baryon_Overdensity, 
          units=r"")

# Weak lensing convergence.
# Eqn 4 of Metzler, White, & Loken (2001, ApJ, 547, 560).
def _convertConvergence(data):
    if not data.pf.parameters.has_key('cosmology_calculator'):
        data.pf.parameters['cosmology_calculator'] = Cosmology(
            HubbleConstantNow=(100.*data.pf.hubble_constant),
            OmegaMatterNow=data.pf.omega_matter, OmegaLambdaNow=data.pf.omega_lambda)
    # observer to lens
    DL = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
        data.pf.parameters['observer_redshift'], data.pf.current_redshift)
    # observer to source
    DS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
        data.pf.parameters['observer_redshift'], data.pf.parameters['lensing_source_redshift'])
    # lens to source
    DLS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
        data.pf.current_redshift, data.pf.parameters['lensing_source_redshift'])
    return (((DL * DLS) / DS) * (1.5e14 * data.pf.omega_matter * 
                                (data.pf.hubble_constant / speed_of_light_cgs)**2 *
                                (1 + data.pf.current_redshift)))
add_field("WeakLensingConvergence", function=_DensityPerturbation, 
          convert_function=_convertConvergence, 
          projection_conversion='mpccm')

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
          validators=[ValidateSpatial(1, ["Density"])])

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
    ds = div_fac * data['dx'].flat[0]
    f  = data["x-velocity"][sl_right,1:-1,1:-1]/ds
    f -= data["x-velocity"][sl_left ,1:-1,1:-1]/ds
    if data.pf.dimensionality > 1:
        ds = div_fac * data['dy'].flat[0]
        f += data["y-velocity"][1:-1,sl_right,1:-1]/ds
        f -= data["y-velocity"][1:-1,sl_left ,1:-1]/ds
    if data.pf.dimensionality > 2:
        ds = div_fac * data['dz'].flat[0]
        f += data["z-velocity"][1:-1,1:-1,sl_right]/ds
        f -= data["z-velocity"][1:-1,1:-1,sl_left ]/ds
    new_field = na.zeros(data["x-velocity"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = f
    return new_field
def _convertDivV(data):
    return data.convert("cm")**-1.0
add_field("DivV", function=_DivV,
            validators=[ValidateSpatial(1,
            ["x-velocity","y-velocity","z-velocity"])],
          units=r"\rm{s}^{-1}", take_log=False,
          convert_function=_convertDivV)

def _AbsDivV(field, data):
    return na.abs(data['DivV'])
add_field("AbsDivV", function=_AbsDivV,
          units=r"\rm{s}^{-1}")

def _Contours(field, data):
    return -na.ones_like(data["Ones"])
add_field("Contours", validators=[ValidateSpatial(0)], take_log=False,
          display_field=False, function=_Contours)
add_field("tempContours", function=_Contours,
          validators=[ValidateSpatial(0), ValidateGridType()],
          take_log=False, display_field=False)

def obtain_velocities(data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    xv = data["x-velocity"] - bv[0]
    yv = data["y-velocity"] - bv[1]
    zv = data["z-velocity"] - bv[2]
    return xv, yv, zv

def _convertSpecificAngularMomentum(data):
    return data.convert("cm")
def _convertSpecificAngularMomentumKMSMPC(data):
    return data.convert("mpc")/1e5

def _SpecificAngularMomentumX(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return yv*rv[2,:] - zv*rv[1,:]
def _SpecificAngularMomentumY(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return -(xv*rv[2,:] - zv*rv[0,:])
def _SpecificAngularMomentumZ(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return xv*rv[1,:] - yv*rv[0,:]
for ax in 'XYZ':
    n = "SpecificAngularMomentum%s" % ax
    add_field(n, function=eval("_%s" % n),
              convert_function=_convertSpecificAngularMomentum,
              units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter("center")])

def _AngularMomentumX(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumX"]
add_field("AngularMomentumX", function=_AngularMomentumX,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", vector_field=False,
         validators=[ValidateParameter('center')])
def _AngularMomentumY(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumY"]
add_field("AngularMomentumY", function=_AngularMomentumY,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", vector_field=False,
         validators=[ValidateParameter('center')])
def _AngularMomentumZ(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumZ"]
add_field("AngularMomentumZ", function=_AngularMomentumZ,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", vector_field=False,
         validators=[ValidateParameter('center')])

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
#add_field("ParticleSpecificAngularMomentum",
#          function=_ParticleSpecificAngularMomentum, particle_type=True,
#          convert_function=_convertSpecificAngularMomentum, vector_field=True,
#          units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])
def _convertSpecificAngularMomentumKMSMPC(data):
    return data.convert("mpc")/1e5
#add_field("ParticleSpecificAngularMomentumKMSMPC",
#          function=_ParticleSpecificAngularMomentum, particle_type=True,
#          convert_function=_convertSpecificAngularMomentumKMSMPC, vector_field=True,
#          units=r"\rm{km}\rm{Mpc}/\rm{s}", validators=[ValidateParameter('center')])

def _ParticleSpecificAngularMomentumX(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    center = data.get_field_parameter('center')
    y = data["particle_position_y"] - center[1]
    z = data["particle_position_z"] - center[2]
    yv = data["particle_velocity_y"] - bv[1]
    zv = data["particle_velocity_z"] - bv[2]
    return yv*z - zv*y
def _ParticleSpecificAngularMomentumY(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    center = data.get_field_parameter('center')
    x = data["particle_position_x"] - center[0]
    z = data["particle_position_z"] - center[2]
    xv = data["particle_velocity_x"] - bv[0]
    zv = data["particle_velocity_z"] - bv[2]
    return -(xv*z - zv*x)
def _ParticleSpecificAngularMomentumZ(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = na.zeros(3, dtype='float64')
    center = data.get_field_parameter('center')
    x = data["particle_position_x"] - center[0]
    y = data["particle_position_y"] - center[1]
    xv = data["particle_velocity_x"] - bv[0]
    yv = data["particle_velocity_y"] - bv[1]
    return xv*y - yv*x
for ax in 'XYZ':
    n = "ParticleSpecificAngularMomentum%s" % ax
    add_field(n, function=eval("_%s" % n), particle_type=True,
              convert_function=_convertSpecificAngularMomentum,
              units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter("center")])
    add_field(n + "KMSMPC", function=eval("_%s" % n), particle_type=True,
              convert_function=_convertSpecificAngularMomentumKMSMPC,
              units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter("center")])

def _ParticleAngularMomentum(field, data):
    return data["ParticleMass"] * data["ParticleSpecificAngularMomentum"]
#add_field("ParticleAngularMomentum",
#          function=_ParticleAngularMomentum, units=r"\rm{g}\/\rm{cm}^2/\rm{s}",
#          particle_type=True, validators=[ValidateParameter('center')])
def _ParticleAngularMomentumMSUNKMSMPC(field, data):
    return data["ParticleMass"] * data["ParticleSpecificAngularMomentumKMSMPC"]
#add_field("ParticleAngularMomentumMSUNKMSMPC",
#          function=_ParticleAngularMomentumMSUNKMSMPC,
#          units=r"M_{\odot}\rm{km}\rm{Mpc}/\rm{s}",
#          particle_type=True, validators=[ValidateParameter('center')])

def _ParticleAngularMomentumX(field, data):
    return data["CellMass"] * data["ParticleSpecificAngularMomentumX"]
add_field("ParticleAngularMomentumX", function=_ParticleAngularMomentumX,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
         validators=[ValidateParameter('center')])
def _ParticleAngularMomentumY(field, data):
    return data["CellMass"] * data["ParticleSpecificAngularMomentumY"]
add_field("ParticleAngularMomentumY", function=_ParticleAngularMomentumY,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
         validators=[ValidateParameter('center')])
def _ParticleAngularMomentumZ(field, data):
    return data["CellMass"] * data["ParticleSpecificAngularMomentumZ"]
add_field("ParticleAngularMomentumZ", function=_ParticleAngularMomentumZ,
         units=r"\rm{g}\/\rm{cm}^2/\rm{s}", particle_type=True,
         validators=[ValidateParameter('center')])


def _ParticleRadius(field, data):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = na.zeros(data["particle_position_x"].shape, dtype='float64')
    for i, ax in enumerate('xyz'):
        r = na.abs(data["particle_position_%s" % ax] - center[i])
        radius += na.minimum(r, na.abs(DW[i]-r))**2.0
    na.sqrt(radius, radius)
    return radius
def _Radius(field, data):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = na.zeros(data["x"].shape, dtype='float64')
    for i, ax in enumerate('xyz'):
        r = na.abs(data[ax] - center[i])
        radius += na.minimum(r, na.abs(DW[i]-r))**2.0
    na.sqrt(radius, radius)
    return radius
def _ConvertRadiusCGS(data):
    return data.convert("cm")
add_field("ParticleRadius", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units=r"\rm{cm}",
          particle_type = True,
          display_name = "Particle Radius")
add_field("Radius", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units=r"\rm{cm}")

def _ConvertRadiusMpc(data):
    return data.convert("mpc")
add_field("RadiusMpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units=r"\rm{Mpc}",
          display_name = "Radius")
add_field("ParticleRadiusMpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units=r"\rm{Mpc}",
          particle_type=True,
          display_name = "Particle Radius")

def _ConvertRadiuskpc(data):
    return data.convert("kpc")
add_field("ParticleRadiuskpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuskpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}",
          display_name = "Radius")

def _ConvertRadiuskpch(data):
    return data.convert("kpch")
add_field("ParticleRadiuskpch", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpch, units=r"\rm{kpc}/\rm{h}",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuskpch", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units=r"\rm{kpc}/\rm{h}",
          display_name = "Radius")

def _ConvertRadiuspc(data):
    return data.convert("pc")
add_field("ParticleRadiuspc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units=r"\rm{pc}",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuspc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units=r"\rm{pc}",
          display_name="Radius")

def _ConvertRadiusAU(data):
    return data.convert("au")
add_field("ParticleRadiusAU", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units=r"\rm{AU}",
          particle_type=True,
          display_name = "Particle Radius")
add_field("RadiusAU", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units=r"\rm{AU}",
          display_name = "Radius")

add_field("ParticleRadiusCode", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          particle_type=True,
          display_name = "Particle Radius (code)")
add_field("RadiusCode", function=_Radius,
          validators=[ValidateParameter("center")],
          display_name = "Radius (code)")

def _RadialVelocity(field, data):
    center = data.get_field_parameter("center")
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = na.zeros(3)
    new_field = ( (data['x']-center[0])*(data["x-velocity"]-bulk_velocity[0])
                + (data['y']-center[1])*(data["y-velocity"]-bulk_velocity[1])
                + (data['z']-center[2])*(data["z-velocity"]-bulk_velocity[2])
                )/data["RadiusCode"]
    if na.any(na.isnan(new_field)): # to fix center = point
        new_field[na.isnan(new_field)] = 0.0
    return new_field
def _RadialVelocityABS(field, data):
    return na.abs(_RadialVelocity(field, data))
def _ConvertRadialVelocityKMS(data):
    return 1e-5
add_field("RadialVelocity", function=_RadialVelocity,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("RadialVelocityABS", function=_RadialVelocityABS,
          units=r"\rm{cm}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("RadialVelocityKMS", function=_RadialVelocity,
          convert_function=_ConvertRadialVelocityKMS, units=r"\rm{km}/\rm{s}",
          validators=[ValidateParameter("center")])
add_field("RadialVelocityKMSABS", function=_RadialVelocityABS,
          convert_function=_ConvertRadialVelocityKMS, units=r"\rm{km}/\rm{s}",
          validators=[ValidateParameter("center")])

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

def _convertDensity(data):
    return data.convert("Density")
def _pdensity(field, data):
    blank = na.zeros(data.ActiveDimensions, dtype='float32')
    if data.NumberOfParticles == 0: return blank
    CICDeposit_3(data["particle_position_x"].astype(na.float64),
                 data["particle_position_y"].astype(na.float64),
                 data["particle_position_z"].astype(na.float64),
                 data["particle_mass"].astype(na.float32),
                 na.int64(data.NumberOfParticles),
                 blank, na.array(data.LeftEdge).astype(na.float64),
                 na.array(data.ActiveDimensions).astype(na.int32),
                 na.float64(data['dx']))
    return blank
add_field("particle_density", function=_pdensity,
          validators=[ValidateGridType()], convert_function=_convertDensity,
          display_name=r"\mathrm{Particle}\/\mathrm{Density})")

def _MagneticEnergy(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    MagneticEnergy field to deal with non-unitary \mu_0.
    """
    return (data["Bx"]**2 + data["By"]**2 + data["Bz"]**2)/2.
add_field("MagneticEnergy",function=_MagneticEnergy,
          units=r"",
          validators = [ValidateDataField("Bx"),
                        ValidateDataField("By"),
                        ValidateDataField("Bz")])

def _VorticitySquared(field, data):
    mylog.debug("Generating vorticity on %s", data)
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["x-velocity"].shape)
    dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
             data["z-velocity"][1:-1,sl_left,1:-1]) \
             / (div_fac*data["dy"].flat[0])
    dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
             data["y-velocity"][1:-1,1:-1,sl_left]) \
             / (div_fac*data["dz"].flat[0])
    new_field[1:-1,1:-1,1:-1] += (dvzdy - dvydz)**2.0
    del dvzdy, dvydz
    dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
             data["x-velocity"][1:-1,1:-1,sl_left]) \
             / (div_fac*data["dz"].flat[0])
    dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
             data["z-velocity"][sl_left,1:-1,1:-1]) \
             / (div_fac*data["dx"].flat[0])
    new_field[1:-1,1:-1,1:-1] += (dvxdz - dvzdx)**2.0
    del dvxdz, dvzdx
    dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
             data["y-velocity"][sl_left,1:-1,1:-1]) \
             / (div_fac*data["dx"].flat[0])
    dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
             data["x-velocity"][1:-1,sl_left,1:-1]) \
             / (div_fac*data["dy"].flat[0])
    new_field[1:-1,1:-1,1:-1] += (dvydx - dvxdy)**2.0
    del dvydx, dvxdy
    new_field = na.abs(new_field)
    return new_field
def _convertVorticitySquared(data):
    return data.convert("cm")**-2.0
add_field("VorticitySquared", function=_VorticitySquared,
          validators=[ValidateSpatial(1,
              ["x-velocity","y-velocity","z-velocity"])],
          units=r"\rm{s}^{-2}",
          convert_function=_convertVorticitySquared)

def _gradPressureX(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Pressure"].shape, dtype='float64')
    ds = div_fac * data['dx'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Pressure"][sl_right,1:-1,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Pressure"][sl_left ,1:-1,1:-1]/ds
    return new_field
def _gradPressureY(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Pressure"].shape, dtype='float64')
    ds = div_fac * data['dy'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Pressure"][1:-1,sl_right,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Pressure"][1:-1,sl_left ,1:-1]/ds
    return new_field
def _gradPressureZ(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Pressure"].shape, dtype='float64')
    ds = div_fac * data['dz'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Pressure"][1:-1,1:-1,sl_right]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Pressure"][1:-1,1:-1,sl_left ]/ds
    return new_field
def _convertgradPressure(data):
    return 1.0/data.convert("cm")
for ax in 'XYZ':
    n = "gradPressure%s" % ax
    add_field(n, function=eval("_%s" % n),
              convert_function=_convertgradPressure,
              validators=[ValidateSpatial(1, ["Pressure"])],
              units=r"\rm{dyne}/\rm{cm}^{3}")

def _gradPressureMagnitude(field, data):
    return na.sqrt(data["gradPressureX"]**2 +
                   data["gradPressureY"]**2 +
                   data["gradPressureZ"]**2)
add_field("gradPressureMagnitude", function=_gradPressureMagnitude,
          validators=[ValidateSpatial(1, ["Pressure"])],
          units=r"\rm{dyne}/\rm{cm}^{3}")

def _gradDensityX(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Density"].shape, dtype='float64')
    ds = div_fac * data['dx'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Density"][sl_right,1:-1,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Density"][sl_left ,1:-1,1:-1]/ds
    return new_field
def _gradDensityY(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Density"].shape, dtype='float64')
    ds = div_fac * data['dy'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Density"][1:-1,sl_right,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Density"][1:-1,sl_left ,1:-1]/ds
    return new_field
def _gradDensityZ(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["Density"].shape, dtype='float64')
    ds = div_fac * data['dz'].flat[0]
    new_field[1:-1,1:-1,1:-1]  = data["Density"][1:-1,1:-1,sl_right]/ds
    new_field[1:-1,1:-1,1:-1] -= data["Density"][1:-1,1:-1,sl_left ]/ds
    return new_field
def _convertgradDensity(data):
    return 1.0/data.convert("cm")
for ax in 'XYZ':
    n = "gradDensity%s" % ax
    add_field(n, function=eval("_%s" % n),
              convert_function=_convertgradDensity,
              validators=[ValidateSpatial(1, ["Density"])],
              units=r"\rm{g}/\rm{cm}^{4}")

def _gradDensityMagnitude(field, data):
    return na.sqrt(data["gradDensityX"]**2 +
                   data["gradDensityY"]**2 +
                   data["gradDensityZ"]**2)
add_field("gradDensityMagnitude", function=_gradDensityMagnitude,
          validators=[ValidateSpatial(1, ["Density"])],
          units=r"\rm{g}/\rm{cm}^{4}")

def _BaroclinicVorticityX(field, data):
    rho2 = data["Density"].astype('float64')**2
    return (data["gradPressureY"] * data["gradDensityZ"] -
            data["gradPressureZ"] * data["gradDensityY"]) / rho2
def _BaroclinicVorticityY(field, data):
    rho2 = data["Density"].astype('float64')**2
    return (data["gradPressureZ"] * data["gradDensityX"] -
            data["gradPressureX"] * data["gradDensityZ"]) / rho2
def _BaroclinicVorticityZ(field, data):
    rho2 = data["Density"].astype('float64')**2
    return (data["gradPressureX"] * data["gradDensityY"] -
            data["gradPressureY"] * data["gradDensityX"]) / rho2
for ax in 'XYZ':
    n = "BaroclinicVorticity%s" % ax
    add_field(n, function=eval("_%s" % n),
          validators=[ValidateSpatial(1, ["Density", "Pressure"])],
          units=r"\rm{s}^{-1}")

def _BaroclinicVorticityMagnitude(field, data):
    return na.sqrt(data["BaroclinicVorticityX"]**2 +
                   data["BaroclinicVorticityY"]**2 +
                   data["BaroclinicVorticityZ"]**2)
add_field("BaroclinicVorticityMagnitude",
          function=_BaroclinicVorticityMagnitude,
          validators=[ValidateSpatial(1, ["Density", "Pressure"])],
          units=r"\rm{s}^{-1}")

def _VorticityX(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["z-velocity"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = (data["z-velocity"][1:-1,sl_right,1:-1] -
                                 data["z-velocity"][1:-1,sl_left,1:-1]) \
                                 / (div_fac*data["dy"].flat[0])
    new_field[1:-1,1:-1,1:-1] -= (data["y-velocity"][1:-1,1:-1,sl_right] -
                                  data["y-velocity"][1:-1,1:-1,sl_left]) \
                                  / (div_fac*data["dz"].flat[0])
    return new_field
def _VorticityY(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["z-velocity"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = (data["x-velocity"][1:-1,1:-1,sl_right] -
                                 data["x-velocity"][1:-1,1:-1,sl_left]) \
                                 / (div_fac*data["dz"].flat[0])
    new_field[1:-1,1:-1,1:-1] -= (data["z-velocity"][sl_right,1:-1,1:-1] -
                                  data["z-velocity"][sl_left,1:-1,1:-1]) \
                                  / (div_fac*data["dx"].flat[0])
    return new_field
def _VorticityZ(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = na.zeros(data["x-velocity"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = (data["y-velocity"][sl_right,1:-1,1:-1] -
                                 data["y-velocity"][sl_left,1:-1,1:-1]) \
                                 / (div_fac*data["dx"].flat[0])
    new_field[1:-1,1:-1,1:-1] -= (data["x-velocity"][1:-1,sl_right,1:-1] -
                                  data["x-velocity"][1:-1,sl_left,1:-1]) \
                                  / (div_fac*data["dy"].flat[0])
    return new_field
def _convertVorticity(data):
    return 1.0/data.convert("cm")
for ax in 'XYZ':
    n = "Vorticity%s" % ax
    add_field(n, function=eval("_%s" % n),
              convert_function=_convertVorticity,
              validators=[ValidateSpatial(1, 
                          ["x-velocity", "y-velocity", "z-velocity"])],
              units=r"\rm{s}^{-1}")

def _VorticityMagnitude(field, data):
    return na.sqrt(data["VorticityX"]**2 +
                   data["VorticityY"]**2 +
                   data["VorticityZ"]**2)
add_field("VorticityMagnitude", function=_VorticityMagnitude,
          validators=[ValidateSpatial(1, 
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units=r"\rm{s}^{-1}")

def _VorticityStretchingX(field, data):
    return data["DivV"] * data["VorticityX"]
def _VorticityStretchingY(field, data):
    return data["DivV"] * data["VorticityY"]
def _VorticityStretchingZ(field, data):
    return data["DivV"] * data["VorticityZ"]
for ax in 'XYZ':
    n = "VorticityStretching%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(0)])
def _VorticityStretchingMagnitude(field, data):
    return na.sqrt(data["VorticityStretchingX"]**2 +
                   data["VorticityStretchingY"]**2 +
                   data["VorticityStretchingZ"]**2)
add_field("VorticityStretchingMagnitude", 
          function=_VorticityStretchingMagnitude,
          validators=[ValidateSpatial(1, 
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units=r"\rm{s}^{-1}")

def _VorticityGrowthX(field, data):
    return -data["VorticityStretchingX"] - data["BaroclinicVorticityX"]
def _VorticityGrowthY(field, data):
    return -data["VorticityStretchingY"] - data["BaroclinicVorticityY"]
def _VorticityGrowthZ(field, data):
    return -data["VorticityStretchingZ"] - data["BaroclinicVorticityZ"]
for ax in 'XYZ':
    n = "VorticityGrowth%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1, 
                          ["x-velocity", "y-velocity", "z-velocity"])],
              units=r"\rm{s}^{-2}")
def _VorticityGrowthMagnitude(field, data):
    result = na.sqrt(data["VorticityGrowthX"]**2 +
                     data["VorticityGrowthY"]**2 +
                     data["VorticityGrowthZ"]**2)
    dot = na.zeros(result.shape)
    for ax in "XYZ":
        dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
    result = na.sign(dot) * result
    return result
add_field("VorticityGrowthMagnitude", function=_VorticityGrowthMagnitude,
          validators=[ValidateSpatial(1, 
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units=r"\rm{s}^{-1}",
          take_log=False)
def _VorticityGrowthMagnitudeABS(field, data):
    return na.sqrt(data["VorticityGrowthX"]**2 +
                   data["VorticityGrowthY"]**2 +
                   data["VorticityGrowthZ"]**2)
add_field("VorticityGrowthMagnitudeABS", function=_VorticityGrowthMagnitudeABS,
          validators=[ValidateSpatial(1, 
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units=r"\rm{s}^{-1}")

def _VorticityGrowthTimescale(field, data):
    domegax_dt = data["VorticityX"] / data["VorticityGrowthX"]
    domegay_dt = data["VorticityY"] / data["VorticityGrowthY"]
    domegaz_dt = data["VorticityZ"] / data["VorticityGrowthZ"]
    return na.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt)
add_field("VorticityGrowthTimescale", function=_VorticityGrowthTimescale,
          validators=[ValidateSpatial(1, 
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units=r"\rm{s}")

########################################################################
# With radiation pressure
########################################################################

def _VorticityRadPressureX(field, data):
    rho = data["Density"].astype('float64')
    return (data["RadAccel2"] * data["gradDensityZ"] -
            data["RadAccel3"] * data["gradDensityY"]) / rho
def _VorticityRadPressureY(field, data):
    rho = data["Density"].astype('float64')
    return (data["RadAccel3"] * data["gradDensityX"] -
            data["RadAccel1"] * data["gradDensityZ"]) / rho
def _VorticityRadPressureZ(field, data):
    rho = data["Density"].astype('float64')
    return (data["RadAccel1"] * data["gradDensityY"] -
            data["RadAccel2"] * data["gradDensityX"]) / rho
def _convertRadAccel(data):
    return data.convert("x-velocity")/data.convert("Time")
for ax in 'XYZ':
    n = "VorticityRadPressure%s" % ax
    add_field(n, function=eval("_%s" % n),
              convert_function=_convertRadAccel,
              validators=[ValidateSpatial(1, 
                   ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
              units=r"\rm{s}^{-1}")

def _VorticityRadPressureMagnitude(field, data):
    return na.sqrt(data["VorticityRadPressureX"]**2 +
                   data["VorticityRadPressureY"]**2 +
                   data["VorticityRadPressureZ"]**2)
add_field("VorticityRadPressureMagnitude",
          function=_VorticityRadPressureMagnitude,
          validators=[ValidateSpatial(1, 
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units=r"\rm{s}^{-1}")

def _VorticityRPGrowthX(field, data):
    return -data["VorticityStretchingX"] - data["BaroclinicVorticityX"] \
           -data["VorticityRadPressureX"]
def _VorticityRPGrowthY(field, data):
    return -data["VorticityStretchingY"] - data["BaroclinicVorticityY"] \
           -data["VorticityRadPressureY"]
def _VorticityRPGrowthZ(field, data):
    return -data["VorticityStretchingZ"] - data["BaroclinicVorticityZ"] \
           -data["VorticityRadPressureZ"]
for ax in 'XYZ':
    n = "VorticityRPGrowth%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1, 
                       ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
              units=r"\rm{s}^{-1}")
def _VorticityRPGrowthMagnitude(field, data):
    result = na.sqrt(data["VorticityRPGrowthX"]**2 +
                     data["VorticityRPGrowthY"]**2 +
                     data["VorticityRPGrowthZ"]**2)
    dot = na.zeros(result.shape)
    for ax in "XYZ":
        dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
    result = na.sign(dot) * result
    return result
add_field("VorticityRPGrowthMagnitude", function=_VorticityGrowthMagnitude,
          validators=[ValidateSpatial(1, 
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units=r"\rm{s}^{-1}",
          take_log=False)
def _VorticityRPGrowthMagnitudeABS(field, data):
    return na.sqrt(data["VorticityRPGrowthX"]**2 +
                   data["VorticityRPGrowthY"]**2 +
                   data["VorticityRPGrowthZ"]**2)
add_field("VorticityRPGrowthMagnitudeABS", 
          function=_VorticityRPGrowthMagnitudeABS,
          validators=[ValidateSpatial(1, 
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units=r"\rm{s}^{-1}")

def _VorticityRPGrowthTimescale(field, data):
    domegax_dt = data["VorticityX"] / data["VorticityRPGrowthX"]
    domegay_dt = data["VorticityY"] / data["VorticityRPGrowthY"]
    domegaz_dt = data["VorticityZ"] / data["VorticityRPGrowthZ"]
    return na.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
add_field("VorticityRPGrowthTimescale", function=_VorticityRPGrowthTimescale,
          validators=[ValidateSpatial(1, 
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units=r"\rm{s}^{-1}")
