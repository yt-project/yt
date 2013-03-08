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
import numpy as np
import inspect
import copy

from yt.funcs import *

from yt.utilities.lib import CICDeposit_3, obtain_rvec, obtain_rv_vec
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
    speed_of_light_cgs, \
    km_per_cm, keV_per_K

from yt.utilities.math_utils import \
    get_sph_r_component, \
    get_sph_theta_component, \
    get_sph_phi_component, \
    get_cyl_r_component, \
    get_cyl_z_component, \
    get_cyl_theta_component, \
    get_cyl_r, get_cyl_theta, \
    get_cyl_z, get_sph_r, \
    get_sph_theta, get_sph_phi, \
    periodic_dist, euclidean_dist

# Note that, despite my newfound efforts to comply with PEP-8,
# I violate it here in order to keep the name/func_name relationship

def _dx(field, data):
    return np.ones(data.ActiveDimensions, dtype=np.float64) * data.dds[0]

add_field('dx', function=_dx, display_field=False,
          validators=[ValidateSpatial(0)])

def _dy(field, data):
    return np.ones(data.ActiveDimensions, dtype=np.float64) * data.dds[1]

add_field('dy', function=_dy, display_field=False,
          validators=[ValidateSpatial(0)])

def _dz(field, data):
    return np.ones(data.ActiveDimensions, dtype=np.float64) * data.dds[2]

add_field('dz', function=_dz,
          display_field=False, validators=[ValidateSpatial(0)])

def _coord_x(field, data):
    dim = data.ActiveDimensions[0]
    return ( ( np.ones(data.ActiveDimensions, dtype=np.float64)
               * np.arange(data.ActiveDimensions[0])[:, None, None] + 0.5 )
             * data['dx'] + data.LeftEdge[0] )

add_field('x', function=_coord_x, display_field=False,
          validators=[ValidateSpatial(0)])

def _coord_y(field, data):
    dim = data.ActiveDimensions[1]
    return ( ( np.ones(data.ActiveDimensions, dtype=np.float64)
               * np.arange(data.ActiveDimensions[1])[None, :, None] + 0.5 )
             * data['dy'] + data.LeftEdge[1] )

add_field('y', function=_coord_y, display_field=False,
          validators=[ValidateSpatial(0)])

def _coord_z(field, data):
    dim = data.ActiveDimensions[2]
    return ( ( np.ones(data.ActiveDimensions, dtype=np.float64)
               * np.arange(data.ActiveDimensions[2])[None, None, :] + 0.5 )
             * data['dz'] + data.LeftEdge[2] )

add_field('z', function=_coord_z, display_field=False,
          validators=[ValidateSpatial(0)])

def _grid_level(field, data):
    return np.ones(data.ActiveDimensions) * data.Level

add_field("grid_level", function=_grid_level,
          validators=[ValidateGridType(), ValidateSpatial(0)])

def _grid_indices(field, data):
    return np.ones(data["ones"].shape) * (data.id - data._id_offset)

add_field("grid_indices", function=_grid_indices, take_log=False,
          validators=[ValidateGridType(), ValidateSpatial(0)])

def _ones_over_dx(field, data):
    return np.ones(data["ones"].shape, dtype=data["density"].dtype) / data['dx']

add_field("ones_over_dx", function=_ones_over_dx, display_field=False)

def _ones(field, data):
    return np.ones(data.shape, dtype=np.float64)

add_field("ones", function=_ones, projection_conversion="unitary",
          display_field=False)
add_field("cells_per_bin", function=_ones, display_field=False)


def _sound_speed(field, data):
    if data.pf["eos_type"] == 1:
        return ( np.ones(data["density"].shape, dtype=np.float64)
                 * data.pf["eos_sound_speed"] )
    return np.sqrt( data.pf.gamma * data["pressure"] / data["density"] )

add_field("sound_speed", function=_sound_speed, units="cm/s")

def _radial_mach_number(field, data):
    """ M{|v|/t_sound} """
    return np.abs(data["radial_velocity"]) / data["sound_speed"]

add_field("radial_mach_number", function=_radial_mach_number)

def _mach_number(field, data):
    """ M{|v|/t_sound} """
    return data["velocity_magnitude"] / data["sound_speed"]

add_field("mach_number", function=_mach_number)

def _courant_time_step(field, data):
    t1 = data["dx"] / (data["sound_speed"] + np.abs(data["velocity_x"]))
    t2 = data["dy"] / (data["sound_speed"] + np.abs(data["velocity_y"]))
    t3 = data["dz"] / (data["sound_speed"] + np.abs(data["velocity_z"]))
    return np.minimum(np.minimum(t1, t2), t3)
def _convert_courant_time_step(data):
    # sound speed and z-velocity are in cm/s, dx is in code
    return data.convert("cm")

add_field("courant_time_step", function=_courant_time_step,
          convert_function=_convert_courant_time_step, units="s")

def _particle_velocity_magnitude(field, data):
    """ M{|v|} """
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity is None:
        bulk_velocity = np.zeros(3)
    return np.sqrt( (data["particle_velocity_x"] - bulk_velocity[0])**2
                    + (data["particle_velocity_y"] - bulk_velocity[1])**2
                    + (data["particle_velocity_z"] - bulk_velocity[2])**2 )

add_field("particle_velocity_magnitude", function=_particle_velocity_magnitude,
          particle_type=True, take_log=False, units="cm/s")

def _velocity_magnitude(field, data):
    """ M{|v|} """
    velocities = obtain_rv_vec(data)
    return np.sqrt(np.sum(velocities**2, axis=0))

add_field("velocity_magnitude", function=_velocity_magnitude,
          take_log=False, units="cm/s")

def _tangential_over_velocity_magnitude(field, data):
    # @todo: can velocity_magnitude be negative?
    return np.abs(data["tangential_velocity"]) / np.abs(data["velocity_magnitude"])

add_field("tangential_over_velocity_magnitude",
          function=_tangential_over_velocity_magnitude, take_log=False)

def _pressure(field, data):
    """ M{(Gamma-1.0)*rho*E} """
    return (data.pf.gamma - 1.0) * data["density"] * data["thermal_energy"]

add_field("pressure", function=_pressure, units="dyne/cm**2")

def _TempkeV(field, data):
    return data["Temperature"] * keV_per_K
add_field("TempkeV", function=_TempkeV, units="keV",
          display_name="Temperature")

def _entropy(field, data):
    if data.has_field_parameter("mu"):
        mw = mh * data.get_field_parameter("mu")
    else:
        mw = mh
    try:
        gammam1 = data.pf["Gamma"] - 1.0
    except:
        gammam1 = 5./3. - 1.0
    return kboltz * data["Temperature"] / \
           ((data["Density"]/mw)**gammam1)
add_field("entropy", units="erg/K", function=_entropy)

### spherical coordinates: r (radius)
def _spherical_r(field, data):
    center = data.get_field_parameter("center")
    coords = obtain_rvec(data).transpose()
    return get_sph_r(vectors, center)

def _convert_spherical_r_cgs(data):
   return data.convert("cm")

add_field("spherical_r", function=_spherical_r,
         validators=[ValidateParameter("center")],
         convert_function=_convert_spherical_r_cgs, units="cm")

### spherical coordinates: theta (angle with respect to normal)
def _spherical_theta(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data).transpose()
    return get_sph_theta(coords, normal)

add_field("spherical_theta", function=_spherical_theta,
         validators=[ValidateParameter("center"), ValidateParameter("normal")])

### spherical coordinates: phi (angle in the plane perpendicular to the normal)
def _spherical_phi(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data).transpose()
    return get_sph_phi(coords, normal)

add_field("spherical_phi", function=_spherical_phi,
         validators=[ValidateParameter("center"), ValidateParameter("normal")])

### cylindrical coordinates: R (radius in the cylinder's plane)
def _cylindrical_r(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data).transpose()
    return get_cyl_r(coords, normal)

def _convert_cylindrical_r_cgs(data):
   return data.convert("cm")

add_field("cylindrical_r", function=_cylindrical_r,
         validators=[ValidateParameter("center"), ValidateParameter("normal")],
         convert_function=_convert_cylindrical_r_cgs, units="cm")
add_field("cylindrical_r_code", function=_cylindrical_r,
          validators=[ValidateParameter("center"), ValidateParameter("normal")])

### cylindrical coordinates: z (height above the cylinder's plane)
def _cylindrical_z(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data).transpose()
    return get_cyl_z(coords, normal)

def _convert_cylindrical_z_cgs(data):
   return data.convert("cm")

add_field("cylindrical_z", function=_cylindrical_z,
          validators=[ValidateParameter("center"), ValidateParameter("normal")],
          convert_function=_convert_cylindrical_z_cgs, units="cm")

### cylindrical coordinates: theta (angle in the cylinder's plane)
def _cylindrical_theta(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data).transpose()
    return get_cyl_theta(coords, normal)

add_field("cylindrical_theta", function=_cylindrical_theta,
          validators=[ValidateParameter("center"), ValidateParameter("normal")])

### The old field DiskAngle is the same as the spherical coordinates'
### 'theta' angle. I'm keeping DiskAngle for backwards compatibility.
# @todo: remove in 3.0?
def _disk_angle(field, data):
    return data["spherical_theta"]

add_field("disk_angle", function=_disk_angle, take_log=False,
          validators=[ValidateParameter("center"), ValidateParameter("normal")],
          display_field=False)

### The old field Height is the same as the cylindrical coordinates' z
### field. I'm keeping Height for backwards compatibility.
# @todo: remove in 3.0?
def _height(field, data):
    return data["cylindrical_z"]

def _convert_height(data):
    return data.convert("cm")

def _convert_height_au(data):
    return data.convert("au")

add_field("height", function=_height, convert_function=_convert_height,
          validators=[ValidateParameter("center"), ValidateParameter("normal")],
          units="cm", display_field=False)
add_field("height_au", function=_height, convert_function=_convert_height_au,
          validators=[ValidateParameter("center"),
                      ValidateParameter("normal")],
          units="AU", display_field=False)

def _cylindrical_radial_velocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data).transpose()
    theta = np.tile(data['cylindrical_theta'], (3, 1)).transpose()
    return get_cyl_r_component(velocities, theta, normal)

def _cylindrical_radial_velocity_absolute(field, data):
    return np.abs(_cylindrical_radial_velocity(field, data))

def _convert_cylindrical_radial_velocity_kms(data):
    return km_per_cm

add_field("cylindrical_radial_velocity", function=_cylindrical_radial_velocity,
          units="cm/s", validators=[ValidateParameter("normal")])
add_field("cylindrical_radial_velocity_absolute",
          function=_cylindrical_radial_velocity_absolute,
          units="cm/s", validators=[ValidateParameter("normal")])
add_field("cylindrical_radial_velocity_kms",
          function=_cylindrical_radial_velocity,
          convert_function=_convert_cylindrical_radial_velocity_kms,
          units="km/s", validators=[ValidateParameter("normal")])
add_field("cylindrical_radial_velocity_kms_absolute",
          function=_cylindrical_radial_velocity_absolute,
          convert_function=_convert_cylindrical_radial_velocity_kms,
          units="km/s", validators=[ValidateParameter("normal")])

def _cylindrical_tangential_velocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data).transpose()
    theta = np.tile(data["cylindrical_theta"], (3, 1)).transpose()
    return get_cyl_theta_component(velocities, theta, normal)

def _cylindrical_tangential_velocity_absolute(field, data):
    return np.abs(_cylindrical_tangential_velocity(field, data))

def _convert_cylindrical_tangential_velocity_kms(data):
    return km_per_cm

add_field("cylindrical_tangential_velocity",
          function=_cylindrical_tangential_velocity,
          units="cm/s", validators=[ValidateParameter("normal")])
add_field("cylindrical_tangential_velocity_absolute",
          function=_cylindrical_tangential_velocity_absolute,
          units="cm/s", validators=[ValidateParameter("normal")])

def _dynamical_time(field, data):
    """
    sqrt(3 pi / (16 G rho))
    """
    return np.sqrt(3.0 * np.pi / (16.0 * G * data["density"]))

add_field("dynamical_time", function=_dynamical_time, units="s")

def jeans_mass(field, data):
    return ( MJ_constant
             * ((data["temperature"] / data["mean_molecular_weight"])**(1.5))
             * (data["density"]**(-0.5)) )

add_field("jeans_mass", function=jeans_mass, units="g")

def _cell_mass(field, data):
    return data["density"] * data["cell_volume"]

add_field("cell_mass", function=_cell_mass, units="g")

def _total_mass(field, data):
    return (data["density"] + data["dark_matter_density"]) * data["cell_volume"]

add_field("total_mass", function=_total_mass, units="g")

def _star_mass(field, data):
    return data["star_density"] * data["cell_volume"]

add_field("star_mass", units="g", function=_star_mass)

def _matter_density(field, data):
    return (data["density"] + data["dark_matter_density"])

add_field("matter_density", function=_matter_density, units="g/cm**3")

def _comoving_density(field, data):
    z = data.pf.current_redshift
    return data["density"] / (1.0 + z)**3

add_field("comoving_density", function=_comoving_density, units="g/cm**3")

# This is rho_total / rho_cr(z).
def _convert_overdensity(data):
    return 1 / (rho_crit_now * data.pf.hubble_constant**2 *
                (1+data.pf.current_redshift)**3)
add_field("Overdensity",function=_matter_density,
          convert_function=_convert_overdensity, units="")

# This is rho_matter / <rho_matter> - 1.0
def _overdensity(field, data):
    omega_m = data.pf.omega_matter
    h = data.pf.hubble_constant
    z = data.pf.current_redshift
    rho_m = rho_crit_now * h**2 * omega_m * (1.0 + z)**3
    return data["matter_density"] / rho_m - 1.0

add_field("overdensity", function=_overdensity)

# This is rho_baryon / <rho_baryon> - 1.0.
def _baryon_overdensity(field, data):
    # @todo: should we provide this field if the dataset doesn't have omega_b?
    if data.pf.has_key('omega_baryon_now'):
        omega_baryon_now = data.pf['omega_baryon_now']
    else:
        omega_baryon_now = 0.0441

    return data["density"] / (omega_baryon_now * rho_crit_now *
                              (data.pf["CosmologyHubbleConstantNow"]**2) *
                              ((1.0 + data.pf["CosmologyCurrentRedshift"])**3))

add_field("baryon_overdensity", function=_baryon_overdensity)

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
add_field("WeakLensingConvergence", function=_overdensity,
          convert_function=_convertConvergence,
          projection_conversion='mpccm')

def _cell_volume(field, data):
    if data["dx"].size == 1:
        try:
            return ( data["dx"] * data["dy"] * data["dx"]
                     * np.ones(data.ActiveDimensions, dtype=np.float64) )
        except AttributeError:
            return data["dx"] * data["dy"] * data["dx"]
    return data["dx"] * data["dy"] * data["dz"]

add_field("cell_volume", units="cm**3", function=_cell_volume)

def _chandra_emissivity(field, data):
    logT0 = np.log10(data["temperature"]) - 7
    return ( data["number_density"].astype(np.float64)**2
             * ( 10**(-0.0103 * logT0**8
                      +0.0417 * logT0**7
                      -0.0636 * logT0**6
                      +0.1149 * logT0**5
                      -0.3151 * logT0**4
                      +0.6655 * logT0**3
                      -1.1256 * logT0**2
                      +1.0026 * logT0**1
                      -0.6984 * logT0)
                 + data["metallicity"] * 10**(0.0305 * logT0**11
                                              -0.0045 * logT0**10
                                              -0.3620 * logT0**9
                                              +0.0513 * logT0**8
                                              +1.6669 * logT0**7
                                              -0.3854 * logT0**6
                                              -3.3604 * logT0**5
                                              +0.4728 * logT0**4
                                              +4.5774 * logT0**3
                                              -2.3661 * logT0**2
                                              -1.6667 * logT0**1
                                              -0.2193 * logT0) ) )

def _convert_chandra_emissivity(data):
    return 1.0  # 1.0e-23*0.76**2

add_field("chandra_emissivity", function=_chandra_emissivity,
          convert_function=_convert_chandra_emissivity,
          projection_conversion="1")

def _xray_emissivity(field, data):
    return ( data["density"].astype(np.float64)**2
             * data["temperature"]**0.5 )

def _convert_xray_emissivity(data):
    return 2.168e60

add_field("xray_emissivity", function=_xray_emissivity,
          convert_function=_convert_xray_emissivity,
          projection_conversion="1")

def _sz_kinetic(field, data):
    vel_axis = data.get_field_parameter("axis")
    if vel_axis > 2:
        raise NeedsParameter(["axis"])
    vel = data["velocity_%s" % ({0: "x", 1: "y", 2: "z"}[vel_axis])]
    return (vel * data["density"])

def _convert_sz_kinetic(data):
    return 0.88 * sigma_thompson / mh / clight

add_field("sz_kinetic", function=_sz_kinetic,
          convert_function=_convert_sz_kinetic,
          validators=[ValidateParameter("axis")])

def _szy(field, data):
    return data["density"] * data["temperature"]

def _convert_szy(data):
    conv = 0.88 / mh * kboltz / (me * clight*clight) * sigma_thompson
    return conv

add_field("szy", function=_szy, convert_function=_convert_szy)

def _averaged_density(field, data):
    nx, ny, nz = data["density"].shape
    new_field = np.zeros((nx-2, ny-2, nz-2), dtype=np.float64)
    weight_field = np.zeros((nx-2, ny-2, nz-2), dtype=np.float64)
    i_i, j_i, k_i = np.mgrid[0:3, 0:3, 0:3]

    for i, j, k in zip(i_i.ravel(), j_i.ravel(), k_i.ravel()):
        sl = [slice(i, nx-(2-i)), slice(j, ny-(2-j)), slice(k, nz-(2-k))]
        new_field += data["density"][sl] * data["cell_mass"][sl]
        weight_field += data["cell_mass"][sl]

    # Now some fancy footwork
    new_field2 = np.zeros((nx, ny, nz))
    new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
    return new_field2

add_field("averaged_density", function=_averaged_density,
          validators=[ValidateSpatial(1, ["density"])])

def _div_v(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    ds = div_fac * data["dx"].flat[0]
    f  = data["x-velocity"][sl_right,1:-1,1:-1]/ds
    f -= data["x-velocity"][sl_left ,1:-1,1:-1]/ds
    if data.pf.dimensionality > 1:
        ds = div_fac * data["dy"].flat[0]
        f += data["y-velocity"][1:-1,sl_right,1:-1]/ds
        f -= data["y-velocity"][1:-1,sl_left ,1:-1]/ds
    if data.pf.dimensionality > 2:
        ds = div_fac * data["dz"].flat[0]
        f += data["z-velocity"][1:-1,1:-1,sl_right]/ds
        f -= data["z-velocity"][1:-1,1:-1,sl_left ]/ds
    new_field = np.zeros(data["x-velocity"].shape, dtype=np.float64)
    new_field[1:-1,1:-1,1:-1] = f
    return new_field

def _convert_div_v(data):
    return data.convert("cm")**-1.0

add_field("div_v", function=_div_v,
          validators=[ValidateSpatial(1, ["velocity_x", "velocity_y",
                                          "velocity_z"])],
          units="1/s", take_log=False, convert_function=_convert_div_v)

def _absolute_div_v(field, data):
    return np.abs(data["div_v"])

add_field("absolute_div_v", function=_absolute_div_v, units="1/s")

def _contours(field, data):
    return -np.ones_like(data["ones"])

add_field("contours", validators=[ValidateSpatial(0)], take_log=False,
          display_field=False, function=_contours)
add_field("temp_contours", function=_contours,
          validators=[ValidateSpatial(0), ValidateGridType()],
          take_log=False, display_field=False)

def obtain_velocities(data):
    return obtain_rv_vec(data)

def _specific_angular_momentum_x(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return yv * rv[2, :] - zv * rv[1, :]

def _specific_angular_momentum_y(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return - (xv * rv[2, :] - zv * rv[0, :])

def _specific_angular_momentum_z(field, data):
    xv, yv, zv = obtain_velocities(data)
    rv = obtain_rvec(data)
    return xv * rv[1, :] - yv * rv[0, :]

add_field("specific_angular_momentum_x", function=_specific_angular_momentum_x,
          units="cm**2/s", validators=[ValidateParameter("center")])
add_field("specific_angular_momentum_y", function=_specific_angular_momentum_y,
          units="cm**2/s", validators=[ValidateParameter("center")])
add_field("specific_angular_momentum_z", function=_specific_angular_momentum_z,
          units="cm**2/s", validators=[ValidateParameter("center")])

def _angular_momentum_x(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumX"]
add_field("AngularMomentumX", function=_angular_momentum_x,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])
def _angular_momentum_y(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumY"]
add_field("AngularMomentumY", function=_angular_momentum_y,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])
def _angular_momentum_z(field, data):
    return data["CellMass"] * data["SpecificAngularMomentumZ"]
add_field("AngularMomentumZ", function=_angular_momentum_z,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])

def _ParticleSpecificAngularMomentum(field, data):
    """
    Calculate the angular of a particle velocity.  Returns a vector for each
    particle.
    """
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = np.zeros(3, dtype=np.float64)
    xv = data["particle_velocity_x"] - bv[0]
    yv = data["particle_velocity_y"] - bv[1]
    zv = data["particle_velocity_z"] - bv[2]
    center = data.get_field_parameter('center')
    coords = np.array([data['particle_position_x'],
                       data['particle_position_y'],
                       data['particle_position_z']], dtype=np.float64)
    new_shape = tuple([3] + [1]*(len(coords.shape)-1))
    r_vec = coords - np.reshape(center,new_shape)
    v_vec = np.array([xv,yv,zv], dtype=np.float64)
    return np.cross(r_vec, v_vec, axis=0)
#add_field("ParticleSpecificAngularMomentum",
#          function=_ParticleSpecificAngularMomentum, particle_type=True,
#          convert_function=_convertSpecificAngularMomentum, vector_field=True,
#          units=r"\rm{cm}^2/\rm{s}", validators=[ValidateParameter('center')])
def _convertSpecificAngularMomentumKMSMPC(data):
    return km_per_cm*data.convert("mpc")
#add_field("ParticleSpecificAngularMomentumKMSMPC",
#          function=_ParticleSpecificAngularMomentum, particle_type=True,
#          convert_function=_convertSpecificAngularMomentumKMSMPC, vector_field=True,
#          units=r"\rm{km}\rm{Mpc}/\rm{s}", validators=[ValidateParameter('center')])

def _ParticleSpecificAngularMomentumX(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = np.zeros(3, dtype=np.float64)
    center = data.get_field_parameter('center')
    y = data["particle_position_y"] - center[1]
    z = data["particle_position_z"] - center[2]
    yv = data["particle_velocity_y"] - bv[1]
    zv = data["particle_velocity_z"] - bv[2]
    return yv*z - zv*y
def _ParticleSpecificAngularMomentumY(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = np.zeros(3, dtype=np.float64)
    center = data.get_field_parameter('center')
    x = data["particle_position_x"] - center[0]
    z = data["particle_position_z"] - center[2]
    xv = data["particle_velocity_x"] - bv[0]
    zv = data["particle_velocity_z"] - bv[2]
    return -(xv*z - zv*x)
def _ParticleSpecificAngularMomentumZ(field, data):
    if data.has_field_parameter("bulk_velocity"):
        bv = data.get_field_parameter("bulk_velocity")
    else: bv = np.zeros(3, dtype=np.float64)
    center = data.get_field_parameter('center')
    x = data["particle_position_x"] - center[0]
    y = data["particle_position_y"] - center[1]
    xv = data["particle_velocity_x"] - bv[0]
    yv = data["particle_velocity_y"] - bv[1]
    return xv*y - yv*x
for ax in 'XYZ':
    n = "ParticleSpecificAngularMomentum%s" % ax
    add_field(n, function=eval("_%s" % n), particle_type=True,
              units="cm**2/s", validators=[ValidateParameter("center")])
    add_field(n + "KMSMPC", function=eval("_%s" % n), particle_type=True,
              convert_function=_convertSpecificAngularMomentumKMSMPC,
              units="cm**2/s", validators=[ValidateParameter("center")])

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
         units="g*cm**2/s", particle_type=True,
         validators=[ValidateParameter('center')])
def _ParticleAngularMomentumY(field, data):
    return data["CellMass"] * data["ParticleSpecificAngularMomentumY"]
add_field("ParticleAngularMomentumY", function=_ParticleAngularMomentumY,
         units="g*cm**2/s", particle_type=True,
         validators=[ValidateParameter('center')])
def _ParticleAngularMomentumZ(field, data):
    return data["CellMass"] * data["ParticleSpecificAngularMomentumZ"]
add_field("ParticleAngularMomentumZ", function=_ParticleAngularMomentumZ,
         units="g*cm**2/s", particle_type=True,
         validators=[ValidateParameter('center')])

def get_radius(data, field_prefix):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = np.zeros(data[field_prefix+"x"].shape, dtype='float64')
    r = radius.copy()
    if any(data.pf.periodicity):
        rdw = radius.copy()
    for i, ax in enumerate('xyz'):
        np.subtract(data["%s%s" % (field_prefix, ax)], center[i], r)
        if data.pf.periodicity[i] == True:
            np.subtract(DW[i], r, rdw)
            np.abs(r, r)
            np.minimum(r, rdw, r)
        np.power(r, 2.0, r)
        np.add(radius, r, radius)
    np.sqrt(radius, radius)
    return radius

def _ParticleRadius(field, data):
    return get_radius(data, "particle_position_")
def _Radius(field, data):
    return get_radius(data, "")

def _ConvertRadiusCGS(data):
    return data.convert("cm")
add_field("ParticleRadius", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units="cm",
          particle_type = True,
          display_name = "Particle Radius")
add_field("Radius", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusCGS, units="cm")

def _ConvertRadiusMpc(data):
    return data.convert("mpc")
add_field("RadiusMpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units="Mpc",
          display_name = "Radius")
add_field("ParticleRadiusMpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusMpc, units="Mpc",
          particle_type=True,
          display_name = "Particle Radius")

def _ConvertRadiuskpc(data):
    return data.convert("kpc")
add_field("ParticleRadiuskpc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units="kpc",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuskpc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units="kpc",
          display_name = "Radius")

def _ConvertRadiuskpch(data):
    return data.convert("kpch")
add_field("ParticleRadiuskpch", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpch, units="kpc/h",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuskpch", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuskpc, units="kpc/h",
          display_name = "Radius")

def _ConvertRadiuspc(data):
    return data.convert("pc")
add_field("ParticleRadiuspc", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units="pc",
          particle_type=True,
          display_name = "Particle Radius")
add_field("Radiuspc", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiuspc, units="pc",
          display_name="Radius")

def _ConvertRadiusAU(data):
    return data.convert("au")
add_field("ParticleRadiusAU", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units="AU",
          particle_type=True,
          display_name = "Particle Radius")
add_field("RadiusAU", function=_Radius,
          validators=[ValidateParameter("center")],
          convert_function = _ConvertRadiusAU, units="AU",
          display_name = "Radius")

add_field("ParticleRadiusCode", function=_ParticleRadius,
          validators=[ValidateParameter("center")],
          particle_type=True,
          display_name = "Particle Radius (code)")
add_field("RadiusCode", function=_Radius,
          validators=[ValidateParameter("center")],
          display_name = "Radius (code)")

def _RadialVelocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data).transpose()
    theta = np.tile(data['sph_theta'], (3, 1)).transpose()
    phi   = np.tile(data['sph_phi'], (3, 1)).transpose()

    return get_sph_r_component(velocities, theta, phi, normal)

def _RadialVelocityABS(field, data):
    return np.abs(_RadialVelocity(field, data))
def _ConvertRadialVelocityKMS(data):
    return km_per_cm
add_field("RadialVelocity", function=_RadialVelocity,
          units="cm/s")
add_field("RadialVelocityABS", function=_RadialVelocityABS,
          units="cm/s")
add_field("RadialVelocityKMS", function=_RadialVelocity,
          convert_function=_ConvertRadialVelocityKMS, units="km/s")
add_field("RadialVelocityKMSABS", function=_RadialVelocityABS,
          convert_function=_ConvertRadialVelocityKMS, units="km/s")

def _TangentialVelocity(field, data):
    return np.sqrt(data["VelocityMagnitude"]**2.0
                 - data["RadialVelocity"]**2.0)
add_field("TangentialVelocity",
          function=_TangentialVelocity,
          take_log=False, units="cm/s")

def _CuttingPlaneVelocityX(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = np.zeros(3)
    v_vec = np.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                - bulk_velocity[...,np.newaxis]
    return np.dot(x_vec, v_vec)
add_field("CuttingPlaneVelocityX",
          function=_CuttingPlaneVelocityX,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="km/s")
def _CuttingPlaneVelocityY(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = np.zeros(3)
    v_vec = np.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                - bulk_velocity[...,np.newaxis]
    return np.dot(y_vec, v_vec)
add_field("CuttingPlaneVelocityY",
          function=_CuttingPlaneVelocityY,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="km/s")

def _CuttingPlaneBx(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    b_vec = np.array([data["B%s" % ax] for ax in 'xyz'])
    return np.dot(x_vec, b_vec)
add_field("CuttingPlaneBx",
          function=_CuttingPlaneBx,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="gauss")
def _CuttingPlaneBy(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    b_vec = np.array([data["B%s" % ax] for ax in 'xyz'])
    return np.dot(y_vec, b_vec)
add_field("CuttingPlaneBy",
          function=_CuttingPlaneBy,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="gauss")

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
          units="Msun")

def _convertDensity(data):
    return data.convert("Density")
def _pdensity(field, data):
    blank = np.zeros(data.ActiveDimensions, dtype='float32')
    if data["particle_position_x"].size == 0: return blank
    CICDeposit_3(data["particle_position_x"].astype(np.float64),
                 data["particle_position_y"].astype(np.float64),
                 data["particle_position_z"].astype(np.float64),
                 data["particle_mass"].astype(np.float32),
                 data["particle_position_x"].size,
                 blank, np.array(data.LeftEdge).astype(np.float64),
                 np.array(data.ActiveDimensions).astype(np.int32),
                 np.float64(data['dx']))
    return blank
add_field("particle_density", function=_pdensity,
          validators=[ValidateGridType()], convert_function=_convertDensity,
          display_name=r"\mathrm{Particle}\/\mathrm{Density})")

def _MagneticEnergy(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    MagneticEnergy field to deal with non-unitary \mu_0.
    """
    return (data["Bx"]**2 + data["By"]**2 + data["Bz"]**2)/(8*np.pi)
add_field("MagneticEnergy",function=_MagneticEnergy,
          units="erg / cm**3",
          display_name=r"\rm{Magnetic}\/\rm{Energy}")

def _BMagnitude(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    BMagnitude field to deal with non-unitary \mu_0.
    """
    return np.sqrt((data["Bx"]**2 + data["By"]**2 + data["Bz"]**2))
add_field("BMagnitude",
          function=_BMagnitude,
          display_name=r"|B|", units="gauss")

def _PlasmaBeta(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    PlasmaBeta field to deal with non-unitary \mu_0.
    """
    return data['Pressure']/data['MagneticEnergy']
add_field("PlasmaBeta",
          function=_PlasmaBeta,
          display_name=r"\rm{Plasma}\/\beta", units="")

def _MagneticPressure(field,data):
    return data['MagneticEnergy']
add_field("MagneticPressure",
          function=_MagneticPressure,
          display_name=r"\rm{Magnetic}\/\rm{Energy}",
          units="erg / cm**3")

def _BPoloidal(field,data):
    normal = data.get_field_parameter("normal")

    Bfields = np.array([data['Bx'], data['By'], data['Bz']])

    theta = data['sph_theta']
    phi   = data['sph_phi']

    return get_sph_theta_component(Bfields, theta, phi, normal)

add_field("BPoloidal", function=_BPoloidal,
          units="gauss",
          validators=[ValidateParameter("normal")])

def _BToroidal(field,data):
    normal = data.get_field_parameter("normal")

    Bfields = np.array([data['Bx'], data['By'], data['Bz']])

    phi   = data['sph_phi']

    return get_sph_phi_component(Bfields, phi, normal)

add_field("BToroidal", function=_BToroidal,
          units="gauss",
          validators=[ValidateParameter("normal")])

def _BRadial(field,data):
    normal = data.get_field_parameter("normal")

    Bfields = np.array([data['Bx'], data['By'], data['Bz']])

    theta = data['sph_theta']
    phi   = data['sph_phi']

    return get_sph_r_component(Bfields, theta, phi, normal)

add_field("BRadial", function=_BPoloidal,
          units="gauss",
          validators=[ValidateParameter("normal")])

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
    new_field = np.zeros(data["x-velocity"].shape)
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
    new_field = np.abs(new_field)
    return new_field
def _convertVorticitySquared(data):
    return data.convert("cm")**-2.0
add_field("VorticitySquared", function=_VorticitySquared,
          validators=[ValidateSpatial(1,
              ["x-velocity","y-velocity","z-velocity"])],
          units="s**-2",
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
    new_field = np.zeros(data["Pressure"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["Pressure"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["Pressure"].shape, dtype=np.float64)
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
              units="dyne/cm**3")

def _gradPressureMagnitude(field, data):
    return np.sqrt(data["gradPressureX"]**2 +
                   data["gradPressureY"]**2 +
                   data["gradPressureZ"]**2)
add_field("gradPressureMagnitude", function=_gradPressureMagnitude,
          validators=[ValidateSpatial(1, ["Pressure"])],
          units="dyne/cm**3")

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
    new_field = np.zeros(data["Density"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["Density"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["Density"].shape, dtype=np.float64)
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
              units="g/cm**4")

def _gradDensityMagnitude(field, data):
    return np.sqrt(data["gradDensityX"]**2 +
                   data["gradDensityY"]**2 +
                   data["gradDensityZ"]**2)
add_field("gradDensityMagnitude", function=_gradDensityMagnitude,
          validators=[ValidateSpatial(1, ["Density"])],
          units="g/cm**4")

def _BaroclinicVorticityX(field, data):
    rho2 = data["Density"].astype(np.float64)**2
    return (data["gradPressureY"] * data["gradDensityZ"] -
            data["gradPressureZ"] * data["gradDensityY"]) / rho2
def _BaroclinicVorticityY(field, data):
    rho2 = data["Density"].astype(np.float64)**2
    return (data["gradPressureZ"] * data["gradDensityX"] -
            data["gradPressureX"] * data["gradDensityZ"]) / rho2
def _BaroclinicVorticityZ(field, data):
    rho2 = data["Density"].astype(np.float64)**2
    return (data["gradPressureX"] * data["gradDensityY"] -
            data["gradPressureY"] * data["gradDensityX"]) / rho2
for ax in 'XYZ':
    n = "BaroclinicVorticity%s" % ax
    add_field(n, function=eval("_%s" % n),
          validators=[ValidateSpatial(1, ["Density", "Pressure"])],
          units="1/s")

def _BaroclinicVorticityMagnitude(field, data):
    return np.sqrt(data["BaroclinicVorticityX"]**2 +
                   data["BaroclinicVorticityY"]**2 +
                   data["BaroclinicVorticityZ"]**2)
add_field("BaroclinicVorticityMagnitude",
          function=_BaroclinicVorticityMagnitude,
          validators=[ValidateSpatial(1, ["Density", "Pressure"])],
          units="1/s")

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
    new_field = np.zeros(data["z-velocity"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["z-velocity"].shape, dtype=np.float64)
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
    new_field = np.zeros(data["x-velocity"].shape, dtype=np.float64)
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
              units="1/s")

def _VorticityMagnitude(field, data):
    return np.sqrt(data["VorticityX"]**2 +
                   data["VorticityY"]**2 +
                   data["VorticityZ"]**2)
add_field("VorticityMagnitude", function=_VorticityMagnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

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
    return np.sqrt(data["VorticityStretchingX"]**2 +
                   data["VorticityStretchingY"]**2 +
                   data["VorticityStretchingZ"]**2)
add_field("VorticityStretchingMagnitude",
          function=_VorticityStretchingMagnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

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
              units="1/s")
def _VorticityGrowthMagnitude(field, data):
    result = np.sqrt(data["VorticityGrowthX"]**2 +
                     data["VorticityGrowthY"]**2 +
                     data["VorticityGrowthZ"]**2)
    dot = np.zeros(result.shape)
    for ax in "XYZ":
        dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
    result = np.sign(dot) * result
    return result
add_field("VorticityGrowthMagnitude", function=_VorticityGrowthMagnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s",
          take_log=False)
def _VorticityGrowthMagnitudeABS(field, data):
    return np.sqrt(data["VorticityGrowthX"]**2 +
                   data["VorticityGrowthY"]**2 +
                   data["VorticityGrowthZ"]**2)
add_field("VorticityGrowthMagnitudeABS", function=_VorticityGrowthMagnitudeABS,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

def _VorticityGrowthTimescale(field, data):
    domegax_dt = data["VorticityX"] / data["VorticityGrowthX"]
    domegay_dt = data["VorticityY"] / data["VorticityGrowthY"]
    domegaz_dt = data["VorticityZ"] / data["VorticityGrowthZ"]
    return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt)
add_field("VorticityGrowthTimescale", function=_VorticityGrowthTimescale,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="s")

########################################################################
# With radiation pressure
########################################################################

def _VorticityRadPressureX(field, data):
    rho = data["Density"].astype(np.float64)
    return (data["RadAccel2"] * data["gradDensityZ"] -
            data["RadAccel3"] * data["gradDensityY"]) / rho
def _VorticityRadPressureY(field, data):
    rho = data["Density"].astype(np.float64)
    return (data["RadAccel3"] * data["gradDensityX"] -
            data["RadAccel1"] * data["gradDensityZ"]) / rho
def _VorticityRadPressureZ(field, data):
    rho = data["Density"].astype(np.float64)
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
              units="1/s")

def _VorticityRadPressureMagnitude(field, data):
    return np.sqrt(data["VorticityRadPressureX"]**2 +
                   data["VorticityRadPressureY"]**2 +
                   data["VorticityRadPressureZ"]**2)
add_field("VorticityRadPressureMagnitude",
          function=_VorticityRadPressureMagnitude,
          validators=[ValidateSpatial(1,
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units="1/s")

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
              units="1/s")
def _VorticityRPGrowthMagnitude(field, data):
    result = np.sqrt(data["VorticityRPGrowthX"]**2 +
                     data["VorticityRPGrowthY"]**2 +
                     data["VorticityRPGrowthZ"]**2)
    dot = np.zeros(result.shape)
    for ax in "XYZ":
        dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
    result = np.sign(dot) * result
    return result
add_field("VorticityRPGrowthMagnitude", function=_VorticityGrowthMagnitude,
          validators=[ValidateSpatial(1,
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units="1/s",
          take_log=False)
def _VorticityRPGrowthMagnitudeABS(field, data):
    return np.sqrt(data["VorticityRPGrowthX"]**2 +
                   data["VorticityRPGrowthY"]**2 +
                   data["VorticityRPGrowthZ"]**2)
add_field("VorticityRPGrowthMagnitudeABS",
          function=_VorticityRPGrowthMagnitudeABS,
          validators=[ValidateSpatial(1,
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units="1/s")

def _VorticityRPGrowthTimescale(field, data):
    domegax_dt = data["VorticityX"] / data["VorticityRPGrowthX"]
    domegay_dt = data["VorticityY"] / data["VorticityRPGrowthY"]
    domegaz_dt = data["VorticityZ"] / data["VorticityRPGrowthZ"]
    return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
add_field("VorticityRPGrowthTimescale", function=_VorticityRPGrowthTimescale,
          validators=[ValidateSpatial(1,
                      ["Density", "RadAccel1", "RadAccel2", "RadAccel3"])],
          units="1/s")
