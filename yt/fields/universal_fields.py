"""
The basic field info container resides here.  These classes, code specific and
universal, are the means by which we access fields across YT, both derived and
native.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import types
import numpy as np
import inspect
import copy

from yt.funcs import *

from yt.data_objects.yt_array import YTArray
from yt.utilities.lib import obtain_rvec, obtain_rv_vec
from yt.utilities.math_utils import resize_vector
from yt.utilities.cosmology import Cosmology
from yt.data_objects.field_info_container import \
    add_field, \
    ValidateDataField, \
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsGridType, \
    NeedsOriginalGrid, \
    NeedsDataField, \
    NeedsProperty, \
    NeedsParameter, \
    NullFunc

from yt.utilities.physical_constants import \
     mass_sun_cgs, \
    mh, \
    me, \
    sigma_thompson, \
    clight, \
    kboltz, \
    G, \
    rho_crit_now, \
    speed_of_light_cgs, \
    km_per_cm

def _get_conv(unit):
    def _conv(data):
        return data.convert(unit)
    return _conv

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

def _grid_level(field, data):
    return np.ones(data.ActiveDimensions)*(data.Level)
add_field("grid_level", function=_grid_level,
          validators=[ValidateGridType(),
                      ValidateSpatial(0)])

def _grid_indices(field, data):
    return np.ones(data["ones"].shape)*(data.id-data._id_offset)
add_field("grid_indices", function=_grid_indices,
          validators=[ValidateGridType(),
                      ValidateSpatial(0)], take_log=False)

def _ones_over_dx(field, data):
    return np.ones(data["ones"].shape,
                   dtype=data["density"].dtype)/data['dx']
add_field("ones_over_dx", function=_ones_over_dx,
          units = "1 / cm",
          display_field=False)

def _zeros(field, data):
    arr = np.zeros(data["ones"].shape, dtype='float64')
    return field.apply_units(arr)

add_field("zeros", function=_zeros,
          units = "",
          projection_conversion="unitary",
          display_field=False)

def _ones(field, data):
    arr = np.ones(data.ires.shape, dtype="float64")
    if data._spatial:
        return data._reshape_vals(arr)
    return field.apply_units(arr)

add_field("ones", function=_ones,
          projection_conversion="unitary",
          units = "",
          display_field=False)

add_field("cells_per_bin", function=_ones,
          display_field = False)

def _sound_speed(field, data):
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
    t1 = data["dx"] / (data["sound_speed"] + np.abs(data["x-velocity"]))
    t2 = data["dy"] / (data["sound_speed"] + np.abs(data["y-velocity"]))
    t3 = data["dz"] / (data["sound_speed"] + np.abs(data["z-velocity"]))
    tr = np.minimum(np.minimum(t1, t2), t3)
    return field.apply_units(tr)

add_field("courant_time_step", function=_courant_time_step,
          units="s")

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

def _entropy(field, data):
    mw = data.get_field_parameter("mu")
    if mw is None:
        mw = 1.0
    mw = mh
    try:
        gammam1 = data.pf.gamma - 1.0
    except:
        gammam1 = 5./3. - 1.0
    return kboltz * data["temperature"] / \
           ((data["density"]/mw)**gammam1)
add_field("entropy", units="erg/K", function=_entropy)

### spherical coordinates: r (radius)
def _spherical_r(field, data):
    center = data.get_field_parameter("center")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
    return get_sph_r(coords)

add_field("spherical_r", function=_spherical_r,
         validators=[ValidateParameter("center")],
         units="cm")

### spherical coordinates: theta (angle with respect to normal)
def _spherical_theta(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
    return get_sph_theta(coords, normal)

add_field("spherical_theta", function=_spherical_theta,
         validators=[ValidateParameter("center"), ValidateParameter("normal")])

### spherical coordinates: phi (angle in the plane perpendicular to the normal)
def _spherical_phi(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
    return get_sph_phi(coords, normal)

add_field("spherical_phi", function=_spherical_phi,
         validators=[ValidateParameter("center"), ValidateParameter("normal")])

### cylindrical coordinates: R (radius in the cylinder's plane)
def _cylindrical_r(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
    return get_cyl_r(coords, normal)

add_field("cylindrical_r", function=_cylindrical_r,
         validators=[ValidateParameter("center"), ValidateParameter("normal")],
         units="cm")

### cylindrical coordinates: z (height above the cylinder's plane)
def _cylindrical_z(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
    return get_cyl_z(coords, normal)

add_field("cylindrical_z", function=_cylindrical_z,
          validators=[ValidateParameter("center"), ValidateParameter("normal")],
          units="cm")

### cylindrical coordinates: theta (angle in the cylinder's plane)
def _cylindrical_theta(field, data):
    center = data.get_field_parameter("center")
    normal = data.get_field_parameter("normal")
    coords = obtain_rvec(data)
    coords[0,...] -= center[0]
    coords[1,...] -= center[1]
    coords[2,...] -= center[2]
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

add_field("height", function=_height, 
          validators=[ValidateParameter("center"), ValidateParameter("normal")],
          units="cm", display_field=False)

def _cylindrical_radial_velocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data)
    theta = resize_vector(data['cylindrical_theta'], velocities)
    return get_cyl_r_component(velocities, theta, normal)

def _cylindrical_radial_velocity_absolute(field, data):
    return np.abs(_cylindrical_radial_velocity(field, data))

add_field("cylindrical_radial_velocity", function=_cylindrical_radial_velocity,
          units="cm/s", validators=[ValidateParameter("normal")])
add_field("cylindrical_radial_velocity_absolute",
          function=_cylindrical_radial_velocity_absolute,
          units="cm/s", validators=[ValidateParameter("normal")])

def _cylindrical_tangential_velocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data)
    theta = data['cylindrical_theta'].copy()
    theta = np.tile(theta, (3,) + (1,)*len(theta.shape))
    return get_cyl_theta_component(velocities, theta, normal)

def _cylindrical_tangential_velocity_absolute(field, data):
    return np.abs(_cylindrical_tangential_velocity(field, data))

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
    MJ_constant = (((5.0 * kboltz) / (G * mh)) ** (1.5)) * \
      (3.0 / (4.0 * np.pi)) ** (0.5) / mass_sun_cgs
    return ( MJ_constant
             * ((data["temperature"] / data["mean_molecular_weight"])**(1.5))
             * (data["density"]**(-0.5)) )

add_field("jeans_mass", function=jeans_mass, units="g")

def _cell_mass(field, data):
    return data["density"] * data["cell_volume"]

add_field("cell_mass", function=_cell_mass, units="g")

def _total_mass(field, data):
    return data["cell_mass"] + data["deposit","all_mass"]

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
def _overdensity(field, data):
    return data["matter_density"] / (rho_crit_now * data.pf.hubble_constant**2 *
                (1+data.pf.current_redshift)**3)
add_field("overdensity", function=_overdensity)

# This is rho_matter / <rho_matter> - 1.0
def _density_perturbation(field, data):
    omega_m = data.pf.omega_matter
    h = data.pf.hubble_constant
    z = data.pf.current_redshift
    rho_m = rho_crit_now * h**2 * omega_m * (1.0 + z)**3
    return data["matter_density"] / rho_m - 1.0

add_field("density_perturbation", function=_overdensity)

# This is rho_baryon / <rho_baryon>
def _baryon_overdensity(field, data):
    # @todo: should we provide this field if the dataset doesn't have omega_b?
    if data.pf.has_key('omega_baryon_now'):
        omega_baryon_now = data.pf['omega_baryon_now']
    else:
        omega_baryon_now = 0.0441

    # These are enzo parameters and should be changed.
    return data["density"] / (omega_baryon_now * rho_crit_now *
                              (data.pf["CosmologyHubbleConstantNow"]**2) *
                              ((1.0 + data.pf["CosmologyCurrentRedshift"])**3))

add_field("baryon_overdensity", function=_baryon_overdensity)

#FIXME

# Weak lensing convergence.
# Eqn 4 of Metzler, White, & Loken (2001, ApJ, 547, 560).
#def _convertConvergence(data):
#    if not data.pf.parameters.has_key('cosmology_calculator'):
#        data.pf.parameters['cosmology_calculator'] = Cosmology(
#            HubbleConstantNow=(100.*data.pf.hubble_constant),
#            OmegaMatterNow=data.pf.omega_matter, OmegaLambdaNow=data.pf.omega_lambda)
#    # observer to lens
#    DL = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
#        data.pf.parameters['observer_redshift'], data.pf.current_redshift)
#    # observer to source
#    DS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
#        data.pf.parameters['observer_redshift'], data.pf.parameters['lensing_source_redshift'])
#    # lens to source
#    DLS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
#        data.pf.current_redshift, data.pf.parameters['lensing_source_redshift'])
#    # TODO: convert 1.5e14 to constants
#    return (((DL * DLS) / DS) * (1.5e14 * data.pf.omega_matter *
#                                (data.pf.hubble_constant / speed_of_light_cgs)**2 *
#                                (1 + data.pf.current_redshift)))
#add_field("weak_lensing_convergence", function=_overdensity,
#          convert_function=_convertConvergence,
#          projection_conversion='mpccm')

def _cell_volume(field, data):
    if data["dx"].size == 1:
        try:
            return ( data["dx"] * data["dy"] * data["dx"]
                     * np.ones(data.ActiveDimensions, dtype=np.float64) )
        except AttributeError:
            return data["dx"] * data["dy"] * data["dx"]
    vol = (data["dx"] * data["dy"] * data["dz"])
    vol.convert_to_cgs()
    return vol

add_field("cell_volume", units="cm**3", function=_cell_volume)

### Begin block that should probably be in an analysis module ###

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

#add_field("chandra_emissivity", function=_chandra_emissivity,
#          convert_function=_convert_chandra_emissivity,
#          projection_conversion="1")

def _xray_emissivity(field, data):
    return ( data["density"].astype(np.float64)**2
             * data["temperature"]**0.5 )

def _convert_xray_emissivity(data):
    return 2.168e60  #TODO: cnvert me to constants

#add_field("xray_emissivity", function=_xray_emissivity,
#          convert_function=_convert_xray_emissivity,
#          projection_conversion="1")

def _sz_kinetic(field, data):
    vel_axis = data.get_field_parameter("axis")
    if vel_axis > 2:
        raise NeedsParameter(["axis"])
    vel = data["velocity_%s" % ({0: "x", 1: "y", 2: "z"}[vel_axis])]
    return (vel * data["density"])

def _convert_sz_kinetic(data):
    return 0.88 * sigma_thompson / mh / clight

#add_field("sz_kinetic", function=_sz_kinetic,
#          convert_function=_convert_sz_kinetic,
#          validators=[ValidateParameter("axis")])

def _szy(field, data):
    return data["density"] * data["temperature"]

def _convert_szy(data):
    conv = 0.88 / mh * kboltz / (me * clight*clight) * sigma_thompson
    return conv

#add_field("szy", function=_szy, convert_function=_convert_szy)

### End block that should probably be in an analysis module ###

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
    new_field2 = YTArray(np.zeros((nx, ny, nz)), 'g/cm**3')
    new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
    return new_field2

add_field("averaged_density", function=_averaged_density,
          validators=[ValidateSpatial(1, ["density"])])

def _velocity_divergence(field, data):
    # We need to set up stencils.
    # This is based on enzo parameters and should probably be changed.
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    ds = div_fac * just_one(data["dx"])
    f  = data["x-velocity"][sl_right,1:-1,1:-1]/ds
    f -= data["x-velocity"][sl_left ,1:-1,1:-1]/ds
    if data.pf.dimensionality > 1:
        ds = div_fac * just_one(data["dy"])
        f += data["y-velocity"][1:-1,sl_right,1:-1]/ds
        f -= data["y-velocity"][1:-1,sl_left ,1:-1]/ds
    if data.pf.dimensionality > 2:
        ds = div_fac * just_one(data["dz"])
        f += data["z-velocity"][1:-1,1:-1,sl_right]/ds
        f -= data["z-velocity"][1:-1,1:-1,sl_left ]/ds
    new_field = YTArray(np.zeros(data["x-velocity"].shape,
                                 dtype=np.float64), '1/s')
    new_field[1:-1,1:-1,1:-1] = f
    return new_field

add_field("velocity_divergence", function=_velocity_divergence,
          validators=[ValidateSpatial(1, ["x-velocity", "y-velocity",
                                          "z-velocity"])],
          units="1/s", take_log=False)

def _velocity_divergence_absolute(field, data):
    return np.abs(data["velocity_divergence"])

add_field("velocity_divergence_absolute", function=_velocity_divergence_absolute, units="1/s")

def _contours(field, data):
    return -YTArray(np.ones_like(data["ones"]))

add_field("contours", validators=[ValidateSpatial(0)], take_log=False,
          display_field=False, function=_contours)
add_field("temp_contours", function=_contours,
          validators=[ValidateSpatial(0), ValidateGridType()],
          take_log=False, display_field=False)

def obtain_velocities(data):
    return obtain_rv_vec(data)

def _specific_angular_momentum_x(field, data):
    xv, yv, zv = obtain_velocities(data)
    center = data.get_field_parameter('center')
    v_vec = obtain_rvec(data)
    v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
    rv = v_vec - center
    return yv * rv[...,2] - zv * rv[...,1]

def _specific_angular_momentum_y(field, data):
    xv, yv, zv = obtain_velocities(data)
    center = data.get_field_parameter('center')
    v_vec = obtain_rvec(data)
    v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
    rv = v_vec - center
    return - (xv * rv[...,2] - zv * rv[...,0])

def _specific_angular_momentum_z(field, data):
    xv, yv, zv = obtain_velocities(data)
    center = data.get_field_parameter('center')
    v_vec = obtain_rvec(data)
    v_vec = np.rollaxis(v_vec, 0, len(v_vec.shape))
    rv = v_vec - center
    return xv * rv[...,1] - yv * rv[...,0]

add_field("specific_angular_momentum_x", function=_specific_angular_momentum_x,
          units="cm**2/s", validators=[ValidateParameter("center")])
add_field("specific_angular_momentum_y", function=_specific_angular_momentum_y,
          units="cm**2/s", validators=[ValidateParameter("center")])
add_field("specific_angular_momentum_z", function=_specific_angular_momentum_z,
          units="cm**2/s", validators=[ValidateParameter("center")])

def _angular_momentum_x(field, data):
    return data["cell_mass"] * data["specific_angular_momentum_x"]
add_field("angular_momentum_x", function=_angular_momentum_x,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])
def _angular_momentum_y(field, data):
    return data["cell_mass"] * data["specific_angular_momentum_y"]
add_field("angular_momentum_y", function=_angular_momentum_y,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])
def _angular_momentum_z(field, data):
    return data["cell_mass"] * data["specific_angular_momentum_z"]
add_field("angular_momentum_z", function=_angular_momentum_z,
         units="g * cm**2 / s", vector_field=False,
         validators=[ValidateParameter('center')])

def get_radius(data, field_prefix):
    center = data.get_field_parameter("center")
    DW = data.pf.domain_right_edge - data.pf.domain_left_edge
    radius = YTArray(np.zeros(data[field_prefix+"x"].shape, dtype='float64'),
                     'cm')
    r = radius.copy()
    if any(data.pf.periodicity):
        rdw = radius.copy()
    for i, ax in enumerate('xyz'):
        np.subtract(data["%s%s" % (field_prefix, ax)],
                    YTArray(center[i], center.units), r)
        if data.pf.periodicity[i] == True:
            np.abs(r, r)
            np.subtract(r, DW[i], rdw)
            np.abs(rdw, rdw)
            np.minimum(r, rdw, r)
        np.power(r, 2.0, r)
        np.add(radius, r, radius)
        if data.pf.dimensionality < i+1:
            break
    np.sqrt(radius, radius)
    return radius

def _radius(field, data):
    return get_radius(data, "")

add_field("radius", function=_radius,
          validators=[ValidateParameter("center")],
          units="cm")

def _radial_velocity(field, data):
    normal = data.get_field_parameter("normal")
    velocities = obtain_rv_vec(data)
    theta = data['spherical_theta'].copy()
    phi = data['spherical_phi'].copy()
    return get_sph_r_component(velocities, theta, phi, normal)

def _radial_velocity_absolute(field, data):
    return np.abs(_radial_velocity(field, data))
add_field("radial_velocity", function=_radial_velocity,
          units="cm/s")
add_field("radial_velocity_absolute", function=_radial_velocity_absolute,
          units="cm/s")

def _tangential_velocity(field, data):
    return np.sqrt(data["velocity_magnitude"]**2.0
                 - data["radial_velocity"]**2.0)
add_field("tangential_velocity",
          function=_tangential_velocity,
          take_log=False, units="cm/s")

def _cutting_plane_velocity_x(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bv = data.get_field_parameter("bulk_velocity")
    if bv == None: bv = np.zeros(3)
    tr  = (data["x-velocity"] - bv[0]) * x_vec[0]
    tr += (data["y-velocity"] - bv[1]) * x_vec[1]
    tr += (data["z-velocity"] - bv[2]) * x_vec[2]
    return tr
add_field("cutting_plane_velocity_x",
          function=_cutting_plane_velocity_x,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="km/s")
def _cutting_plane_velocity_y(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    bv = data.get_field_parameter("bulk_velocity")
    if bv == None: bv = np.zeros(3)
    tr  = (data["x-velocity"] - bv[0]) * y_vec[0]
    tr += (data["y-velocity"] - bv[1]) * y_vec[1]
    tr += (data["z-velocity"] - bv[2]) * y_vec[2]
    return tr
add_field("cutting_plane_velocity_y",
          function=_cutting_plane_velocity_y,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="km/s")

def _cutting_plane_magnetic_field_x(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    tr  = data["magnetic_field_x"] * x_vec[0]
    tr += data["magnetic_field_y"] * x_vec[1]
    tr += data["magnetic_field_z"] * x_vec[2]
    return tr
add_field("cutting_plane_magnetic_field_x",
          function=_cutting_plane_magnetic_field_x,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="gauss")
def _cutting_plane_magnetic_field_y(field, data):
    x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                           for ax in 'xyz']
    tr  = data["magnetic_field_x"] * y_vec[0]
    tr += data["magnetic_field_y"] * y_vec[1]
    tr += data["magnetic_field_z"] * y_vec[2]
    return tr
add_field("cutting_plane_magnetic_field_y",
          function=_cutting_plane_magnetic_field_y,
          validators=[ValidateParameter("cp_%s_vec" % ax)
                      for ax in 'xyz'], units="gauss")

def _mean_molecular_weight(field,data):
    return (data["density"] / (mh *data["number_density"]))
add_field("mean_molecular_weight", function=_mean_molecular_weight, units=r"")

def _pdensity(field, data):
    pmass = data[('deposit','all_mass')]
    pmass /= data["cell_volume"]
    return field.apply_units(pmass)
add_field("particle_density", function=_pdensity,
          validators=[ValidateGridType()],
          units = "g / cm**3",
          display_name=r"\mathrm{Particle}\/\mathrm{Density}")

def _magnetic_energy(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    magnetic_energy field to deal with non-unitary \mu_0.
    """
    return (data["magnetic_field_x"]**2 +
            data["magnetic_field_y"]**2 +
            data["magnetic_field_z"]**2)/(8*np.pi)
add_field("magnetic_energy",function=_magnetic_energy,
          units="erg / cm**3",
          display_name=r"\rm{Magnetic}\/\rm{Energy}")

def _magnetic_field_magnitude(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    magnetic_field_magnitude field to deal with non-unitary \mu_0.
    """
    return np.sqrt((data["magnetic_field_x"]**2 +
                    data["magnetic_field_y"]**2 +
                    data["magnetic_field_z"]**2))
add_field("magnetic_field_magnitude",
          function=_magnetic_field_magnitude,
          display_name=r"|B|", units="gauss")

def _plasma_beta(field,data):
    """This assumes that your front end has provided Bx, By, Bz in
    units of Gauss. If you use MKS, make sure to write your own
    PlasmaBeta field to deal with non-unitary \mu_0.
    """
    return data['pressure']/data['magnetic_energy']
add_field("plasma_beta",
          function=_plasma_beta,
          display_name=r"\rm{Plasma}\/\beta", units="")

def _magnetic_pressure(field,data):
    return data['magnetic_energy']
add_field("magnetic_pressure",
          function=_magnetic_pressure,
          display_name=r"\rm{Magnetic}\/\rm{Pressure}",
          units="erg / cm**3")

def _magnetic_field_poloidal(field,data):
    normal = data.get_field_parameter("normal")
    d = data['magnetic_field_x']

    Bfields = YTArray([data['magnetic_field_x'],
                       data['magnetic_field_y'],
                       data['magnetic_field_z']],
                       d.units)
    
    theta = data['spherical_theta']
    phi   = data['spherical_phi']
    
    return get_sph_theta_component(Bfields, theta, phi, normal)

add_field("magnetic_field_poloidal", function=_magnetic_field_poloidal,
          units="gauss",
          validators=[ValidateParameter("normal")])

def _magnetic_field_toroidal(field,data):
    normal = data.get_field_parameter("normal")
    d = data['magnetic_field_x']

    Bfields = YTArray([data['magnetic_field_x'],
                       data['magnetic_field_y'],
                       data['magnetic_field_z']],
                       d.units)
    
    phi = data['spherical_phi']
    
    return get_sph_phi_component(Bfields, phi, normal)

add_field("magnetic_field_toroidal", function=_magnetic_field_toroidal,
          units="gauss",
          validators=[ValidateParameter("normal")])

def _magnetic_field_radial(field,data):
    normal = data.get_field_parameter("normal")
    d = data['magnetic_field_x']

    Bfields = YTArray([data['magnetic_field_x'],
                       data['magnetic_field_y'],
                       data['magnetic_field_z']],
                       d.units)
    
    theta = data['spherical_theta']
    phi   = data['spherical_phi']
    
    return get_sph_r_component(Bfields, theta, phi, normal)

add_field("magnetic_field_radial", function=_magnetic_field_radial,
          units="gauss",
          validators=[ValidateParameter("normal")])

def _vorticity_squared(field, data):
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
    new_field = YTArray(np.zeros(data["x-velocity"].shape), 'cm/s')
    dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
             data["z-velocity"][1:-1,sl_left,1:-1]) \
             / (div_fac*just_one(data["dy"]))
    dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
             data["y-velocity"][1:-1,1:-1,sl_left]) \
             / (div_fac*just_one(data["dz"]))
    new_field[1:-1,1:-1,1:-1] += (dvzdy - dvydz)**2.0
    del dvzdy, dvydz
    dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
             data["x-velocity"][1:-1,1:-1,sl_left]) \
             / (div_fac*just_one(data["dz"]))
    dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
             data["z-velocity"][sl_left,1:-1,1:-1]) \
             / (div_fac*just_one(data["dx"]))
    new_field[1:-1,1:-1,1:-1] += (dvxdz - dvzdx)**2.0
    del dvxdz, dvzdx
    dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
             data["y-velocity"][sl_left,1:-1,1:-1]) \
             / (div_fac*just_one(data["dx"]))
    dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
             data["x-velocity"][1:-1,sl_left,1:-1]) \
             / (div_fac*just_one(data["dy"]))
    new_field[1:-1,1:-1,1:-1] += (dvydx - dvxdy)**2.0
    del dvydx, dvxdy
    new_field = np.abs(new_field)
    return new_field

add_field("vorticity_squared", function=_vorticity_squared,
          validators=[ValidateSpatial(1,
              ["x-velocity","y-velocity","z-velocity"])],
          units="s**-2")

def _pressure_gradient_x(field, data):
    # We need to set up stencils
    # This is based on enzo parameters and should probably be changed.    
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["pressure"].shape, dtype=np.float64),
                        'dyne/cm**3')
    ds = div_fac * just_one(data['dx'])
    new_field[1:-1,1:-1,1:-1]  = data["pressure"][sl_right,1:-1,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["pressure"][sl_left ,1:-1,1:-1]/ds
    return new_field
def _pressure_gradient_y(field, data):
    # We need to set up stencils
    # This is based on enzo parameters and should probably be changed.    
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["pressure"].shape, dtype=np.float64),
                        'dyne/cm**3')
    ds = div_fac * just_one(data['dy'])
    new_field[1:-1,1:-1,1:-1]  = data["pressure"][1:-1,sl_right,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["pressure"][1:-1,sl_left ,1:-1]/ds
    return new_field
def _pressure_gradient_z(field, data):
    # We need to set up stencils
    # This is based on enzo parameters and should probably be changed.
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["pressure"].shape, dtype=np.float64),
                        'dyne/cm**3')
    ds = div_fac * just_one(data['dz'])
    new_field[1:-1,1:-1,1:-1]  = data["pressure"][1:-1,1:-1,sl_right]/ds
    new_field[1:-1,1:-1,1:-1] -= data["pressure"][1:-1,1:-1,sl_left ]/ds
    return new_field

for ax in 'xyz':
    n = "pressure_gradient_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1, ["pressure"])],
              units="dyne/cm**3")

def _pressure_gradient_magnitude(field, data):
    return np.sqrt(data["pressure_gradient_x"]**2 +
                   data["pressure_gradient_y"]**2 +
                   data["pressure_gradient_z"]**2)
add_field("pressure_gradient_magnitude", function=_pressure_gradient_magnitude,
          validators=[ValidateSpatial(1, ["pressure"])],
          units="dyne/cm**3")

def _density_gradient_x(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["density"].shape, dtype=np.float64),
                        'g/cm**4')
    ds = div_fac * just_one(data['dx'])
    new_field[1:-1,1:-1,1:-1]  = data["density"][sl_right,1:-1,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["density"][sl_left ,1:-1,1:-1]/ds
    return new_field
def _density_gradient_y(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["density"].shape, dtype=np.float64),
                        'g/cm**4')
    ds = div_fac * just_one(data['dy'])
    new_field[1:-1,1:-1,1:-1]  = data["density"][1:-1,sl_right,1:-1]/ds
    new_field[1:-1,1:-1,1:-1] -= data["density"][1:-1,sl_left ,1:-1]/ds
    return new_field
def _density_gradient_z(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = YTArray(np.zeros(data["density"].shape, dtype=np.float64),
                        'g/cm**4')
    ds = div_fac * just_one(data['dz'])
    new_field[1:-1,1:-1,1:-1]  = data["density"][1:-1,1:-1,sl_right]/ds
    new_field[1:-1,1:-1,1:-1] -= data["density"][1:-1,1:-1,sl_left ]/ds
    return new_field

for ax in 'xyz':
    n = "density_gradient_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1, ["density"])],
              units="g/cm**4")

def _density_gradient_magnitude(field, data):
    return np.sqrt(data["density_gradient_x"]**2 +
                   data["density_gradient_y"]**2 +
                   data["density_gradient_z"]**2)
add_field("density_gradient_magnitude", function=_density_gradient_magnitude,
          validators=[ValidateSpatial(1, ["density"])],
          units="g/cm**4")

def _baroclinic_vorticity_x(field, data):
    rho2 = data["density"].astype(np.float64)**2
    return (data["pressure_gradient_y"] * data["density_gradient_z"] -
            data["pressure_gradient_z"] * data["density_gradient_z"]) / rho2
def _baroclinic_vorticity_y(field, data):
    rho2 = data["density"].astype(np.float64)**2
    return (data["pressure_gradient_z"] * data["density_gradient_x"] -
            data["pressure_gradient_x"] * data["density_gradient_z"]) / rho2
def _baroclinic_vorticity_z(field, data):
    rho2 = data["density"].astype(np.float64)**2
    return (data["pressure_gradient_x"] * data["density_gradient_y"] -
            data["pressure_gradient_y"] * data["density_gradient_x"]) / rho2
for ax in 'xyz':
    n = "baroclinic_vorticity_%s" % ax
    add_field(n, function=eval("_%s" % n),
          validators=[ValidateSpatial(1, ["density", "pressure"])],
          units="1/s")

def _baroclinic_vorticity_Magnitude(field, data):
    return np.sqrt(data["baroclinic_vorticity_x"]**2 +
                   data["baroclinic_vorticity_y"]**2 +
                   data["baroclinic_vorticity_z"]**2)
add_field("baroclinic_vorticity_Magnitude",
          function=_baroclinic_vorticity_Magnitude,
          validators=[ValidateSpatial(1, ["density", "pressure"])],
          units="1/s")

def _vorticity_x(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = \
      YTArray(np.zeros(data["z-velocity"].shape, dtype=np.float64), '1/s')
    new_field[1:-1,1:-1,1:-1] = (data["z-velocity"][1:-1,sl_right,1:-1] -
                                 data["z-velocity"][1:-1,sl_left,1:-1]) \
                                 / (div_fac*just_one(data["dy"]))
    new_field[1:-1,1:-1,1:-1] -= (data["y-velocity"][1:-1,1:-1,sl_right] -
                                  data["y-velocity"][1:-1,1:-1,sl_left]) \
                                  / (div_fac*just_one(data["dz"]))
    return new_field
def _vorticity_y(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = \
      YTArray(np.zeros(data["z-velocity"].shape, dtype=np.float64), '1/s')
    new_field[1:-1,1:-1,1:-1] = (data["x-velocity"][1:-1,1:-1,sl_right] -
                                 data["x-velocity"][1:-1,1:-1,sl_left]) \
                                 / (div_fac*just_one(data["dz"]))
    new_field[1:-1,1:-1,1:-1] -= (data["z-velocity"][sl_right,1:-1,1:-1] -
                                  data["z-velocity"][sl_left,1:-1,1:-1]) \
                                  / (div_fac*just_one(data["dx"]))
    return new_field
def _vorticity_z(field, data):
    # We need to set up stencils
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    new_field = \
      YTArray(np.zeros(data["z-velocity"].shape, dtype=np.float64), '1/s')
    new_field[1:-1,1:-1,1:-1] = (data["y-velocity"][sl_right,1:-1,1:-1] -
                                 data["y-velocity"][sl_left,1:-1,1:-1]) \
                                 / (div_fac*just_one(data["dx"]))
    new_field[1:-1,1:-1,1:-1] -= (data["x-velocity"][1:-1,sl_right,1:-1] -
                                  data["x-velocity"][1:-1,sl_left,1:-1]) \
                                  / (div_fac*just_one(data["dy"]))
    return new_field

for ax in 'xyz':
    n = "vorticity_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1,
                          ["x-velocity", "y-velocity", "z-velocity"])],
              units="1/s")

def _vorticity_magnitude(field, data):
    return np.sqrt(data["vorticity_x"]**2 +
                   data["vorticity_y"]**2 +
                   data["vorticity_z"]**2)
add_field("vorticity_magnitude", function=_vorticity_magnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

def _vorticity_stretching_x(field, data):
    return data["velocity_divergence"] * data["vorticity_x"]
def _vorticity_stretching_y(field, data):
    return data["velocity_divergence"] * data["vorticity_y"]
def _vorticity_stretching_z(field, data):
    return data["velocity_divergence"] * data["vorticity_z"]
for ax in 'xyz':
    n = "vorticity_stretching_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(0)])
def _vorticity_stretching_magnitude(field, data):
    return np.sqrt(data["vorticity_stretching_x"]**2 +
                   data["vorticity_stretching_y"]**2 +
                   data["vorticity_stretching_z"]**2)
add_field("vorticity_stretching_magnitude",
          function=_vorticity_stretching_magnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

def _vorticity_growth_x(field, data):
    return -data["vorticity_stretching_x"] - data["baroclinic_vorticity_x"]
def _vorticity_growth_y(field, data):
    return -data["vorticity_stretching_y"] - data["baroclinic_vorticity_y"]
def _vorticity_growth_z(field, data):
    return -data["vorticity_stretching_z"] - data["baroclinic_vorticity_z"]
for ax in 'xyz':
    n = "vorticity_growth_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1,
                          ["x-velocity", "y-velocity", "z-velocity"])],
              units="1/s")
def _vorticity_growth_magnitude(field, data):
    result = np.sqrt(data["vorticity_growth_x"]**2 +
                     data["vorticity_growth_y"]**2 +
                     data["vorticity_growth_z"]**2)
    dot = YTArray(np.zeros(result.shape), '1/s')
    for ax in "xyz":
        dot += data["vorticity_%s" % ax] * data["vorticity_growth_%s" % ax]
    result = np.sign(dot) * result
    return result
add_field("vorticity_growth_magnitude", function=_vorticity_growth_magnitude,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s",
          take_log=False)
def _vorticity_growth_magnitude_absolute(field, data):
    return np.sqrt(data["vorticity_growth_x"]**2 +
                   data["vorticity_growth_y"]**2 +
                   data["vorticity_growth_z"]**2)
add_field("vorticity_growth_magnitude_absolute", function=_vorticity_growth_magnitude_absolute,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="1/s")

def _vorticity_growth_timescale(field, data):
    domegax_dt = data["vorticity_x"] / data["vorticity_growth_x"]
    domegay_dt = data["vorticity_y"] / data["vorticity_growth_y"]
    domegaz_dt = data["vorticity_z"] / data["vorticity_growth_z"]
    return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
add_field("vorticity_growth_timescale", function=_vorticity_growth_timescale,
          validators=[ValidateSpatial(1,
                      ["x-velocity", "y-velocity", "z-velocity"])],
          units="s")

########################################################################
# With radiation pressure
########################################################################

def _vorticity_radiation_pressure_x(field, data):
    rho = data["density"].astype(np.float64)
    return (data["radiation_acceleration_y"] * data["density_gradient_z"] -
            data["radiation_acceleration_z"] * data["density_gradient_y"]) / rho
def _vorticity_radiation_pressure_y(field, data):
    rho = data["density"].astype(np.float64)
    return (data["radiation_acceleration_z"] * data["density_gradient_x"] -
            data["radiation_acceleration_x"] * data["density_gradient_z"]) / rho
def _vorticity_radiation_pressure_z(field, data):
    rho = data["density"].astype(np.float64)
    return (data["radiation_acceleration_x"] * data["density_gradient_y"] -
            data["radiation_acceleration_y"] * data["density_gradient_x"]) / rho

for ax in 'xyz':
    n = "vorticity_radiation_pressure_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1,
                   ["density",
                    "radiation_acceleration_x",
                    "radiation_acceleration_y",
                    "radiation_acceleration_z"])],
              units="1/s")

def _vorticity_radiation_pressure_magnitude(field, data):
    return np.sqrt(data["vorticity_radiation_pressure_x"]**2 +
                   data["vorticity_radiation_pressure_y"]**2 +
                   data["vorticity_radiation_pressure_z"]**2)
add_field("vorticity_radiation_pressure_magnitude",
          function=_vorticity_radiation_pressure_magnitude,
          validators=[ValidateSpatial(1,
                      ["density",
                       "radiation_acceleration_x",
                       "radiation_acceleration_y",
                       "radiation_acceleration_z"])],
          units="1/s")

def _vorticity_radiation_pressure_growth_x(field, data):
    return -data["vorticity_stretching_x"] - data["baroclinic_vorticity_x"] \
           -data["vorticity_radiation_pressure_x"]
def _vorticity_radiation_pressure_growth_y(field, data):
    return -data["vorticity_stretching_y"] - data["baroclinic_vorticity_y"] \
           -data["vorticity_radiation_pressure_y"]
def _vorticity_radiation_pressure_growth_z(field, data):
    return -data["vorticity_stretching_z"] - data["baroclinic_vorticity_z"] \
           -data["vorticity_radiation_pressure_z"]
for ax in 'xyz':
    n = "vorticity_radiation_pressure_growth_%s" % ax
    add_field(n, function=eval("_%s" % n),
              validators=[ValidateSpatial(1,
                       ["density",
                        "radiation_acceleration_x",
                        "radiation_acceleration_y",
                        "radiation_acceleration_z"])],
              units="1/s")
def _vorticity_radiation_pressure_growth_magnitude(field, data):
    result = np.sqrt(data["vorticity_radiation_pressure_growth_x"]**2 +
                     data["vorticity_radiation_pressure_growth_y"]**2 +
                     data["vorticity_radiation_pressure_growth_z"]**2)
    dot = YTArray(np.zeros(result.shape), '1/s')
    for ax in "xyz":
        dot += data["Vorticity%s" % ax] * data["vorticity_growth_%s" % ax]
    result = np.sign(dot) * result
    return result
add_field("vorticity_radiation_pressure_growth_magnitude", function=_vorticity_growth_magnitude,
          validators=[ValidateSpatial(1,
                      ["density", "radiation_acceleration_x", "radiation_acceleration_y", "radiation_acceleration_z"])],
          units="1/s",
          take_log=False)
def _vorticity_radiation_pressure_growth_magnitude_absolute(field, data):
    return np.sqrt(data["vorticity_radiation_pressure_growth_x"]**2 +
                   data["vorticity_radiation_pressure_growth_y"]**2 +
                   data["vorticity_radiation_pressure_growth_z"]**2)
add_field("vorticity_radiation_pressure_growth_magnitude_absolute",
          function=_vorticity_radiation_pressure_growth_magnitude_absolute,
          validators=[ValidateSpatial(1,
                      ["density",
                       "radiation_acceleration_x",
                       "radiation_acceleration_y",
                       "radiation_acceleration_z"])],
          units="1/s")

def _vorticity_radiation_pressure_growth_timescale(field, data):
    domegax_dt = data["vorticity_x"] / data["vorticity_radiation_pressure_growth_x"]
    domegay_dt = data["vorticity_y"] / data["vorticity_radiation_pressure_growth_y"]
    domegaz_dt = data["vorticity_z"] / data["vorticity_radiation_pressure_growth_z"]
    return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
add_field("vorticity_radiation_pressure_growth_timescale", function=_vorticity_radiation_pressure_growth_timescale,
          validators=[ValidateSpatial(1,
                      ["density",
                       "radiation_acceleration_x",
                       "radiation_acceleration_y",
                       "radiation_acceleration_z"])],
          units="1/s")
