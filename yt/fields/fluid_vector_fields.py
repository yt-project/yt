"""
Complex fluid fields.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.fields.derived_field import \
    ValidateSpatial

from .field_plugin_registry import \
    register_field_plugin

from yt.funcs import \
    just_one

from .vector_operations import \
    create_magnitude_field, \
    create_squared_field

@register_field_plugin
def setup_fluid_vector_fields(registry, ftype = "gas", slice_info = None):
    unit_system = registry.ds.unit_system
    # slice_info would be the left, the right, and the factor.
    # For example, with the old Enzo-ZEUS fields, this would be:
    # slice(None, -2, None)
    # slice(1, -1, None)
    # 1.0
    # Otherwise, we default to a centered difference.
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info
    sl_center = slice(1, -1, None)

    def _baroclinic_vorticity_x(field, data):
        rho2 = data[ftype, "density"].astype(np.float64)**2
        return (data[ftype, "pressure_gradient_y"] *
                data[ftype, "density_gradient_z"] -
                data[ftype, "pressure_gradient_z"] *
                data[ftype, "density_gradient_z"]) / rho2
    def _baroclinic_vorticity_y(field, data):
        rho2 = data[ftype, "density"].astype(np.float64)**2
        return (data[ftype, "pressure_gradient_z"] *
                data[ftype, "density_gradient_x"] -
                data[ftype, "pressure_gradient_x"] *
                data[ftype, "density_gradient_z"]) / rho2
    def _baroclinic_vorticity_z(field, data):
        rho2 = data[ftype, "density"].astype(np.float64)**2
        return (data[ftype, "pressure_gradient_x"] *
                data[ftype, "density_gradient_y"] -
                data[ftype, "pressure_gradient_y"] *
                data[ftype, "density_gradient_x"]) / rho2

    bv_validators = [ValidateSpatial(1, [(ftype, "density"), (ftype, "pressure")])]
    for ax in 'xyz':
        n = "baroclinic_vorticity_%s" % ax
        registry.add_field((ftype, n), function=eval("_%s" % n),
                           validators=bv_validators,
                           units=unit_system["frequency"]**2)

    create_magnitude_field(registry, "baroclinic_vorticity", 
                           unit_system["frequency"]**2,
                           ftype=ftype, slice_info=slice_info,
                           validators=bv_validators)

    def _vorticity_x(field, data):
        f  = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
              data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
              / (div_fac*just_one(data["index", "dy"]))
        f -= (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
              data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
              / (div_fac*just_one(data["index", "dz"]))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _vorticity_y(field, data):
        f  = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
              data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
              / (div_fac*just_one(data["index", "dz"]))
        f -= (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
              data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
              / (div_fac*just_one(data["index", "dx"]))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _vorticity_z(field, data):
        f  = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
              data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
              / (div_fac*just_one(data["index", "dx"]))
        f -= (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
              data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
              / (div_fac*just_one(data["index", "dy"]))
        new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                              dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field

    vort_validators = [ValidateSpatial(1,
                            [(ftype, "velocity_x"),
                             (ftype, "velocity_y"),
                             (ftype, "velocity_z")])]
    for ax in 'xyz':
        n = "vorticity_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units=unit_system["frequency"],
                           validators=vort_validators)
    create_magnitude_field(registry, "vorticity", unit_system["frequency"],
                           ftype=ftype, slice_info=slice_info,
                           validators=vort_validators)
    create_squared_field(registry, "vorticity", unit_system["frequency"]**2,
                         ftype=ftype, slice_info=slice_info,
                         validators=vort_validators)

    def _vorticity_stretching_x(field, data):
        return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_x"]
    def _vorticity_stretching_y(field, data):
        return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_y"]
    def _vorticity_stretching_z(field, data):
        return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_z"]
    for ax in 'xyz':
        n = "vorticity_stretching_%s" % ax
        registry.add_field((ftype, n), 
                           function=eval("_%s" % n),
                           units = unit_system["frequency"]**2,
                           validators=vort_validators)

    create_magnitude_field(registry, "vorticity_stretching", 
                           unit_system["frequency"]**2,
                           ftype=ftype, slice_info=slice_info,
                           validators=vort_validators)

    def _vorticity_growth_x(field, data):
        return -data[ftype, "vorticity_stretching_x"] - \
          data[ftype, "baroclinic_vorticity_x"]
    def _vorticity_growth_y(field, data):
        return -data[ftype, "vorticity_stretching_y"] - \
          data[ftype, "baroclinic_vorticity_y"]
    def _vorticity_growth_z(field, data):
        return -data[ftype, "vorticity_stretching_z"] - \
          data[ftype, "baroclinic_vorticity_z"]
    for ax in 'xyz':
        n = "vorticity_growth_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units=unit_system["frequency"]**2,
                           validators=vort_validators)

    def _vorticity_growth_magnitude(field, data):
        result = np.sqrt(data[ftype, "vorticity_growth_x"]**2 +
                         data[ftype, "vorticity_growth_y"]**2 +
                         data[ftype, "vorticity_growth_z"]**2)
        dot = data.ds.arr(np.zeros(result.shape), "")
        for ax in "xyz":
            dot += (data[ftype, "vorticity_%s" % ax] *
                    data[ftype, "vorticity_growth_%s" % ax]).to_ndarray()
        result = np.sign(dot) * result
        return result
    
    registry.add_field((ftype, "vorticity_growth_magnitude"),
              function=_vorticity_growth_magnitude,
              units=unit_system["frequency"]**2,
              validators=vort_validators,
              take_log=False)

    def _vorticity_growth_magnitude_absolute(field, data):
        return np.sqrt(data[ftype, "vorticity_growth_x"]**2 +
                       data[ftype, "vorticity_growth_y"]**2 +
                       data[ftype, "vorticity_growth_z"]**2)
    
    registry.add_field((ftype, "vorticity_growth_magnitude_absolute"),
                       function=_vorticity_growth_magnitude_absolute,
                       units=unit_system["frequency"]**2,
                       validators=vort_validators)

    def _vorticity_growth_timescale(field, data):
        domegax_dt = data[ftype, "vorticity_x"] / data[ftype, "vorticity_growth_x"]
        domegay_dt = data[ftype, "vorticity_y"] / data[ftype, "vorticity_growth_y"]
        domegaz_dt = data[ftype, "vorticity_z"] / data[ftype, "vorticity_growth_z"]
        return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
    
    registry.add_field((ftype, "vorticity_growth_timescale"),
                       function=_vorticity_growth_timescale,
                       units=unit_system["time"],
                       validators=vort_validators)

    ########################################################################
    # With radiation pressure
    ########################################################################

    def _vorticity_radiation_pressure_x(field, data):
        rho = data[ftype, "density"].astype(np.float64)
        return (data[ftype, "radiation_acceleration_y"] *
                data[ftype, "density_gradient_z"] -
                data[ftype, "radiation_acceleration_z"] *
                data[ftype, "density_gradient_y"]) / rho
    def _vorticity_radiation_pressure_y(field, data):
        rho = data[ftype, "density"].astype(np.float64)
        return (data[ftype, "radiation_acceleration_z"] *
                data[ftype, "density_gradient_x"] -
                data[ftype, "radiation_acceleration_x"] *
                data[ftype, "density_gradient_z"]) / rho
    def _vorticity_radiation_pressure_z(field, data):
        rho = data[ftype, "density"].astype(np.float64)
        return (data[ftype, "radiation_acceleration_x"] *
                data[ftype, "density_gradient_y"] -
                data[ftype, "radiation_acceleration_y"] *
                data[ftype, "density_gradient_x"]) / rho

    vrp_validators = [ValidateSpatial(1,
                            [(ftype, "density"),
                             (ftype, "radiation_acceleration_x"),
                             (ftype, "radiation_acceleration_y"),
                             (ftype, "radiation_acceleration_z")])]
    for ax in 'xyz':
        n = "vorticity_radiation_pressure_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units=unit_system["frequency"]**2,
                           validators=vrp_validators)

    create_magnitude_field(registry, "vorticity_radiation_pressure",
                           unit_system["frequency"]**2,
                           ftype=ftype, slice_info=slice_info,
                           validators=vrp_validators)

    def _vorticity_radiation_pressure_growth_x(field, data):
        return -data[ftype, "vorticity_stretching_x"] - \
          data[ftype, "baroclinic_vorticity_x"] \
          -data[ftype, "vorticity_radiation_pressure_x"]
    def _vorticity_radiation_pressure_growth_y(field, data):
        return -data[ftype, "vorticity_stretching_y"] - \
          data[ftype, "baroclinic_vorticity_y"] \
          -data[ftype, "vorticity_radiation_pressure_y"]
    def _vorticity_radiation_pressure_growth_z(field, data):
        return -data[ftype, "vorticity_stretching_z"] - \
          data[ftype, "baroclinic_vorticity_z"] \
          -data[ftype, "vorticity_radiation_pressure_z"]
               
    for ax in 'xyz':
        n = "vorticity_radiation_pressure_growth_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units=unit_system["frequency"]**2,
                           validators=vrp_validators)

    def _vorticity_radiation_pressure_growth_magnitude(field, data):
        result = np.sqrt(data[ftype, "vorticity_radiation_pressure_growth_x"]**2 +
                         data[ftype, "vorticity_radiation_pressure_growth_y"]**2 +
                         data[ftype, "vorticity_radiation_pressure_growth_z"]**2)
        dot = data.ds.arr(np.zeros(result.shape), "")
        for ax in "xyz":
            dot += (data[ftype, "vorticity_%s" % ax] *
                    data[ftype, "vorticity_growth_%s" % ax]).to_ndarray()
        result = np.sign(dot) * result
        return result
    
    registry.add_field((ftype, "vorticity_radiation_pressure_growth_magnitude"),
                       function=_vorticity_radiation_pressure_growth_magnitude,
                       units=unit_system["frequency"]**2,
                       validators=vrp_validators,
                       take_log=False)

    def _vorticity_radiation_pressure_growth_magnitude_absolute(field, data):
        return np.sqrt(data[ftype, "vorticity_radiation_pressure_growth_x"]**2 +
                       data[ftype, "vorticity_radiation_pressure_growth_y"]**2 +
                       data[ftype, "vorticity_radiation_pressure_growth_z"]**2)
    
    registry.add_field((ftype, "vorticity_radiation_pressure_growth_magnitude_absolute"),
                       function=_vorticity_radiation_pressure_growth_magnitude_absolute,
                       units="s**(-2)",
                       validators=vrp_validators)

    def _vorticity_radiation_pressure_growth_timescale(field, data):
        domegax_dt = data[ftype, "vorticity_x"] / \
          data[ftype, "vorticity_radiation_pressure_growth_x"]
        domegay_dt = data[ftype, "vorticity_y"] / \
          data[ftype, "vorticity_radiation_pressure_growth_y"]
        domegaz_dt = data[ftype, "vorticity_z"] / \
          data[ftype, "vorticity_radiation_pressure_growth_z"]
        return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
    
    registry.add_field((ftype, "vorticity_radiation_pressure_growth_timescale"),
                       function=_vorticity_radiation_pressure_growth_timescale,
                       units=unit_system["time"],
                       validators=vrp_validators)

    def _shear(field, data):
        """
        Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
                             (dvx/dz + dvz/dx)^2 ]^(0.5)
        where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
        and is in units of s^(-1)
        (it's just like vorticity except add the derivative pairs instead
         of subtracting them)
        """
        
        if data.ds.dimensionality > 1:
            dvydx = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
                    data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
                    / (div_fac*just_one(data["index", "dx"]))
            dvxdy = (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
                    data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
                    / (div_fac*just_one(data["index", "dy"]))
            f  = (dvydx + dvxdy)**2.0
            del dvydx, dvxdy
        if data.ds.dimensionality > 2:
            dvzdy = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
                    data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
                    / (div_fac*just_one(data["index", "dy"]))
            dvydz = (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
                    data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
                    / (div_fac*just_one(data["index", "dz"]))
            f += (dvzdy + dvydz)**2.0
            del dvzdy, dvydz
            dvxdz = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
                    data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
                    / (div_fac*just_one(data["index", "dz"]))
            dvzdx = (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
                    data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
                    / (div_fac*just_one(data["index", "dx"]))
            f += (dvxdz + dvzdx)**2.0
            del dvxdz, dvzdx
        np.sqrt(f, out=f)
        new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_x"]), f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    
    registry.add_field((ftype, "shear"),
                       function=_shear,
                       validators=[ValidateSpatial(1,
                        [(ftype, "velocity_x"),
                         (ftype, "velocity_y"),
                         (ftype, "velocity_z")])],
                        units=unit_system["frequency"])

    def _shear_criterion(field, data):
        """
        Divide by c_s to leave shear in units of length**-1, which 
        can be compared against the inverse of the local cell size (1/dx) 
        to determine if refinement should occur.
        """
        
        return data[ftype, "shear"] / data[ftype, "sound_speed"]

    registry.add_field((ftype, "shear_criterion"),
                       function=_shear_criterion,
                       units=unit_system["length"]**-1,
                       validators=[ValidateSpatial(1,
                        [(ftype, "sound_speed"),
                         (ftype, "velocity_x"),
                         (ftype, "velocity_y"),
                         (ftype, "velocity_z")])])

    def _shear_mach(field, data):
        """
        Dimensionless shear (shear_mach) is defined nearly the same as shear, 
        except that it is scaled by the local dx/dy/dz and the local sound speed.
        So it results in a unitless quantity that is effectively measuring 
        shear in mach number.  

        In order to avoid discontinuities created by multiplying by dx/dy/dz at
        grid refinement boundaries, we also multiply by 2**GridLevel.

        Shear (Mach) = [(dvx + dvy)^2 + (dvz + dvy)^2 +
                        (dvx + dvz)^2  ]^(0.5) / c_sound
        """
        
        if data.ds.dimensionality > 1:
            dvydx = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
                     data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
                    / div_fac
            dvxdy = (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
                     data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
                    / div_fac
            f  = (dvydx + dvxdy)**2.0
            del dvydx, dvxdy
        if data.ds.dimensionality > 2:
            dvzdy = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
                     data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
                    / div_fac
            dvydz = (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
                     data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
                    / div_fac
            f += (dvzdy + dvydz)**2.0
            del dvzdy, dvydz
            dvxdz = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
                     data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
                    / div_fac
            dvzdx = (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
                     data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
                    / div_fac
            f += (dvxdz + dvzdx)**2.0
            del dvxdz, dvzdx
        f *= (2.0**data["index", "grid_level"][sl_center, sl_center, sl_center] /
              data[ftype, "sound_speed"][sl_center, sl_center, sl_center])**2.0
        np.sqrt(f, out=f)
        new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_x"]), f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    
    registry.add_field((ftype, "shear_mach"),
                       function=_shear_mach,
                       units="",
                       validators=[ValidateSpatial(1,
                        [(ftype, "sound_speed"),
                         (ftype, "velocity_x"),
                         (ftype, "velocity_y"),
                         (ftype, "velocity_z")])])
