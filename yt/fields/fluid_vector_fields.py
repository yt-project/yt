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
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsParameter

from .field_plugin_registry import \
    register_field_plugin

from yt.funcs import \
    just_one

from .vector_operations import \
     create_magnitude_field
    
@register_field_plugin
def setup_fluid_vector_fields(registry, ftype = "gas", slice_info = None):
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
        registry.add_field((ftype, n), function=eval("_%s" % n),
                           validators=[ValidateSpatial(1, ["density", "pressure"])],
                           units="s**(-2)")

    create_magnitude_field(registry, "baroclinic_vorticity", "s**(-2)",
                           ftype=ftype, slice_info=slice_info,
                           validators=[ValidateSpatial(1,
                            ["density", "pressure"])])

    def _vorticity_x(field, data):
        f  = (data["velocity_z"][sl_center,sl_right,sl_center] -
              data["velocity_z"][sl_center,sl_left,sl_center]) \
              / (div_fac*just_one(data["dy"]).in_cgs())
        f -= (data["velocity_y"][sl_center,sl_center,sl_right] -
              data["velocity_y"][sl_center,sl_center,sl_left]) \
              / (div_fac*just_one(data["dz"].in_cgs()))
        new_field = data.pf.arr(np.zeros_like(data["velocity_z"], dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _vorticity_y(field, data):
        f  = (data["velocity_x"][sl_center,sl_center,sl_right] -
              data["velocity_x"][sl_center,sl_center,sl_left]) \
              / (div_fac*just_one(data["dz"]))
        f -= (data["velocity_z"][sl_right,sl_center,sl_center] -
              data["velocity_z"][sl_left,sl_center,sl_center]) \
              / (div_fac*just_one(data["dx"]))
        new_field = data.pf.arr(np.zeros_like(data["velocity_z"], dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field
    def _vorticity_z(field, data):
        f  = (data["velocity_y"][sl_right,sl_center,sl_center] -
              data["velocity_y"][sl_left,sl_center,sl_center]) \
              / (div_fac*just_one(data["dx"]))
        f -= (data["velocity_x"][sl_center,sl_right,sl_center] -
              data["velocity_x"][sl_center,sl_left,sl_center]) \
              / (div_fac*just_one(data["dy"]))
        new_field = data.pf.arr(np.zeros_like(data["velocity_z"], dtype=np.float64),
                                f.units)
        new_field[sl_center, sl_center, sl_center] = f
        return new_field

    for ax in 'xyz':
        n = "vorticity_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units="1/s",
                           validators=[ValidateSpatial(1,
                            ["velocity_x", "velocity_y", "velocity_z"])])
    create_magnitude_field(registry, "vorticity", "1/s",
                           ftype=ftype, slice_info=slice_info,
                           validators=[ValidateSpatial(1,
                            ["velocity_x", "velocity_y", "velocity_z"])])

    def _vorticity_stretching_x(field, data):
        return data["velocity_divergence"] * data["vorticity_x"]
    def _vorticity_stretching_y(field, data):
        return data["velocity_divergence"] * data["vorticity_y"]
    def _vorticity_stretching_z(field, data):
        return data["velocity_divergence"] * data["vorticity_z"]
    for ax in 'xyz':
        n = "vorticity_stretching_%s" % ax
        registry.add_field((ftype, n), 
                           function=eval("_%s" % n),
                           units = "s**(-2)",
                           validators=[ValidateSpatial(1,
                            ["velocity_x", "velocity_y", "velocity_z"])])

    create_magnitude_field(registry, "vorticity_stretching", "s**(-2)",
                           ftype=ftype, slice_info=slice_info,
                           validators=[ValidateSpatial(1,
                            ["velocity_x", "velocity_y", "velocity_z"])])

    def _vorticity_growth_x(field, data):
        return -data["vorticity_stretching_x"] - data["baroclinic_vorticity_x"]
    def _vorticity_growth_y(field, data):
        return -data["vorticity_stretching_y"] - data["baroclinic_vorticity_y"]
    def _vorticity_growth_z(field, data):
        return -data["vorticity_stretching_z"] - data["baroclinic_vorticity_z"]
    for ax in 'xyz':
        n = "vorticity_growth_%s" % ax
        registry.add_field((ftype, n),
                           function=eval("_%s" % n),
                           units="s**(-2)",
                           validators=[ValidateSpatial(1,
                            ["velocity_x", "velocity_y", "velocity_z"])])

    def _vorticity_growth_magnitude(field, data):
        result = np.sqrt(data["vorticity_growth_x"]**2 +
                         data["vorticity_growth_y"]**2 +
                         data["vorticity_growth_z"]**2)
        dot = data.pf.arr(np.zeros(result.shape), "")
        for ax in "xyz":
            dot += (data["vorticity_%s" % ax] *
                    data["vorticity_growth_%s" % ax]).to_ndarray()
        result = np.sign(dot) * result
        return result
    
    registry.add_field((ftype, "vorticity_growth_magnitude"),
              function=_vorticity_growth_magnitude,
              units="s**(-2)",
              validators=[ValidateSpatial(1,
               ["velocity_x", "velocity_y", "velocity_z"])],
              take_log=False)

    def _vorticity_growth_magnitude_absolute(field, data):
        return np.sqrt(data["vorticity_growth_x"]**2 +
                       data["vorticity_growth_y"]**2 +
                       data["vorticity_growth_z"]**2)
    
    registry.add_field((ftype, "vorticity_growth_magnitude_absolute"),
                       function=_vorticity_growth_magnitude_absolute,
                       units="s**(-2)",
                       validators=[ValidateSpatial(1,
                        ["velocity_x", "velocity_y", "velocity_z"])])

    def _vorticity_growth_timescale(field, data):
        domegax_dt = data["vorticity_x"] / data["vorticity_growth_x"]
        domegay_dt = data["vorticity_y"] / data["vorticity_growth_y"]
        domegaz_dt = data["vorticity_z"] / data["vorticity_growth_z"]
        return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
    
    registry.add_field((ftype, "vorticity_growth_timescale"),
                       function=_vorticity_growth_timescale,
                       units="s",
                       validators=[ValidateSpatial(1,
                        ["velocity_x", "velocity_y", "velocity_z"])])

    # ########################################################################
    # # With radiation pressure
    # ########################################################################

    # def _vorticity_radiation_pressure_x(field, data):
    #     rho = data["density"].astype(np.float64)
    #     return (data["radiation_acceleration_y"] * data["density_gradient_z"] -
    #             data["radiation_acceleration_z"] * data["density_gradient_y"]) / rho
    # def _vorticity_radiation_pressure_y(field, data):
    #     rho = data["density"].astype(np.float64)
    #     return (data["radiation_acceleration_z"] * data["density_gradient_x"] -
    #             data["radiation_acceleration_x"] * data["density_gradient_z"]) / rho
    # def _vorticity_radiation_pressure_z(field, data):
    #     rho = data["density"].astype(np.float64)
    #     return (data["radiation_acceleration_x"] * data["density_gradient_y"] -
    #             data["radiation_acceleration_y"] * data["density_gradient_x"]) / rho

    # for ax in 'xyz':
    #     n = "vorticity_radiation_pressure_%s" % ax
    #     add_field(n, function=eval("_%s" % n),
    #               validators=[ValidateSpatial(1,
    #                    ["density",
    #                     "radiation_acceleration_x",
    #                     "radiation_acceleration_y",
    #                     "radiation_acceleration_z"])],
    #               units="1/s")

    # def _vorticity_radiation_pressure_magnitude(field, data):
    #     return np.sqrt(data["vorticity_radiation_pressure_x"]**2 +
    #                    data["vorticity_radiation_pressure_y"]**2 +
    #                    data["vorticity_radiation_pressure_z"]**2)
    # add_field("vorticity_radiation_pressure_magnitude",
    #           function=_vorticity_radiation_pressure_magnitude,
    #           validators=[ValidateSpatial(1,
    #                       ["density",
    #                        "radiation_acceleration_x",
    #                        "radiation_acceleration_y",
    #                        "radiation_acceleration_z"])],
    #           units="1/s")

    # def _vorticity_radiation_pressure_growth_x(field, data):
    #     return -data["vorticity_stretching_x"] - data["baroclinic_vorticity_x"] \
    #            -data["vorticity_radiation_pressure_x"]
    # def _vorticity_radiation_pressure_growth_y(field, data):
    #     return -data["vorticity_stretching_y"] - data["baroclinic_vorticity_y"] \
    #            -data["vorticity_radiation_pressure_y"]
    # def _vorticity_radiation_pressure_growth_z(field, data):
    #     return -data["vorticity_stretching_z"] - data["baroclinic_vorticity_z"] \
    #            -data["vorticity_radiation_pressure_z"]
    # for ax in 'xyz':
    #     n = "vorticity_radiation_pressure_growth_%s" % ax
    #     add_field(n, function=eval("_%s" % n),
    #               validators=[ValidateSpatial(1,
    #                        ["density",
    #                         "radiation_acceleration_x",
    #                         "radiation_acceleration_y",
    #                         "radiation_acceleration_z"])],
    #               units="1/s")
    # def _vorticity_radiation_pressure_growth_magnitude(field, data):
    #     result = np.sqrt(data["vorticity_radiation_pressure_growth_x"]**2 +
    #                      data["vorticity_radiation_pressure_growth_y"]**2 +
    #                      data["vorticity_radiation_pressure_growth_z"]**2)
    #     dot = data.pf.arr(np.zeros(result.shape), '1/s')
    #     for ax in "xyz":
    #         dot += data["Vorticity%s" % ax] * data["vorticity_growth_%s" % ax]
    #     result = np.sign(dot) * result
    #     return result
    # add_field("vorticity_radiation_pressure_growth_magnitude", function=_vorticity_growth_magnitude,
    #           validators=[ValidateSpatial(1,
    #                       ["density", "radiation_acceleration_x", "radiation_acceleration_y", "radiation_acceleration_z"])],
    #           units="1/s",
    #           take_log=False)
    # def _vorticity_radiation_pressure_growth_magnitude_absolute(field, data):
    #     return np.sqrt(data["vorticity_radiation_pressure_growth_x"]**2 +
    #                    data["vorticity_radiation_pressure_growth_y"]**2 +
    #                    data["vorticity_radiation_pressure_growth_z"]**2)
    # add_field("vorticity_radiation_pressure_growth_magnitude_absolute",
    #           function=_vorticity_radiation_pressure_growth_magnitude_absolute,
    #           validators=[ValidateSpatial(1,
    #                       ["density",
    #                        "radiation_acceleration_x",
    #                        "radiation_acceleration_y",
    #                        "radiation_acceleration_z"])],
    #           units="1/s")

    # def _vorticity_radiation_pressure_growth_timescale(field, data):
    #     domegax_dt = data["vorticity_x"] / data["vorticity_radiation_pressure_growth_x"]
    #     domegay_dt = data["vorticity_y"] / data["vorticity_radiation_pressure_growth_y"]
    #     domegaz_dt = data["vorticity_z"] / data["vorticity_radiation_pressure_growth_z"]
    #     return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
    # add_field("vorticity_radiation_pressure_growth_timescale", function=_vorticity_radiation_pressure_growth_timescale,
    #           validators=[ValidateSpatial(1,
    #                       ["density",
    #                        "radiation_acceleration_x",
    #                        "radiation_acceleration_y",
    #                        "radiation_acceleration_z"])],
    #           units="1/s")


    # def _Shear(field, data):
    #     """
    #     Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
    #                          (dvx/dz + dvz/dx)^2 ]^(0.5)
    #     where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
    #     and is in units of s^(-1)
    #     (it's just like vorticity except add the derivative pairs instead
    #      of subtracting them)
    #     """
    #     # We need to set up stencils
    #     if data.pf["HydroMethod"] == 2:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(1,-1,None)
    #         div_fac = 1.0
    #     else:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(2,None,None)
    #         div_fac = 2.0
    #     new_field = np.zeros(data["x-velocity"].shape)
    #     if data.pf.dimensionality > 1:
    #         dvydx = (data["y-velocity"][sl_right,sl_center,sl_center] -
    #                 data["y-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac*data["dx"].flat[0])
    #         dvxdy = (data["x-velocity"][sl_center,sl_right,sl_center] -
    #                 data["x-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac*data["dy"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvydx + dvxdy)**2.0
    #         del dvydx, dvxdy
    #     if data.pf.dimensionality > 2:
    #         dvzdy = (data["z-velocity"][sl_center,sl_right,sl_center] -
    #                 data["z-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac*data["dy"].flat[0])
    #         dvydz = (data["y-velocity"][sl_center,sl_center,sl_right] -
    #                 data["y-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac*data["dz"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvzdy + dvydz)**2.0
    #         del dvzdy, dvydz
    #         dvxdz = (data["x-velocity"][sl_center,sl_center,sl_right] -
    #                 data["x-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac*data["dz"].flat[0])
    #         dvzdx = (data["z-velocity"][sl_right,sl_center,sl_center] -
    #                 data["z-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac*data["dx"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvxdz + dvzdx)**2.0
    #         del dvxdz, dvzdx
    #     new_field = new_field**0.5
    #     new_field = np.abs(new_field)
    #     return new_field
    # def _convertShear(data):
    #     return data.convert("cm")**-1.0
    # #add_field("Shear", function=_Shear,
    # #          validators=[ValidateSpatial(1,
    # #              ["x-velocity","y-velocity","z-velocity"])],
    # #          units=r"\rm{s}^{-1}",
    # #          convert_function=_convertShear, take_log=False)

    # def _ShearCriterion(field, data):
    #     """
    #     Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
    #                          (dvx/dz + dvz/dx)^2 ]^(0.5)
    #     where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
    #     and is in units of s^(-1)
    #     (it's just like vorticity except add the derivative pairs instead
    #      of subtracting them)

    #     Divide by c_s to leave Shear in units of cm**-1, which 
    #     can be compared against the inverse of the local cell size (1/dx) 
    #     to determine if refinement should occur.
    #     """
    #     # We need to set up stencils
    #     if data.pf["HydroMethod"] == 2:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(1,-1,None)
    #         div_fac = 1.0
    #     else:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(2,None,None)
    #         div_fac = 2.0
    #     new_field = np.zeros(data["x-velocity"].shape)
    #     if data.pf.dimensionality > 1:
    #         dvydx = (data["y-velocity"][sl_right,sl_center,sl_center] -
    #                 data["y-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac*data["dx"].flat[0])
    #         dvxdy = (data["x-velocity"][sl_center,sl_right,sl_center] -
    #                 data["x-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac*data["dy"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvydx + dvxdy)**2.0
    #         del dvydx, dvxdy
    #     if data.pf.dimensionality > 2:
    #         dvzdy = (data["z-velocity"][sl_center,sl_right,sl_center] -
    #                 data["z-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac*data["dy"].flat[0])
    #         dvydz = (data["y-velocity"][sl_center,sl_center,sl_right] -
    #                 data["y-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac*data["dz"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvzdy + dvydz)**2.0
    #         del dvzdy, dvydz
    #         dvxdz = (data["x-velocity"][sl_center,sl_center,sl_right] -
    #                 data["x-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac*data["dz"].flat[0])
    #         dvzdx = (data["z-velocity"][sl_right,sl_center,sl_center] -
    #                 data["z-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac*data["dx"].flat[0])
    #         new_field[sl_center,sl_center,sl_center] += (dvxdz + dvzdx)**2.0
    #         del dvxdz, dvzdx
    #     new_field /= data["SoundSpeed"]**2.0
    #     new_field = new_field**(0.5)
    #     new_field = np.abs(new_field)
    #     return new_field

    # def _convertShearCriterion(data):
    #     return data.convert("cm")**-1.0
    # #add_field("ShearCriterion", function=_ShearCriterion,
    # #          validators=[ValidateSpatial(1,
    # #              ["x-velocity","y-velocity","z-velocity", "SoundSpeed"])],
    # #          units=r"\rm{cm}^{-1}",
    # #          convert_function=_convertShearCriterion, take_log=False)

    # def _ShearMach(field, data):
    #     """
    #     Dimensionless Shear (ShearMach) is defined nearly the same as shear, 
    #     except that it is scaled by the local dx/dy/dz and the local sound speed.
    #     So it results in a unitless quantity that is effectively measuring 
    #     shear in mach number.  

    #     In order to avoid discontinuities created by multiplying by dx/dy/dz at
    #     grid refinement boundaries, we also multiply by 2**GridLevel.

    #     Shear (Mach) = [(dvx + dvy)^2 + (dvz + dvy)^2 +
    #                     (dvx + dvz)^2  ]^(0.5) / c_sound
    #     """
    #     # We need to set up stencils
    #     if data.pf["HydroMethod"] == 2:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(1,-1,None)
    #         div_fac = 1.0
    #     else:
    #         sl_left = slice(None,-2,None)
    #         sl_right = slice(2,None,None)
    #         div_fac = 2.0
    #     new_field = np.zeros(data["x-velocity"].shape)
    #     if data.pf.dimensionality > 1:
    #         dvydx = (data["y-velocity"][sl_right,sl_center,sl_center] -
    #                 data["y-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac)
    #         dvxdy = (data["x-velocity"][sl_center,sl_right,sl_center] -
    #                 data["x-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac)
    #         new_field[sl_center,sl_center,sl_center] += (dvydx + dvxdy)**2.0
    #         del dvydx, dvxdy
    #     if data.pf.dimensionality > 2:
    #         dvzdy = (data["z-velocity"][sl_center,sl_right,sl_center] -
    #                 data["z-velocity"][sl_center,sl_left,sl_center]) \
    #                 / (div_fac)
    #         dvydz = (data["y-velocity"][sl_center,sl_center,sl_right] -
    #                 data["y-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac)
    #         new_field[sl_center,sl_center,sl_center] += (dvzdy + dvydz)**2.0
    #         del dvzdy, dvydz
    #         dvxdz = (data["x-velocity"][sl_center,sl_center,sl_right] -
    #                 data["x-velocity"][sl_center,sl_center,sl_left]) \
    #                 / (div_fac)
    #         dvzdx = (data["z-velocity"][sl_right,sl_center,sl_center] -
    #                 data["z-velocity"][sl_left,sl_center,sl_center]) \
    #                 / (div_fac)
    #         new_field[sl_center,sl_center,sl_center] += (dvxdz + dvzdx)**2.0
    #         del dvxdz, dvzdx
    #     new_field *= ((2.0**data.level)/data["SoundSpeed"])**2.0
    #     new_field = new_field**0.5
    #     new_field = np.abs(new_field)
    #     return new_field
    # #add_field("ShearMach", function=_ShearMach,
    # #          validators=[ValidateSpatial(1,
    # #              ["x-velocity","y-velocity","z-velocity","SoundSpeed"])],
    # #          units=r"\rm{Mach}",take_log=False)

