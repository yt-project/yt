"""
Here are some fields that are specific to fluids.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import \
    just_one

from .derived_field import \
    ValidateParameter, \
    ValidateSpatial

from .field_plugin_registry import \
    register_field_plugin

from .vector_operations import \
     create_averaged_field, \
     create_magnitude_field, \
     create_vector_fields

from yt.utilities.physical_constants import \
    mh, \
    kboltz

from yt.utilities.physical_ratios import \
    metallicity_sun

from yt.units.yt_array import \
    YTArray

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

@register_field_plugin
def setup_fluid_fields(registry, ftype = "gas", slice_info = None):
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

    create_vector_fields(registry, "velocity", "cm / s", ftype, slice_info)

    def _cell_mass(field, data):
        return data[ftype, "density"] * data["index", "cell_volume"]

    registry.add_field((ftype, "cell_mass"),
        function=_cell_mass,
        units="g")

    def _sound_speed(field, data):
        tr = data.ds.gamma * data[ftype, "pressure"] / data[ftype, "density"]
        return np.sqrt(tr)
    registry.add_field((ftype, "sound_speed"),
             function=_sound_speed,
             units="cm/s")

    def _radial_mach_number(field, data):
        """ Radial component of M{|v|/c_sound} """
        tr = data[ftype, "radial_velocity"] / data[ftype, "sound_speed"]
        return np.abs(tr)
    registry.add_field((ftype, "radial_mach_number"),
             function=_radial_mach_number,
             units = "")

    def _kin_energy(field, data):
        return 0.5*data[ftype, "density"] * ( data[ftype, "velocity_x"]**2.0
                                              + data[ftype, "velocity_y"]**2.0
                                              + data[ftype, "velocity_z"]**2.0 )
    registry.add_field(("gas", "kinetic_energy"),
                       function = _kin_energy,
                       units = "erg / cm**3")

    def _mach_number(field, data):
        """ M{|v|/c_sound} """
        return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]
    registry.add_field((ftype, "mach_number"),
            function=_mach_number,
            units = "")

    def _courant_time_step(field, data):
        t1 = data["index", "dx"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_x"]))
        t2 = data["index", "dy"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_y"]))
        t3 = data["index", "dz"] / (data[ftype, "sound_speed"]
                        + np.abs(data[ftype, "velocity_z"]))
        tr = np.minimum(np.minimum(t1, t2), t3)
        return tr

    registry.add_field((ftype, "courant_time_step"),
             function=_courant_time_step,
             units="s")

    def _pressure(field, data):
        """ M{(Gamma-1.0)*rho*E} """
        tr = (data.ds.gamma - 1.0) \
           * (data[ftype, "density"] * data[ftype, "thermal_energy"])
        return tr

    registry.add_field((ftype, "pressure"),
             function=_pressure,
             units="dyne/cm**2")

    def _kT(field, data):
        return (kboltz*data[ftype, "temperature"]).in_units("keV")
    registry.add_field((ftype, "kT"),
                       function=_kT,
                       units="keV",
                       display_name="Temperature")

    def _entropy(field, data):
        mw = data.get_field_parameter("mu")
        if mw is None:
            mw = 1.0
        mw *= mh
        gammam1 = 2./3.
        tr = data[ftype,"kT"] / ((data[ftype, "density"]/mw)**gammam1)
        return data.apply_units(tr, field.units)
    registry.add_field((ftype, "entropy"),
             units="keV*cm**2",
             function=_entropy)

    def _metallicity(field, data):
        tr = data[ftype, "metal_density"] / data[ftype, "density"]
        tr /= metallicity_sun
        return data.apply_units(tr, "Zsun")
    registry.add_field((ftype, "metallicity"),
             function=_metallicity,
             units="Zsun")

    def _metal_mass(field, data):
        return data[ftype, "metal_density"] * data["index", "cell_volume"]
    registry.add_field((ftype, "metal_mass"),
                       function=_metal_mass,
                       units="g")

    def _number_density(field, data):
        field_data = np.zeros_like(data["gas", "%s_number_density" % \
                                        data.ds.field_info.species_names[0]])
        for species in data.ds.field_info.species_names:
            field_data += data["gas", "%s_number_density" % species]
        return field_data
    registry.add_field((ftype, "number_density"),
                       function = _number_density,
                       units="cm**-3")
    
    def _mean_molecular_weight(field, data):
        return (data[ftype, "density"] / (mh * data[ftype, "number_density"]))
    registry.add_field((ftype, "mean_molecular_weight"),
              function=_mean_molecular_weight,
              units="")

    setup_gradient_fields(registry, (ftype, "pressure"), "dyne/cm**2",
                          slice_info)

    setup_gradient_fields(registry, (ftype, "density"), "g / cm**3",
                          slice_info)

    create_averaged_field(registry, "density", "g/cm**3",
                          ftype=ftype, slice_info=slice_info,
                          weight="cell_mass")

def setup_gradient_fields(registry, grad_field, field_units, slice_info = None):
    assert(isinstance(grad_field, tuple))
    ftype, fname = grad_field
    if slice_info is None:
        sl_left = slice(None, -2, None)
        sl_right = slice(2, None, None)
        div_fac = 2.0
    else:
        sl_left, sl_right, div_fac = slice_info
    slice_3d = [slice(1, -1), slice(1, -1), slice(1, -1)]

    def grad_func(axi, ax):
        slice_3dl = slice_3d[:]
        slice_3dr = slice_3d[:]
        slice_3dl[axi] = sl_left
        slice_3dr[axi] = sl_right
        def func(field, data):
            ds = div_fac * data["index", "dx"]
            f  = data[grad_field][slice_3dr]/ds[slice_3d]
            f -= data[grad_field][slice_3dl]/ds[slice_3d]
            new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                    f.units)
            new_field[slice_3d] = f
            return new_field
        return func

    grad_units = "(%s) / cm" % field_units
    for axi, ax in enumerate('xyz'):
        f = grad_func(axi, ax)
        registry.add_field((ftype, "%s_gradient_%s" % (fname, ax)),
                           function = f,
                           validators = [ValidateSpatial(1, [grad_field])],
                           units = grad_units)
    create_magnitude_field(registry, "%s_gradient" % fname,
                           grad_units, ftype=ftype,
                           validators = [ValidateSpatial(1, [grad_field])])
