"""
These are common particle deposition fields.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2013 Matthew Turk.  All Rights Reserved.

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

import numpy as np

from yt.funcs import *
from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType, \
    NullFunc, \
    TranslationFunc
from yt.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_sun_cgs, \
    mh

def particle_deposition_functions(ptype, coord_name, mass_name, registry):
    orig = set(registry.keys())
    def particle_count(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, method = "count")
        return d

    registry.add_field(("deposit", "%s_count" % ptype),
             function = particle_count,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Count}" % ptype,
             projection_conversion = '1')

    def particle_mass(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        return d

    registry.add_field(("deposit", "%s_mass" % ptype),
             function = particle_mass,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Mass}" % ptype,
             units = r"\mathrm{g}",
             projected_units = r"\mathrm{g}\/\mathrm{cm}",
             projection_conversion = 'cm')

    def particle_density(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
        d /= data["gas","CellVolume"]
        return d

    registry.add_field(("deposit", "%s_density" % ptype),
             function = particle_density,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    def particle_cic(field, data):
        pos = data[ptype, coord_name]
        d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
        d /= data["gas","CellVolume"]
        return d

    registry.add_field(("deposit", "%s_cic" % ptype),
             function = particle_cic,
             validators = [ValidateSpatial()],
             display_name = "\\mathrm{%s CIC Density}" % ptype,
             units = r"\mathrm{g}/\mathrm{cm}^{3}",
             projected_units = r"\mathrm{g}/\mathrm{cm}^{2}",
             projection_conversion = 'cm')

    # Now some translation functions.

    registry.add_field((ptype, "ParticleMass"),
            function = TranslationFunc((ptype, mass_name)),
            particle_type = True,
            units = r"\mathrm{g}")

    def _ParticleMassMsun(field, data):
        return data[ptype, mass_name].copy()
    def _conv_Msun(data):
        return 1.0/mass_sun_cgs

    registry.add_field((ptype, "ParticleMassMsun"),
            function = _ParticleMassMsun,
            convert_function = _conv_Msun,
            particle_type = True,
            units = r"\mathrm{M}_\odot")

    def particle_mesh_ids(field, data):
        pos = data[ptype, coord_name]
        ids = np.zeros(pos.shape[0], dtype="float64") - 1
        # This is float64 in name only.  It will be properly cast inside the
        # deposit operation.
        #_ids = ids.view("float64")
        data.deposit(pos, [ids], method = "mesh_id")
        return ids
    registry.add_field((ptype, "mesh_id"),
            function = particle_mesh_ids,
            validators = [ValidateSpatial()],
            particle_type = True)

    return list(set(registry.keys()).difference(orig))


def particle_scalar_functions(ptype, coord_name, vel_name, registry):

    # Now we have to set up the various velocity and coordinate things.  In the
    # future, we'll actually invert this and use the 3-component items
    # elsewhere, and stop using these.
    
    # Note that we pass in _ptype here so that it's defined inside the closure.
    orig = set(registry.keys())

    def _get_coord_funcs(axi, _ptype):
        def _particle_velocity(field, data):
            return data[_ptype, vel_name][:,axi]
        def _particle_position(field, data):
            return data[_ptype, coord_name][:,axi]
        return _particle_velocity, _particle_position
    for axi, ax in enumerate("xyz"):
        v, p = _get_coord_funcs(axi, ptype)
        registry.add_field((ptype, "particle_velocity_%s" % ax),
            particle_type = True, function = v)
        registry.add_field((ptype, "particle_position_%s" % ax),
            particle_type = True, function = p)

    return list(set(registry.keys()).difference(orig))

def particle_vector_functions(ptype, coord_names, vel_names, registry):

    # This will column_stack a set of scalars to create vector fields.
    orig = set(registry.keys())

    def _get_vec_func(_ptype, names):
        def particle_vectors(field, data):
            return np.column_stack([data[_ptype, name] for name in names])
        return particle_vectors
    registry.add_field((ptype, "Coordinates"),
                       function=_get_vec_func(ptype, coord_names),
                       particle_type=True)
    registry.add_field((ptype, "Velocities"),
                       function=_get_vec_func(ptype, vel_names),
                       particle_type=True)

    return list(set(registry.keys()).difference(orig))
