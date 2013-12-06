"""
Halo Catalog-specific fields




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.funcs import mylog
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.data_objects.yt_array import \
    YTArray

from yt.utilities.physical_constants import \
    mh, \
    mass_sun_cgs

import yt.utilities.lib as amr_utils

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_mass / code_length**3"
vel_units = "code_length / code_time"

class RockstarFieldInfo(FieldInfoContainer):
    known_other_fields = (
    )

    known_particle_fields = (
        ("particle_identifier", ("", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_mposition_x", ("code_length", [], None)),
        ("particle_mposition_y", ("code_length", [], None)),
        ("particle_mposition_z", ("code_length", [], None)),
        ("particle_velocity_x", ("km / s", [], None)),
        ("particle_velocity_y", ("km / s", [], None)),
        ("particle_velocity_z", ("km / s", [], None)),
        ("particle_bvelocity_x", ("km / s", [], None)),
        ("particle_bvelocity_y", ("km / s", [], None)),
        ("particle_bvelocity_z", ("km / s", [], None)),
        ("particle_mass", ("Msun / h", [], None)),
        ("particle_radius", ("code_length", [], None)),
        ("child_radius", ("code_length", [], None)),
        ("vmax_r", ("km / s", [], None)),
    # These fields I don't have good definitions for yet.
    ('mgrav', ("", [], None)),
    ('vmax', ("", [], None)),
    ('rvmax', ("", [], None)),
    ('rs', ("", [], None)),
    ('klypin_rs', ("", [], None)),
    ('vrms', ("", [], None)),
    ('JX', ("", [], None)),
    ('JY', ("", [], None)),
    ('JZ', ("", [], None)),
    ('energy', ("", [], None)),
    ('spin', ("", [], None)),
    ('alt_m1', ("", [], None)),
    ('alt_m2', ("", [], None)),
    ('alt_m3', ("", [], None)),
    ('alt_m4', ("", [], None)),
    ('Xoff', ("", [], None)),
    ('Voff', ("", [], None)),
    ('b_to_a', ("", [], None)),
    ('c_to_a', ("", [], None)),
    ('Ax', ("", [], None)),
    ('Ay', ("", [], None)),
    ('Az', ("", [], None)),
    ('bullock_spin', ("", [], None)),
    ('kin_to_pot', ("", [], None)),
    ('num_p', ("", [], None)),
    ('num_child_particles', ("", [], None)),
    ('p_start', ("", [], None)),
    ('desc', ("", [], None)),
    ('flags', ("", [], None)),
    ('n_core', ("", [], None)),
    ('min_pos_err', ("", [], None)),
    ('min_vel_err', ("", [], None)),
    ('min_bulkvel_err', ("", [], None)),
)
