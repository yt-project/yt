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

from yt.units.yt_array import YTArray
from yt.utilities.lib import obtain_rvec, obtain_rv_vec
from yt.utilities.math_utils import resize_vector
from yt.utilities.cosmology import Cosmology
from yt.fields.derived_field import \
    ValidateGridType, \
    ValidateParameter, \
    ValidateSpatial, \
    NeedsParameter
from yt.fields.field_info_container import \
    FieldInfoContainer

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

UniversalFields = FieldInfoContainer(None, [])
add_field = UniversalFields.add_field

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
    new_field2 = data.pf.arr(np.zeros((nx, ny, nz)), 'g/cm**3')
    new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
    return new_field2

add_field("averaged_density", function=_averaged_density,
          validators=[ValidateSpatial(1, ["density"])])

