#-----------------------------------------------------------------------------
# Copyright (c) 2016, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.testing import \
    assert_allclose_units, \
    fake_random_ds

def test_vector_component_conversions():
    ds = fake_random_ds(16)

    ad = ds.all_data()

    vmag = ad['velocity_magnitude']
    vmag_cart = np.sqrt(
        ad['velocity_x']**2 +
        ad['velocity_y']**2 +
        ad['velocity_z']**2
    )

    assert_allclose_units(vmag, vmag_cart)

    vmag_cyl = np.sqrt(
        ad['velocity_cylindrical_radius']**2 +
        ad['velocity_cylindrical_theta']**2 +
        ad['velocity_cylindrical_z']**2
    )

    assert_allclose_units(vmag, vmag_cyl)

    vmag_sph = np.sqrt(
        ad['velocity_spherical_radius']**2 +
        ad['velocity_spherical_theta']**2 +
        ad['velocity_spherical_phi']**2
    )

    assert_allclose_units(vmag, vmag_sph)
