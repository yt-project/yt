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
from yt.units import km, s

def random_unit_vector(prng):
    v = prng.random_sample(3)
    while (v == 0).all():
        v = prng.random_sample(3)
    return v / np.sqrt((v**2).sum())

def test_vector_component_conversions():
    prng = np.random.RandomState(8675309)
    v_rand = np.power(10, prng.random_sample(3) * 5)
    normals = [[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1],
               random_unit_vector(prng),
               random_unit_vector(prng)]

    for bv in [[0, 0, 0], [1e5, 1e5, 1e5], v_rand]:
        for normal in normals:
            ds = fake_random_ds(16)

            ad = ds.all_data()

            bulk = bv*km/s

            ad.set_field_parameter('bulk_velocity', bulk)
            ad.set_field_parameter('normal', normal)

            vmag = ad['velocity_magnitude']
            vmag_cart = np.sqrt(
                (ad['velocity_x'] - bulk[0])**2 +
                (ad['velocity_y'] - bulk[1])**2 +
                (ad['velocity_z'] - bulk[2])**2
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

            for i, d in enumerate('xyz'):
                assert_allclose_units(ad['velocity_%s' % d] - bulk[i],
                                      ad['relative_velocity_%s' % d])
