"""
halo_catalog frontend tests



"""

#-----------------------------------------------------------------------------
# Copyright (c) yt Development Team. All rights reserved.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.convenience import \
    load as yt_load
from yt.frontends.halo_catalog.data_structures import \
    HaloCatalogDataset
from yt.frontends.ytdata.utilities import \
    save_as_dataset
from yt.testing import \
    assert_array_equal, \
    requires_module, \
    TempDirTest
from yt.units.yt_array import \
    YTArray, \
    YTQuantity

def fake_halo_catalog(data):
    filename = "catalog.0.h5"

    ftypes = dict((field, '.') for field in data)
    extra_attrs = {"data_type": "halo_catalog",
                   "num_halos": data['particle_mass'].size}

    ds = {'cosmological_simulation': 1,
          'omega_lambda': 0.7,
          'omega_matter': 0.3,
          'hubble_constant': 0.7,
          'current_redshift': 0,
          'current_time': YTQuantity(1, 'yr'),
          'domain_left_edge': YTArray(np.zeros(3), 'cm'),
          'domain_right_edge': YTArray(np.ones(3), 'cm')}
    save_as_dataset(
        ds, filename, data,
        field_types=ftypes, extra_attrs=extra_attrs)
    return filename

class HaloCatalogTest(TempDirTest):
    @requires_module('h5py')
    def test_halo_catalog(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = ['particle_%s' % name for name in
                  ['mass'] + ['position_%s' % ax for ax in 'xyz']]
        units = ['g'] + ['cm']*3
        data = dict((field, YTArray(rs.random_sample(n_halos), unit))
                    for field, unit in zip(fields, units))

        fn = fake_halo_catalog(data)
        ds = yt_load(fn)

        assert isinstance(ds, HaloCatalogDataset)

        for field in fields:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)

    @requires_module('h5py')
    def test_halo_catalog_boundary_particles(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = ['particle_%s' % name for name in
                  ['mass'] + ['position_%s' % ax for ax in 'xyz']]
        units = ['g'] + ['cm']*3
        data = dict((field, YTArray(rs.random_sample(n_halos), unit))
                    for field, unit in zip(fields, units))

        data['particle_position_x'][0] = 1.0
        data['particle_position_x'][1] = 0.0
        data['particle_position_y'][2] = 1.0
        data['particle_position_y'][3] = 0.0
        data['particle_position_z'][4] = 1.0
        data['particle_position_z'][5] = 0.0

        fn = fake_halo_catalog(data)
        ds = yt_load(fn)

        assert isinstance(ds, HaloCatalogDataset)

        for field in ['particle_mass']:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)
