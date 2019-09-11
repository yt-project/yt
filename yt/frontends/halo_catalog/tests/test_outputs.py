"""
Title: test_halo_catalog.py
Purpose: halo_catalog frontend tests
Notes:
    Copyright (c) yt Development Team. All rights reserved.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import numpy as np
import pytest

from yt.convenience import load
from yt.frontends.halo_catalog.data_structures import HaloCatalogDataset
from yt.frontends.ytdata.utilities import save_as_dataset
from yt.testing import \
    assert_array_equal, \
    requires_module
from yt.units.yt_array import YTArray, YTQuantity
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


#============================================
#              TestHaloCatalog
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.usefixtures('temp_dir')
class TestHaloCatalog(fw.AnswerTest):

    #-----
    # test_halo_caatalog
    #-----
    @requires_module('h5py')
    def test_halo_catalog(self):
        rs = np.random.RandomState(3670474)
        n_halos = 100
        fields = ['particle_%s' % name for name in
                  ['mass'] + ['position_%s' % ax for ax in 'xyz']]
        units = ['g'] + ['cm']*3
        data = dict((field, YTArray(rs.random_sample(n_halos), unit))
                    for field, unit in zip(fields, units))
        fn = utils.fake_halo_catalog(data)
        ds = load(fn)
        assert isinstance(ds, HaloCatalogDataset)
        for field in fields:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)

    #-----
    # test_halo_catalog_boundary_particles
    #-----
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
        fn = utils.fake_halo_catalog(data)
        ds = load(fn)
        assert isinstance(ds, HaloCatalogDataset)
        for field in ['particle_mass']:
            f1 = data[field].in_base()
            f1.sort()
            f2 = ds.r[field].in_base()
            f2.sort()
            assert_array_equal(f1, f2)
