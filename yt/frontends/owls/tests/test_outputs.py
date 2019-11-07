"""
Title: test_owls.py
Purpose: OWLS frontend tests using the snapshot_033 dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
import pytest

from yt.testing import requires_file
from yt.frontends.owls.api import OWLSDataset
from yt.data_objects.particle_filters import add_particle_filter
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
os33 = "snapshot_033/snap_033.0.hdf5"




#============================================
#                 TestOwls
#============================================
@pytest.mark.answer_test
@pytest.mark.usefixtures('answer_file')
class TestOwls(fw.AnswerTest):
    #-----
    # test_snapshot_033
    #-----
    @pytest.mark.big_data
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(os33)
    def test_snapshot_033(self, f, w, d, a, ds_os33):
        self.hashes.update(self.sph_answer(ds_os33, 'snap_033', 2*128**3, f, w, d, a))

    #-----
    # test_OWLSDataset
    #-----
    @requires_file(os33)
    def test_OWLSDataset(self, ds_os33):
        assert isinstance(ds_os33, OWLSDataset)
    
    #-----
    # test_OWLS_particlefilter
    #----- 
    @utils.requires_ds(os33)
    def test_OWLS_particlefilter(self, ds_os33):
        ds = ds_os33
        ad = ds.all_data()
        def cold_gas(pfilter, data):
            temperature = data[pfilter.filtered_type, "Temperature"]
            filter = (temperature.in_units('K') <= 1e5)
            return filter
        add_particle_filter("gas_cold", function=cold_gas, filtered_type='PartType0', requires=["Temperature"])
        ds.add_particle_filter('gas_cold')
        mask = (ad['PartType0','Temperature'] <= 1e5)
        assert ad['PartType0','Temperature'][mask].shape == ad['gas_cold','Temperature'].shape
