"""
Title: test_owls.py
Purpose: OWLS frontend tests using the snapshot_033 dataset
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.
"""
from collections import OrderedDict

import pytest

from yt.testing import requires_file
from yt.frontends.owls.api import OWLSDataset
from yt.data_objects.particle_filters import add_particle_filter
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


# Test data
os33 = "snapshot_033/snap_033.0.hdf5"


# Globals
# This maps from field names to weight field names to use for projections
_fields = OrderedDict(
    [
        (("gas", "density"), None),
        (("gas", "temperature"), None),
        (("gas", "temperature"), ("gas", "density")),
        (('gas', 'He_p0_number_density'), None),
        (('gas', 'velocity_magnitude'), None),
        (("deposit", "all_density"), None),
        (("deposit", "all_count"), None),
        (("deposit", "all_cic"), None),
        (("deposit", "PartType0_density"), None),
        (("deposit", "PartType4_density"), None),
    ]
)


# Answer file
answer_file = 'owls_answers.yaml'


#============================================
#                 TestOwls
#============================================
@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestOwls(fw.AnswerTest):
    #-----
    # test_snapshot_033
    #-----
    @pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
        reason="Skipping test_jet because --answer-big-data was not set."
    )
    @utils.requires_ds(os33)
    def test_snapshot_033(self, ds_os33):
        hashes = self.sph_answer(ds_os33, 'snap_033', 2*128**3, _fields)
        hashes = {'snapshot_033' : hashes}
        utils.handle_hashes(self.save_dir, answer_file, hashes, self.answer_store)

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
