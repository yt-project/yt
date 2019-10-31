from collections import OrderedDict

import pytest

from yt.fields.xray_emission_fields import \
    add_xray_emissivity_field
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


@pytest.mark.answer_test
@pytest.mark.big_data
@pytest.mark.usefixtures('answer_file')
class TestXRayFields(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(sloshing)
    def test_sloshing_apec(self, field, dobj_name, axis):
        ds = utils.data_dir_load(sloshing)
        pv_hd = self.projection_values_test(ds, axis, field, None, dobj_name)
        self.hashes.update({'projection_values' : pv_hd})
        fv_hd = self.field_values_test(ds, field, dobj_name)
        self.hashes.update({'field_values' : fv_hd})

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(d9p)
    def test_d9p_cloudy(self, field, dobj_name, axis):
        ds = utils.data_dir_load(d9p)
        pv_hd = self.projection_values_test(ds, axis, field, None, dobj_name)
        self.hashes.update({'projection_values' : pv_hd})
        fv_hd = self.field_values_test(ds, field, dobj_name)
        self.hashes.update({'field_values' : fv_hd}
