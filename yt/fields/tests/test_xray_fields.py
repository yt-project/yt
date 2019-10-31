from collections import OrderedDict

import pytest

from yt.fields.xray_emission_fields import \
    add_xray_emissivity_field
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils


sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


def setup():
    from yt.config import ytcfg
    ytcfg["yt","__withintesting"] = "True"



@pytest.mark.answer_test
@pytest.mark.big_data
@pytest.mark.usefixtures('answer_file')
class TestXRayFields(fw.AnswerTest):
    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(sloshing)
    def test_sloshing_apec(self):
        ds = utils.data_dir_load(sloshing)
        fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="apec", 
                                           metallicity=0.3)
        dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
        self.hashes['projection_values'] = OrderedDict()
        self.hashes['field_values'] = OrderedDict()
        for field in fields:
            self.hashes['projection_values'][field] = OrderedDict()
            self.hashes['field_values'][field] = OrderedDict()
            for dobj_name in dso:
                self.hashes['projection_values'][field][dobj_name] = OrderedDict()
                for axis in [0, 1, 2]:
                    pv_hd = self.projection_values_test(ds, axis, field, None, dobj_name)
                    self.hashes['projection_values'][field][dobj_name][axis] = pv_hd
                fv_hd = self.field_values_test(ds, field, dobj_name)
                self.hashes['field_values'][field][dobj_name] = fv_hd

    @pytest.mark.usefixtures('hashing')
    @utils.requires_ds(d9p)
    def test_d9p_cloudy(self):
        ds = utils.data_dir_load(d9p)
        fields = add_xray_emissivity_field(ds, 0.5, 2.0, redshift=ds.current_redshift,
                                           table_type="cloudy", cosmology=ds.cosmology,
                                           metallicity=("gas", "metallicity"))
        dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
        self.hashes['projection_values'] = OrderedDict()
        self.hashes['field_values'] = OrderedDict()
        for field in fields:
            self.hashes['projection_values'][field] = OrderedDict()
            self.hashes['field_values'][field] = OrderedDict()
            for dobj_name in dso:
                self.hashes['projection_values'][field][dobj_name] = OrderedDict()
                for axis in [0, 1, 2]:
                    pv_hd = self.projection_values_test(ds, axis, field, None, dobj_name)
                    self.hashes['projection_values'][field][dobj_name][axis] = pv_hd
                fv_hd = self.field_values_test(ds, field, dobj_name)
                self.hashes['field_values'][field][dobj_name] = fv_hd
