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



@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
@pytest.mark.skipif(not pytest.config.getvalue('--answer-big-data'),
    reason="--answer-big-data not set.")
@pytest.mark.usefixtures('answer_file')
class TestXRayFields(fw.AnswerTest):
    @utils.requires_ds(sloshing)
    def test_sloshing_apec(self):
        ds = utils.data_dir_load(sloshing)
        fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="apec", 
                                           metallicity=0.3)
        dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
        hd = OrderedDict()
        hd['projection_values'] = OrderedDict()
        hd['field_values'] = OrderedDict()
        for field in fields:
            hd['projection_values'][field] = OrderedDict()
            hd['field_values'][field] = OrderedDict()
            for dobj_name in dso:
                hd['projection_values'][field][dobj_name] = OrderedDict()
                for axis in [0, 1, 2]:
                    pv_hd = utils.generate_hash(
                        self.projection_values_test(ds, axis, field, 
                                               None, dobj_name)
                    )
                    hd['projection_values'][field][dobj_name][axis] = pv_hd
                fv_hd = utils.generate_hash(
                    self.field_values_test(ds, field, dobj_name)
                )
                hd['field_values'][field][dobj_name] = fv_hd
        hashes = {'sloshing_apec' : hd}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)

    @utils.requires_ds(d9p)
    def test_d9p_cloudy(self):
        ds = utils.data_dir_load(d9p)
        fields = add_xray_emissivity_field(ds, 0.5, 2.0, redshift=ds.current_redshift,
                                           table_type="cloudy", cosmology=ds.cosmology,
                                           metallicity=("gas", "metallicity"))
        dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
        hd = OrderedDict()
        hd['projection_values'] = OrderedDict()
        hd['field_values'] = OrderedDict()
        for field in fields:
            hd['projection_values'][field] = OrderedDict()
            hd['field_values'][field] = OrderedDict()
            for dobj_name in dso:
                hd['projection_values'][field][dobj_name] = OrderedDict()
                for axis in [0, 1, 2]:
                    pv_hd = utils.generate_hash(
                        self.projection_values_test(ds, axis, field, 
                                               None, dobj_name)
                    )
                    hd['projection_values'][field][dobj_name][axis] = pv_hd
                fv_hd = utils.generate_hash(
                    self.field_values_test(ds, field, dobj_name)
                )
                hd['field_values'][field][dobj_name] = fv_hd
        hashes = {'d9p_cloudy' : hd}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)
