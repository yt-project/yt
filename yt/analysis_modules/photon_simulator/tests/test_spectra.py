from collections import OrderedDict

import pytest

from yt.analysis_modules.photon_simulator.api import \
    TableApecModel, XSpecThermalModel
from yt.config import ytcfg
from yt.testing import requires_module, fake_random_ds
import yt.utilities.answer_testing.framework as fw
from yt.utilities.answer_testing import utils

def setup():
    ytcfg["yt", "__withintesting"] = "True"

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

ds = fake_random_ds(64)


@pytest.mark.skipif(not pytest.config.getvalue('--with-answer-testing'),
    reason="--with-answer-testing not set.")
class TestSpectra(fw.AnswerTest):
    @requires_module("xspec")
    @requires_module("astropy")
    def test_apec(self):
        settings = {"APECROOT":xray_data_dir+"/apec_v2.0.2"}
        xmod = XSpecThermalModel("apec", 0.1, 10.0, 10000, thermal_broad=True,
                                 settings=settings)
        xmod.prepare_spectrum(0.2)
        xcspec, xmspec = xmod.get_spectrum(6.0)
        spec1 = xcspec+0.3*xmspec
        amod = TableApecModel(xray_data_dir, 0.1, 10.0, 
                              10000, thermal_broad=True)
        amod.prepare_spectrum(0.2)
        acspec, amspec = amod.get_spectrum(6.0)
        spec2 = acspec+0.3*amspec
        def spec1_test():
            return spec1.v
        def spec2_test():
            return spec2.v
        ga_hd1 = utils.generate_hash(self.generic_array_test(ds, spec1_test))
        ga_hd2 = utils.generate_hash(self.generic_array_test(ds, spec2_test))
        hashes = OrderedDict()
        hashes['generic_array1'] = ga_hd1
        hashes['generic_array2'] = ga_hd2
        hashes = {'apec' : hashes}
        utils.handle_hashes(self.save_dir, self.answer_file, hashes, self.answer_store)
        xmod.cleanup_spectrum()
