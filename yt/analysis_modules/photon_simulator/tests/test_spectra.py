from yt.analysis_modules.photon_simulator.api import \
    TableApecModel, XSpecThermalModel
from yt.testing import requires_module, fake_random_ds
from yt.utilities.answer_testing.framework import \
    GenericArrayTest
from yt.config import ytcfg

def setup():
    ytcfg["yt", "__withintesting"] = "True"

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

ds = fake_random_ds(64)

@requires_module("xspec")
@requires_module("astropy")
def test_apec():

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

    for test in [GenericArrayTest(ds, spec1_test),
                 GenericArrayTest(ds, spec2_test)]:
        test_apec.__name__ = test.description
        yield test

    xmod.cleanup_spectrum()
