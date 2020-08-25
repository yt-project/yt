import pytest

from yt.utilities.answer_testing.answer_tests import field_values, projection_values

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


@pytest.mark.answer_test
@pytest.mark.big_data
class TestXRayFields:
    answer_file = None
    saved_hashes = None

    @pytest.mark.usefixtures("hashing")
    def test_sloshing_apec(self, ds, f, d, a):
        pv = projection_values(ds, a, f, None, d)
        self.hashes.update({"projection_values": pv})
        fv = field_values(ds, f, d)
        self.hashes.update({"field_values": fv})

    @pytest.mark.usefixtures("hashing")
    def test_d9p_cloudy(self, ds, f, d, a):
        pv = projection_values(ds, a, f, None, d)
        self.hashes.update({"projection_values": pv})
        fv = field_values(ds, f, d)
        self.hashes.update({"field_values": fv})
