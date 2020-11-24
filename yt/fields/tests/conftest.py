import pytest

from yt.fields.xray_emission_fields import add_xray_emissivity_field
from yt.utilities.answer_testing.testing_utilities import data_dir_load

sloshing = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0300"
d9p = "D9p_500/10MpcBox_HartGal_csf_a0.500.d"


dso = [None, ("sphere", ("m", (0.1, "unitary")))]
axes = [0, 1, 2]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == "test_sloshing_apec":
        try:
            ds = data_dir_load(sloshing)
        except FileNotFoundError:
            pytest.skip("Data not found.")
        fields = add_xray_emissivity_field(
            ds, 0.5, 7.0, table_type="apec", metallicity=0.3
        )
    if metafunc.function.__name__ == "test_d9p_cloudy":
        try:
            ds = data_dir_load(d9p)
        except FileNotFoundError:
            pytest.skip("Data not found.")
        fields = add_xray_emissivity_field(
            ds,
            0.5,
            2.0,
            redshift=ds.current_redshift,
            table_type="cloudy",
            cosmology=ds.cosmology,
            metallicity=("gas", "metallicity"),
        )
    if metafunc.function.__name__ in ["test_sloshing_apec", "test_d9p_cloudy"]:
        metafunc.parametrize(
            "ds, f",
            [(ds, f) for f in fields],
            ids=[f.__repr__() for f in fields],
            indirect=True,
        )
        metafunc.parametrize("d", dso, ids=["None", "sphere"], indirect=True)
        metafunc.parametrize("a", axes, ids=["0", "1", "2"], indirect=True)
