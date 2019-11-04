import pytest

from yt.utilities.answer_testing import utils


dso = [ None, ("sphere", ("m", (0.1, 'unitary')))]
axes = [0, 1, 2]


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == "test_sloshing_apec":
        ds = utils.data_dir_load(sloshing)
        fields = add_xray_emissivity_field(ds, 0.5, 7.0, table_type="apec", 
                                           metallicity=0.3)
        metafunc.parametrize('field', fields, ids=[f.__repr__ for f in fields])
    if metafunc.function.__name__ == "test_d9p_cloudy":
        ds = utils.data_dir_load(d9p)
        fields = add_xray_emissivity_field(ds, 0.5, 2.0, redshift=ds.current_redshift,
                                           table_type="cloudy", cosmology=ds.cosmology,
                                           metallicity=("gas", "metallicity"))
    if metafunc.function.__name__ in ["test_sloshing_apec", "test_d9p_cloudy"]:
        metafunc.parametrize('field', fields, ids=[f.__repr__ for f in fields])
        metafunc.parametrize('dobj_name', dso, ids=['None', 'sphere'])
        metafunc.parametrize('axis', axes, ids=['0', '1', '2'])
