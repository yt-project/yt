"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
isothermal_h5 = "IsothermalCollapse/snap_505.hdf5"

iso_kwargs = dict(bounding_box=[[-3, 3], [-3, 3], [-3, 3]])

fields = {
        ("gas", "density"): None,
        ("gas", "temperature"): None,
        ("gas", "temperature"): ('gas', 'density'),
        ('gas', 'velocity_magnitude'): None,
        ("deposit", "all_density"): None,
        ("deposit", "all_count"): None,
        ("deposit", "all_cic"): None,
        ("deposit", "PartType0_density"): None
    }


def pytest_generate_tests(metafunc):
    if metafunc.function.__name__ == 'test_iso_collapse':
        metafunc.parametrize('f, w', [(k, v) for k, v in fields.items()],
            ids=['dens', 'temp-None', 'temp-dens', 'velocity_magnitude', 'all_density',
            'all_count', 'all_cic', 'PartType0Den'])


@pytest.fixture(scope='class')
def ds_isothermal_h5():
    ds = utils.data_dir_load(isothermal_h5, kwargs=iso_kwargs)
    return ds
