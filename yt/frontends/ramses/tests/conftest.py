"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
output_00080 = "output_00080/info_00080.txt"
ramsesNonCosmo = 'DICEGalaxyDisk_nonCosmological/output_00002/info_00002.txt'
ramsesExtraFieldsSmall = 'ramses_extra_fields_small/output_00001'
ramses_rt = "ramses_rt_00088/output_00088/info_00088.txt"
ramses_sink = "ramses_sink_00016/output_00016/info_00016.txt"
ramses_new_format = "ramses_new_format/output_00002/info_00002.txt"
ramses_empty_record = "ramses_empty_record/output_00003/info_00003.txt"


# Test parameters. Format:
# {test1: {param1 : [(val1, val2,...), (id1, id2,...)], param2 : ...}, test2: ...}
test_params = {
    'test_output_00080' : {
        'a' : [(0, 1, 2,), ('0', '1', '2')],
        'd' : [(None, ('sphere', ('max', (0.1, 'unitary')))), ('None', 'sphere')],
        'w' : [(None, 'density'), ('None', 'dneisty')],
        'f' : [('temperature', 'density', 'velocity_magnitude',
                ('deposit', 'all_density'), ('deposit', 'all_count')),
                ('temperature', 'density', 'velocity_magnitude', 'dep_all_density',
                'dep_all_count')]
    }
}


def pytest_generate_tests(metafunc):
    # Loop over each test in test_params
    for test_name, params in test_params.items():
        if metafunc.function.__name__ == test_name:
            # Parametrize
            for param_name, param_vals in params.items():
                metafunc.parametrize(param_name, param_vals[0], ids=param_vals[1])


@pytest.fixture(scope='class')
def ds_output_00080():
    ds = utils.data_dir_load(output_00080)
    assert str(ds) == "info_00080"
    return ds

@pytest.fixture(scope='class')
def ds_ramses_rt():
    ds = utils.data_dir_load(ramses_rt)
    return ds

@pytest.fixture(scope='class')
def ds_ramses_sink():
    ds = utils.data_dir_load(ramses_sink)
    return ds

@pytest.fixture(scope='class')
def ds_ramses_new_format():
    ds = utils.data_dir_load(ramses_new_format)
    return ds
