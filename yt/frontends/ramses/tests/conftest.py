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
