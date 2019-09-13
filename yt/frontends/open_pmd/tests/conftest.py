"""
Title: conftest.py
Purpose: Contains fixtures for loading data.
"""
import pytest

from yt.utilities.answer_testing import utils


# Test data
twoD = "example-2d/hdf5/data00000100.h5"
threeD = "example-3d/hdf5/data00000100.h5"
noFields = "no_fields/data00000400.h5"
noParticles = "no_particles/data00000400.h5"


@pytest.fixture(scope='class')
def ds_threeD():
    ds = utils.data_dir_load(threeD)
    assert str(ds) == "data00000100.h5"
    return ds

@pytest.fixture(scope='class')
def ds_twoD():
    ds = utils.data_dir_load(twoD)
    assert str(ds) == "data00000100.h5"
    return ds

@pytest.fixture(scope='class')
def ds_noFields():
    ds = utils.data_dir_load(noFields)
    assert str(ds) == "data00000400.h5"
    return ds

@pytest.fixture(scope='class')
def ds_noParticles():
    ds = utils.data_dir_load(noParticles)
    assert str(ds) == "data00000400.h5"
    return ds
